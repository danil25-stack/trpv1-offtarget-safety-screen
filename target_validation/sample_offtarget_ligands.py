"""Sample a property-matched active/inactive ligand set per off-target
(RARG, MMP3) for docking-based discrimination validation.

Input: data/processed/skin_offtarget_screen/{gene}_chembl_actives_inactives.csv
(from fetch_chembl_bioactivity.py). Output:
data/processed/skin_offtarget_screen/{gene}_docking_validation_set.csv

Why matched, not just "top N most potent actives + N random inactives":
if actives and inactives differ systematically in MW/logP (plausible --
mature SAR series cluster tighter in chemical space than a target's
historical inactive/discontinued compounds), a later docking-score
discrimination check could just be reflecting that property bias rather
than real pocket complementarity. See the caveat already flagged in
chembl_bioactivity_offtarget.md.

Method:
1. Pull full_mwt + alogp per candidate molecule from ChEMBL's /molecule
   endpoint (these are ChEMBL-computed properties, not locally computed --
   no cheminformatics tooling runs on this machine, per project convention
   that all such compute belongs on the vast.ai instance).
2. Actives: pick N spread across potency deciles (not just the N most
   potent) so the sample isn't one tight chemotype cluster.
3. Inactives: greedy nearest-neighbor match to the selected actives in
   standardized (MW, logP) space, without replacement.
4. Crude salt stripping (split canonical_smiles on '.', keep the longest
   fragment by character count) -- string-level, not a chemistry
   computation, since the real embedding happens downstream with RDKit on
   vast.ai.

Usage:
    python sample_offtarget_ligands.py [--n 30]
"""
import argparse
import csv
import math
import time
from pathlib import Path

import requests

DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "processed" / "skin_offtarget_screen"
MOLECULE_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule"
CHUNK = 10
INACTIVE_LABELS = {"inactive", "inactive_no_binding_detected"}


def get_json(url, params, retries=12, backoff=3):
    last_exc = None
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=90)
            if r.status_code == 200 and r.content:
                return r.json()
            last_exc = RuntimeError(f"HTTP {r.status_code}, {len(r.content)} bytes")
        except requests.RequestException as e:
            last_exc = e
        wait = min(backoff * (attempt + 1), 20)
        print(f"    retry {attempt + 1}/{retries} after {last_exc} (sleeping {wait}s)")
        time.sleep(wait)
    raise RuntimeError(f"giving up after {retries} retries: {last_exc}")


def fetch_properties(molecule_ids):
    """molecule_chembl_id -> (full_mwt, alogp), skipping molecules with either missing."""
    props = {}
    ids = sorted(set(molecule_ids))
    for start in range(0, len(ids), CHUNK):
        chunk = ids[start:start + CHUNK]
        print(f"  molecule properties {start}-{start + len(chunk)}/{len(ids)}")
        data = get_json(MOLECULE_URL, {
            "molecule_chembl_id__in": ",".join(chunk),
            "limit": CHUNK,
            "format": "json",
        })
        for m in data["molecules"]:
            mp = m.get("molecule_properties") or {}
            mwt, alogp = mp.get("full_mwt"), mp.get("alogp")
            if mwt is None or alogp is None:
                continue
            try:
                props[m["molecule_chembl_id"]] = (float(mwt), float(alogp))
            except (TypeError, ValueError):
                continue
        time.sleep(0.3)
    return props


def strip_salt(smiles):
    parts = smiles.split(".")
    return max(parts, key=len)


def load_rows(gene):
    path = DATA_DIR / f"{gene}_chembl_actives_inactives.csv"
    with open(path) as f:
        return list(csv.DictReader(f))


def pick_actives(active_rows, n):
    """Spread the sample across potency deciles of the active class rather
    than just taking the N most potent (which would all be one tight
    chemotype cluster in a mature SAR series)."""
    active_rows = sorted(active_rows, key=lambda r: float(r["standard_value_nM"]))
    if len(active_rows) <= n:
        return active_rows
    picked = []
    for i in range(n):
        idx = round(i * (len(active_rows) - 1) / (n - 1))
        picked.append(active_rows[idx])
    # dedupe (rounding can collide at the ends for small pools)
    seen = set()
    out = []
    for r in picked:
        if r["molecule_chembl_id"] not in seen:
            seen.add(r["molecule_chembl_id"])
            out.append(r)
    return out


def match_inactives(actives, inactive_rows, props):
    """Greedy nearest-neighbor match (standardized MW, logP), no replacement."""
    mwts = [props[r["molecule_chembl_id"]][0] for r in inactive_rows]
    alogps = [props[r["molecule_chembl_id"]][1] for r in inactive_rows]
    mwt_mean, mwt_sd = sum(mwts) / len(mwts), (sum((x - sum(mwts) / len(mwts)) ** 2 for x in mwts) / len(mwts)) ** 0.5 or 1.0
    alogp_mean, alogp_sd = sum(alogps) / len(alogps), (sum((x - sum(alogps) / len(alogps)) ** 2 for x in alogps) / len(alogps)) ** 0.5 or 1.0

    def z(mid):
        mwt, alogp = props[mid]
        return (mwt - mwt_mean) / mwt_sd, (alogp - alogp_mean) / alogp_sd

    pool = list(inactive_rows)
    matched = []
    for a in actives:
        az_mwt, az_alogp = z(a["molecule_chembl_id"])
        best_idx, best_dist = None, math.inf
        for i, cand in enumerate(pool):
            cz_mwt, cz_alogp = z(cand["molecule_chembl_id"])
            dist = (az_mwt - cz_mwt) ** 2 + (az_alogp - cz_alogp) ** 2
            if dist < best_dist:
                best_dist, best_idx = dist, i
        if best_idx is None:
            break
        matched.append({**pool.pop(best_idx), "matched_active_id": a["molecule_chembl_id"]})
    return matched


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=30, help="actives per target (inactives matched 1:1)")
    args = ap.parse_args()

    for gene in ("RARG", "MMP3"):
        print(f"{gene}:")
        rows = load_rows(gene)
        active_rows = [r for r in rows if r["label"] == "active"]
        inactive_rows = [r for r in rows if r["label"] in INACTIVE_LABELS]
        print(f"  pool: {len(active_rows)} active, {len(inactive_rows)} inactive")

        # Only fetch MW/logP for the N actives we'd actually pick (potency-
        # spread selection needs standard_value_nM only) plus the full
        # inactive pool (needed for nearest-neighbor matching) -- fetching
        # properties for the *entire* active pool first (up to 1387 for
        # MMP3) was the original approach and made this step needlessly
        # slow against a flaky API for no benefit, since >95% of those
        # properties would never be used.
        active_candidates = pick_actives(active_rows, args.n)
        candidate_ids = [r["molecule_chembl_id"] for r in active_candidates] + \
                        [r["molecule_chembl_id"] for r in inactive_rows]
        props = fetch_properties(candidate_ids)
        active_candidates = [r for r in active_candidates if r["molecule_chembl_id"] in props]
        inactive_rows = [r for r in inactive_rows if r["molecule_chembl_id"] in props]
        print(f"  with valid MW/logP: {len(active_candidates)} active (of {args.n} picked), "
              f"{len(inactive_rows)} inactive pool")

        actives = active_candidates
        inactives = match_inactives(actives, inactive_rows, props)
        print(f"  sampled: {len(actives)} active, {len(inactives)} inactive (property-matched)")

        out_rows = []
        for r in actives:
            mwt, alogp = props[r["molecule_chembl_id"]]
            out_rows.append({
                "molecule_chembl_id": r["molecule_chembl_id"],
                "smiles": strip_salt(r["canonical_smiles"]),
                "label": "active",
                "standard_value_nM": r["standard_value_nM"],
                "full_mwt": mwt, "alogp": alogp,
                "matched_active_id": "",
            })
        for r in inactives:
            mwt, alogp = props[r["molecule_chembl_id"]]
            out_rows.append({
                "molecule_chembl_id": r["molecule_chembl_id"],
                "smiles": strip_salt(r["canonical_smiles"]),
                "label": "inactive",
                "standard_value_nM": r["standard_value_nM"],
                "full_mwt": mwt, "alogp": alogp,
                "matched_active_id": r["matched_active_id"],
            })

        out_path = DATA_DIR / f"{gene}_docking_validation_set.csv"
        with open(out_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=[
                "molecule_chembl_id", "smiles", "label", "standard_value_nM",
                "full_mwt", "alogp", "matched_active_id",
            ])
            writer.writeheader()
            writer.writerows(out_rows)
        print(f"  wrote {len(out_rows)} rows -> {out_path}")


if __name__ == "__main__":
    main()
