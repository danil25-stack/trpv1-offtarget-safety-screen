"""Fetch known active/inactive small molecules for RARG and MMP3 from ChEMBL --
the missing ground-truth check for skin_offtarget_pocket_comparison.md's
cross-docking results: without this, nothing confirms Vina's scoring function
discriminates real binders from non-binders on these two off-target folds at
all, which is exactly the "cross-receptor score comparisons carry additional
uncertainty" limitation flagged in OFFTARGET_SAFETY_ARTICLE_EN.md Discussion.

BindingDB's documented REST endpoints (getLigandsByUniprot[s]) were tried
first and return HTTP 200 with an empty body for every UniProt tested,
including well-annotated targets unrelated to this project (e.g. carbonic
anhydrase, P00918) -- the service looks dead, not merely lacking RARG/MMP3
data. ChEMBL is used as the sole source here; revisit BindingDB via its bulk
TSV dump if independent corroboration is later needed.

Target ChEMBL IDs resolved via /target?target_components__accession=<UniProt>
and confirmed SINGLE PROTEIN, Homo sapiens:
  RARG (UniProt P13631) -> CHEMBL2003
  MMP3 (UniProt P08254) -> CHEMBL283

The EBI ChEMBL API is intermittently flaky (transient 500s on otherwise-valid
queries, confirmed by immediate retry succeeding) -- every request goes
through get_json()'s retry-with-backoff.

Usage:
    python fetch_chembl_bioactivity.py
"""
import csv
import time
from pathlib import Path

import requests

RAW_DIR = Path(__file__).resolve().parents[2] / "data" / "raw" / "chembl"
OUT_DIR = Path(__file__).resolve().parents[2] / "data" / "processed" / "skin_offtarget_screen"
BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/activity"

TARGETS = {
    "RARG": "CHEMBL2003",
    "MMP3": "CHEMBL283",
}

STANDARD_TYPES = "IC50,EC50,Ki,Kd"
PAGE_LIMIT = 1000

# nM thresholds for the active/inactive call -- shared across both targets
# rather than target-tuned, matching the common ChEMBL-derived-benchmark
# convention (e.g. ExCAPE-DB): <=1000 nM active, >10000 nM inactive, the
# 1000-10000 nM band left as ambiguous rather than forced into either class.
ACTIVE_NM = 1000.0
INACTIVE_NM = 10000.0

ACTIVITY_FIELDS = [
    "molecule_chembl_id", "canonical_smiles", "standard_type",
    "standard_relation", "standard_value", "standard_units",
    "pchembl_value", "activity_comment", "data_validity_comment",
    "potential_duplicate", "assay_type", "assay_description",
    "document_chembl_id", "document_year",
]


def get_json(url, params, retries=6, backoff=3):
    last_exc = None
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=60)
            if r.status_code == 200 and r.content:
                return r.json()
            last_exc = RuntimeError(f"HTTP {r.status_code}, {len(r.content)} bytes")
        except requests.RequestException as e:
            last_exc = e
        wait = backoff * (attempt + 1)
        print(f"    retry {attempt + 1}/{retries} after {last_exc} (sleeping {wait}s)")
        time.sleep(wait)
    raise RuntimeError(f"giving up after {retries} retries: {last_exc}")


def fetch_activities(target_chembl_id):
    rows = []
    params = {
        "target_chembl_id": target_chembl_id,
        "standard_type__in": STANDARD_TYPES,
        "standard_units": "nM",
        "limit": PAGE_LIMIT,
        "offset": 0,
        "format": "json",
    }
    url = BASE_URL
    total = None
    while True:
        data = get_json(url, params)
        page = data["page_meta"]
        if total is None:
            total = page["total_count"]
            print(f"  {total} activities total")
        for a in data["activities"]:
            rows.append({k: a.get(k) for k in ACTIVITY_FIELDS})
        print(f"  fetched {len(rows)}/{total}")
        if not page["next"]:
            break
        # 'next' is a path with its own querystring already encoded (offset
        # advanced server-side) -- follow it directly rather than
        # reconstructing params, and stop passing 'params' on subsequent
        # requests since they're baked into the path.
        url = "https://www.ebi.ac.uk" + page["next"]
        params = None
        time.sleep(0.3)
    return rows


def classify(rows):
    """Dedupe by molecule_chembl_id and label active/inactive/ambiguous.

    Two distinct evidence types, kept separate rather than merged into one
    "best value" per molecule:
    - '=' records: a measured potency. Best (lowest) value per molecule
      drives active/inactive/ambiguous the usual way.
    - '>' records for a molecule with NO '=' record at all: "no activity
      detected up to X nM" -- a cleaner negative than a weak '=' value, but
      only usable for molecules that were never also measured with '='
      (otherwise the '=' record is the real potency and the '>' record from
      a different, less sensitive assay would just be noise).
    """
    valid = [r for r in rows if not r["data_validity_comment"] and r["standard_value"] is not None]
    by_molecule = {}
    for r in valid:
        by_molecule.setdefault(r["molecule_chembl_id"], []).append(r)

    out = []
    for mid, recs in by_molecule.items():
        eq_recs = []
        for r in recs:
            if r["standard_relation"] != "=":
                continue
            try:
                eq_recs.append((float(r["standard_value"]), r))
            except (TypeError, ValueError):
                continue
        if eq_recs:
            value, r = min(eq_recs, key=lambda t: t[0])
            label = "active" if value <= ACTIVE_NM else "inactive" if value > INACTIVE_NM else "ambiguous"
            out.append({**r, "standard_value_nM": value, "label": label})
            continue

        gt_recs = []
        for r in recs:
            if r["standard_relation"] != ">":
                continue
            try:
                gt_recs.append((float(r["standard_value"]), r))
            except (TypeError, ValueError):
                continue
        if gt_recs:
            value, r = max(gt_recs, key=lambda t: t[0])  # strongest ">" threshold = strongest negative evidence
            if value >= INACTIVE_NM:
                out.append({**r, "standard_value_nM": value, "label": "inactive_no_binding_detected"})

    out.sort(key=lambda r: (r["label"] != "active", r["standard_value_nM"]))
    return out


def main():
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    for gene, target_chembl_id in TARGETS.items():
        print(f"{gene} ({target_chembl_id}): fetching activities...")
        rows = fetch_activities(target_chembl_id)

        raw_csv = RAW_DIR / f"{gene}_activities_raw.csv"
        with open(raw_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=ACTIVITY_FIELDS)
            writer.writeheader()
            writer.writerows(rows)
        print(f"  wrote {len(rows)} raw rows -> {raw_csv}")

        classified = classify(rows)
        n_active = sum(1 for r in classified if r["label"] == "active")
        n_inactive = sum(1 for r in classified if r["label"] == "inactive")
        n_no_binding = sum(1 for r in classified if r["label"] == "inactive_no_binding_detected")
        n_ambiguous = sum(1 for r in classified if r["label"] == "ambiguous")
        print(f"  {len(classified)} unique molecules classified: "
              f"{n_active} active (<={ACTIVE_NM:.0f} nM), "
              f"{n_inactive} inactive (measured '=' >{INACTIVE_NM:.0f} nM), "
              f"{n_no_binding} inactive (no binding detected, '>' only), "
              f"{n_ambiguous} ambiguous")

        out_csv = OUT_DIR / f"{gene}_chembl_actives_inactives.csv"
        fieldnames = ACTIVITY_FIELDS + ["standard_value_nM", "label"]
        with open(out_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(classified)
        print(f"  wrote {len(classified)} classified rows -> {out_csv}")


if __name__ == "__main__":
    main()
