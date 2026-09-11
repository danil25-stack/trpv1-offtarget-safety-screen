"""Bootstrap CI, Mann-Whitney significance, and PR-AUC for the ground-truth
discrimination check (OFFTARGET_SAFETY_ARTICLE Section 3.4 / Table 3).

ROC-AUC alone doesn't convey how much a 30-active/30-inactive panel could
have produced that number by chance -- this reports the same three docking
results (RARG/Vina, MMP3/Vina, MMP3/AD4Zn) with a 95% bootstrap CI, a
one-sided Mann-Whitney U test against AUC=0.5, and PR-AUC as a check against
ROC-AUC alone looking optimistic on a small, balanced panel.

Also reports the EC50-mixing sensitivity check: RARG's panel is ~60%
functional-assay (EC50) rather than direct-binding (Ki/Kd) or enzymatic
(IC50) records; this recomputes AUC restricted to the binding/enzymatic-only
subset to rule out EC50/binding-assay mixing as the source of the RARG
result.

Usage:
    python discrimination_stats.py
"""
import csv
from pathlib import Path

import numpy as np
from scipy import stats
from sklearn.metrics import average_precision_score, roc_auc_score

DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "processed" / "skin_offtarget_screen"
ACTIVE_NM = 1000.0
INACTIVE_NM = 10000.0
N_BOOTSTRAP = 10000
SEED = 42


def load(path):
    rows = list(csv.DictReader(open(path)))
    scores = np.array([-float(r["vina_score"]) for r in rows])  # flip sign: higher = more active
    labels = np.array([1 if r["label"] == "active" else 0 for r in rows])
    return scores, labels


def bootstrap_auc_ci(scores, labels, n=N_BOOTSTRAP, seed=SEED):
    rng = np.random.default_rng(seed)
    idx_pos = np.where(labels == 1)[0]
    idx_neg = np.where(labels == 0)[0]
    aucs = []
    for _ in range(n):
        bp = rng.choice(idx_pos, len(idx_pos), replace=True)
        bn = rng.choice(idx_neg, len(idx_neg), replace=True)
        bidx = np.concatenate([bp, bn])
        try:
            aucs.append(roc_auc_score(labels[bidx], scores[bidx]))
        except ValueError:
            continue
    return np.percentile(aucs, [2.5, 97.5])


def discrimination_report(name, path):
    scores, labels = load(path)
    auc = roc_auc_score(labels, scores)
    ci_lo, ci_hi = bootstrap_auc_ci(scores, labels)
    pos, neg = scores[labels == 1], scores[labels == 0]
    _, p_value = stats.mannwhitneyu(pos, neg, alternative="greater")
    pr_auc = average_precision_score(labels, scores)
    print(f"{name}: AUC={auc:.3f} [95% CI {ci_lo:.3f}-{ci_hi:.3f}], "
          f"Mann-Whitney p={p_value:.4f}, PR-AUC={pr_auc:.3f}")


def classify_with_type(raw_path):
    """Recover the ChEMBL standard_type of the record that actually set
    each molecule's active/inactive label (mirrors fetch_chembl_bioactivity.
    classify(), but keeps the deciding record's standard_type)."""
    rows = list(csv.DictReader(open(raw_path)))
    by_mol = {}
    for r in rows:
        by_mol.setdefault(r["molecule_chembl_id"], []).append(r)
    result = {}
    for mid, recs in by_mol.items():
        eq = []
        for r in recs:
            if r["standard_relation"] != "=" or r["data_validity_comment"]:
                continue
            try:
                eq.append((float(r["standard_value"]), r["standard_type"]))
            except (TypeError, ValueError):
                continue
        if eq:
            value, st = min(eq, key=lambda t: t[0])
            result[mid] = st
            continue
        gt = []
        for r in recs:
            if r["standard_relation"] != ">":
                continue
            try:
                gt.append((float(r["standard_value"]), r["standard_type"]))
            except (TypeError, ValueError):
                continue
        if gt:
            value, st = max(gt, key=lambda t: t[0])
            if value >= INACTIVE_NM:
                result[mid] = st
    return result


def ec50_sensitivity_check(target, docked_csv, raw_csv):
    type_map = classify_with_type(raw_csv)
    docked = list(csv.DictReader(open(docked_csv)))
    scores, labels, deciding_types = [], [], []
    for r in docked:
        scores.append(-float(r["vina_score"]))
        labels.append(1 if r["label"] == "active" else 0)
        deciding_types.append(type_map.get(r["molecule_chembl_id"], "UNKNOWN"))
    scores, labels = np.array(scores), np.array(labels)
    deciding_types = np.array(deciding_types)

    full_auc = roc_auc_score(labels, scores)
    mask = deciding_types != "EC50"
    n_ec50 = (~mask).sum()
    print(f"{target}: {n_ec50}/{len(labels)} labels set by an EC50 record; full AUC={full_auc:.3f}")
    if n_ec50 == 0:
        print(f"  no EC50-mixing concern -- panel is entirely IC50/Ki/Kd")
        return
    if mask.sum() >= 6 and len(set(labels[mask])) == 2:
        sub_auc = roc_auc_score(labels[mask], scores[mask])
        print(f"  binding/enzymatic-only subset: n={mask.sum()} AUC={sub_auc:.3f}")


if __name__ == "__main__":
    print("=== Ground-truth discrimination: bootstrap CI, significance, PR-AUC ===")
    discrimination_report("RARG (Vina)", DATA_DIR / "RARG_docking_validation_results.csv")
    discrimination_report("MMP3 (Vina)", DATA_DIR / "MMP3_docking_validation_results.csv")
    discrimination_report("MMP3 (AD4Zn)", DATA_DIR / "MMP3_ad4zn_docking_validation_results.csv")

    print("\n=== EC50/binding-assay mixing sensitivity check ===")
    ec50_sensitivity_check("RARG", DATA_DIR / "RARG_docking_validation_results.csv",
                            DATA_DIR / "RARG_chembl_actives_inactives.csv")
    ec50_sensitivity_check("MMP3", DATA_DIR / "MMP3_docking_validation_results.csv",
                            DATA_DIR / "MMP3_chembl_actives_inactives.csv")
