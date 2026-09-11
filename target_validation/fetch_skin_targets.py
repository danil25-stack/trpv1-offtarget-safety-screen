"""Build a candidate list of druggable skin-expressed off-target proteins for
pocket-similarity screening against TRPV1's vanilloid pocket.

Source list: Human Protein Atlas RNA tissue-specificity API, filtered to genes
classified skin-enriched / group-enriched / tissue-enhanced. Each gene is then
cross-referenced against Open Targets target class + small-molecule tractability
(pocket-quality evidence) to separate real druggable pockets from structural
proteins (keratins, collagens, etc.) that dominate the raw HPA skin list.

Usage:
    python fetch_skin_targets.py            # re-fetch HPA list + re-annotate
    python fetch_skin_targets.py --annotate-only  # reuse cached HPA TSV
"""
import csv
import sys
import time
from pathlib import Path

import requests

RAW_DIR = Path(__file__).resolve().parents[2] / "data" / "raw" / "hpa"
OUT_DIR = Path(__file__).resolve().parents[2] / "data" / "processed" / "skin_offtarget_screen"
HPA_TSV = RAW_DIR / "skin_tissue_specific_genes.tsv"
OUT_CSV = OUT_DIR / "skin_targets_annotated.csv"

HPA_API_URL = "https://www.proteinatlas.org/api/search_download.php"
HPA_SEARCH = "tissue_category_rna:skin;Tissue enriched,Group enriched,Tissue enhanced"
OT_API_URL = "https://api.platform.opentargets.org/api/v4/graphql"
CHUNK = 15


def fetch_hpa_list():
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    params = {
        "format": "tsv",
        "columns": "g,eg,rnatsm,rnats",
        "search": HPA_SEARCH,
    }
    r = requests.get(HPA_API_URL, params=params, timeout=60)
    r.raise_for_status()
    HPA_TSV.write_bytes(r.content)


def load_hpa():
    rows = []
    with open(HPA_TSV) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            rows.append({
                "gene": r["Gene"],
                "ensembl": r["Ensembl"],
                "specificity": r["RNA tissue specificity"],
                "detail": r["RNA tissue specific nTPM"],
            })
    return rows


def batch_query(chunk):
    fields = []
    variables = {}
    var_defs = []
    for i, row in enumerate(chunk):
        alias = f"t{i}"
        var_defs.append(f"$id{i}: String!")
        variables[f"id{i}"] = row["ensembl"]
        fields.append(f'''
          {alias}: target(ensemblId: $id{i}) {{
            approvedSymbol
            biotype
            targetClass {{ label }}
            tractability {{ modality value label }}
          }}
        ''')
    query = "query Batch(" + ", ".join(var_defs) + ") {" + "".join(fields) + "}"
    r = requests.post(OT_API_URL, json={"query": query, "variables": variables}, timeout=60)
    r.raise_for_status()
    data = r.json()
    if "errors" in data:
        print("  errors:", data["errors"][:2])
    return data.get("data", {})


def annotate(hpa_rows):
    results = []
    for start in range(0, len(hpa_rows), CHUNK):
        chunk = hpa_rows[start:start + CHUNK]
        print(f"  querying {start}-{start + len(chunk)} / {len(hpa_rows)}")
        try:
            data = batch_query(chunk)
        except Exception as e:
            print(f"    FAILED chunk at {start}: {e}")
            time.sleep(2)
            continue
        for i, row in enumerate(chunk):
            t = data.get(f"t{i}")
            if not t:
                results.append({**row, "biotype": "NOT_FOUND", "target_class": "", "pocket_quality": "", "sm_approved": False})
                continue
            classes = "; ".join(c["label"] for c in (t.get("targetClass") or []))
            tract = t.get("tractability") or []
            sm = [x for x in tract if x["modality"] == "SM" and x["value"]]
            sm_labels = {x["label"] for x in sm}
            pocket_quality = (
                "High-Quality Pocket" if "High-Quality Pocket" in sm_labels else
                "Med-Quality Pocket" if "Med-Quality Pocket" in sm_labels else
                "Ligand evidence only" if sm_labels & {"High-Quality Ligand", "Structure with Ligand"} else
                ""
            )
            results.append({
                **row,
                "biotype": t.get("biotype", ""),
                "target_class": classes,
                "pocket_quality": pocket_quality,
                "sm_approved": "Approved Drug" in sm_labels,
            })
        time.sleep(0.3)
    return results


def main():
    if "--annotate-only" not in sys.argv or not HPA_TSV.exists():
        print("Fetching HPA skin-specificity gene list...")
        fetch_hpa_list()
    hpa_rows = load_hpa()
    print(f"Loaded {len(hpa_rows)} HPA skin-specific genes")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results = annotate(hpa_rows)
    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["gene", "ensembl", "specificity", "detail", "biotype", "target_class", "pocket_quality", "sm_approved"])
        writer.writeheader()
        writer.writerows(results)
    print(f"Wrote {len(results)} rows -> {OUT_CSV}")


if __name__ == "__main__":
    main()
