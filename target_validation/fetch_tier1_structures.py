"""Fetch best-available structure per Tier1 skin off-target gene (see
skin_offtarget_targets.md): prefer experimental PDB (largest chain span, best
resolution), fall back to AlphaFold model. Writes protein-only, single-chain
PDBs to data/raw/structures/tier1_skin_targets/{GENE}.pdb

Requires: requests, biopython (available in the `pocket-compare` conda env
or docker/environment-docking.yml once fpocket is added there).
"""
import csv
import re
import sys
import time
from pathlib import Path

import requests
from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select

OUT_DIR = Path(__file__).resolve().parents[2] / "data" / "raw" / "structures" / "tier1_skin_targets"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TIER1_GENES = [
    "ADRB2", "CA12", "RARG", "BLMH", "MMP3", "RORA", "CD207", "FCER1A",
    "GSTA3", "SULT1E1", "SULT2B1", "FABP9", "ANXA8", "DSP", "ITGB4", "SFN",
    "CHP2", "FBXW7", "GAN", "GLTP", "NTF4", "RAB3D", "SERPINB2", "SERPINB5",
    "SRY", "TREX2",
]


class ProteinChainOnly(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_residue(self, residue):
        return residue.id[0] == " " and residue.get_parent().id == self.chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id


def uniprot_lookup(gene):
    # NB: UniProt's `gene:X` query is a free-text match, not an exact-symbol
    # filter -- it can return an unrelated gene whose synonym/alias happens to
    # match (e.g. querying "FBXW7" returned B0L3A2 / FBXW7-AS1, a completely
    # different antisense-derived receptor gene, not FBXW7 itself). Always
    # verify the returned primary gene symbol equals the query before using it.
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"gene:{gene} AND organism_id:9606 AND reviewed:true",
        "fields": "accession,gene_names,xref_pdb",
        "format": "json",
    }
    r = requests.get(url, params=params, timeout=30)
    r.raise_for_status()
    data = r.json()
    entry = None
    for candidate in data.get("results", []):
        primary_gene = candidate["genes"][0]["geneName"]["value"] if candidate.get("genes") else None
        if primary_gene == gene:
            entry = candidate
            break
    if entry is None:
        print(f"  WARNING: no UniProt entry with exact primary gene symbol '{gene}' "
              f"(got {[c['genes'][0]['geneName']['value'] for c in data.get('results', []) if c.get('genes')]})")
        return None, []
    accession = entry["primaryAccession"]
    pdb_hits = []
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref["database"] != "PDB":
            continue
        props = {p["key"]: p["value"] for p in xref.get("properties", [])}
        method = props.get("Method", "")
        resolution_str = props.get("Resolution", "")
        chains_str = props.get("Chains", "")
        m = re.match(r"([A-Za-z0-9]+)(?:/[A-Za-z0-9]+)*=(\d+)-(\d+)", chains_str)
        if not m:
            continue
        chain_id, start, end = m.group(1), int(m.group(2)), int(m.group(3))
        span = end - start
        try:
            resolution = float(resolution_str.replace(" A", ""))
        except ValueError:
            resolution = 99.0
        if method not in ("X-ray", "EM"):
            continue
        pdb_hits.append({
            "pdb_id": xref["id"], "chain": chain_id, "span": span, "resolution": resolution,
        })
    pdb_hits.sort(key=lambda h: (-h["span"], h["resolution"]))
    return accession, pdb_hits


def fetch_pdb_chain(pdb_id, chain_id, out_path):
    r = requests.get(f"https://files.rcsb.org/download/{pdb_id}.cif", timeout=60)
    if r.status_code != 200:
        return False
    cif_path = out_path.with_suffix(".cif.tmp")
    cif_path.write_bytes(r.content)
    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure(pdb_id, str(cif_path))
        model = structure[0]
        if chain_id not in model:
            # some mmCIF author chain IDs differ from UniProt's; try first chain as fallback
            available = [c.id for c in model]
            if not available:
                return False
            chain_id_use = available[0]
        else:
            chain_id_use = chain_id
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(out_path), ProteinChainOnly(chain_id_use))
        n_atoms = sum(1 for _ in open(out_path) if _.startswith("ATOM"))
        return n_atoms > 200
    except Exception as e:
        print(f"    parse error {pdb_id}: {e}")
        return False
    finally:
        cif_path.unlink(missing_ok=True)


def fetch_alphafold(accession, out_path):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v4.pdb"
    r = requests.get(url, timeout=60)
    if r.status_code != 200:
        return False
    out_path.write_bytes(r.content)
    return True


def main():
    log = []
    for gene in TIER1_GENES:
        out_path = OUT_DIR / f"{gene}.pdb"
        if out_path.exists():
            print(f"{gene}: already fetched, skipping")
            continue
        print(f"{gene}: looking up...")
        try:
            accession, pdb_hits = uniprot_lookup(gene)
        except Exception as e:
            print(f"  FAILED uniprot lookup: {e}")
            log.append((gene, "FAILED", "", ""))
            continue
        if accession is None:
            print("  no uniprot entry found")
            log.append((gene, "FAILED", "no uniprot entry", ""))
            continue

        source = None
        used = ""
        for hit in pdb_hits[:5]:
            print(f"  trying PDB {hit['pdb_id']} chain {hit['chain']} (span={hit['span']}, res={hit['resolution']})")
            if fetch_pdb_chain(hit["pdb_id"], hit["chain"], out_path):
                source = "PDB"
                used = f"{hit['pdb_id']}:{hit['chain']}"
                break
            time.sleep(0.3)
        if source is None:
            print(f"  falling back to AlphaFold ({accession})")
            if fetch_alphafold(accession, out_path):
                source = "AlphaFold"
                used = accession
            else:
                print("  FAILED alphafold fetch too")
        log.append((gene, source or "FAILED", used, accession))
        time.sleep(0.3)

    print("\n=== SUMMARY ===")
    for gene, source, used, accession in log:
        print(f"{gene:10s} {source:10s} {used:15s} {accession}")


if __name__ == "__main__":
    main()
