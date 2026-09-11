"""Validation ligand set: 7 known TRPV1 agonists spanning ~5 orders of
magnitude in EC50, all sharing the capsaicinoid/vanillamide pharmacophore
so they dock into the same (agonist-locked, 11CK) pocket conformation.

Isomeric (stereo-accurate) SMILES from PubChem -- ConnectivitySMILES lacks
E/Z double-bond geometry, which matters for the shape of these long
unsaturated acyl chains.
"""

# name -> (isomeric SMILES, EC50 in nM, PubChem CID)
LIGANDS = {
    "resiniferatoxin": (
        "C[C@@H]1C[C@]2([C@H]3[C@H]4[C@]1([C@@H]5C=C(C(=O)[C@]5(CC(=C4)COC(=O)CC6=CC(=C(C=C6)O)OC)O)C)"
        "O[C@](O3)(O2)CC7=CC=CC=C7)C(=C)C",
        0.011, 5702546,
    ),
    "phenylacetylrinvanil": (
        "CCCCCC[C@H](C/C=C\\CCCCCCCC(=O)NCC1=CC(=C(C=C1)O)OC)OC(=O)CC2=CC=CC=C2",
        0.090, 25171466,
    ),
    "olvanil": (
        "CCCCCCCC/C=C\\CCCCCCCC(=O)NCC1=CC(=C(C=C1)O)OC",
        0.6, 5311093,
    ),
    "arvanil": (
        "CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCC(=O)NCC1=CC(=C(C=C1)O)OC",
        0.5, 6449767,
    ),
    "rinvanil": (
        "CCCCCC[C@H](C/C=C\\CCCCCCCC(=O)NCC1=CC(=C(C=C1)O)OC)O",
        6.0, 10252323,
    ),
    "capsaicin": (
        "CC(C)/C=C/CCCCC(=O)NCC1=CC(=C(C=C1)O)OC",
        712.0, 1548943,
    ),
    "nonivamide": (
        "CCCCCCCCC(=O)NCC1=CC(=C(C=C1)O)OC",
        1400.0, 2998,
    ),
}

# Competitive antagonist(s) for the full-tetramer, unrestrained gating pilot
# (does the pore relax toward closed with an antagonist bound in the same
# agonist-locked pocket, vs. stay open with an agonist -- see STATUS.md).
# name -> (isomeric SMILES, PubChem CID) -- no EC50 here, antagonists are
# characterized by IC50/Ki, not relevant to this qualitative drift check.
# SMILES/CID verified directly against PubChem's REST API (2026-08-06), not
# just a search-engine summary, after an earlier literature-citation error
# this same session.
ANTAGONISTS = {
    "capsazepine": (
        "C1CC2=CC(=C(C=C2CN(C1)C(=S)NCCC3=CC=C(C=C3)Cl)O)O",
        2733484,
    ),
}
