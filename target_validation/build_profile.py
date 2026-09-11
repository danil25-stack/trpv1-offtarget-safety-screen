"""Build a documented target profile for a Track A (topical analgesic) candidate:
genetic validation, safety liabilities, paralog/off-target family, expression &
localization, tractability. Pulls from the Open Targets Platform GraphQL API
(public, no license gate) and layers in curated literature caveats that the API
doesn't capture.

Usage:
    python build_profile.py SCN9A TRPV1
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common.opentargets_client import resolve_ensembl_id, fetch_target_profile

RAW_DIR = Path(__file__).resolve().parents[2] / "data" / "processed" / "target_profiles"
REPORT_DIR = Path(__file__).resolve().parent / "profiles"

# Literature findings that materially affect interpretation of a target's genetic
# validation but aren't in the Open Targets API — kept separate and explicit so
# they don't get silently lost between curated-DB evidence and population evidence.
MANUAL_CAVEATS = {
    "SCN9A": [
        "UK Biobank burden analysis (medRxiv 2025.03.12.25323817, ~148k carriers of "
        "59 known pathogenic gain-of-function variants) found NO increase in chronic/"
        "neuropathic pain prevalence or analgesic/opioid prescriptions at the population "
        "level, despite well-established monogenic channelopathies (erythromelalgia, "
        "PEPD). Rare high-penetrance Mendelian phenotypes do not reliably show up in "
        "biobank-scale burden PheWAS — treat curated OMIM/ClinVar phenotype evidence and "
        "population PheWAS burden evidence as separate, non-interchangeable evidence types.",
    ],
    "TRPV1": [
        "No monogenic human channelopathy comparable to SCN9A is established for TRPV1 — "
        "genetic validation here rests more on GWAS/expression/pharmacology than on rare "
        "high-penetrance human variants. Clinical precedent instead comes from topical "
        "capsaicin desensitization (Qutenza, an approved TRPV1 agonist patch) and from the "
        "failure of systemic TRPV1 antagonists due to on-target hyperthermia — an argument "
        "specifically in favor of a topical/local route for this target, not systemic.",
        "Off-target paralog literature pass (2026-08-10): of the five sequence paralogs, "
        "TRPV6 is the most systemically consequential to hit (KO mice: impaired intestinal "
        "Ca2+ absorption, reduced bone mineral density, alopecia/dermatitis, severe male "
        "infertility) and TRPV5 is a plausible renal Ca2+-handling risk (rate-limiting "
        "channel for distal-tubule Ca2+ reabsorption); neither had liability text in the "
        "original profile. TRPV4's 'pulmonary edema' liability is direction-dependent: the "
        "clinical antagonist GSK2798745 was well-tolerated in heart-failure patients with no "
        "cardiopulmonary safety signal, so off-target *inhibition* by a TRPV1 antagonist is "
        "not expected to reproduce the edema risk (that's an activation-direction liability). "
        "TRPV3 has direct human proof of consequence at the skin (Olmsted syndrome, from "
        "gain-of-function mutations) — relevant given TRPV1 and TRPV3 share the topical site "
        "of action, though GoF disease phenotype doesn't necessarily predict what pharmacological "
        "antagonism would do. TRPV2 looks acutely low-risk (adult whole-body KO and "
        "neutralizing-antibody dosing show no overt phenotype) but has an uncharacterized "
        "developmental cardiac role. TRPA1 (mechanistically adjacent, not a sequence paralog) "
        "has an acceptable Phase 1 safety record in two independent clinical antagonist "
        "programs (Lilly's LY3526318, Roche/Genentech), both of which nonetheless failed on "
        "efficacy — so an off-target TRPA1 hit looks safety-neutral, not liability-free by "
        "default, more 'clinically unproven for this indication.' Net: no paralog in this "
        "family blocks the topical program on safety grounds alone, but TRPV5/TRPV6 warrant "
        "an explicit selectivity check in any hit-to-lead off-target panel, which wasn't "
        "previously called out.",
    ],
    "SCN10A": [
        "The SCN10A locus carries one of the strongest GWAS signals for cardiac conduction "
        "(PR interval) and atrial fibrillation risk in the human genome -- reflected here in "
        "'atrial fibrillation' and 'cardiac arrhythmia' scoring almost as high as the pain "
        "phenotype itself. There is a long-running debate over whether this is a true SCN10A "
        "(Nav1.8) effect or a regulatory effect of the locus on the adjacent SCN5A (Nav1.5, "
        "cardiac) promoter -- either way, cardiac conduction must be tracked as a first-order "
        "risk, not a paralog-selectivity afterthought. Suzetrigine (Journavx, Vertex), a "
        "Nav1.8-selective blocker, was FDA-approved for acute pain in 2025, so the liability "
        "appears manageable at achieved clinical exposures -- but this was a systemic, "
        "peripherally-restricted-by-expression strategy, not a topical one; the safety margin "
        "for a topical formulation still needs its own exposure/PK justification.",
    ],
    "CNR2": [
        "Open Targets' direct genetic association score for 'Pain' is weak (~0.1, not in the "
        "top-5 disease list) despite CB2 having a well-established pharmacological rationale "
        "as an analgesic/anti-inflammatory target -- the case for CB2 here is a mechanistic/"
        "pharmacology case, not a human-genetics case, and should not be scored the same way "
        "as SCN9A's or SCN10A's genetic evidence. Curated safety liabilities are also "
        "bidirectional (both increased and decreased inflammation reported), consistent with "
        "known context/dose-dependence of cannabinoid receptor signaling.",
    ],
    "P2RX3": [
        "The strongest validated human indication for P2X3 antagonism is chronic cough "
        "(gefapixant, approved), not pain -- 'Cough' dominates the associated-disease list "
        "and a pain association isn't in the top 5. The clinically dose-limiting off-target "
        "effect (taste disturbance) comes from the P2X2/P2X3 heteromer, i.e. from P2RX2, not "
        "from sequence-paralog identity alone -- make sure P2RX2 co-expression/heteromer "
        "pharmacology is tracked explicitly rather than assumed to show up in a generic "
        "paralog-identity ranking.",
    ],
}

# Paralogs with a known mechanism-based liability if hit off-target, keyed by
# the queried gene's symbol. Supplements raw sequence-identity ranking with why
# a given paralog actually matters.
KEY_OFF_TARGET_FAMILY = {
    "SCN9A": {
        "SCN5A": "cardiac (Nav1.5) — arrhythmia risk",
        "SCN4A": "skeletal muscle (Nav1.4) — myotonia/paralysis risk",
        "SCN1A": "CNS (Nav1.1) — seizure risk (Dravet syndrome gene)",
        "SCN2A": "CNS (Nav1.2) — seizure risk",
        "SCN8A": "CNS (Nav1.6) — seizure risk, essential for action potential propagation",
        "SCN10A": "DRG/pain (Nav1.8) — co-expressed nociceptor channel, mechanistically adjacent",
        "SCN11A": "DRG/pain (Nav1.9) — co-expressed nociceptor channel, mechanistically adjacent",
    },
    "TRPV1": {
        "TRPV2": "developmental cardiac contractility/Ca2+ handling (KO cardiomyopathy in growing hearts); "
                 "adult whole-body KO and neutralizing-antibody dosing show no overt phenotype -- acute "
                 "pharmacological inhibition looks low-risk, developmental/chronic exposure less characterized",
        "TRPV3": "keratinocytes/skin barrier & thermosensation (shares topical site of action); human GoF "
                 "mutations cause Olmsted syndrome (painful keratoderma, pruritus, skin barrier failure) -- "
                 "direct human proof this channel matters at the skin, though that's a GoF disease phenotype, "
                 "not necessarily what antagonism would produce",
        "TRPV4": "osmo/mechanosensation; activation implicated in cardiogenic pulmonary edema, BUT clinical "
                 "TRPV4 *antagonist* GSK2798745 was well-tolerated through repeat dosing in heart-failure "
                 "patients with no significant cardiovascular/pulmonary safety signal -- direction matters: "
                 "TRPV1-antagonist off-target TRPV4 *inhibition* trends protective for edema, not causative",
        "TRPV5": "kidney distal-tubule Ca2+ reabsorption, rate-limiting step of renal Ca2+ handling -- "
                 "off-target inhibition is a plausible hypercalciuria/renal Ca2+ handling risk",
        "TRPV6": "intestinal Ca2+ absorption; KO mice show impaired intestinal Ca2+ uptake, reduced bone "
                 "mineral density/osteoporosis, alopecia/dermatitis, and severely impaired male fertility -- "
                 "the most systemically consequential TRPV paralog to hit off-target if these phenotypes "
                 "translate to pharmacological inhibition",
        "TRPA1": "co-expressed nociceptor channel, mechanistically adjacent (not a sequence paralog); "
                 "clinical TRPA1 antagonists (LY3526318, Roche/Genentech compound) had acceptable Phase 1 "
                 "safety but failed efficacy in OA/CLBP/DPNP and chronic cough -- off-target TRPA1 hit looks "
                 "safety-neutral but track record of clinical failure is itself informative",
    },
    "SCN10A": {
        "SCN5A": "cardiac (Nav1.5) — arrhythmia risk; SCN10A locus itself is GWAS-linked to AF/conduction",
        "SCN4A": "skeletal muscle (Nav1.4) — myotonia/paralysis risk",
        "SCN9A": "DRG/pain (Nav1.7) — co-expressed nociceptor channel, mechanistically adjacent",
        "SCN11A": "DRG/pain (Nav1.9) — co-expressed nociceptor channel, mechanistically adjacent",
    },
    "TRPA1": {
        "TRPV1": "co-expressed nociceptor channel, mechanistically adjacent (not a sequence paralog)",
    },
    "CNR2": {
        "CNR1": "CNS/psychoactive (CB1) — the entire rationale for CB2 selectivity is avoiding this receptor",
    },
    "P2RX3": {
        "P2RX2": "sensory neurons/taste buds — forms the P2X2/3 heteromer responsible for the "
                 "clinically dose-limiting taste-disturbance liability seen with P2X3 antagonists",
    },
}


def render_markdown(profile: dict, gene_symbol: str) -> str:
    lines = []
    lines.append(f"# Target Profile: {profile['approvedSymbol']} ({profile['approvedName']})")
    lines.append("")
    lines.append(f"- Ensembl ID: `{profile['id']}`")
    lines.append(f"- Biotype: {profile['biotype']}")
    if profile.get("targetClass"):
        classes = ", ".join(c["label"] for c in profile["targetClass"])
        lines.append(f"- Target class: {classes}")
    lines.append(f"- Essential (cell-fitness, DepMap-derived): {profile['isEssential']}")
    lines.append("")

    if profile.get("functionDescriptions"):
        lines.append("## Function")
        lines.append(profile["functionDescriptions"][0])
        lines.append("")

    lines.append("## Genetic constraint (gnomAD)")
    lines.append("")
    lines.append("| Constraint | Observed | Expected | o/e | Score |")
    lines.append("|---|---|---|---|---|")
    for c in profile.get("geneticConstraint") or []:
        lines.append(
            f"| {c['constraintType']} | {c['obs']} | {c['exp']:.1f} | "
            f"{c['oe']:.2f} ({c['oeLower']:.2f}-{c['oeUpper']:.2f}) | {c['score']:.2e} |"
        )
    lines.append("")
    lines.append(
        "_Low o/e for `lof` = fewer loss-of-function variants tolerated in the general "
        "population than expected by chance -> gene under purifying selection against LoF._"
    )
    lines.append("")

    lines.append("## Tractability")
    lines.append("")
    tractable = [t for t in profile.get("tractability", []) if t["value"]]
    if tractable:
        lines.append("| Modality | Evidence |")
        lines.append("|---|---|")
        for t in tractable:
            lines.append(f"| {t['modality']} | {t['label']} |")
    else:
        lines.append("_No tractability evidence returned._")
    lines.append("")

    lines.append("## Safety -- curated liabilities (Open Targets)")
    lines.append("")
    liabilities = profile.get("safetyLiabilities") or []
    if liabilities:
        for s in liabilities:
            lines.append(f"- **{s['event']}** (source: {s['datasource']})")
            for e in s.get("effects") or []:
                lines.append(f"  - effect: {e.get('direction')} / dosing: {e.get('dosing')}")
    else:
        lines.append(
            "_No curated safety liabilities in Open Targets for this target -- "
            "does not mean none exist, only that none are database-curated yet._"
        )
    lines.append("")

    lines.append("## Safety -- off-target family (paralog selectivity risk)")
    lines.append("")
    lines.append(
        "Human paralogs ranked by sequence identity. High identity = harder to achieve "
        "selectivity; annotated rows are paralogs with a known mechanism-based liability "
        "if inhibited/activated off-target."
    )
    lines.append("")
    lines.append("| Paralog | % identity | Known liability if hit |")
    lines.append("|---|---|---|")
    known = KEY_OFF_TARGET_FAMILY.get(gene_symbol, {})
    # NB: `isHighConfidence` is not used as a filter here -- for many targets the API
    # returns the literal string "NULL" rather than a boolean/null, which is truthy in
    # Python and silently defeats a naive confidence filter. homologyType is the
    # reliable signal for "this is a real paralog".
    paralogs = [
        h for h in profile.get("homologues", [])
        if h["speciesName"] == "Human" and h["homologyType"] in ("other_paralog", "within_species_paralog")
    ]
    paralogs.sort(key=lambda h: h["targetPercentageIdentity"], reverse=True)
    for h in paralogs:
        liability = known.get(h["targetGeneSymbol"], "")
        lines.append(f"| {h['targetGeneSymbol']} | {h['targetPercentageIdentity']:.1f}% | {liability} |")
    for extra_gene, liability in known.items():
        if not any(h["targetGeneSymbol"] == extra_gene for h in paralogs):
            lines.append(f"| {extra_gene} | n/a (not a sequence paralog) | {liability} |")
    lines.append("")

    lines.append("## Expression / localization")
    lines.append("")
    if profile.get("subcellularLocations"):
        locs = ", ".join(sorted({l["location"] for l in profile["subcellularLocations"]}))
        lines.append(f"- Subcellular location: {locs}")
        lines.append("")
    lines.append("Top tissues by baseline (bulk RNA) expression:")
    lines.append("")
    lines.append("| Tissue | Median expression | Unit |")
    lines.append("|---|---|---|")
    rows = [r for r in profile.get("baselineExpression", {}).get("rows", []) if r.get("median") is not None]
    seen = set()
    dedup = []
    for r in sorted(rows, key=lambda r: r["median"], reverse=True):
        name = r["tissueBiosample"]["biosampleName"] if r["tissueBiosample"] else "unknown"
        if name in seen:
            continue
        seen.add(name)
        dedup.append((name, r["median"], r["unit"]))
        if len(dedup) >= 10:
            break
    for name, median, unit in dedup:
        lines.append(f"| {name} | {median:.2f} | {unit} |")
    lines.append("")
    lines.append(
        "_Caveat: this is bulk GTEx-derived expression and does not include dorsal root "
        "ganglion (DRG) or peripheral sensory neuron / skin nerve-ending samples "
        "specifically -- the tissue most relevant to a topical analgesic's site of action. "
        "Cross-check against DRG-specific expression atlases from the literature before "
        "concluding local target abundance at the intended application site._"
    )
    lines.append("")

    if MANUAL_CAVEATS.get(gene_symbol):
        lines.append("## Additional caveats (literature, not in Open Targets)")
        lines.append("")
        for note in MANUAL_CAVEATS[gene_symbol]:
            lines.append(f"- {note}")
        lines.append("")

    lines.append("## Associated diseases (top by overall association score)")
    lines.append("")
    lines.append("| Disease | Score |")
    lines.append("|---|---|")
    diseases = sorted(
        profile.get("associatedDiseases", {}).get("rows", []), key=lambda r: r["score"], reverse=True
    )[:10]
    for d in diseases:
        lines.append(f"| {d['disease']['name']} | {d['score']:.3f} |")
    lines.append("")

    return "\n".join(lines)


def main(gene_symbols):
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    for gene in gene_symbols:
        print(f"Fetching {gene}...")
        ensembl_id = resolve_ensembl_id(gene)
        profile = fetch_target_profile(ensembl_id)
        (RAW_DIR / f"{gene}.json").write_text(json.dumps(profile, indent=2))
        report = render_markdown(profile, gene)
        (REPORT_DIR / f"{gene}.md").write_text(report)
        print(f"  -> {REPORT_DIR / f'{gene}.md'}")


if __name__ == "__main__":
    genes = sys.argv[1:] or ["SCN9A", "TRPV1"]
    main(genes)
