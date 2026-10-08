# Off-target safety screen for a topical TRPV1 agonist — short summary

*Companion summary to the full manuscript
(`OFFTARGET_SAFETY_ARTICLE_EN.md`), prepared for co-author review.
Code and supporting materials:
https://github.com/danil25-stack/trpv1-offtarget-safety-screen*

## The question

Off-target safety review for a small molecule normally stops at the
target's sequence-paralog family. For a topical TRPV1 agonist that leaves
two gaps. Whole-protein sequence identity is a poor proxy for similarity at
the binding site that actually matters pharmacologically; and sequence
relatedness says nothing about *unrelated* proteins that happen to combine
a geometrically similar pocket with co-expression at the drug's local site
of action. This screen asks both questions computationally, cheaply enough
to run before a lead compound exists.

## What was done

1. **Paralog family (TRPV2–6, plus TRPA1).** Curated identity, liabilities
   and tractability from Open Targets, with a literature review separating
   whole-protein identity from binding-site-level conservation.
2. **Outside the family.** Human Protein Atlas skin-expression filter (604
   candidate genes) → Open Targets druggability tiering (486 excluded as
   structural/barrier proteins with no pharmacological pocket; 26 of the
   remaining 118 had confirmed pocket evidence) → fpocket geometric pocket
   detection and an eight-descriptor pocket-fingerprint comparison against
   TRPV1's vanilloid pocket. The fpocket rule was validated first on TRPV1
   itself, where its own top-ranked cavity independently recovered all four
   literature-established vanilloid-pocket residues (Y511, S512, T550,
   E570).
3. **Cross-docking.** The programme's seven-ligand TRPV1 reference series
   (resiniferatoxin → nonivamide, ~5 orders of magnitude in EC50) docked
   with AutoDock Vina into the flagged off-targets, with per-target
   redocking validation of the protocol.
4. **Ground-truth discrimination check.** Property-matched panels of 30
   known actives and 30 known inactives per off-target, pulled from ChEMBL
   and docked into the same receptor and box, to ask whether the scoring
   function can separate real binders from non-binders *at all* on that
   pocket — a precondition for trusting any cross-receptor comparison.
5. **ESM-2 pocket comparison.** Local protein-language-model
   representations restricted to the four pocket-lining residues, mapped
   across paralogs by alignment rather than by residue number, as a
   quantitative test of the qualitative paralog claim.

## Key results

Within the family, identity-ranking misorders the real binding-site risk:
the vanilloid-pocket region is reported conserved across TRPV1–3 regardless
of overall identity (and is experimentally transplantable), so TRPV2/TRPV3
pocket risk may exceed what 45.0% and 39.0% identity suggest. TRPV6 carries
the most systemically consequential profile if hit (knockout: impaired
intestinal Ca²⁺ absorption, reduced bone mineral density,
alopecia/dermatitis, severe male infertility).

Outside the family, **RARG — retinoic acid receptor gamma — is the finding.**
It has no meaningful sequence identity or evolutionary relationship to
TRPV1, so a paralog-only review would never have considered it. Six of the
seven reference ligands dock into RARG's retinoid pocket with scores
numerically comparable to, or slightly more favourable than, their scores
against TRPV1 itself; resiniferatoxin, the bulkiest ligand, is the
consistent outlier.

| Target | Scoring | n (active/inactive) | ROC-AUC | 95% CI | *p* vs 0.5 |
|---|---|---|---|---|---|
| RARG | Vina | 30/30 | **0.811** | 0.686–0.914 | <0.0001 |
| MMP3 | Vina | 30/30 | 0.503 | 0.351–0.654 | 0.485 |
| MMP3 | AD4Zn | 30/30 | 0.551 | 0.394–0.694 | 0.251 |

The discrimination check is what separates the two candidates. RARG passes
it convincingly, on a compound set and a question entirely unrelated to the
TRPV1 series — so the RARG cross-docking signal reflects the scoring
function genuinely working on that pocket. MMP3 fails it outright, and
stays near chance even under AD4Zn, a force field built specifically for
zinc coordination. That rules out "the catalytic zinc is being scored as an
inert rigid body" as the sole explanation and means the MMP3 cross-docking
numbers cannot be read as evidence in either direction.

RARG's clinical relevance is anchored rather than hypothetical: the
structure used for the comparison (PDB 6FX0) comes from the medicinal
chemistry programme that produced trifarotene, an approved RARγ-selective
topical retinoid, and tazarotene is an older topical precedent. So
small-molecule occupancy of this pocket in this anatomical site is the
mechanism of at least one approved topical drug. The signal also spans most
of the series rather than concentrating in the smallest ligands, so the
MW/logP window favoured for topical delivery (200–450 Da, logP 1.0–4.0)
would not filter this risk out.

## What changed since the August draft

- The ground-truth discrimination check and the AD4Zn re-run are new, and
  they **downgrade MMP3** from a co-priority target to a methodologically
  unresolved watch item.
- The MMP3 receptor was corrected to a holo structure (PDB 1HFS) retaining
  all resolved Zn²⁺ and Ca²⁺ ions; the earlier apo run understated the
  MMP3 signal.
- The ESM-2 comparison was added. It supports the TRPV5/TRPV6 part of the
  paralog claim (0.93–0.94 versus 0.965–0.969 for TRPV2–4) but **not** the
  TRPV4-specific part, and is presented with that limitation stated.
- All references were verified against Crossref/PubMed; three were wrong
  (notably the TRPV4/GSK2798745 paper, which is Goyal et al., *Am J
  Cardiovasc Drugs* 2019;19(3):335–342, not Goldsmith in *Clin
  Pharmacokinetics*), eleven placeholders were filled, and two statements
  that overstated their sources were corrected — the TRPA1 entry (one
  programme with a hepatotoxicity signal, not two with acceptable safety)
  and the BLMH/flagellate-dermatitis link (a hypothesis, not an established
  mechanism).

## Where a second opinion would help most

1. Is the recommendation pitched correctly — RARG as a counter-screen worth
   running early, rather than as a demonstrated liability?
2. Is reporting MMP3 as an unresolved watch item the right call, or should
   it be dropped from the manuscript altogether?
3. Journal fit for a hypothesis-generating computational triage paper of
   this size.

## Limitations in brief

This is a hypothesis-generating computational screen, not a confirmatory
safety study: no experimental cross-reactivity data were generated.
Cross-receptor Vina score comparison is inherently less reliable than
single-target ranking; both receptor structures are single static
snapshots; the eight-descriptor fingerprint is a coarse proxy for true
pocket shape and electrostatics; bulk tissue RNA expression does not
resolve dorsal-root-ganglion or nerve-terminal expression specifically; and
only the highest-confidence tier (26 of 118 genes) was screened, leaving
several plausible pain/inflammation-relevant candidates (PTGS1, LTB4R/
LTB4R2, HCAR2/HCAR3) unevaluated.
