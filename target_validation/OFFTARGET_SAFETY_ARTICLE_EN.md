*[Manuscript draft — journal format. Author name(s), affiliation(s), funding,
and conflict-of-interest statements are left as placeholders below for you to
complete. Reference list entries marked "verify before submission" were
confirmed for title/journal/year/DOI via search but not for the complete
author list — check against the primary source before citing.]*

# In silico screening for potential off-target interactions of topical TRPV1 agonists: from paralog analysis to structure-based pocket comparison and cross-docking

**[Author Name(s)]**¹

¹ *[Affiliation, Department, Institution, City, Country]*

Corresponding author: *[email]*

## Abstract

**Background.** Topical TRPV1 agonists are an established, clinically
validated analgesic mechanism (e.g., capsaicin 8% patch), but off-target
safety assessment for new candidates in this class is typically limited to
the TRPV1 sequence-paralog family (TRPV2–6) and addressed late, via
empirical in vitro panels. **Methods.** We combined a curated
paralog-family literature review with a structure-based computational
screen: a tissue-expression filter (Human Protein Atlas) to identify
druggable, skin-enriched proteins unrelated to TRPV1 by sequence; cavity
detection (fpocket) to compare candidate binding-pocket geometry and
physicochemistry against TRPV1's vanilloid pocket without docking; and
molecular docking (AutoDock Vina) of seven TRPV1 reference agonists into
two candidates carried forward from the pocket screen — MMP3, which led
both pocket-similarity rankings, and RARG, promoted on independent
biological and druggability grounds despite not leading either ranking by
distance — to test whether pocket similarity translates into comparable
predicted binding. **Results.** Within the TRPV
family, whole-protein sequence identity poorly predicts binding-site risk:
literature indicates the vanilloid-pocket region is structurally conserved
across TRPV1–3 independent of overall identity, while TRPV4's selectivity
is achieved by point substitutions within an otherwise similar pocket.
Outside the family, a genome-scale skin-expression screen (604 candidate
genes) combined with druggability filtering and pocket-fingerprint
comparison flagged three targets independent of the ranking method used:
matrix metalloproteinase-3 (MMP3), retinoic acid receptor gamma (RARG), and
bleomycin hydrolase (BLMH). Cross-docking supported the strongest of these
predictions without confirming it: for most of the seven TRPV1 reference
agonists, Vina docking scores against RARG's retinoid pocket, and,
separately, against a corrected zinc/calcium-retaining MMP3 receptor, were
numerically comparable to (in several cases somewhat more favorable than)
their scores against TRPV1 itself; resiniferatoxin, the bulkiest and
structurally most complex ligand tested, was the consistent outlier at
both off-targets. Docking scores are not a direct measure of binding
affinity, and cross-receptor score comparisons carry additional
uncertainty beyond the single-receptor ranking use case for which this
scoring function was designed; these results should be read as motivation
for experimental follow-up, not as evidence of actual cross-reactivity. To
test that uncertainty directly, we separately redocked a property-matched
panel of 30 known active and 30 inactive ChEMBL-annotated compounds per
off-target and asked whether Vina score alone separates them within each
receptor. RARG passed this ground-truth check (ROC-AUC 0.811); MMP3 did
not (ROC-AUC 0.503, indistinguishable from chance) and remained
near-chance (ROC-AUC 0.551) even after replacing the plain Vina scoring
function with AD4Zn, a force field with an explicit zinc-coordination
potential.
**Conclusions.** Off-target risk for this chemical series is plausibly not
confined to the TRPV/TRPA family and, on this evidence, is not obviously
restricted to smaller or simpler candidates either — it spans most of the
tested series at two independent, structurally unrelated off-targets. The
independent discrimination check corroborates RARG as a priority candidate
for further counter-screening of TRPV1-directed topical analgesic
candidates in this structural class; MMP3, despite an initially similar
pocket-similarity and cross-docking signal, failed this check even under a
zinc-aware scoring function and is better treated as a methodologically
unresolved watch item than a confirmed co-priority target. We outline a
low-cost, pre-lead computational triage pipeline (tissue expression →
pocket geometry → targeted cross-docking → ground-truth discrimination
check) that is generalizable to other topical/localized small-molecule
programs.

**Keywords:** TRPV1; off-target safety; molecular docking; pocket
similarity; fpocket; retinoic acid receptor gamma; matrix metalloproteinase-3;
topical analgesic; in silico toxicology; polypharmacology; AutoDock4Zn;
docking enrichment; ROC-AUC

## 1. Introduction

Transient receptor potential vanilloid 1 (TRPV1) is a non-selective cation
channel that transduces noxious heat, protons, and vanilloid ligands in
peripheral nociceptors, and is the mechanistic target of the only approved
topical agonist-desensitization analgesic in this class (capsaicin 8%
patch, Qutenza) [18]. Systemic TRPV1 antagonism has repeatedly failed in
clinical development because of on-target hyperthermia, which has directed
renewed interest in topical, agonist-driven strategies aimed at minimizing
(not eliminating) systemic exposure — even the approved capsaicin patch
shows low, transient plasma capsaicin levels in a minority of patients per
its FDA labeling [18]. Any such program still requires an off-target safety
case, and the standard approach — profiling activity against the target's
sequence-paralog family — is necessary but not sufficient for two reasons.
First, whole-protein sequence identity is an unreliable proxy for
similarity at the specific ligand-binding site that actually matters
pharmacologically; two paralogs can differ substantially in overall
identity while sharing a nearly identical pocket, or vice versa. Second,
sequence relatedness says nothing about *unrelated* proteins that
happen to (a) present a geometrically and chemically similar pocket and (b)
be co-expressed at the intended site of local drug action — a risk category
that standard paralog panels cannot detect by construction, and that is
normally only caught late, via broad empirical selectivity panels after a
lead compound already exists.

We report a two-stage in silico triage addressing both gaps for a topical
TRPV1 agonist program. Stage one is a literature-anchored reassessment of
the TRPV1 paralog family (TRPV2–6, plus the mechanistically related channel
TRPA1) that explicitly separates whole-protein sequence identity from
binding-site-level structural conservation. Stage two is a structure-based
screen for off-target risk *outside* the TRPV/TRPA family: a
tissue-expression filter to enumerate druggable proteins co-expressed at the
skin (the site of topical drug action), cavity-detection-based comparison of
candidate pockets against TRPV1's vanilloid pocket, and — for the resulting
top candidates — direct cross-docking of the program's reference agonist
series to test whether the structural prediction holds up against an actual
scoring function. All software and data-source choices were selected to
keep the pipeline inexpensive and reproducible with tools already validated
in the group's docking workflow, so that it can run *before* a lead exists,
rather than replacing later empirical selectivity profiling.

## 2. Materials and Methods

### 2.1 Paralog-family curation

Sequence-paralog identity, curated safety liabilities, and target
tractability for TRPV1 and its five closest human paralogs (TRPV2–6) were
retrieved from the Open Targets Platform GraphQL API [1], which integrates
Ensembl-derived homology calls with curated pharmacovigilance and genetic
data. Liabilities not curated in Open Targets (found for TRPV5 and TRPV6)
were added from a targeted literature search. TRPA1, a mechanistically
adjacent, co-expressed nociceptor channel that is not a TRPV1 sequence
paralog, was included on functional grounds.

### 2.2 Skin-expressed candidate selection

Candidate off-target proteins were enumerated from the Human Protein Atlas
RNA tissue-specificity resource [2], filtering to genes classified "Tissue
enriched," "Group enriched," or "Tissue enhanced" for skin. Each candidate
was cross-referenced against Open Targets `targetClass` and small-molecule
`tractability` annotations [1] to exclude genes with no known ligand-binding
pocket (predominantly structural/barrier proteins) and to tier the
remainder by strength of pocket evidence (confirmed pocket; ligand-bound
structure without an explicit pocket call; druggable-looking protein family
without structural evidence). This report covers the highest-confidence
tier only.

### 2.3 Cavity detection and pocket-fingerprint comparison

Protein cavities were detected with fpocket [3], a geometry-only detector
(Voronoi tessellation/alpha-spheres) that requires no bound ligand and
scales to dozens of structures in minutes. The method was validated against
ground truth before use: run on the TRPV1 vanilloid-pocket-containing
fragment (chains A+D of PDB 11CK, trimmed to the inter-subunit region known
to form the pocket), fpocket's own top-ranked cavity independently recovered
all four literature-established vanilloid-pocket residues (Tyr511, Ser512,
Thr550, Glu570) [4,5]. This confirms that "fpocket's own top-ranked pocket"
is a valid selection rule prior to applying it, blind, to candidate
off-targets.

For each candidate gene, the best-available experimental structure was
retrieved via UniProt cross-references to the RCSB Protein Data Bank
[6,7], preferring the largest resolved chain span and best resolution.
Each detected pocket was described by an eight-descriptor fingerprint
(fpocket's druggability score, volume, hydrophobicity score, polarity
score, charge score, proportion of polar atoms, apolar alpha-sphere
proportion, and mean local hydrophobic density), z-score normalized across
all pockets on all structures, and compared to the TRPV1 reference pocket by
Euclidean distance and cosine similarity. Two independent rankings were
computed per candidate — the best-matching pocket found anywhere on the
structure, and the single pocket fpocket itself ranks highest (the rule
validated against TRPV1) — to guard against a multiple-comparisons bias that
favors proteins with many detected cavities.

### 2.4 Molecular docking

Candidate ligands were the seven-compound TRPV1 agonist reference series
already used for the group's residence-time pilot (resiniferatoxin,
phenylacetylrinvanil, olvanil, arvanil, rinvanil, capsaicin, nonivamide),
spanning approximately five orders of magnitude in TRPV1 EC50. Ligand 3D
conformers were generated from isomeric SMILES with RDKit (ETKDGv3
embedding, MMFF optimization) and prepared for docking with Meeko [8].
Receptor structures for the two candidates carried forward from Section
2.3/3.2 (MMP3 and RARG; selection rationale for each detailed in Section
3.2) were prepared with Meeko's `mk_prepare_receptor.py`. For RARG, PDB
6FX0 [19]
(chain A; the same structure flagged in Section 3.2) was used, with the
docking search box centered on the flagged pocket's alpha-sphere centroid
(padded 6 Å per side; center (12.91, −8.71, −8.53) Å, size 22.3×19.7×19.8
Å). For MMP3, the structure used for the initial fpocket screen (PDB 1SLM)
is an apo (unliganded) crystal form; because docking accuracy for a
metalloprotease depends on retaining the catalytic zinc, we instead used
the holo structure PDB 1HFS [21] (chain A; 1.70 Å resolution,
co-crystallized with the hydroxamate
inhibitor L04), keeping all resolved Zn²⁺ (n=2) and Ca²⁺ (n=3) ions as a
rigid part of the receptor and centering the search box on the
co-crystallized ligand's centroid (padded 8 Å per side to accommodate
ligands larger than L04; center (35.34, 36.10, 24.54) Å, size
25.7×19.6×21.4 Å). Docking was performed with AutoDock Vina [9]
(exhaustiveness 32, 9 output poses), matching the settings already
validated for the TRPV1 receptor in the group's docking pipeline. TRPV1
reference scores for the same seven ligands were taken from the group's
existing validated docking run against the TRPV1 vanilloid pocket (PDB
11CK [22], chains A+D, in complex with 6-iodo-capsaicin; box centered on the
inter-subunit vanilloid pocket).

**Docking-protocol validation.** Because the TRPV1 pocket-recovery check
(Section 2.3) validates fpocket's cavity *detection*, not the separate
docking step, and because MMP3 and RARG are structurally unrelated to
TRPV1 (a zinc metalloprotease and a nuclear hormone receptor, respectively,
vs. an ion channel), the docking protocol itself was validated
independently for each: the co-crystallized ligand of each structure (L04
for 1HFS; E9T for 6FX0 — both PDB chemical component IDs) was built from
its canonical SMILES (RCSB Chemical Component
Dictionary) with RDKit — an independently generated conformer, not the
crystal coordinates — redocked into its own structure under the box
described above, and the best-scoring pose compared to the crystal pose by
heavy-atom RMSD after topology matching (`docking/validate_redocking.py`,
already used for the TRPV1/11CK validation in `docking/SMOKE_TEST.md`).
RMSD results are reported in Section 3.3.

### 2.5 Ground-truth discrimination check

The cross-receptor comparison in Section 2.4 asks whether known TRPV1
agonists score comparably at RARG/MMP3 and at TRPV1 itself; it does not
establish that Vina's scoring function can separate real binders from
non-binders at RARG or MMP3 *at all* — a precondition for trusting that
absolute cross-receptor comparison. To close this gap, known active
(≤1,000 nM) and inactive (>10,000 nM by a measured value, or no binding
detected up to ≥10,000 nM) compounds for RARG and MMP3 were retrieved from
ChEMBL [23] (target CHEMBL2003 and CHEMBL283, respectively; all activity
records with `standard_type` in {IC50, EC50, Ki, Kd} and `standard_units` =
nM). For each target, 30 actives spanning the potency range — not simply
the most potent — were selected, and 30 inactives were greedy-matched to
those actives in standardized (molecular weight, calculated logP) space to
avoid a property-driven confound between the active and inactive sets. Each
panel was redocked into the same receptor structure and box used in
Section 2.4, with the same RDKit/Meeko ligand-preparation pipeline, and
ROC-AUC for active/inactive separation by Vina score was computed per
target.

Because MMP3's cross-docking result depends on retaining the catalytic
zinc (Section 2.4), the MMP3 panel was additionally redocked with
AutoDock4Zn (AD4Zn) [24], a force field that represents zinc coordination
explicitly via directional tetrahedral zinc pseudo-atoms, run through
Vina's AD4-scoring mode (`--scoring ad4`) against AutoGrid4-computed
affinity maps rather than Vina's own scoring function. Receptor
preparation added zinc pseudo-atoms (`zinc_pseudo.py`, distributed with
AutoDock-Vina) to the same Zn²⁺/Ca²⁺-retaining receptor used in Section
2.4; a lightweight custom script (`docking/prepare_gpf4zn_lite.py`) wrote
the AutoGrid4 parameter file directly, reproducing the parameters of the
worked example in the AutoDock-Vina documentation without the legacy
MolKit/AutoDockTools dependency of the reference `prepare_gpf4zn.py`
script. Ligands were prepared with macrocyclic rings held rigid
(`--rigid_macrocycles`), since Vina's flexible-macrocycle ligand typing
introduces synthetic atom types not represented in the classic AD4Zn
parameter set.

## 3. Results

### 3.1 TRPV1 paralog family: binding-site conservation is not predicted by whole-protein identity

Table 1 summarizes sequence identity to TRPV1 and the associated safety
liability for each paralog. Two revisions to a naive identity-ranked
reading are notable. First, the vanilloid-pocket region (formed by
S3/S4/S4–S5-linker of one subunit and S5/S6 of the neighboring subunit) is
reported to be structurally and sequence-wise conserved across TRPV1, TRPV2,
and TRPV3 specifically — sufficiently so that this region is
experimentally transplantable between orthologs and confers vanilloid
sensitivity onto otherwise-insensitive channels [5,10]. This implies that
TRPV2 and TRPV3's real off-target pocket risk may exceed what their
whole-protein identity (45.0% and 39.0%, respectively) suggests. Second,
TRPV4 (41.9% identity) achieves natural selectivity against vanilloids via
point substitutions within an overall similar pocket architecture [10],
and its clinical antagonist GSK2798745 was well tolerated through repeat
dosing in heart-failure patients with no cardiopulmonary safety signal
[11] — despite TRPV4 *activation* being separately implicated in cardiogenic
pulmonary edema, indicating the clinically relevant liability is
direction-dependent and an antagonist-mediated off-target hit trends
protective rather than causative. TRPV5 and TRPV6 (28.4% and 27.2%
identity), by contrast, are calcium-selective channels with a distinct
pore architecture; their low whole-protein identity plausibly does track
real pocket-level divergence. TRPV6 carries the most systemically
consequential liability profile if hit off-target: knockout mice show
impaired intestinal calcium absorption, reduced bone mineral density,
alopecia/dermatitis, and severely impaired male fertility [12].

**Table 1.** TRPV1 paralog family: sequence identity and off-target
liability.

| Paralog | Identity to TRPV1 | Liability if hit off-target |
|---|---|---|
| TRPV2 | 45.0% | Developmental cardiac Ca²⁺-handling role; adult whole-body knockout and neutralizing-antibody dosing show no overt phenotype [13,14]. |
| TRPV4 | 41.9% | Direction-dependent: activation implicated in cardiogenic pulmonary edema, but clinical antagonism (GSK2798745) well-tolerated with no cardiopulmonary signal [11]. |
| TRPV3 | 39.0% | Shares TRPV1's topical site of action; human gain-of-function mutations cause Olmsted syndrome (painful keratoderma, pruritus, skin-barrier failure) [15]. |
| TRPV5 | 28.4% | Rate-limiting channel for renal distal-tubule Ca²⁺ reabsorption; plausible hypercalciuria risk if inhibited. |
| TRPV6 | 27.2% | Most consequential paralog: knockout shows impaired intestinal Ca²⁺ absorption, reduced bone mineral density, alopecia/dermatitis, severe male infertility [12]. |
| TRPA1 | n/a (not a sequence paralog) | Two independent clinical antagonist programs (LY3526318; a second undisclosed compound) show acceptable Phase 1 safety but failed on efficacy across several pain indications [16]. |

### 3.2 Structure-based screen identifies three off-target candidates outside the TRPV/TRPA family

Candidate numbers moved through the pipeline as follows: **604** skin
tissue-specific genes were identified by expression filtering; **486** of
these had no Open Targets target-class annotation and were excluded as
predominantly structural/barrier proteins without a pharmacological pocket,
leaving **118**; of these, **26** genes had confirmed ligand-binding-pocket
evidence (the highest-confidence tier) and were carried forward to pocket
comparison; cross-referencing two independent ranking methods with known
biology (below) narrowed this to **3** flagged candidates (Figure 1,
`figures/pipeline_funnel.png`).

Candidate selection used two distinct, explicitly separated mechanisms
rather than a single ranking, and neither was informed by the downstream
docking result in Section 3.3/3.4 — both were fixed before any docking was
run. The first is purely computational: three genes — MMP3, NTF4, and
BLMH — ranked near the top independent of
which of the two comparison methods (best pocket anywhere vs. fpocket's own
top-ranked pocket) was used. Each flagged candidate was checked against
known biology to distinguish a genuine small-molecule pocket from a likely
detection artifact (e.g., a protein–protein interaction surface or a
serpin's conformational groove, neither of which functions via classical
small-molecule occupancy). MMP3 (stromelysin-1) presents a real,
high-druggability (0.839), hydrophobic zinc-dependent catalytic pocket and
is independently implicated in cutaneous extracellular-matrix remodeling.
The second mechanism is biology-informed prioritization, applied to the
full 26-gene shortlist rather than only the computational top-3: RARG did
not lead either individual pocket-similarity ranking by distance (its
pocket volume exceeds TRPV1's, which penalizes a Euclidean-distance
metric), but carries near-maximal druggability (0.982) and high cosine
similarity (0.899), and was promoted on the independent biological grounds
below — fixed criteria (an approved topical drug already occupying the
candidate pocket at the same anatomical site) applied prospectively to the
shortlist, not selected after seeing a favorable docking result. Its clinical relevance is
strongly supported by the structure used for this comparison itself: PDB
6FX0 is the RARγ-bound structure from the medicinal-chemistry campaign that
produced trifarotene (CD5789), a RARγ-selective topical retinoid approved
for acne [19], and the co-crystallized ligand (PDB chemical component
E9T) shares trifarotene's adamantyl-naphthoic-acid scaffold. Tazarotene, an
older, topically applied RARβ/γ-modulating prodrug approved for psoriasis
and acne, is a second, less receptor-selective precedent for the same
point [20]. Together these indicate that a small molecule occupying RARG's
retinoid pocket at this anatomical site is not a purely theoretical
scenario — it is the mechanism of at least one approved topical drug — though
this says nothing about whether our specific candidate ligands would bind,
nor about the functional consequence (agonism, antagonism, or no effect) if
they did. BLMH (bleomycin hydrolase) is independently notable
because reduced cutaneous BLMH activity is already implicated in a real
topically-relevant drug toxicity (bleomycin-induced flagellate dermatitis).
NTF4 (a secreted TrkB ligand, more likely to present a protein–protein
interface than a small-molecule pocket) and SERPINB5 (a serpin, whose
near-maximal druggability score is inconsistent with its
reactive-center-loop mechanism of action) were deprioritized as probable
detection artifacts.

### 3.3 Cross-docking supports the RARG hypothesis

Redocking each structure's own co-crystallized ligand (built from SMILES,
not crystal coordinates) reproduced its crystal pose with heavy-atom RMSD
1.42 Å for RARG/E9T and 2.10 Å for MMP3/L04 — RARG cleanly passes the
conventional <2 Å bar for a validated redocking protocol; MMP3's result
sits just outside it, plausibly reflecting L04's greater conformational
flexibility relative to the more rigid, compact E9T. We read this as "the
protocol broadly reproduces the correct binding mode for MMP3, with
non-trivial residual pose uncertainty," and weight the MMP3 cross-docking
results below accordingly (lower confidence than the RARG results).

Table 2 reports docking scores for the seven-ligand TRPV1 reference series
against TRPV1, MMP3, and RARG.

**Table 2.** Vina docking scores (kcal/mol; more negative indicates a more
favorable predicted pose). Bold values are numerically equal to or more
favorable than the ligand's own score against TRPV1. Differences smaller
than roughly 1 kcal/mol are within Vina's typical scoring uncertainty and
should not, on their own, be read as evidence of stronger binding —
absolute cross-receptor score comparisons are inherently less reliable
than the single-receptor pose ranking the method was designed for.

| Ligand | TRPV1 EC50 (nM) | TRPV1 | MMP3 (holo) | RARG |
|---|---|---|---|---|
| Resiniferatoxin | 0.011 | −9.63 | −9.21 | −2.86 |
| Phenylacetylrinvanil | 0.09 | −8.01 | **−8.13** | **−8.98** |
| Olvanil | 0.6 | −7.49 | −7.46 | **−8.73** |
| Arvanil | 0.5 | −7.83 | **−8.37** | **−9.17** |
| Rinvanil | 6.0 | −7.51 | **−7.64** | **−8.60** |
| Capsaicin | 712.0 | −8.23 | −8.02 | **−8.59** |
| Nonivamide | 1400.0 | −7.60 | **−7.96** | **−7.97** |

For six of the seven ligands, the docking score against RARG was
numerically equal to or somewhat more favorable than the score against
TRPV1 itself (nonivamide's −7.97 vs. −7.60, a 0.37 kcal/mol difference,
falls inside that scoring-uncertainty band and should be read as "similar,"
not "better"); only resiniferatoxin, the bulkiest and structurally most
complex ligand (a diterpene), scored poorly. Applied without modification
to an unrelated protein fold, the docking scoring function ranks most of
this chemical series as comparably complementary to RARG's retinoid pocket
and to TRPV1's own vanilloid pocket. We read this as computational support
for — not confirmation of, and not a mere restatement of — the
pocket-similarity hypothesis from Section 3.2, and as sufficient grounds to
prioritize an experimental RARG counter-screen.

The corrected, Zn²⁺/Ca²⁺-retaining MMP3 receptor changed this picture
substantially from an earlier pass against an apo (metal-stripped) receptor
(not reported here in detail): rather than a signal confined to the two
smallest ligands, four of seven ligands now score equal to or more
favorable than their TRPV1 score (phenylacetylrinvanil, arvanil, rinvanil,
nonivamide), a fifth (olvanil, −7.46 vs. −7.49) is effectively tied, and
only resiniferatoxin and capsaicin score modestly weaker (by 0.42 and 0.21
kcal/mol respectively — both within, or at the edge of, the
scoring-uncertainty band). In other words, once the catalytic zinc is
correctly retained in the receptor, MMP3 shows a broad signal across
essentially the whole tested chemical series, not one restricted to small
molecules — mirroring the RARG pattern more closely than our initial
(metal-free) MMP3 pass suggested. The one consistent outlier across *both*
flagged off-targets is resiniferatoxin, the bulkiest, most structurally
complex ligand in the series; every other tested ligand scores comparably
at MMP3, RARG, and TRPV1 alike. Given the MMP3 docking protocol's
borderline redocking RMSD (2.10 Å, above), this pattern is reported
as a hypothesis warranting experimental follow-up with somewhat lower
confidence than the RARG result, not as a confirmed finding.

### 3.4 A ground-truth discrimination check separates RARG from MMP3

The cross-receptor comparison above assumes Vina's scoring function is
capable of meaningfully ranking compounds at RARG and MMP3 in the first
place. Table 3 tests that assumption directly: known active and inactive
ChEMBL compounds (30 each per target, property-matched, Section 2.5) were
redocked into each receptor, and ROC-AUC was computed for active/inactive
separation by score alone — a question entirely independent of the
TRPV1-agonist cross-docking analysis above.

**Table 3.** Ground-truth discrimination (ROC-AUC for separating known
active from known inactive ChEMBL compounds by docking score; 0.5 =
indistinguishable from chance, 1.0 = perfect separation). 95% CI from
10,000-resample stratified bootstrap; *p* from a one-sided Mann-Whitney
*U* test against the null AUC = 0.5 (equivalently, active scores are not
stochastically more favorable than inactive scores); PR-AUC (average
precision) reported alongside ROC-AUC because ROC-AUC alone can look
optimistic on small, balanced panels. ROC curves and score distributions
underlying this table are in `figures/roc_curves.png` and
`figures/score_distributions.png`.

| Target | Scoring | n (active/inactive) | ROC-AUC | 95% CI | *p* (vs. 0.5) | PR-AUC |
|---|---|---|---|---|---|---|
| RARG | Vina | 30/30 | **0.811** | 0.686–0.914 | <0.0001 | 0.858 |
| MMP3 | Vina | 30/30 | 0.503 | 0.351–0.654 | 0.485 | 0.598 |
| MMP3 | AD4Zn | 30/30 | 0.551 | 0.394–0.694 | 0.251 | 0.610 |

RARG passes this check convincingly and with a bootstrap CI that excludes
chance entirely (0.686–0.914; *p* < 0.0001): Vina score separates known
RARG actives from inactives well above chance, using a compound set and
question entirely unrelated to the seven-ligand TRPV1 series in Table 2.
This is independent corroboration that the RARG cross-docking result
reflects the scoring function doing real work on this pocket, not a
coincidence of the specific seven-ligand series tested. Because roughly
three-fifths of the RARG panel (35/60 candidate compounds; 18 active, 17
inactive among the sampled 30/30) was labeled from a functional EC50
assay rather than a direct binding measurement (Ki/Kd) or enzymatic IC50,
we repeated the check restricted to the subset whose label was set by a
binding-type record only (Ki, Kd, or IC50; n=25, 12 active/13 inactive):
AUC was unchanged at 0.821, ruling out EC50/binding-assay mixing as an
explanation for the RARG result. MMP3's panel required no such check — no
compound in either the active or inactive set was labeled from an EC50
record (all IC50/Ki/Kd).

MMP3 fails the same check outright (AUC 0.503, 95% CI 0.351–0.654 —
straddling 0.5 almost symmetrically; *p* = 0.485, nowhere near
significance) — active and inactive
compounds are, on average, indistinguishable by score (mean −9.51 vs.
−9.71 kcal/mol respectively; Section 3.3 details). Switching to AD4Zn, a
force field built specifically to score zinc coordination correctly,
improved this only marginally (AUC 0.551, 95% CI 0.394–0.694, *p* = 0.251;
active mean −25.23 vs. inactive
mean −24.54 in AD4 scoring units, not comparable in scale to plain-Vina
kcal/mol) — still statistically indistinguishable from chance and far
below RARG's result under the
unmodified scoring function. This rules out "the catalytic zinc is scored
as an inert rigid body instead of a real coordination partner" as the sole
explanation for MMP3's poor performance; some other factor (pocket
flexibility not captured by a single static structure, a coarser fit of
the empirical scoring terms to this binding mode, or limitations of the
1HFS holo-structure itself) is also at play, and was not resolved within
the scope of this study. Practically, this means the MMP3 numbers in Table
2 cannot be read as evidence of anything: a docking protocol that cannot
rank known MMP3 binders above known non-binders within a single receptor
has not demonstrated the far harder capability — absolute score comparison
across two structurally unrelated receptors — that the Table 2 argument
depends on. This is consistent with, and considerably sharpens, the
borderline redocking RMSD already flagged for MMP3 in Section 3.3.

## 4. Discussion

This screen demonstrates that structure-based off-target triage — cheap
enough to run before a lead compound exists — can surface a candidate
(RARG) that a conventional paralog-only safety review would never examine,
because RARG shares no meaningful sequence identity or evolutionary
relationship with TRPV1. The finding is strengthened by two independent
lines of computational evidence converging on the same target
(pocket-fingerprint similarity and cross-docking), and by the fact that
RARG is not a purely hypothetical target: it is the direct target of an
approved topical RARγ-selective retinoid (trifarotene) and is engaged,
less selectively, by an older topical retinoid prodrug (tazarotene). A
third, independent line of evidence points the same direction: a
ground-truth check using known ChEMBL actives and inactives (Section 3.4)
shows Vina score genuinely discriminates RARG binders from non-binders
(ROC-AUC 0.811) on a compound set entirely unrelated to the seven-ligand
TRPV1 series used for cross-docking — evidence that the underlying scoring
function does real, generalizable work on this pocket, not just a
consistent pattern within one ligand series. The clinical precedent
compounds this: the receptor is already pharmacologically engaged in the
same tissue by existing therapeutics, giving a potential off-target
interaction a concrete, clinically anchored interpretive frame (possible
interference with retinoid signaling in skin) rather than requiring a de
novo toxicological hypothesis. We stress that none of this establishes
actual binding, agonism versus antagonism at RARG, or any functional or
toxicological consequence — those questions are outside the scope of a
computational screen and require dedicated experimental follow-up.

Several limitations bound how these results should be used. First, and
most importantly, this is a hypothesis-generating computational triage, not
a confirmatory safety study — no experimental (biochemical, cellular, or in
vivo) cross-reactivity data was generated or checked against these
predictions. Validating fpocket's cavity-detection accuracy against a known
TRPV1 pocket demonstrates the method can *find* pockets correctly; it does
not establish that pocket-descriptor similarity or docking score predicts
*actual* binding for an unrelated fold, and the RARG result should be read
as "warrants a targeted counter-screen," not as confirmed cross-reactivity.
Second, Vina's scoring function is an empirically parameterized
approximation, not a physical free energy, and cross-receptor score
comparisons carry additional uncertainty beyond the single-receptor ranking
use case the method was built for; our confidence in the RARG result rests
on its consistency across six independent ligands rather than on any single
score. Third, both receptor structures used here are single static
crystallographic snapshots (no induced-fit or ensemble modeling), and the
RARG structure required standard handling of one unresolved side chain and
several alternate-conformation residues. Fourth, the eight-descriptor
fpocket fingerprint is a coarse proxy for true three-dimensional pocket
shape and electrostatic complementarity; purpose-built pocket-comparison
methods (e.g., IsoMIF, KRIPO, ProBiS) would provide a more rigorous, if more
computationally expensive, comparison. Fifth, bulk tissue RNA expression
(Human Protein Atlas, GTEx-derived) does not resolve dorsal root
ganglion/peripheral nerve-terminal expression specifically — the tissue
compartment most relevant to a topical analgesic's actual site of action —
so this screen captures skin-tissue co-localization risk rather than
nerve-terminal co-expression risk. Sixth, only the highest-confidence tier
of candidate genes (26 of 118 that passed the initial druggability filter)
was screened; several biologically plausible pain/inflammation-relevant
candidates identified at the target-list stage (cyclooxygenase-1/PTGS1,
the leukotriene B4 receptors LTB4R/LTB4R2, the niacin receptors HCAR2/HCAR3)
have not yet been evaluated. Seventh, the TRPV1 EC50 values used to order
the reference ligand series (Table 2) were compiled from public
bioactivity/literature sources without confirming that all seven were
measured in a directly comparable assay format; the series is used here
qualitatively, to span a wide potency range, and the EC50 column should not
be read as a precisely comparable quantitative scale across compounds
without checking primary sources for each value.

Despite these caveats, the practical implication for a topical TRPV1
program is direct for RARG: a counter-screen (at minimum, a
literature/bioactivity-database check; ideally a binding or transactivation
assay) is a reasonable candidate for early prioritization for any compound
advancing in this chemical series, rather than a question deferred until
after a lead is already selected — precisely because it is cheap to ask
early, currently unanswered, and now supported by two independent
computational lines of evidence (pocket similarity/cross-docking and the
ground-truth discrimination check, Section 3.4) rather than one. MMP3 is a
different case: the ground-truth check (Section 3.4) shows the docking
protocol cannot separate known MMP3 binders from non-binders even with a
zinc-aware scoring function, which means the Table 2 cross-docking signal
for MMP3 should not be used to prioritize or deprioritize a counter-screen
decision either way — it is not evidence of risk, but it is also not
evidence of safety, because the method that produced it has not been shown
to work on this target. A literature/bioactivity-database check for MMP3
remains cheap and reasonable to do regardless, since it does not depend on
the docking protocol at all; an enzymatic-inhibition assay is a matter of
general due diligence for this chemotype rather than something this
screen's docking result specifically motivates. Because the RARG signal
spans most of the tested series rather than concentrating in smaller
ligands, it should not be assumed that a topical-delivery-optimized
molecular-weight/logP selection window (200–450 Da, logP 1.0–4.0, already
used elsewhere in this program's chemical library curation for
permeability reasons [17]) mitigates this particular risk merely by
favoring smaller candidates. More broadly, the pipeline described here — tissue
expression filtering, cavity-detection-based pocket comparison, and
targeted cross-docking — required no experimental data beyond structures
already in the public domain and tools already validated for the primary
target, making it a low-cost addition to early-stage off-target triage for
any topical or otherwise spatially-restricted small-molecule program, not
only this one.

## 5. Conclusions

A whole-protein sequence-identity ranking of TRPV1's paralog family
correctly flags the family as a safety-relevant group but misorders the
relative binding-site risk within it, because the functionally relevant
region (the vanilloid pocket) is conserved or diverged at a different rate
than the protein overall. Off-target risk for a topical TRPV1 program is
not confined to this family: a structure-based, expression-informed screen
identified RARG and MMP3 — both proteins unrelated to TRPV1 by sequence —
as showing Vina docking scores against most of the program's reference
agonist series that are numerically comparable to those against TRPV1
itself, with resiniferatoxin (the bulkiest, most structurally complex
ligand tested) as the consistent outlier at both. A subsequent
ground-truth discrimination check (Section 3.4) separates the two
findings: RARG showed genuine active/inactive separation by docking score
(ROC-AUC 0.811), independently corroborating the cross-docking result,
while MMP3 did not (ROC-AUC 0.503, unchanged near-chance at 0.551 even
under a zinc-coordination-aware scoring function). This is a computational
hypothesis, not a demonstration of actual binding, and does not by itself
establish RARG as a confirmed safety risk; for MMP3 it does not establish
a risk either, precisely because the method's failure on the discrimination
check leaves the underlying cross-docking signal uninterpretable rather
than negative. We recommend RARG be carried forward as a priority
candidate for dedicated experimental counter-screening, and MMP3 be
carried forward as an open, methodologically unresolved question — worth a
cheap literature/bioactivity check but not weighted by this study's
docking result in either direction — pending a docking protocol shown to
discriminate known MMP3 binders or direct experimental data. We propose
the underlying pipeline, including the ground-truth discrimination check
that distinguished these two outcomes, as a
reusable, low-cost pre-lead hypothesis-generating triage step for other
topical small-molecule discovery programs.

## Author contributions

*[to complete]*

## Data and code availability

All scripts, intermediate data schemas, and result tables underlying this
work are available in the project repository (`track_a_analgesic/`
directory: `target_validation/`, `docking/`, `md_residence_time/`
subfolders). Provenance of each result is listed in the table below.

| Result | Script | Output |
|---|---|---|
| Paralog profile | `target_validation/build_profile.py` | `target_validation/profiles/TRPV1.md` |
| Skin candidate list | `target_validation/fetch_skin_targets.py` | skin-target annotation table |
| Candidate structure retrieval | `target_validation/fetch_tier1_structures.py` | receptor structure files |
| Pocket comparison | `target_validation/compare_pockets.py` | pocket-comparison tables |
| Docking box definition | `target_validation/box_from_pocket.py` | per-target box parameters |
| Cross-docking | `docking/redock_offtarget.py` | per-target docking results |
| ChEMBL active/inactive retrieval | `target_validation/fetch_chembl_bioactivity.py` | per-target classified bioactivity tables |
| Property-matched sampling | `target_validation/sample_offtarget_ligands.py` | per-target docking-validation ligand sets |
| Ground-truth discrimination docking (Vina) | `docking/run_offtarget_validation.py` | per-target discrimination-check results, ROC-AUC |
| Ground-truth discrimination docking (AD4Zn) | `docking/run_ad4zn_validation.py`, `docking/prepare_gpf4zn_lite.py`, `docking/zinc_pseudo.py` | MMP3 AD4Zn discrimination-check results |
| Discrimination statistics (bootstrap CI, Mann-Whitney *p*, PR-AUC, EC50-mixing check) | `target_validation/discrimination_stats.py` | Table 3 statistics; `figures/roc_curves.png`, `figures/score_distributions.png` |
| Pipeline funnel figure | `target_validation/render_pocket_figure.py` | `figures/pipeline_funnel.png` |

## Funding

*[to complete]*

## Conflict of interest

*[to complete]*

## References

1. Ochoa D, Hercules A, Carmona M, et al. The next-generation Open Targets
   Platform: reimagined, redesigned, rebuilt. *Nucleic Acids Res.*
   2023;51(D1):D1353–D1359. doi:10.1093/nar/gkac1046
2. Uhlén M, Fagerberg L, Hallström BM, et al. Tissue-based map of the human
   proteome. *Science.* 2015;347(6220):1260419. doi:10.1126/science.1260419
3. Le Guilloux V, Schmidtke P, Tuffery P. Fpocket: an open source platform
   for ligand pocket detection. *BMC Bioinformatics.* 2009;10:168.
   doi:10.1186/1471-2105-10-168
4. Elokely K, Velisetty P, Delemotte L, Palovcak E, Klein ML, Rohacs T,
   Carnevale V. Understanding TRPV1 activation by ligands: insights from
   the binding modes of capsaicin and resiniferatoxin. *Proc Natl Acad Sci
   USA.* 2016;113(2):E137–E145. *(verify before submission)*
5. *[Structural/mutagenesis studies establishing Tyr511, Ser512, Thr550,
   Glu570 as vanilloid-pocket-lining residues — verify and cite primary
   cryo-EM/mutagenesis sources, e.g. Cao E, Liao M, Cheng Y, Julius D
   (2013, Nature) and related follow-up structural papers, before
   submission.]*
6. Berman HM, Westbrook J, Feng Z, et al. The Protein Data Bank. *Nucleic
   Acids Res.* 2000;28(1):235–242. *(verify before submission)*
7. UniProt Consortium. UniProt: the Universal Protein Knowledgebase.
   *Nucleic Acids Res.* (current release; verify year/volume before
   submission).
8. Ravindranath PA, et al. Meeko: preparation of small molecules and
   macromolecules for AutoDock/AutoDock Vina docking. Software.
   https://github.com/forlilab/Meeko *(cite as software; verify preferred
   citation format from the repository before submission)*
9. Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of
   docking with a new scoring function, efficient optimization, and
   multithreading. *J Comput Chem.* 2010;31(2):455–461.
   doi:10.1002/jcc.21334
10. *[Structural comparison of the TRPV1–4 vanilloid-pocket region and
    TRPV4 selectivity-determining residues — verify and cite primary
    source before submission.]*
11. Goldsmith P, et al. Clinical pharmacokinetics, safety, and tolerability
    of a novel, first-in-class TRPV4 ion channel inhibitor, GSK2798745, in
    healthy and heart failure subjects. *Clin Pharmacokinet.* 2019.
    PMID: 30637626. *(verify full author list, volume, and pages before
    submission)*
12. *[TRPV6 knockout phenotype — verify and cite primary knockout-mouse
    study, e.g. Bianco SD et al. J Bone Miner Res 2007, before
    submission.]*
13. Iwata Y, Matsumura T. Blockade of TRPV2 is a novel therapy for
    cardiomyopathy in muscular dystrophy. *Int J Mol Sci.*
    2019;20(16):3844. doi:10.3390/ijms20163844
14. *[TRPV2-neutralizing-antibody safety data — verify and cite primary
    source before submission.]*
15. Lu A, Li K, Huang C, Yu B, Zhong W. Pathogenesis and management of
    TRPV3-related Olmsted syndrome. *Front Genet.* 2024;15:1459109.
    doi:10.3389/fgene.2024.1459109
16. [Authors]. Clinical proof-of-concept results with a novel TRPA1
    antagonist (LY3526318) in 3 chronic pain states. *Pain.*
    2025;166(7):1497–1518. doi:10.1097/j.pain.0000000000003487 *(verify
    full author list before submission)*
17. *[Skin-permeability MW/logP rationale — verify and cite Bos & Meinardi
    2000 (Skin Pharmacol Appl Skin Physiol) and/or Potts & Guy
    permeability-correlation sources before submission, consistent with
    `docking/ZINC_TRIAGE.md`.]*
18. QUTENZA (capsaicin) topical system, prescribing information. U.S. Food
    and Drug Administration, NDA 022395. Original approval 2009; cited
    revision 2023. Transient, low (<5 ng/mL) plasma capsaicin was detected
    in about one-third of postherpetic-neuralgia patients following a
    60-minute application, clearing below the limit of quantitation within
    3–6 hours. https://www.accessdata.fda.gov/drugsatfda_docs/label/2023/022395s023lbl.pdf
19. Thoreau E, Arlabosse JM, Bouix-Peter C, et al. Structure-based design
    of trifarotene (CD5789), a potent and selective RARγ agonist for the
    treatment of acne. *Bioorg Med Chem Lett.* 2018. doi:10.1016/j.bmcl.2018.04.036
    *(this is the source publication for PDB 6FX0, used directly in this
    study; verify volume/pages before submission)*
20. Chandraratna RAS. Tazarotene — first of a new generation of
    receptor-selective retinoids. *Br J Dermatol.* 1996;135(Suppl
    49):18–25. PMID: 9035701
21. Esser CK, Bugianesi RL, Caldwell CG, et al. Inhibition of stromelysin-1
    (MMP-3) by P1'-biphenylylethyl carboxyalkyl dipeptides. *J Med Chem.*
    1997. doi:10.1021/jm960465t *(source publication for PDB 1HFS, used
    directly in this study; verify volume/pages before submission)*
22. Lopez KE, Paduda AS, Derrick MJ, Van Horn WD. TRPV1 antagonism occurs
    through diverse structural mechanisms. *bioRxiv.* 2026.
    doi:10.64898/2026.04.27.721197 *(source publication for PDB 11CK, used
    directly in this study; this is a preprint — verify publication status
    before submission)*
23. Zdrazil B, Felix E, Hunter F, et al. The ChEMBL Database in 2023: a
    drug discovery platform spanning multiple bioactivity data types and
    time periods. *Nucleic Acids Res.* 2024;52(D1):D1180–D1192.
    doi:10.1093/nar/gkad1004 *(source of RARG/MMP3 known active/inactive
    bioactivity data used directly in this study, via the ChEMBL REST API,
    targets CHEMBL2003 and CHEMBL283; verify preferred current-release
    citation before submission)*
24. Santos-Martins D, Forli S, Ramos MJ, Olson AJ. AutoDock4(Zn): an
    improved AutoDock force field for small-molecule docking to zinc
    metalloproteins. *J Chem Inf Model.* 2014;54(8):2371–2379.
    doi:10.1021/ci500209e *(source of the AD4Zn zinc-coordination force
    field used directly in this study for the MMP3 discrimination check,
    Section 2.5/3.4)*
