# Off-target safety triage for a topical TRPV1-directed analgesic (Track A)

**Status: draft for review.** Written up from the working sessions of
2026-08-10/11; every number here traces back to a script and a data file in
this repo (paths given throughout) rather than to hand recollection.

## Abstract

Before committing further medicinal-chemistry effort to Track A's
TRPV1-targeted topical analgesic, we asked which off-target proteins pose a
plausible safety risk and why. We addressed this in two stages. First, the
classical approach: the TRPV1 sequence-paralog family (TRPV2-6) plus the
mechanistically adjacent channel TRPA1, ranked by whole-protein sequence
identity and annotated against the literature. Second, because whole-protein
identity is a poor proxy for whether an off-target protein's *ligand-binding
pocket* actually resembles TRPV1's vanilloid pocket, we built a
structure-based screen: a genome-scale skin-expression filter to find
candidate off-targets outside the TRPV/TRPA family, a cavity-detection
pipeline (fpocket) to compare pocket geometry/chemistry against TRPV1's
vanilloid pocket without running any docking, and -- for the two top hits --
actual cross-docking of our seven validation ligands to check whether the
structural resemblance translates into comparable predicted binding. The
headline result: **RARG (retinoic acid receptor gamma)**, an existing
approved topical drug target unrelated to TRPV1 by sequence or fold, binds
six of our seven candidate-class ligands as well as or better than TRPV1
itself does, by a docking scoring function applied blind to an unrelated
protein. This is a candidate-class liability worth carrying forward as an
explicit counter-screen, not a footnote.

## 1. Motivation

Track A's mechanism (TRPV1 agonist-driven desensitization, Qutenza-style, via
a topical/local route) makes the standard analgesic off-target panel
incomplete on its own. Two separate risks need separate methods:

1. **Family-level risk**: does the drug hit a *sequence relative* of TRPV1
   (TRPV2-6) via the same or a homologous binding site, given they're
   evolutionarily and structurally related?
2. **Tissue-coincidence risk**: does the drug hit an *unrelated* protein that
   happens to (a) share a similarly-shaped binding pocket and (b) be
   co-expressed at the same site of action (skin, for a topical product),
   regardless of sequence relationship?

Standard practice covers (1) reasonably well (paralog panels are routine).
(2) is normally addressed empirically, late, and expensively (a broad in
vitro panel after a lead exists). We asked whether a structural pre-screen
could narrow that panel *before* lead optimization, using tools already in
the project's docking stack.

## 2. Family-level analysis: TRPV1 sequence paralogs

### 2.1 Method

Paralog identity and curated safety liabilities pulled from the Open Targets
Platform GraphQL API (`target_validation/build_profile.py`, reusing the same
pipeline built for the six-target comparison in
[comparison_all_candidates.md](profiles/comparison_all_candidates.md)).
Supplemented with a targeted literature pass per paralog (2026-08-10 session)
specifically to fill gaps the API doesn't curate (TRPV5/TRPV6 had no
liability text in the original pull).

### 2.2 Results

| Paralog | Sequence identity to TRPV1 | Liability if hit off-target |
|---|---|---|
| TRPV2 | 45.0% | Developmental cardiac Ca2+-handling role (cardiomyopathy in TRPV2-deficient growing hearts); adult whole-body knockout and neutralizing-antibody dosing show no overt phenotype, so acute pharmacological inhibition looks lower-risk than the developmental data alone suggests. |
| TRPV4 | 41.9% | Direction-dependent. TRPV4 *activation* is implicated in cardiogenic pulmonary edema, but the clinical TRPV4 *antagonist* GSK2798745 was well-tolerated through repeat dosing in heart-failure patients with no cardiopulmonary safety signal (Cheung et al., 2019, PubMed:30637626) -- an off-target TRPV1-antagonist hit on TRPV4 trends protective, not causative, for this specific liability. |
| TRPV3 | 39.0% | Shares TRPV1's topical site of action (keratinocytes/skin barrier). Human gain-of-function mutations cause Olmsted syndrome (painful palmoplantar keratoderma, pruritus, skin-barrier failure) -- direct proof this channel matters clinically at the skin, though that is a GoF disease phenotype and doesn't necessarily predict the consequence of antagonism. |
| TRPV5 | 28.4% | Rate-limiting channel for renal distal-tubule Ca2+ reabsorption; off-target inhibition is a plausible hypercalciuria risk. Not curated in Open Targets; added from the primary literature. |
| TRPV6 | 27.2% | Most systemically consequential paralog if hit: knockout mice show impaired intestinal Ca2+ absorption, reduced bone mineral density (osteoporosis), alopecia/dermatitis, and severely impaired male fertility. Not curated in Open Targets; added from the primary literature. |
| TRPA1 | n/a (not a sequence paralog; mechanistically adjacent, co-expressed nociceptor channel) | Two independent clinical antagonist programs (Lilly's LY3526318; a Roche/Genentech compound) had an acceptable Phase 1 safety profile but both failed on efficacy (knee OA, chronic low back pain, diabetic peripheral neuropathy for LY3526318; chronic cough for the Roche/Genentech compound) -- an off-target TRPA1 hit looks safety-neutral, but "clinically unproven for this indication," not "clinically validated as safe to hit." |

Full detail and sourcing: [profiles/TRPV1.md](profiles/TRPV1.md) (see "Safety
-- off-target family" table and the appended literature caveat).

### 2.3 A methodological correction made mid-analysis

The initial paralog ranking (by % whole-protein sequence identity) implicitly
treats identity as a proxy for *binding-site* risk. It isn't one. TRPV1's
vanilloid pocket is formed by specific residues (Y511, S512, T550, E570,
identified from capsaicin/RTX/capsazepine-bound cryo-EM structures) in the
S3/S4/S4-S5-linker region of one subunit plus S5/S6 of the neighboring
subunit. The literature indicates this specific region is **structurally and
sequence-wise highly conserved between TRPV1, TRPV2, and TRPV3** -- the
TM3/4 region is experimentally transplantable between orthologs and confers
vanilloid sensitivity onto otherwise-insensitive channels. TRPV4's
selectivity against vanilloids, despite similar pocket architecture, is
attributable to specific residue substitutions *within* that conserved
scaffold. This means TRPV2 and TRPV3's real off-target pocket risk may be
*higher* than their whole-protein identity (45%, 39%) suggests, while TRPV4's
may be somewhat *lower* than 41.9% implies, because nature has already
engineered selectivity into an otherwise-similar pocket. TRPV5/TRPV6, by
contrast, are Ca2+-selective channels with a different pore/gating
architecture altogether -- their low whole-protein identity plausibly *does*
track real pocket-level divergence, unlike TRPV2-4.

This observation is what motivated Section 3: if whole-protein identity is
an unreliable proxy for pocket risk even *within* the TRPV family, it says
nothing at all about off-target risk *outside* the family -- which requires
comparing actual pocket structure, not gene relatedness.

## 3. Structure-based screen for off-target proteins outside the TRPV/TRPA family

### 3.1 Candidate list: skin-expressed, druggable, non-TRPV/TRPA

**Source list.** Human Protein Atlas RNA tissue-specificity API, filtered to
genes classified `Tissue enriched` / `Group enriched` / `Tissue enhanced` for
skin: 604 genes (2026-08-10 pull, `data/raw/hpa/skin_tissue_specific_genes.tsv`).

**Druggability filter.** Each gene cross-referenced against Open Targets
`targetClass` and small-molecule `tractability` (same fields as Section 2),
via a batched GraphQL query (`target_validation/fetch_skin_targets.py`). 486
of 604 genes returned no Open Targets target class at all -- these are
overwhelmingly structural/barrier proteins (keratins, small proline-rich
proteins, late cornified envelope genes, filaggrin-family genes) with no
pharmacological pocket, and were dropped regardless of how strongly
skin-specific their expression is.

**Tiering** of the remaining 118 genes by strength of pocket evidence
(`data/processed/skin_offtarget_screen/skin_targets_annotated.csv`; full
writeup in
[skin_offtarget_targets.md](skin_offtarget_targets.md)):

- **Tier 1** (n=26): Open Targets flags a confirmed ligand-binding pocket
  (High/Med-Quality Pocket in its structure-based druggability prediction).
- **Tier 2** (n=59, TRPV3 excluded as already covered): a real
  ligand-bound structure or high-quality chemical probe exists, but the
  automated pocket classifier didn't tag it explicitly.
- **Tier 3** (n=22): classified into a druggable-looking family
  (GPCR/kinase/protease/transporter) by sequence, but no structural evidence
  yet.

This writeup covers the follow-up done on **Tier 1**; Tiers 2-3 (which
include some biologically obvious candidates for a topical/inflammatory
context -- PTGS1/COX-1, the leukotriene B4 receptors LTB4R/LTB4R2,
HCAR2/HCAR3, several kallikreins) are flagged but not yet screened.

### 3.2 Pocket-vs-pocket structural comparison

**Tool**: fpocket, a geometry-only cavity detector (Voronoi
tessellation/alpha-spheres) -- no ligand is posed, so this is much cheaper
than docking and scales to dozens of targets in minutes.

**Method validation.** Run first on TRPV1 itself (`receptor_AD_raw.pdb`, the
inter-subunit chain A+D fragment from PDB 11CK already validated in the MD
residence-time pilot as containing the true vanilloid pocket -- see
`md_residence_time/STATUS.md`). fpocket's own top-ranked cavity (by its
internal scoring, not by druggability alone) independently recovered all
four literature residues (Y511, S512, T550, E570), with a druggability score
of 0.783 and volume 758 A^3. This validates using "fpocket's own top-ranked
pocket" as the comparison rule for the 26 candidates -- it isn't a free
parameter, it was checked against ground truth first.

**Structure retrieval.** For each of the 26 Tier 1 genes: UniProt lookup for
the best-resolution, largest-span experimental PDB structure (falling back to
AlphaFold if none existed; in practice all 26 had an experimental structure).
One bug caught and fixed during this step: UniProt's free-text `gene:X`
query is not an exact-symbol filter -- querying "FBXW7" silently returned an
unrelated antisense-derived receptor gene (FBXW7-AS1/DEspR, accession
B0L3A2) instead of FBXW7 itself (Q969H0). Fixed by verifying the returned
primary gene symbol matches the query exactly before accepting a hit
(`target_validation/fetch_tier1_structures.py`); confirmed the other 25
genes were unaffected.

**Fingerprint comparison.** An 8-descriptor vector per detected pocket
(druggability score, volume, hydrophobicity score, polarity score, charge
score, proportion of polar atoms, apolar alpha-sphere proportion, mean local
hydrophobic density), z-score normalized across every pocket found on every
structure, compared to the TRPV1 reference by Euclidean distance and cosine
similarity (`target_validation/compare_pockets.py`).

Two rankings were computed and cross-checked against each other rather than
trusting either alone: the best-matching pocket found *anywhere* on a given
structure (which favors proteins fpocket found many candidate cavities on --
a multiple-comparisons problem), and the single pocket fpocket itself ranks
highest (the same selection rule validated against TRPV1's known pocket,
so a fairer one-shot comparison). Full data and both rankings:
[skin_offtarget_pocket_comparison.md](skin_offtarget_pocket_comparison.md).

### 3.3 Results: which targets to flag

Three genes land near the top of *both* independent rankings: **MMP3, NTF4,
BLMH**. Each flagged hit was then checked for biological plausibility --
does the detected "pocket" correspond to a real, known small-molecule
binding site, or is it a probable detection artifact (e.g. a
protein-protein interaction surface, a serpin's conformational groove)?

- **MMP3 (stromelysin-1)**: real, deep, hydrophobic zinc-dependent catalytic
  pocket, high druggability (0.839). Mechanistically plausible cross-reactor
  for a hydrophobic, elongated vanilloid-shaped ligand. Also independently
  relevant to skin biology (extracellular-matrix remodeling, wound healing).
  **Flagged, high confidence.**
- **RARG**: did not lead either individual ranking by distance (6th in the
  fairer one-shot ranking), but has near-maximal druggability (0.982) and
  high cosine similarity (0.899) -- its rank was suppressed mainly by
  having a larger pocket volume than TRPV1's (1032 vs 758 A^3), not by
  chemical dissimilarity. Promoted to a flag on biological-plausibility
  grounds specifically: **RARG is already an approved topical skin drug
  target** (tazarotene, for psoriasis/acne) -- an off-target hit here
  competes with an existing therapeutic mechanism in the same tissue, not a
  hypothetical one. **Flagged, high confidence** (and subsequently the
  strongest-confirmed hit -- see Section 4).
- **BLMH (bleomycin hydrolase)**: moderate druggability (0.404), but
  clinically notable because *low* BLMH activity in skin is already
  implicated in a real topically-relevant drug toxicity (bleomycin-induced
  flagellate dermatitis). **Flagged, moderate confidence.**
- **NTF4, CHP2, SERPINB5**: down-weighted. NTF4 is a secreted growth factor
  (TrkB ligand) -- its detected "pocket" is more likely a protein-protein
  interaction surface than a small-molecule site. SERPINB5's near-maximal
  druggability score (0.995) is suspicious for a serpin, whose function is
  mediated by reactive-center-loop insertion rather than small-molecule
  pocket occupancy -- probable detection artifact. CHP2 has no strong prior
  either way and wasn't pursued further.

Caveats specific to this stage (multiple-comparisons risk, a likely
fusion-protein artifact discounting the ADRB2 result, the coarseness of an
8-descriptor fingerprint vs. a true 3D shape/electrostatics comparison) are
detailed in [skin_offtarget_pocket_comparison.md](skin_offtarget_pocket_comparison.md#caveats)
and carry over to this document's limitations (Section 5).

## 4. Cross-docking validation: does the structural resemblance predict actual binding?

Pocket-descriptor similarity is, at best, a hypothesis generator. The
decisive check is whether real candidate-class ligands actually dock
favorably. We redocked the seven TRPV1 validation ligands already used to
calibrate the MD residence-time pilot (`md_residence_time/ligands.py`:
resiniferatoxin, phenylacetylrinvanil, olvanil, arvanil, rinvanil, capsaicin,
nonivamide -- an agonist series spanning ~5 orders of magnitude in TRPV1
EC50) into MMP3 and RARG, using the same Vina/Meeko pipeline already
validated for TRPV1 in `docking/SMOKE_TEST.md`. Receptors: PDB 1SLM (MMP3)
and PDB 6FX0 (RARG), prepared with `mk_prepare_receptor.py`; docking box
centered on the fpocket-flagged pocket's alpha-sphere centroid (+6 A
padding). Run on two vast.ai instances in parallel, one target per instance.

| Ligand | TRPV1 EC50 (nM) | Vina score, TRPV1 | Vina score, MMP3 | Vina score, RARG |
|---|---|---|---|---|
| resiniferatoxin | 0.011 | -9.63 | -5.18 | -2.86 |
| phenylacetylrinvanil | 0.09 | -8.01 | -6.70 | **-8.98** |
| olvanil | 0.6 | -7.49 | -6.98 | **-8.73** |
| arvanil | 0.5 | -7.83 | -7.51 | **-9.17** |
| rinvanil | 6.0 | -7.51 | -7.28 | **-8.60** |
| capsaicin | 712.0 | -8.23 | **-8.43** | **-8.59** |
| nonivamide | 1400.0 | -7.60 | -7.70 | -7.97 |

(Vina score in kcal/mol; more negative = stronger predicted binding.
Bold = matches or exceeds the ligand's own TRPV1 score.)

**RARG: six of seven ligands dock as well as or better than at TRPV1
itself.** Only resiniferatoxin -- the bulkiest, most structurally complex
ligand in the set, a diterpene -- fits poorly. Applied blind to an unrelated
protein fold, Vina's scoring function ranks most of this chemical series as
comparably or more complementary to RARG's retinoid pocket than to TRPV1's
own vanilloid pocket. This is an independent confirmation of the
pocket-similarity flag, not a restatement of it.

**MMP3: size-dependent.** Capsaicin and nonivamide -- the smallest,
simplest molecules in the set, and the closest analogs to what Track A's
ZINC diversity-triage MW/logP window (200-450 Da, logP 1.0-4.0, anchored
explicitly on capsaicin -- `docking/ZINC_TRIAGE.md`) would actually select
-- match or slightly exceed their TRPV1 score at MMP3. The bulkier,
longer-chain ligands score consistently worse at MMP3 than at TRPV1.

**Implication for the ZINC triage window.** The off-target risk at both
targets, and especially RARG, is concentrated in exactly the small,
capsaicin-like chemical space the topical MW/logP filter already selects
for. This is not a risk that a "prefer smaller, more permeable molecules"
heuristic screens out on its own -- if anything, smaller molecules are the
ones binding RARG and MMP3 comparably well.

## 5. Limitations

- **Hypothesis-generating, not confirmatory.** No wet-lab or literature
  cross-reactivity data was checked against these predictions. The TRPV1
  self-validation (Section 3.2) proves fpocket recovers the *known* pocket;
  it does not prove descriptor similarity or Vina score predicts *actual*
  binding for an unrelated fold. RARG's result should be read as "worth an
  explicit counter-screen," not as "confirmed cross-reactivity."
- **Vina score is an approximate, empirically-fit scoring function**, not a
  free energy. Comparing absolute scores across different receptors/pocket
  sizes has its own caveats beyond the usual per-target ranking use case Vina
  is built for; the direction and magnitude of the RARG result (consistent
  across 6/7 ligands, not a single outlier) is what makes it credible, not
  any individual score in isolation.
- **Single static receptor conformation per target**, no induced-fit or
  ensemble docking. RARG's structure (6FX0) required
  `--allow_bad_res --default_altloc A` to handle one unresolved Tyr180
  side chain and seven altloc residues -- a normal crystallographic
  wrinkle, but a reminder these are single experimental snapshots.
- **8-descriptor fpocket fingerprint is coarse**: no 3D pocket-shape
  alignment, no explicit electrostatics, no knowledge of which residues
  would sterically clash with a specific ligand. A purpose-built pocket
  comparison tool (IsoMIF, KRIPO, ProBiS) would be more rigorous; fpocket
  was chosen specifically because it needed no new tooling and runs in
  minutes across dozens of targets.
- **HPA/GTEx bulk tissue expression doesn't resolve dorsal root
  ganglion/peripheral nerve-terminal expression** -- the tissue most
  relevant to a topical analgesic's actual site of action. A target could
  matter at the nerve terminal without being "skin-enriched" in bulk RNA,
  or vice versa. This caveat also applies to the original TRPV1 profile
  (Section 2) and is inherited here, not new to this stage.
- **GPCR structures may be compromised by crystallization fusion
  constructs.** ADRB2's top fpocket-ranked pocket had implausibly low
  druggability (0.017) for a validated orthosteric GPCR site, suggesting
  fpocket's #1 pick was a fusion-protein (T4-lysozyme/BRIL) artifact rather
  than the real pocket -- ADRB2's negative result should be discounted, not
  read as genuine dissimilarity.
- **Only Tier 1 (26/118 candidate genes) has been screened.** Tier 2's
  biologically-obvious pain/inflammation-adjacent candidates (PTGS1,
  LTB4R/LTB4R2, HCAR2/3) have not yet been run through either the pocket
  comparison or redocking.

## 6. Recommendations / next steps

1. **RARG counter-screen becomes a standard step for any Track A
   hit-to-lead candidate** in this chemical series, not an optional
   follow-up -- at minimum a ChEMBL/literature bioactivity check, ideally a
   real RARG binding or transactivation assay before advancing a compound.
2. **MMP3 stays on the watchlist**, weighted toward smaller candidates
   specifically (the capsaicin/nonivamide-sized end of the series).
3. **BLMH stays a lower-priority watch item**; not redocked in this pass.
4. Extend the same pocket-comparison + redocking pipeline to **Tier 2's
   pain/inflammation-relevant candidates** (PTGS1, LTB4R/LTB4R2, HCAR2/3,
   GPRC5D) if screening bandwidth allows -- these were flagged as
   biologically notable at the target-list stage (Section 3.1) but never
   structurally screened.
5. If RARG cross-reactivity is corroborated further, consider whether it
   changes the chemotype-selection strategy for the ZINC triage pool
   (Section 4's implication) rather than treating it as a downstream
   counter-screen only.

## Provenance / how to reproduce

| Stage | Script | Output |
|---|---|---|
| TRPV1 paralog profile | `target_validation/build_profile.py` | `target_validation/profiles/TRPV1.md` |
| Skin candidate list | `target_validation/fetch_skin_targets.py` | `data/processed/skin_offtarget_screen/skin_targets_annotated.csv` |
| Tier 1 structure retrieval | `target_validation/fetch_tier1_structures.py` | `data/raw/structures/tier1_skin_targets/` |
| Pocket comparison | `target_validation/compare_pockets.py` | `data/processed/skin_offtarget_screen/pocket_comparison*.csv` |
| Docking box from pocket | `target_validation/box_from_pocket.py` | (box.txt per target) |
| Cross-docking | `docking/redock_offtarget.py` | `data/processed/skin_offtarget_screen/redocking/{MMP3,RARG}/` |

All raw/processed data is gitignored per repo convention (`data/raw/*`,
`data/processed/*`) -- rerun the scripts above to regenerate.
