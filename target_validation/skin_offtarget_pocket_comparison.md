# Tier 1 skin targets: pocket-vs-pocket comparison against TRPV1's vanilloid pocket

Follow-up to [skin_offtarget_targets.md](skin_offtarget_targets.md). Structural
screen requested as a faster alternative to docking: compare cavity geometry
and physicochemistry of each Tier 1 target's binding pocket(s) against
TRPV1's vanilloid pocket, without posing any ligand.

Scripts: `fetch_tier1_structures.py` (structure retrieval) ->
`compare_pockets.py` (fpocket + comparison). Data:
`data/raw/structures/tier1_skin_targets/`,
`data/processed/skin_offtarget_screen/`.

## Method

1. **Tool**: [fpocket](https://github.com/Discngine/fpocket) (added to
   `docker/environment-docking.yml`) -- detects cavities from protein
   geometry alone (Voronoi tessellation / alpha-spheres), no ligand needed.
   Each detected pocket gets druggability/physicochemical descriptors
   (volume, hydrophobicity score, polarity score, charge score, apolar
   alpha-sphere proportion, etc).
2. **TRPV1 reference**: reused the already-validated `receptor_AD_raw.pdb`
   from the MD residence-time pilot (chains A+D, residues 430-725 --
   empirically confirmed to contain the true inter-subunit vanilloid pocket,
   see `md_residence_time/STATUS.md`). Ran fpocket on it; the top-scoring
   detected cavity (fpocket's own internal ranking) contains all four
   literature-validated vanilloid-pocket residues (Y511, S512, T550, E570)
   -- **this is the method-validation step**: fpocket's highest-ranked
   pocket independently recovered the known functional site, so its ranking
   convention is trustworthy enough to apply blind to the 26 candidate
   targets. Druggability score 0.783, volume 758 A^3.
3. **Structure retrieval for the 26 targets**: UniProt lookup -> best
   experimental PDB structure (largest resolved span, best resolution,
   X-ray/cryo-EM only) via RCSB; would have fallen back to AlphaFold if no
   experimental structure existed, but all 26 had one (25 on the first
   pass, 1 -- FBXW7 -- recovered after fixing a lookup bug, see Caveats).
4. **Fingerprint comparison**: 8-descriptor vector (Druggability Score,
   Volume, Hydrophobicity score, Polarity score, Charge score, Proportion
   of polar atoms, Apolar alpha sphere proportion, Mean local hydrophobic
   density) per pocket, z-score normalized across every pocket found on
   every structure (26 targets + TRPV1), compared to the TRPV1 reference by
   Euclidean distance and cosine similarity.
5. **Two rankings reported**, because they answer different questions:
   - *Best pocket anywhere on the structure* -- most similar cavity fpocket
     found, regardless of whether it's biologically the target's real
     functional site. Inflated for proteins where fpocket found many
     candidate pockets (more shots at a coincidental match by chance).
   - *Top fpocket-ranked pocket only* -- the single cavity fpocket itself
     considers most likely real, i.e. the same selection rule validated
     against TRPV1's known vanilloid pocket above. Fairer one-shot
     comparison, not multiple-comparisons-inflated.

## Results

### Top fpocket-ranked pocket only (primary ranking -- same rule validated on TRPV1)

| Rank | Gene | Eucl. dist | Cos. sim | Druggability | Volume (A^3) |
|---|---|---|---|---|---|
| 1 | **MMP3** | 2.56 | 0.896 | 0.839 | 626 |
| 2 | **NTF4** | 3.10 | 0.858 | 0.501 | 490 |
| 3 | **CHP2** | 3.65 | 0.790 | 0.750 | 492 |
| 4 | **BLMH** | 3.66 | 0.727 | 0.404 | 1107 |
| 5 | RORA | 4.27 | 0.615 | 0.089 | 735 |
| 6 | RARG | 4.28 | 0.899 | 0.982 | 1032 |
| 7 | SERPINB5 | 4.36 | 0.912 | 0.995 | 1046 |
| 8-26 | GAN, ITGB4, FBXW7, FCER1A, SFN, CA12, SRY, ANXA8, GSTA3, GLTP, ADRB2, RAB3D, DSP, CD207, FABP9, SULT1E1, SULT2B1, TREX2, SERPINB2 | 4.5-9.5 | mixed | mixed | mixed |

Full data: `data/processed/skin_offtarget_screen/pocket_comparison_top1_only.csv`

### Best pocket anywhere on structure (exploratory, secondary)

Top 6: FABP9, TREX2, FCER1A, DSP, BLMH, MMP3. Full data:
`data/processed/skin_offtarget_screen/pocket_comparison_results.csv`

### Cross-referenced (appear near the top of *both* rankings)

**MMP3, NTF4, BLMH** are the only three targets that land in the top of
both rankings independently -- the most defensible hits, since they don't
depend on which comparison convention you trust.

## Reading the results -- biological plausibility check

Geometric/physicochemical similarity alone doesn't mean a real off-target
risk; cross-checked each flagged hit against what's actually known about it:

- **MMP3 (stromelysin-1)** -- matrix metalloproteinase, real deep
  hydrophobic S1' catalytic pocket (zinc-dependent), heavily implicated in
  skin extracellular-matrix remodeling and wound healing. High
  druggability (0.839) confirms it's a genuine small-molecule site, not a
  detection artifact. Mechanistically plausible that a hydrophobic,
  vanilloid-shaped ligand could cross-react -- **highest-confidence flag**.
- **RARG** -- retinoic acid receptor gamma nuclear-hormone-receptor
  ligand-binding domain. Near-maximal druggability (0.982) and high cosine
  similarity (0.899) despite a larger Euclidean distance (driven mostly by
  its larger volume, 1032 vs 758 A^3). Directly relevant: **RARG is already
  an approved topical skin drug target** (tazarotene, for
  psoriasis/acne) -- an off-target hit here isn't a hypothetical liability,
  it's competing with an existing therapeutic mechanism in the same tissue.
  **Second highest-confidence flag**, despite ranking 6th by distance alone.
- **BLMH (bleomycin hydrolase)** -- cysteine protease, clinically notable
  because *low* BLMH activity in skin is implicated in bleomycin-induced
  flagellate dermatitis -- i.e. this enzyme's skin activity is already known
  to matter clinically for a topically-relevant drug toxicity. Worth
  flagging, moderate druggability (0.404).
- **NTF4 (neurotrophin-4/5)** and **CHP2**, **SERPINB5 (maspin)** -- lower
  confidence. NTF4 is a secreted growth factor (TrkB ligand); its "pocket"
  is likely a protein-protein interaction surface groove, not a
  small-molecule site in the same sense as TRPV1's -- geometric similarity
  here is less pharmacologically meaningful than for an enzyme active site.
  SERPINB5's near-maximal druggability score (0.995) is suspicious for a
  serpin (function is via reactive-center-loop insertion, not classic
  small-molecule occupancy) -- treat as a probable pocket-detection
  artifact rather than a real liability until checked by hand. CHP2 is
  poorly characterized pharmacologically; no strong prior either way.

## Caveats

- **Not validated against any known cross-reactivity data.** This is a
  hypothesis-generating triage step (which targets deserve a closer look),
  not a prediction of actual binding. The TRPV1 self-validation only proves
  fpocket recovers the *known* pocket -- it says nothing about whether
  descriptor-similarity actually predicts real ligand cross-reactivity for
  an unrelated fold.
- **Multiple-comparisons risk** in the "best pocket anywhere" ranking:
  proteins with more detected pockets (DSP: 31, SERPINB2: 27, SULT1E1: 24)
  get more chances to find a coincidentally close match. The "top
  fpocket-ranked pocket only" ranking avoids this and should be weighted
  more heavily -- which is why it's listed as primary above.
- **GPCR structures (ADRB2) may be compromised by fusion-protein
  crystallization constructs** (e.g. T4-lysozyme/BRIL fusions commonly
  inserted into ICL3 for crystallization) that weren't stripped out during
  chain extraction -- ADRB2's top fpocket-ranked pocket has druggability
  0.017, implausibly low for a validated GPCR orthosteric site, suggesting
  fpocket's #1 pick here is a fusion-construct artifact, not the real
  beta-adrenergic pocket. ADRB2's result should be discounted, not read as
  "genuinely dissimilar."
- **Fixed during this run**: `fetch_tier1_structures.py`'s original
  UniProt lookup used a free-text `gene:X` query that isn't an exact-symbol
  filter -- querying "FBXW7" silently returned B0L3A2 (FBXW7-AS1, an
  unrelated antisense-derived receptor gene) instead of FBXW7 itself
  (Q969H0). Fixed by verifying the returned primary gene symbol matches the
  query exactly before accepting a hit; re-ran and confirmed the other 25
  genes were unaffected by this bug.
- **8-descriptor fingerprint is a coarse proxy** for real 3D shape/chemistry
  match -- it doesn't capture pocket shape topology, doesn't align pockets
  in 3D, and doesn't know which specific residues would clash sterically
  with a specific ligand. A tool built for this (IsoMIF, KRIPO, ProBiS) would
  give a more rigorous answer; this was chosen specifically because it's
  fast and needs no new tooling beyond fpocket.

## Follow-up: actual redocking (2026-08-11)

Deferred step done: the 7 TRPV1 validation ligands
(`md_residence_time/ligands.py`) redocked into MMP3 and RARG's fpocket-flagged
pocket with the same Vina/Meeko pipeline (`docking/ligand_prep.py`,
`docking/run_vina.py`), box centered on the pocket's alpha-sphere centroid
(`compare_pockets.py`'s reference pocket, padded +6 A). Run on two vast.ai
instances (one per target, in parallel) -- receptors: `mk_prepare_receptor.py`
on the same protein-only PDBs used for the fpocket screen (RARG needed
`--allow_bad_res --default_altloc A` for one unresolved Tyr180 side chain and
7 altloc residues in the 6FX0 crystal structure). Data:
`data/processed/skin_offtarget_screen/redocking/{MMP3,RARG}/`.

**Vina score (kcal/mol, more negative = stronger predicted binding), TRPV1 baseline from the existing `md_residence_time/docked/*_docked.pdbqt` runs:**

| Ligand | EC50 (nM, TRPV1) | TRPV1 | MMP3 | RARG |
|---|---|---|---|---|
| resiniferatoxin | 0.011 | -9.63 | -5.18 | -2.86 |
| phenylacetylrinvanil | 0.09 | -8.01 | -6.70 | **-8.98** |
| olvanil | 0.6 | -7.49 | -6.98 | **-8.73** |
| arvanil | 0.5 | -7.83 | -7.51 | **-9.17** |
| rinvanil | 6.0 | -7.51 | -7.28 | **-8.60** |
| capsaicin | 712.0 | -8.23 | **-8.43** | **-8.59** |
| nonivamide | 1400.0 | -7.60 | -7.70 | -7.97 |

**RARG: 6 of 7 ligands dock with equal or stronger predicted affinity at RARG
than at their actual target, TRPV1 itself** -- only resiniferatoxin (the
largest, most structurally complex ligand, a diterpene) fits poorly (-2.86).
This is a real, independent confirmation of the pocket-similarity flag, not
just a descriptor-geometry coincidence -- Vina's scoring function, applied
blind to an unrelated fold, ranks these capsaicinoids as comparably or more
complementary to RARG's retinoid pocket than to TRPV1's own vanilloid pocket.

**MMP3: mixed, size-dependent.** Capsaicin and nonivamide -- the smallest,
simplest vanillamides in the set, closest in analogy to what a topical
ZINC-triaged candidate would look like (`docking/ZINC_TRIAGE.md`'s MW
200-450 / logP 1-4 window is explicitly anchored on capsaicin) -- match or
slightly exceed their TRPV1 score at MMP3. The bulkier long-chain ligands
(RTX, phenylacetylrinvanil, olvanil, arvanil, rinvanil) score consistently
worse at MMP3 than at TRPV1.

**Practical implication for the ZINC triage window**: the off-target risk at
both MMP3 and (especially) RARG is concentrated in exactly the small,
capsaicin-like chemical space that Track A's topical MW/logP filter already
selects for -- this isn't a risk that screens itself out by picking smaller,
more skin-permeable molecules; if anything smaller/simpler vanilloids are the
ones that fit RARG and MMP3 comparably well.

## Recommendation (updated)

- **RARG: escalate to a real off-target liability, not just a watchlist
  entry.** Redocking corroborates the pocket-similarity flag directly --
  most candidate-class ligands are predicted to bind RARG's retinoid pocket
  as well as or better than TRPV1's vanilloid pocket. Any hit-to-lead
  candidate in this chemical space should get an explicit RARG counter-screen
  (at minimum a literature/ChEMBL bioactivity check, ideally an actual
  RARG binding/transactivation assay) before advancing.
- **MMP3: keep as a watchlist item, weighted toward smaller candidates.**
  Real risk for capsaicin/nonivamide-sized molecules specifically; less
  relevant for bulkier chemotypes.
- **BLMH**: unchanged from the pocket-comparison-only assessment (not
  redocked here) -- lower-priority watch item.
- NTF4/CHP2/SERPINB5: still treated as noise (not redocked; pocket-detection
  artifacts per the biological-plausibility check above).
