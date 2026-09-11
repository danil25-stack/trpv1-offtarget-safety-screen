# TRPV1 off-target safety screen

In silico off-target safety triage for a topical TRPV1-directed analgesic,
extracted from a larger multi-target drug discovery project. Full writeup
and manuscript draft: [`target_validation/OFFTARGET_SAFETY_ARTICLE_EN.md`](target_validation/OFFTARGET_SAFETY_ARTICLE_EN.md)
(Ukrainian translation: `OFFTARGET_SAFETY_ARTICLE_UA.md`; `.docx` versions
of both are provided alongside).

## What this is

Standard off-target safety review for a small molecule normally stops at
the target's sequence-paralog family. This asks two questions that a
paralog-only review misses:

1. **Within the TRPV1 paralog family (TRPV2–6, TRPA1)** — does
   whole-protein sequence identity actually predict binding-site risk, or
   can it be misleading?
2. **Outside the family** — can a genome-scale tissue-expression filter +
   structure-based pocket comparison surface an *unrelated* protein that
   happens to share a similar binding pocket and the same site of action
   (skin, for a topical drug), before a lead compound even exists?

Pipeline: HPA skin-expression filter → Open Targets druggability tiering →
fpocket pocket-fingerprint comparison against TRPV1's vanilloid pocket →
cross-docking (AutoDock Vina / AD4Zn) of a 7-ligand TRPV1 reference series
into the top hits → a ChEMBL-derived ground-truth discrimination check
(property-matched active/inactive panels, ROC-AUC) to test whether the
docking scoring function can actually separate real binders from
non-binders at each off-target.

**Headline result:** RARG (retinoic acid receptor gamma — the target of an
approved topical retinoid) binds most of the reference series comparably
to TRPV1 itself and passes the discrimination check (ROC-AUC 0.811); MMP3
shows a similar cross-docking pattern but fails the same check even under
a zinc-aware scoring function (ROC-AUC ~0.50–0.55), so it is reported as a
methodologically unresolved watch item rather than a confirmed risk.

## Layout

```
target_validation/   pocket comparison, ChEMBL bioactivity, stats, figures, writeups, article
docking/              cross-docking (Vina + AD4Zn) against the flagged off-targets
common/               small Open Targets API client used by build_profile.py
docker/               conda environment + Dockerfile for the compute environment
```

Reproduction order and each script's output are listed in the article's
"Data and code availability" table (Section: end of
`OFFTARGET_SAFETY_ARTICLE_EN.md`).

## Notes on this extraction

- Raw/intermediate data (structures, ChEMBL pulls, docking poses) is not
  included — rerun the scripts to regenerate. Only final figures and
  written results are kept.
- A handful of cross-references in the writeups point to files that live
  in the parent (private) project and aren't included here: the TRPV1
  MD residence-time pilot (`md_residence_time/STATUS.md`), the docking
  smoke test (`docking/SMOKE_TEST.md`), and the ZINC diversity-triage
  writeup (`docking/ZINC_TRIAGE.md`). They're cited for context, not
  required to run anything in this repo.
- This is a hypothesis-generating computational screen, not a
  confirmatory safety study — see the article's Limitations section.
