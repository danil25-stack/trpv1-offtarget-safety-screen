# ChEMBL bioactivity pull for RARG and MMP3

Follow-up to [skin_offtarget_pocket_comparison.md](skin_offtarget_pocket_comparison.md).
Closes two gaps flagged against `OFFTARGET_SAFETY_ARTICLE_EN.md`: (1) the
Discussion's promised "literature/bioactivity-database check" for RARG and
MMP3 was never actually run, and (2) the cross-docking results (Section 3.3)
compare Vina scores *across* receptors with no check that Vina's scoring
function discriminates real binders from non-binders *within* either
receptor at all -- which is exactly the "cross-receptor score comparisons
carry additional uncertainty beyond the single-receptor ranking use case"
limitation already flagged in the Discussion.

Script: `fetch_chembl_bioactivity.py`. Raw data:
`data/raw/chembl/{RARG,MMP3}_activities_raw.csv`. Classified data:
`data/processed/skin_offtarget_screen/{RARG,MMP3}_chembl_actives_inactives.csv`.

## Method

1. **Target resolution**: ChEMBL `/target?target_components__accession=<UniProt>`,
   confirmed `SINGLE PROTEIN` / `Homo sapiens` for both. RARG (UniProt
   P13631) -> `CHEMBL2003`; MMP3 (UniProt P08254) -> `CHEMBL283`.
2. **BindingDB tried first, found non-functional.** Both documented REST
   endpoints (`getLigandsByUniprot`, `getLigandsByUniprots`) return HTTP 200
   with an empty body for every UniProt accession tested, including targets
   unrelated to this project with deep BindingDB coverage (e.g. carbonic
   anhydrase II, P00918) -- this is a dead/broken service response, not an
   absence of RARG/MMP3 data specifically. ChEMBL is used as the sole
   bioactivity source below; a BindingDB bulk TSV download (`BindingDB_All`)
   remains a fallback if independent corroboration is needed later, but
   wasn't pulled here (multi-GB file, out of scope for this pass).
3. **Activity retrieval**: all ChEMBL `activity` records for each target
   with `standard_type` in {IC50, EC50, Ki, Kd} and `standard_units = nM`,
   paginated (1000/page). The EBI ChEMBL API returned intermittent,
   non-reproducible HTTP 500s on otherwise-valid queries during this pull
   (confirmed transient -- immediate retry of the identical request often
   succeeded) -- `fetch_chembl_bioactivity.py`'s `get_json()` retries every
   request up to 6x with backoff; re-running the script is expected to need
   a handful of retries and should not be read as a query-correctness
   problem if it happens again.
4. **Classification per molecule** (`classify()` in the fetch script):
   - Records flagged by ChEMBL's own `data_validity_comment` are excluded.
   - If a molecule has at least one `=` (measured) record, its most potent
     value drives the call: **active** if <=1000 nM, **inactive** if
     >10000 nM, **ambiguous** in between (deliberately not forced into
     either class -- this is a shared, untuned cutoff, not
     target-optimized).
   - If a molecule has *no* `=` record but does have a `>` (tested, no
     activity detected) record, and the strongest reported `>` threshold is
     >=10000 nM, it's labeled **inactive_no_binding_detected** -- a cleaner
     negative than a weak measured value, but only used when no measured
     value exists for that molecule (a `>` record alongside a `=` record for
     the same molecule is a less-sensitive assay, not real information, so
     it's dropped in that case).

## Results

| Target | Total activity records | Unique molecules classified | Active (<=1000 nM) | Inactive (measured >10000 nM) | Inactive (no binding, `>`-only) | Ambiguous |
|---|---|---|---|---|---|---|
| RARG | 695 | 409 | 242 | 2 | 96 | 69 |
| MMP3 | 3035 | 2228 | 1387 | 199 | 265 | 377 |

RARG's most potent record (Ki 0.04 nM) and MMP3's activity volume (thousands
of records, consistent with MMP3's long medicinal-chemistry SAR history as a
hydroxamate-inhibitor target) are both consistent with what's known about
these targets qualitatively -- not independently confirmed compound-by-compound
here, but not a data sanity failure either.

## Caveats

- **RARG's negative set is thin even after adding `>`-only inactives** (98
  total vs. MMP3's 464). RARG is a well-studied nuclear receptor with a
  strong medicinal-chemistry bias toward publishing potent agonists;
  genuinely inactive compounds are less likely to be reported at all. Any
  downstream enrichment/discrimination check (the natural next step -- dock
  a sampled active/inactive set into each receptor and check whether Vina
  score actually separates them) will have much less statistical power for
  RARG than for MMP3 as a result. Worth flagging in the article rather than
  silently treating both target's validation as equally powered.
- **No property-matching between actives and inactives yet.** These
  classified sets are raw ChEMBL pulls, not curated decoy sets -- if actives
  and inactives differ systematically in MW/logP (plausible: ChEMBL's active
  compounds for a mature SAR series tend to cluster tighter in chemical
  space than the historical inactive/discontinued-series compounds), a
  future docking-based discrimination check could partly reflect that
  property bias rather than genuine pocket complementarity. Should be
  addressed at the sampling stage before docking, not here.
- **Assay heterogeneity not filtered.** `assay_type`/`assay_description`
  are retained in the raw and classified CSVs but not filtered on here
  (e.g. binding vs. functional/cellular assays are pooled) -- reasonable for
  a first pass at this scale, but a stricter version could restrict to
  binding assays (`assay_type = B`) only for a cleaner potency comparison.
- Neither dataset has been used for docking yet -- this pull only
  establishes ground-truth active/inactive labels. The next step is
  sampling a size-manageable, MW/logP-matched active + inactive subset per
  target and redocking with the same Vina/Meeko pipeline and boxes already
  validated in `skin_offtarget_pocket_comparison.md` and
  `OFFTARGET_SAFETY_ARTICLE_EN.md` Section 2.4, then checking whether Vina
  score actually enriches actives over inactives (e.g. ROC-AUC) at RARG and
  MMP3 specifically -- this is the check that would make the cross-docking
  section's absolute score comparisons defensible.
