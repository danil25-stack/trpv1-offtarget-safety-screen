# Skin off-target candidate list (non-TRPV/TRPA)

Purpose: candidate list of skin-expressed proteins for a future pocket-vs-pocket
structural similarity screen against TRPV1's vanilloid pocket -- i.e. off-target
risk from *binding-site* similarity + *tissue-exposure* overlap, independent of
the TRPV/TRPA family (already covered in
[comparison_all_candidates.md](profiles/comparison_all_candidates.md) and
[TRPV1.md](profiles/TRPV1.md)).

Built by `fetch_skin_targets.py`. Raw data: `data/raw/hpa/skin_tissue_specific_genes.tsv`.
Annotated data: `data/processed/skin_offtarget_screen/skin_targets_annotated.csv`.

## Method

1. **Source list**: Human Protein Atlas RNA tissue-specificity API
   (`proteinatlas.org/api`), filtered to genes classified `Tissue enriched`,
   `Group enriched`, or `Tissue enhanced` for skin -- i.e. genes whose RNA
   expression is markedly higher in skin than most other tissues. 604 genes
   returned (2026-08-10 pull).
2. **Druggability filter**: each gene cross-referenced against Open Targets
   `targetClass` + small-molecule `tractability` (the same fields used to build
   the TRPV1/SCN9A/etc. profiles). Open Targets' SM tractability bucket
   includes a druggability-prediction pocket flag (`High-Quality Pocket` /
   `Med-Quality Pocket`), which is what separates "has a real ligand-binding
   cavity" from the many structural/barrier genes (keratins, small proline-rich
   proteins, late cornified envelope genes, etc.) that dominate the raw skin
   list -- **486 of the 604 genes have no Open Targets target class at all**
   and were dropped.
3. Genes tiered by strength of pocket evidence (below). TRPV3 and TRPM1
   (already-covered TRP family) excluded per scope.

## Caveat carried over from the TRPV1 profile

HPA/GTEx-style bulk "skin" expression does not resolve dorsal root ganglion or
peripheral nerve-ending expression specifically -- a target could be
DRG/nociceptor-terminal-expressed (directly relevant to a topical analgesic's
site of action) without showing up as "skin"-enriched in bulk tissue RNA, and
conversely a keratinocyte/skin-enriched gene here isn't necessarily present at
the nerve terminal itself. This list captures bulk-tissue co-localization risk
at the application site, not nociceptor-specific co-expression.

## Tier 1 -- confirmed ligand-binding pocket (High/Med-Quality Pocket), n=26

Strongest candidates for a pocket-similarity comparison -- Open Targets'
structure-based druggability prediction has already flagged a real cavity.

| Gene | Ensembl | HPA specificity | Target class | Approved SM drug exists |
|---|---|---|---|---|
| ADRB2 | ENSG00000169252 | Tissue enhanced | Family A GPCR, Adrenergic receptor | Yes |
| CA12 | ENSG00000074410 | Tissue enhanced | Enzyme, Lyase (carbonic anhydrase) | Yes |
| RARG | ENSG00000172819 | Tissue enhanced | Nuclear hormone receptor (RAR) | Yes |
| BLMH | ENSG00000108578 | Tissue enhanced | Enzyme, Hydrolase | No |
| MMP3 | ENSG00000149968 | Tissue enhanced | Metallo protease (M10A) | No |
| RORA | ENSG00000069667 | Tissue enhanced | Nuclear hormone receptor (RORα) | No |
| CD207 | ENSG00000116031 | Tissue enriched | Membrane receptor (Langerin) | No |
| FCER1A | ENSG00000179639 | Tissue enhanced | Membrane receptor (IgE Fc receptor) | No |
| GSTA3 | ENSG00000174156 | Group enriched | Enzyme, Transferase | No |
| SULT1E1 | ENSG00000109193 | Tissue enhanced | Enzyme, Transferase | No |
| SULT2B1 | ENSG00000088002 | Tissue enhanced | Enzyme, Transferase | No |
| FABP9 | ENSG00000205186 | Tissue enriched | Fatty acid binding protein | No |
| ANXA8 | ENSG00000265190 | Tissue enhanced | Unclassified protein | No |
| DSP | ENSG00000096696 | Tissue enhanced | Unclassified protein (desmoplakin) | No |
| ITGB4 | ENSG00000132470 | Tissue enhanced | Unclassified protein (integrin β4) | No |
| SFN | ENSG00000175793 | Tissue enhanced | Unclassified protein (14-3-3σ) | No |
| CHP2 | ENSG00000166869 | Group enriched | -- | No |
| FBXW7 | ENSG00000109670 | Tissue enhanced | -- | No |
| GAN | ENSG00000261609 | Tissue enriched | -- | No |
| GLTP | ENSG00000139433 | Tissue enhanced | -- | No |
| NTF4 | ENSG00000225950 | Tissue enhanced | -- | No |
| RAB3D | ENSG00000105514 | Tissue enhanced | -- | No |
| SERPINB2 | ENSG00000197632 | Tissue enhanced | -- | No |
| SERPINB5 | ENSG00000206075 | Tissue enhanced | -- | No |
| SRY | ENSG00000184895 | Tissue enhanced | -- | No |
| TREX2 | ENSG00000183479 | Tissue enhanced | -- | No |

## Tier 2 -- SM ligand evidence but no explicit pocket flag, n=59 (TRPV3 excluded)

Real structures with a bound ligand or high-quality chemical probe exist, but
Open Targets' automated pocket-quality classifier didn't tag it -- still worth
including in a structural screen, just slightly lower prior.

Notable ones for a topical/nociceptor-adjacent context: **PTGS1** (COX-1),
**FGFR3**, **KLK5/KLK7/KLK8/KLK14** (kallikreins, desquamation proteases,
several already drug-discovery targets for Netherton syndrome), **LTB4R /
LTB4R2** (leukotriene B4 receptors -- inflammatory pain-adjacent GPCRs),
**HCAR2/HCAR3** (niacin receptors, flushing-relevant GPCRs), **GPRC5D**
(current multiple myeloma ADC/CAR-T target, skin-restricted expression is the
reason it's used there), **SLC6A2** (norepinephrine transporter), **NOD2**,
**PTK6**, **EPHB3**, **CYP26B1**, **TYR/TYRP1** (melanogenesis enzymes).

Full list in `data/processed/skin_offtarget_screen/skin_targets_annotated.csv`
(filter `pocket_quality == "Ligand evidence only"`).

## Tier 3 -- receptor/enzyme/channel class, no SM tractability evidence yet, n=22

Lowest-confidence tier: Open Targets classifies these into a druggable-looking
family (GPCR, kinase, protease, transporter) by sequence/domain, but has no
structural or ligand evidence yet. Would need an AlphaFold model rather than
an experimental structure for any pocket comparison. Includes **ATP12A**
(gastric-type H+/K+-ATPase, skin-expressed paralog of the drug-targeted gastric
pump), **FZD10**, **GPR12**, **GPR87**, **TNFRSF18**, several kallikreins/CYPs.

Full list in the CSV, filter `pocket_quality == ""` and `target_class != ""`.

## Explicitly excluded

- **TRPV3** (Group enriched, skin) and **TRPM1** (Group enriched, skin) --
  already covered under the TRPV/TRPA family analysis; out of scope per this
  request.
- The **486 genes with no Open Targets target class** -- overwhelmingly
  structural/barrier proteins (keratins `KRT*`, small proline-rich proteins
  `SPRR*`, late cornified envelope genes `LCE*`, corneodesmosin, filaggrin
  family, etc.). No ligand-binding pocket in the pharmacological sense; not
  useful for a docking/pocket-similarity screen regardless of skin-expression
  strength. Full list still in the CSV (`target_class == ""`) if a different
  use case needs it (e.g. literature-based liability check rather than
  structural).

## Next step (not done here per scope: "без докинга")

Tier 1 (and selectively Tier 2, prioritizing the pain/inflammation-adjacent
hits called out above -- PTGS1, LTB4R/LTB4R2, HCAR2/3) are the actual
candidates for the fpocket-based pocket-fingerprint comparison against TRPV1's
vanilloid pocket discussed separately.
