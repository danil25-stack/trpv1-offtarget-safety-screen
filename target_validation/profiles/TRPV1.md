# Target Profile: TRPV1 (transient receptor potential cation channel subfamily V member 1)

- Ensembl ID: `ENSG00000196689`
- Biotype: protein_coding
- Target class: Voltage-gated ion channel, Ion channel, Transient receptor potential channel
- Essential (cell-fitness, DepMap-derived): False

## Function
Non-selective calcium permeant cation channel involved in detection of noxious chemical and thermal stimuli (PubMed:11050376, PubMed:11243859, PubMed:11226139, PubMed:12077606). Seems to mediate proton influx and may be involved in intracellular acidosis in nociceptive neurons. Involved in mediation of inflammatory pain and hyperalgesia. Sensitized by a phosphatidylinositol second messenger system activated by receptor tyrosine kinases, which involves PKC isozymes and PCL. Activated by vanilloids, like capsaicin, and temperatures higher than 42 degrees Celsius (PubMed:37117175). Upon activation, exhibits a time- and Ca(2+)-dependent outward rectification, followed by a long-lasting refractory state. Mild extracellular acidic pH (6.5) potentiates channel activation by noxious heat and vanilloids, whereas acidic conditions (pH <6) directly activate the channel. Can be activated by endogenous compounds, including 12-hydroperoxytetraenoic acid and bradykinin. Acts as ionotropic endocannabinoid receptor with central neuromodulatory effects. Triggers a form of long-term depression (TRPV1-LTD) mediated by the endocannabinoid anandamine in the hippocampus and nucleus accumbens by affecting AMPA receptors endocytosis. {ECO:0000250|UniProtKB:O35433, ECO:0000269|PubMed:11050376, ECO:0000269|PubMed:11226139, ECO:0000269|PubMed:11243859, ECO:0000269|PubMed:12077606, ECO:0000269|PubMed:37117175}.

## Genetic constraint (gnomAD)

| Constraint | Observed | Expected | o/e | Score |
|---|---|---|---|---|
| syn | 514 | 452.6 | 1.14 (1.06-1.22) | -1.57e+00 |
| mis | 1061 | 1064.6 | 1.00 (0.95-1.05) | 4.01e-02 |
| lof | 83 | 79.3 | 1.05 (0.88-1.26) | 6.34e-30 |

_Low o/e for `lof` = fewer loss-of-function variants tolerated in the general population than expected by chance -> gene under purifying selection against LoF._

## Tractability

| Modality | Evidence |
|---|---|
| SM | Approved Drug |
| SM | Structure with Ligand |
| SM | High-Quality Ligand |
| SM | Druggable Family |
| AB | UniProt loc high conf |
| AB | GO CC high conf |
| AB | UniProt loc med conf |
| AB | UniProt SigP or TMHMM |
| PR | Small Molecule Binder |
| OC | Advanced Clinical |

## Safety -- curated liabilities (Open Targets)

- **neurotoxicity** (source: ClinPGx)

## Safety -- off-target family (paralog selectivity risk)

Human paralogs ranked by sequence identity. High identity = harder to achieve selectivity; annotated rows are paralogs with a known mechanism-based liability if inhibited/activated off-target.

| Paralog | % identity | Known liability if hit |
|---|---|---|
| TRPV2 | 45.0% | developmental cardiac contractility/Ca2+ handling (KO cardiomyopathy in growing hearts); adult whole-body KO and neutralizing-antibody dosing show no overt phenotype -- acute pharmacological inhibition looks low-risk, developmental/chronic exposure less characterized |
| TRPV4 | 41.9% | osmo/mechanosensation; activation implicated in cardiogenic pulmonary edema, BUT clinical TRPV4 *antagonist* GSK2798745 was well-tolerated through repeat dosing in heart-failure patients with no significant cardiovascular/pulmonary safety signal -- direction matters: TRPV1-antagonist off-target TRPV4 *inhibition* trends protective for edema, not causative |
| TRPV3 | 39.0% | keratinocytes/skin barrier & thermosensation (shares topical site of action); human GoF mutations cause Olmsted syndrome (painful keratoderma, pruritus, skin barrier failure) -- direct human proof this channel matters at the skin, though that's a GoF disease phenotype, not necessarily what antagonism would produce |
| TRPV5 | 28.4% | kidney distal-tubule Ca2+ reabsorption, rate-limiting step of renal Ca2+ handling -- off-target inhibition is a plausible hypercalciuria/renal Ca2+ handling risk |
| TRPV6 | 27.2% | intestinal Ca2+ absorption; KO mice show impaired intestinal Ca2+ uptake, reduced bone mineral density/osteoporosis, alopecia/dermatitis, and severely impaired male fertility -- the most systemically consequential TRPV paralog to hit off-target if these phenotypes translate to pharmacological inhibition |
| TRPA1 | n/a (not a sequence paralog) | co-expressed nociceptor channel, mechanistically adjacent (not a sequence paralog); clinical TRPA1 antagonists (LY3526318, Roche/Genentech compound) had acceptable Phase 1 safety but failed efficacy in OA/CLBP/DPNP and chronic cough -- off-target TRPA1 hit looks safety-neutral but track record of clinical failure is itself informative |

## Expression / localization

- Subcellular location: Cell membrane ; Multi-pass membrane protein, Cell projection, dendritic spine membrane ; Multi-pass membrane protein, Postsynaptic cell membrane ; Multi-pass membrane protein

Top tissues by baseline (bulk RNA) expression:

| Tissue | Median expression | Unit |
|---|---|---|
| endocrine pancreas | 28.90 | CPM(pseudobulk sum[counts]) |
| exocrine pancreas | 13.66 | CPM(pseudobulk sum[counts]) |
| pituitary gland | 10.44 | TPM |

_Caveat: this is bulk GTEx-derived expression and does not include dorsal root ganglion (DRG) or peripheral sensory neuron / skin nerve-ending samples specifically -- the tissue most relevant to a topical analgesic's site of action. Cross-check against DRG-specific expression atlases from the literature before concluding local target abundance at the intended application site._

## Additional caveats (literature, not in Open Targets)

- No monogenic human channelopathy comparable to SCN9A is established for TRPV1 — genetic validation here rests more on GWAS/expression/pharmacology than on rare high-penetrance human variants. Clinical precedent instead comes from topical capsaicin desensitization (Qutenza, an approved TRPV1 agonist patch) and from the failure of systemic TRPV1 antagonists due to on-target hyperthermia — an argument specifically in favor of a topical/local route for this target, not systemic.
- Off-target paralog literature pass (2026-08-10): of the five sequence paralogs, TRPV6 is the most systemically consequential to hit (KO mice: impaired intestinal Ca2+ absorption, reduced bone mineral density, alopecia/dermatitis, severe male infertility) and TRPV5 is a plausible renal Ca2+-handling risk (rate-limiting channel for distal-tubule Ca2+ reabsorption); neither had liability text in the original profile. TRPV4's 'pulmonary edema' liability is direction-dependent: the clinical antagonist GSK2798745 was well-tolerated in heart-failure patients with no cardiopulmonary safety signal, so off-target *inhibition* by a TRPV1 antagonist is not expected to reproduce the edema risk (that's an activation-direction liability). TRPV3 has direct human proof of consequence at the skin (Olmsted syndrome, from gain-of-function mutations) — relevant given TRPV1 and TRPV3 share the topical site of action, though GoF disease phenotype doesn't necessarily predict what pharmacological antagonism would do. TRPV2 looks acutely low-risk (adult whole-body KO and neutralizing-antibody dosing show no overt phenotype) but has an uncharacterized developmental cardiac role. TRPA1 (mechanistically adjacent, not a sequence paralog) has an acceptable Phase 1 safety record in two independent clinical antagonist programs (Lilly's LY3526318, Roche/Genentech), both of which nonetheless failed on efficacy — so an off-target TRPA1 hit looks safety-neutral, not liability-free by default, more 'clinically unproven for this indication.' Net: no paralog in this family blocks the topical program on safety grounds alone, but TRPV5/TRPV6 warrant an explicit selectivity check in any hit-to-lead off-target panel, which wasn't previously called out.

## Associated diseases (top by overall association score)

| Disease | Score |
|---|---|
| migraine disorder | 0.629 |
| Pain | 0.611 |
| Headache | 0.610 |
| Fever | 0.609 |
| arthritic joint disease | 0.607 |
| Cough | 0.605 |
| Nasal congestion | 0.603 |
| common cold | 0.602 |
| Back pain | 0.598 |
| pharyngitis | 0.598 |
