# Local (pocket-residue) ESM-2 embedding comparison, TRPV1 vs TRPV2-6

Follow-up to `OFFTARGET_SAFETY_WRITEUP.md` Section 2.3's methodological
correction: whole-protein sequence identity is a poor proxy for
binding-site-level risk within the TRPV family, because the
vanilloid-pocket region is reported to be structurally/sequence-wise
conserved across TRPV1-3 independent of overall identity, while TRPV4
achieves selectivity via point substitutions within an otherwise similar
pocket. That claim was previously qualitative (literature-sourced). This
adds a quantitative check using ESM-2 protein-language-model embeddings,
scoped specifically to avoid reproducing the same whole-protein-identity
limitation it's meant to test.

## Method

**Why not whole-sequence embedding similarity.** A single mean-pooled
ESM-2 vector per protein would just re-derive a whole-protein similarity
ranking -- the same limitation already flagged for raw % identity. Instead:

1. **Pocket-residue mapping** (`pocket_residue_alignment.py`, local/CPU,
   Biopython global pairwise alignment, BLOSUM62): TRPV1's four
   vanilloid-pocket-lining residues (Y511, S512, T550, E570 -- confirmed
   against UniProt Q8NER1's canonical sequence, matching
   `md_residence_time/STATUS.md` and the writeup's Section 2.3/3.2) were
   mapped onto the aligned position in each paralog (UniProt Q9Y5S1/TRPV2,
   Q8NET8/TRPV3, Q9HBA0/TRPV4, Q9NQA5/TRPV5, Q9H1D0/TRPV6). All four
   positions sit in a strongly conserved local motif in every paralog
   (e.g. TRPV1's T550 window `LALGW[T]NMLYY` vs. TRPV4's `LVLGW[M]NALYF`)
   -- alignment confidence is high; only TRPV6 has a real gap (no aligned
   residue for S512). Full mapping + physicochemical-class comparison:
   `pocket_residue_alignment.csv`.
2. **ESM-2 embedding** (`esm_pocket_embedding_similarity.py`,
   `esm2_t33_650M_UR50D`, run on a vast.ai GPU instance per project
   convention -- fair-esm, ~2.5GB checkpoint, one forward pass per full
   sequence so each residue's representation carries real sequence
   context, not a stripped-out peptide fragment): per-residue final-layer
   representations extracted for all 6 sequences, then a 5-residue window
   (+/-2) mean-pooled around each of the 4 pocket positions and
   concatenated into one local pocket vector per protein.
3. **Comparison**: TRPV1's local vector vs. each paralog's, by cosine
   similarity, computed **only over positions aligned in both proteins**
   (see bug note below) -- plus a whole-protein mean-pooled cosine
   similarity computed the same way, reported side by side specifically
   to show it does *not* discriminate (predicted going in, per point 1).

**Bug caught and fixed before trusting the result**: the first pass
zero-padded TRPV6's unaligned S512 position into the concatenated vector
instead of excluding it. This mechanically depressed TRPV6's local cosine
similarity (0.7825) purely because concatenating a zero-block into one
operand penalizes cosine similarity regardless of the other 3 positions'
true resemblance -- not a real biology signal. Fixed by restricting the
comparison, per paralog, to only the positions aligned in that paralog
(TRPV6: n=3 of 4; all others: n=4 of 4). TRPV6's corrected score (0.9407)
is in line with TRPV5, not the outlier the buggy version suggested.

## Results

| Paralog | Local pocket-residue cosine sim | Whole-protein mean-pooled cosine sim | Positions compared |
|---|---|---|---|
| TRPV2 | 0.9646 | 0.9933 | 4/4 |
| TRPV3 | 0.9648 | 0.9963 | 4/4 |
| TRPV4 | 0.9689 | 0.9976 | 4/4 |
| TRPV5 | 0.9344 | 0.9929 | 4/4 |
| TRPV6 | 0.9407 | 0.9930 | 3/4 (gap at S512) |

Raw data: `esm_pocket_similarity.csv`.

## Interpretation

- **The whole-protein baseline is exactly as uninformative as predicted**:
  0.993-0.998 across all five paralogs, a ~0.005 spread that tracks
  nothing biologically interesting -- it does not even reproduce the
  whole-protein % identity ranking from Table 1 (TRPV4 has the *highest*
  global cosine here despite being the 2nd-lowest-identity paralog in that
  table). This confirms mean-pooled whole-sequence embeddings are not a
  useful lens for this question, independent of the identity-based
  argument already made in Section 2.3.
- **The local metric does separate a real axis**: TRPV2/TRPV3/TRPV4 cluster
  tightly (0.965-0.969) and TRPV5/TRPV6 sit measurably lower (0.93-0.94).
  This lines up with the literature-based claim that TRPV5/6 (Ca2+-selective,
  distinct pore/gating architecture) diverge more at the pocket than
  TRPV2-4 do -- a second, independent (embedding-based, not
  identity-based) line of evidence for that specific claim.
- **It does *not* corroborate the TRPV4-specific part of the claim.**
  TRPV4 has zero same-physicochemical-class matches at all four pocket
  positions (`pocket_residue_alignment.csv`: Y511->S, S512->F, T550->M,
  E570->Q, all class changes) -- more substitution than TRPV2 or TRPV3 by
  that discrete measure -- yet it has the *highest* local ESM-2 cosine
  similarity of the whole group. Plausible reading: a 5-residue window
  mean-pool is dominated by the strongly conserved flanking motif shared
  across TRPV1-4 (see the alignment windows above), so the embedding
  signal here is closer to "how conserved is the local structural context"
  than "how conservative is the specific substitution" -- those are related
  but not the same question, and this method answers the former more than
  the latter. Should be read as a limitation of the window-pooling
  approach, not as evidence against the literature's TRPV4 point-substitution
  argument (which is about the functional consequence of one specific
  residue change, not local-context similarity).

## Limitations

- n=1 embedding per protein (no ensemble/uncertainty estimate on the
  cosine values); a ~0.03-0.04 spread (TRPV2-4 vs. TRPV5/6) is the
  headline signal and should be treated as suggestive, not as a
  statistically tested separation (no ground-truth labels exist to build
  a discrimination check here, unlike Section 3.4's ChEMBL-based check).
- TRPV6's comparison uses one fewer position (3 vs. 4) than the others,
  because of a real alignment gap, not a filtering choice -- makes it not
  perfectly apples-to-apples with the other four paralogs' scores.
- 5-residue window size (+/-2) was not tuned; a narrower (single-residue)
  or wider window could shift the balance between "substitution-specific"
  and "local-context" signal discussed above -- worth a sensitivity check
  before leaning on this method for a stronger claim.
- Not attempted here: applying this same local-embedding approach to the
  RARG/MMP3 off-target comparison (Section 3 of the main writeup) --
  meaningless there, since ESM per-residue positional comparison requires
  a real sequence alignment between homologous proteins, and RARG/MMP3 are
  unrelated to TRPV1 by sequence (that's the whole reason the fpocket
  geometric screen was used for them instead).

## Provenance

| Stage | Script | Output |
|---|---|---|
| Pocket-residue alignment | `pocket_residue_alignment.py` | `pocket_residue_alignment.csv` |
| ESM-2 local/global embedding similarity | `esm_pocket_embedding_similarity.py` | `esm_pocket_similarity.csv` |
