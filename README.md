# Computational Investigation of HIV-1 Protease

**Sequence Diversity, Conservation and Antiretroviral Drug-Resistance Interpretation**

## Project overview

This undergraduate computational biochemistry project combines HIV-1 protease sequence analysis, multiple-sequence alignment, conservation analysis, structural inspection and Stanford HIVdb resistance interpretation.

**Research question:** How does sequence variation in a small sample of HIV-1 protease sequences relate to conservation patterns, structural context and Stanford HIVdb protease-inhibitor resistance predictions?

## Main findings

- 31,769 raw Los Alamos records were screened.
- 735 records were exactly 297 nt and contained only A/C/G/T.
- Exact duplicate removal left 700 unique sequences.
- A pilot sample of 50 unique sequences was selected (random seed 42).
- MAFFT alignment contained 50 sequences and 297 alignment columns.
- Translation produced 50 sequences of 99 amino acids; no internal stop codons were detected.
- Sample conservation was lowest at position 93 (52%), positions 37 and 63 (54% each), position 77 (62%), position 36 (64%) and position 41 (66%).
- Frequent substitutions against the local P04585-derived reference included S37N (27/50), L63P (27/50), I93L (24/50), V77I (19/50), R41K (17/50), E35D (13/50), M36I (13/50) and I13V (12/50).
- Stanford HIVdb version 10.2 predicted all 50 sequences susceptible to darunavir/r.
- Two sequences had high-level predicted resistance to nelfinavir.

These results describe this pilot sample only. They are not estimates of global prevalence.

## Workflow

1. Retrieve the UniProt P04585 Gag-Pol precursor and extract the 99-aa protease segment.
2. Retrieve and inspect PDB 1HIV in UCSF ChimeraX.
3. Download subtype-B protease nucleotide records from the Los Alamos HIV Sequence Database.
4. Filter to 297-nt A/C/G/T records and remove exact duplicates.
5. Select 50 unique sequences using random seed 42.
6. Align with MAFFT.
7. Translate and analyse conservation/substitutions with Biopython.
8. Submit sequences to Stanford HIVdb and interpret the exported summary tables.

## MAFFT command

```bash
mafft --auto --thread -1 data/sequences/HIV1_protease_sample50.fasta > data/sequences/HIV1_protease_aligned.fasta
```

## Project structure

```text
HIV1_Protease_Project/
├── README.md
├── data/
│   ├── sequences/
│   ├── structures/
│   └── ligands/
├── scripts/
├── analysis/
├── notebooks/
├── results/
│   ├── figures/
│   └── tables/
├── literature/
└── docs/
    └── HIV1_Protease_Computational_Project_Report.docx
```

## Important interpretation caveats

- The sample of 50 is a pilot sample and is not representative of global HIV-1 diversity.
- Treatment histories were not established. The sample was **not** filtered to ensure drug-naive sequences.
- Sequence variation is not automatically drug resistance. Stanford HIVdb distinguishes major/accessory resistance mutations from other substitutions.
- HIVdb predictions are computational interpretations, not phenotypic assays or clinical outcomes.
- ChimeraX reported 14 donor–acceptor assignments in the 1HIV contact analysis. Explicit hydrogen atoms were not present for those assignments, so describe them as geometric proximity observations—not confirmed hydrogen bonds.
- The local P04585-derived 99-aa reference has valine at protease position 3, whereas Stanford HIVdb's subtype-B consensus reference has isoleucine at position 3. This can reverse the notation (local V3I versus HIVdb I3V) for that position. Do not merge local substitution counts with HIVdb mutation categories without harmonizing the reference and numbering.

## Reproducibility

The project used Python/Biopython for sequence processing, MAFFT for alignment, UCSF ChimeraX for structural inspection and Stanford HIVdb algorithm version 10.2 (algorithm date 2026-04-26) for resistance interpretation.

Before public release, record the exact installed versions of Python, Biopython, MAFFT and UCSF ChimeraX. Preserve the scripts and document the source and processing of all input data.

## Resources

- UniProtKB accession P04585
- Protein Data Bank entry 1HIV
- Los Alamos National Laboratory HIV Sequence Database
- MAFFT
- Biopython
- Stanford HIV Drug Resistance Database (HIVdb)
- UCSF ChimeraX

## Report

The full scientific report is in `docs/HIV1_Protease_Computational_Project_Report.docx`.
