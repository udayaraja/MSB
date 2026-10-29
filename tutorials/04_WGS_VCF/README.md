# Module 4 — WGS: from reads to a checked VCF (shell commands)

Part of the Day 1 series (**R basics → Linux → GWAS**; this module applies the Linux shell skills of Module 2 to whole-genome sequence data).

Real whole-genome reads of a Holstein calf and its two parents go through FastQC, fastp, BWA-MEM and bcftools to a quality-checked VCF; then a complete VCF of 132 Holstein is audited (QA/QC). Every step is a plain shell command: each command is printed and its full output is shown.

## Run it (participants)
[**Open in Colab**](https://colab.research.google.com/github/udayaraja/MSB/blob/main/tutorials/04_WGS_VCF/Module4_WGS_VCF_Shell.ipynb)

1. Open the notebook in Google Colab and click **Connect** (keep the default Python runtime; the cells use the shell via `%%bash`).
2. Run the cells top to bottom (or **Runtime → Run all**). Colab shows "This notebook was not authored by Google": click **Run anyway**.
3. Setup (cells 0a–0d: data and tools) takes about 3–4 minutes; Parts A–B about 1 minute; Part C under 1 minute. Total about 9 minutes on a free Colab. The cells can be run again in the same session.
4. The last cell (optional) downloads your VCFs, the MultiQC report and all figures.

Each shell cell prints the commands it runs in green (`$ ...`) followed by their complete output. The only Python cells are the short **Figure** cells that plot the tables made by the shell commands.

## What you do
| Part | Steps |
|---|---|
| **0 Setup** | get the data, install FastQC, MultiQC, fastp, BWA, samtools, bcftools (micromamba, pinned versions) |
| **A Look at the files** | 1 FASTQ (reads) · 2 FASTA (reference) · 3 VCF (what we want to make) |
| **B Reads to VCF** | 4 FastQC · 5 fastp trim · 6 FastQC again · 7 BWA index · 8 BWA-MEM align · 9 call variants · 10 filter · 11 read your VCF · 12 view a variant in the reads · 13 is the call set good (Ts/Tv, Mendelian errors) · 14 what we made |
| **C QA/QC of a complete VCF** | 1 what is in it · 2 Ts/Tv · 3 QC tags · 4 missingness · 5 MAF · 6 Hardy-Weinberg · 7 odd animals (heterozygosity, gtcheck, PCA) · 9 filter · 10 everyday commands · 11 final checks |

Nine figures are drawn along the way (read quality, coverage and insert size, call quality, filter trade-off, allele frequencies, Hardy-Weinberg, heterozygosity per animal, pairs of animals and PCA, variant density).

## Run locally (Linux / WSL)
    git clone https://github.com/udayaraja/MSB.git
    cd MSB/tutorials/04_WGS_VCF
    bash Module4_WGS_VCF_Shell.sh

The script does the same as the notebook (tools are installed to `~/genart-lab` with micromamba; results go to `./Day1_WGS`; figures are saved as PNG files in `Day1_WGS/results/wgs/figures/`).

## Contents
| path | what |
|---|---|
| `Module4_WGS_VCF_Shell.ipynb` | the tutorial notebook (Colab) |
| `Module4_WGS_VCF_Shell.sh` | the same commands as one plain shell script |
| `data/reads/` | real reads of a calf (`dna_A`) and its two parents (`dna_B`, `dna_C`), cut to three 400-kb windows (genes ALB, DGAT1, UMOD) so that every step takes seconds; `animals.tsv`, `SOURCES.md` |
| `data/ref/` | `mini.fa` (the 3-contig reference: windows of cattle genome ARS-UCD2.0) and `regions.tsv` |
| `data/session3/` | `damona132_chr14_0-10Mb.vcf.gz` (+ `.tbi`): published genotypes of 132 Holstein, chromosome 14, first 10 Mb, 41,151 SNPs; `SOURCES.md` |
| `figures/` | 9 PNGs from a full run of the notebook |

## Tools tested
FastQC 0.13.0, MultiQC 1.35, fastp 1.4.0, BWA 0.7.19, samtools 1.24, bcftools 1.24 (full run in Google Colab, 27 Oct 2026; the same numbers on Linux).

## Sources and citation
* Reads: European Nucleotide Archive, project PRJEB53518 (GIGA, University of Liège), public data; Oget-Ebrad C, Kadri NK, Moreira GCM, Karim L, Coppieters W, Georges M, Druet T (2022) Benchmarking phasing software with a whole-genome sequenced cattle pedigree. BMC Genomics 23:130.
* 132-animal VCF: Druet T, Forneris NS (2026) Zenodo doi:10.5281/zenodo.23011638, CC BY 4.0.
* How the data were cut and what was changed: `data/reads/SOURCES.md`, `data/session3/SOURCES.md`.

Author of this module: Pratik Pathade (originally built for the GenArt_workshop series; public copy also at https://github.com/PratikPathade/GenART-Day1-WGS). This copy, hosted in MSB for open access, has Part C (QA/QC) amended by Udayaraja GK: the genome-build cross-check (previously C8) was removed, and the data-retrieval cell clones this repo (`udayaraja/MSB`) directly rather than GenArt_workshop.
