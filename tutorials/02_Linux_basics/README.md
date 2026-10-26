# Module 2 — Linux / Command Line for Bioinformatics

Part of the **R / Linux / Omics** teaching series (Module 2 of 3: R basics → **Linux** → GWAS).

## What this is
A self-contained Google Colab notebook that teaches the command line through real bioinformatics file
types, in the same gentle step-by-step style as Module 1: learning objectives up front, one tool at a
time, small **Exercise** boxes, and a toy dataset (a **VCF**, PLINK **`.bim`/`.fam`**, and a **FASTA**)
built and reused throughout. **Nothing to download.**

## Run it (participants)
1. Open `Module2_Linux_Bioinformatics.ipynb` in Google Colab.
2. **Keep the default runtime** — this module uses the shell, *not* R (no runtime change needed).
3. Run cells top to bottom. **Run the "Setup" cell first** — it creates the toy files.

Each code cell begins with `%%bash`, which runs the whole cell as Linux commands. Only a free Google account is needed.

## Covers
the shell, filesystem & paths (`pwd`, `ls`, `cd`, `mkdir`, absolute vs relative, tab-completion & history) ·
**wildcards** (`*`, `?`) · file management (`cp`, `mv`, `rm`, with a no-undo caution) ·
viewing files (`cat`, `head`, `tail`) · counting (`wc`) · columns (`cut`) · searching (`grep`) ·
tallying (`sort | uniq -c`) · stream editing (`sed`) · column filtering (`awk`) ·
**pipes** (`|`) and **redirection** (`>`) · **shell variables** · **loops** · writing a reusable **`.sh` script** ·
plus a concept note on where it scales (conda/mamba, SLURM/PBS). Draws on Data Carpentry's *Shell Genomics*.

## Notes for the instructor
- Uses `%%bash` cells (not `!`-prefixed lines) so multi-line scripts and `awk`'s `$1`/`$2` work without Colab's `$`-interpolation quirks.
- Every command was executed in a real bash shell and its output verified before release (e.g. 3 variants on chr15, the IGF1-region filter, per-chromosome counts, PASS filtering).
- Toy files only; the CFA15 / IGF1 region references foreshadow the real dog GWAS in Module 3.
