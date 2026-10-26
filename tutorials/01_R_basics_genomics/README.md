# Module 1 — R Basics for Genomics

Part of the **R / Linux / Omics** teaching series (Module 1 of 3: **R basics → Linux → GWAS**).

## What this is
A self-contained Google Colab notebook that teaches R basics through a genomics lens, in a gentle
step-by-step style (modelled on the HBC *Intro to R* course): learning objectives up front, one concept
at a time, a **single dog dataset built up and reused throughout**, and a small **Exercise** box after the
key sections. Mostly base R so it reads simply. All data is inline — **nothing to download**.

## Run it (participants)
1. Upload `Module1_R_Basics_Genomics.ipynb` to Google Colab (or open from GitHub/Drive).
2. **Runtime → Change runtime type → `R` → Save.** (Requires the R kernel, not Python.)
3. Run cells top to bottom with **Shift+Enter**.

Only a free Google account is needed. Runtime ~15–25 min.

## Covers (syllabus order, one small example each)
R & RStudio (vs Colab's R runtime) + assignment · data types · operators · data structures · vectors · factors ·
matrices (genotypes) · data frames (phenotypes) · lists · functions · reading in & inspecting data ·
data wrangling across data frames/matrices/lists (base R, with a light dplyr pointer) ·
writing to file & exporting plots (`write.csv`/`ggsave`, with a Colab download note).

Three in-line **Exercise** boxes plus a closing summary, two "try it yourself" tasks, and further-reading pointers
(including the HBC *Intro to R* course this style follows).

## Notes for the instructor
- The Bioconductor cell (Section 8) compiles packages on first Colab run (~3–6 min) and is clearly marked optional; Modules 2–3 do not depend on it.
- Every code cell was executed and verified under R 4.3.3 before release.
- Toy/synthetic data only, per the project's teaching-vs-production separation. The CFA15/IGF1 references foreshadow the real dog GWAS in Module 3.
