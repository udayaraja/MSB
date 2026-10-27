# iMAP-Dairy GWAS tutorial (subset)

Runnable, Colab-compatible walkthrough of a complete GWAS on a small DKC subset
(588 cows × 12,747 SNPs; BTA1 + BTA14 full density, other autosomes thinned; traits my305, rfi, mastitis).

## Run on Google Colab
1. Open directly: [Open in Colab](https://colab.research.google.com/github/udayaraja/MSB/blob/main/tutorials/03_GWAS_iMAP/GWAS_tutorial_iMAP.ipynb)
2. Runtime → Run all. The setup cell installs PLINK2 + GCTA from bioconda (~1–2 min).
3. The data-retrieval cell shallow-clones this MSB repo automatically and reads `data/` from `tutorials/03_GWAS_iMAP/data/` within it — no Drive mount, no manual upload.
Total run time after install: < 1 min.

## Run locally (Linux / WSL; macOS see note)
    conda env create -f environment.yml && conda activate imap-gwas
    jupyter lab GWAS_tutorial_iMAP.ipynb
macOS Apple silicon: GCTA from bioconda may need an x86 (Rosetta) env: `CONDA_SUBDIR=osx-64 conda env create -f environment.yml`.

## Contents
| path | what |
|---|---|
| `GWAS_tutorial_iMAP.ipynb` | tutorial, executed with all outputs/figures |
| `data/` | PLINK subset + `pheno.txt` (sanitized/synthetic subset, safe for public tutorial use) |
| `figures/` | 13 PNGs: phenotypes, QC, LD decay, PCA, GRM, h², Manhattan, QQ, Miami, model comparison, regional, conditional, effect plots |
| `out/` | `gwas_all_results.tsv.gz` (all SNPs × traits × models), `lead_snps.tsv`, `summary_by_trait.tsv` |
| `environment.yml` | conda environment |

Source: `ml_sets/setB_dense` (PLINK 1.9 QC: --geno 0.05 --maf 0.01, low-density animals removed) + `ml_targets.csv`.
Tools tested: PLINK v2.0.0-a.7.8, GCTA 1.94.1.
