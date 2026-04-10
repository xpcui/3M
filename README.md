# 3M

Title: Identifying the best predictive biomarker in pharmacogenomics through multiple comparisons with the best

Authors: Song Zhai, Judong Shen, Jason C. Hsu, and Xinping Cui

Author of the code: Song Zhai (song.zhai@merck.com)

## R Functions

The R functions can be found in **./R**:
- Functions-Fig.R: Functions to generate plots
- Functions-3M.R: Our proposed 3M method
- Functions-Ranking.R: Ranking-based methods, including p-value ranking, Lasso, and Interaction Trees (IT)
- Functions-IT_p1.R: Supplementary functions for Interaction Trees (IT)
- Functions-IT_p2.R: Supplementary functions for Interaction Trees (IT)

The code was produced with the following versions of R and packages:

**R version 4.4.1**

Platform: x86_64-w64-mingw32/x64 (64-bit)

Running under: Windows 10 x64

**attached base packages:**

stats   graphics   grDevices   utils   datasets   methods   base

**other attached packages:**

ggplot2   dplyr   data.table   glmnet   gridExtra

## Examples

Examples to run functions can be found in **./vignettes**:
- reproduce_results_Fig2.Rmd: Code to reproduce Fig2 in the manuscript
- reproduce_results_Fig3.Rmd: Code to reproduce Fig3 in the manuscript
- reproduce_results_FigS1.Rmd: Code to reproduce FigS1 in the manuscript
- reproduce_results_FigS2.Rmd: Code to reproduce FigS2 in the manuscript
- reproduce_results_FigS4.Rmd: Code to reproduce FigS4 in the manuscript
- reproduce_results_FigS5.Rmd: Code to reproduce FigS5 in the manuscript
- reproduce_results_FigS6.Rmd: Code to reproduce FigS6 in the manuscript
- reproduce_results_FigS7.Rmd: Code to reproduce FigS7 in the manuscript
- reproduce_results_Tab2.Rmd: Code to generate results based on the IMPROVE-IT mockup data

The genotype data and the IMPROVE-IT mockup data can be found in **./Data**:
- geno.Rdata: Genotype matrix for simulations
- ./IMPROVE-IT/mockup.Rdata: IMPROVE-IT mockup data
