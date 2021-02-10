How fast do biological rock crusts grow?
========

[![Twitter Follow](https://img.shields.io/twitter/follow/espadrine.svg?style=social&label=Follow)](https://twitter.com/RoeyAngel)   ![license](https://img.shields.io/github/license/mashape/apistatus.svg?style=flat-square)


Analysis of stable-isotope probing experiment included in the paper: [Aeration and mineral composition of soil determine microbial CUE and OM retention] 


Overview
--------
    ├── DADA2_pseudo # ASV table, taxonomy phylogenetic tree and other data files
    ├── 00_cutadapt.log
    ├── 00_cutadapt_v2.0.sh
    ├── 01_DADA2_16S_merge_V8.6.R
    ├── 01_run_DADA2_16S_V8.6.sh
    ├── 02_Decontamination.html
    ├── 02_Decontamination.md
    ├── 02_Decontamination.Rmd
    ├── 03_Filter_taxa.html
    ├── 03_Filter_taxa.md
    ├── 03_Filter_taxa.Rmd
    ├── 04_calc_tree.log
    ├── 04_calc_tree_V2.0.sh
    ├── 05_calc_tree.log
    ├── 05_Diff_abund.html
    ├── 05_Diff_abund.md
    ├── 05_Diff_abund.Rmd
    ├── Anaerobic_CUE.Rproj # R project file
    ├── AnCUE_Metadata.csv # Metadata file
    ├── AnCUE_Metadata_decontam.csv
    ├── DESeq2_byTime_a-0.05_LFC0-322.txt # DESeq2 results table
    ├── fems-microbiology-ecology.csl
    ├── LICENSE # Copyright information
    ├── README.md # Overview of the repo
    └── references.bib  # Bibtex formatted refereces cited in the RMD file

Reproducing the analysis
--------
The RMD file is best executed using [knitr](https://yihui.name/knitr/) on [RStudio](https://www.rstudio.com/). 
