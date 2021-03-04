Aeration and mineral composition of soil determine microbial CUE and OM retention
========

[![Twitter Follow](https://img.shields.io/twitter/follow/espadrine.svg?style=social&label=Follow)](https://twitter.com/RoeyAngel)   ![license](https://img.shields.io/github/license/mashape/apistatus.svg?style=flat-square)


Analysis of stable-isotope probing experiment included in the paper: [Aeration and mineral composition of soil determine microbial CUE and OM retention] 


Overview
--------
    ├── DADA2_pseudo # ASV table, taxonomy phylogenetic tree and other data files
        ├── 01_run_DADA2_16S_pseudo_V8.6.log
        ├── DADA2.Seqs_decontam.fa
        ├── DADA2.Seqs_decontam_filtered.fa
        ├── DADA2.Seqs.fa
        ├── DADA2.seqtab_nochim_decontam.tsv
        ├── DADA2.seqtab_nochim.RDS
        ├── DADA2.seqtab_nochim.tsv
        ├── DADA2.taxa_silva_decontam.tsv
        ├── DADA2.taxa_silva.tsv
        ├── DADA2.track_each.tsv
        ├── decontam_contaminants.csv
        ├── Ps_obj_decontam_filt3.Rds
        ├── Ps_obj_decontam.Rds
        ├── QualProf.pdf
        ├── Read1Errors.pdf
        ├── Read2Errors.pdf
        └── Tree
            ├── DADA2.filter
            ├── DADA2.Seqs_decontam_filtered.filtered.align
            ├── DADA2.Seqs_decontam_filtered.filtered.align.bionj
            ├── DADA2.Seqs_decontam_filtered.filtered.align.iqtree
            ├── DADA2.Seqs_decontam_filtered.filtered.align.log
            ├── DADA2.Seqs_decontam_filtered.filtered.align.treefile
            └── DADA2.Seqs_decontam_filtered.FT.tree
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

Viewing and reproducing the analysis
--------
The MD files can be used to display the output on github, while HTML files can be used for offline viewing. 
The RMD files are best compiled using [knitr](https://yihui.name/knitr/) on [RStudio](https://www.rstudio.com/). 
