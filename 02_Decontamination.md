---
title: "Anaerobic CUE"
subtitle: "02 Decontaminate dataset"
author: "Roey Angel"
email: "roey.angel@bc.cas.cz"
date: "2021-02-01"
bibliography: references.bib
link-citations: yes
csl: fems-microbiology-ecology.csl
output:
  rmarkdown::html_document:
    toc: true
    toc_float: true
    keep_md: true
    number_sections: false
    highlight: "pygments"
    theme: "flatly"
    dev: "png"
    df_print: "kable"
    fig_caption: true
    code_folding: "show"
---






## Identify and remove contaminant ASVs
Decontamination of sequence library based on [Introduction to decontam](https://benjjneb.github.io/decontam/vignettes/decontam_intro.html) and Davis and colleagues [-@davis_simple_2018]. Decontamination is based on correlating sequence abundance frequencies to initial DNA concentrations used for PCR and also on examining sequence abundance prevalence patterns in the negative control samples.

### Setting general parameters:

```r
set.seed(1000)
samples_prep_path <- "./"
data_path <- "./DADA2_pseudo/"
Metadata_table <- "./AnCUE_Metadata.csv"
Seq_table <- "DADA2.seqtab_nochim.tsv"
Tax_table <- "DADA2.taxa_silva.tsv"
Seq_file <- "DADA2.Seqs.fa"
```

### Reading in raw data and generate phyloseq object

```r
# read OTU mat from data file
read_tsv(paste0(data_path, Seq_table), 
                        trim_ws = TRUE) %>% 
  column_to_rownames("ASV") %>% 
  t() %>% 
  as.data.frame() -> # not tibble because we need row names
  abundance_mat # convert to abundance matrix
```

```
## 
## ── Column specification ────────────────────────────────────────────────────────
## cols(
##   .default = col_double(),
##   ASV = col_character()
## )
## ℹ Use `spec()` for the full column specifications.
```

```r
# get short names of samples
# abundance_mat %>% 
#   rownames() %>% 
#   str_remove("^Roey[0-9]{3,}-?") %>% 
#   str_split("_", simplify = T) %>% 
#   .[, 1] ->
#   short_names

# Read metadata file
read_csv(paste0(samples_prep_path, Metadata_table), 
         trim_ws = TRUE) %>% 
  mutate(`16S copies` = replace(`16S copies`, which(`16S copies` == 0 | is.na(`16S copies`)), 1)) %>%  # add pseudo count
  filter(merged_sample_name %in% str_remove(rownames(abundance_mat), "_L001")) %>% # remove metadata rows if the samples did not go through qual processing
  mutate(`Library size` = rowSums(abundance_mat)) %>% # Add lib size
  mutate(to_names = merged_sample_name) %>% 
  mutate(across(c(
    "Site",
    "Oxygen",
    "Glucose",
    "Label (13C)"), 
    ~factor(.))) %>% 
  mutate(`Density zone` = factor(ifelse(`Density (g ml-1)` > 1.795, "Heavy", "Light"), levels = c("Light", "Heavy"))) %>% # critical for DESeq2 that the reference is the first level
  column_to_rownames("to_names") ->
  Metadata
```

```
## 
## ── Column specification ────────────────────────────────────────────────────────
## cols(
##   .default = col_double(),
##   merged_sample_name = col_character(),
##   Read1_file = col_character(),
##   Sample = col_character(),
##   Site = col_character(),
##   Oxygen = col_character(),
##   Glucose = col_character(),
##   `Label (13C)` = col_character(),
##   `TNA ext. Date` = col_character(),
##   `PCR date` = col_character(),
##   Control = col_logical()
## )
## ℹ Use `spec()` for the full column specifications.
```

```r
# Order abundance_mat samples according to the metadata
sample_order <- match(str_remove(rownames(abundance_mat), "_L001"), rownames(Metadata))
abundance_mat %<>% arrange(sample_order)
rownames(abundance_mat) <- rownames(Metadata) # needed for pyhloseq

# read taxonomy from data file
Raw_tax_data <- read_tsv(paste0(data_path, Tax_table), 
                        trim_ws = TRUE, col_names = TRUE)
```

```
## 
## ── Column specification ────────────────────────────────────────────────────────
## cols(
##   ASV = col_character(),
##   Kingdom = col_character(),
##   Phylum = col_character(),
##   Class = col_character(),
##   Order = col_character(),
##   Family = col_character(),
##   Genus = col_character(),
##   Species = col_character(),
##   `Kingdom (BS)` = col_double(),
##   `Phylum (BS)` = col_double(),
##   `Class (BS)` = col_double(),
##   `Order (BS)` = col_double(),
##   `Family (BS)` = col_double(),
##   `Genus (BS)` = col_double(),
##   `Species (BS)` = col_double()
## )
```

```r
Raw_tax_data %<>%
  mutate_all(~(replace(., is.na(.), "Unclassified"))) # I think mutaute_all is unnecessary here because replace(., is.na(.), "Unclassified") alone should work

Raw_tax_data %>%
  dplyr::select(.,
         `Kingdom (BS)`,
         `Phylum (BS)`,
         `Class (BS)`,
         `Order (BS)`,
         `Family (BS)`,
         `Genus (BS)`) %>%
  cbind(Name = colnames(abundance_mat),. ) ->
  Taxonomy.bs

Raw_tax_data %>%
  dplyr::select(.,
         Kingdom,
         Phylum,
         Class,
         Order,
         Family,
         Genus) %>% 
  # map_dfr(., as_factor) %>% 
  # map_dfr(fct_expand, "Rare")  %>%  
  as.matrix() -> # must be a matrix or phyloseq drops row names and gives and error
  Taxonomy
row.names(Taxonomy) <- colnames(abundance_mat)
# colnames(Taxonomy) <-
#   c("Domain", "Phylum", "Class", "Order", "Family", "Genus")

# read sequence data
ASV_seqs <- readDNAStringSet(
  file = paste0(data_path, Seq_file),
  format = "fasta", 
  nrec = -1L, 
  skip = 0L, 
  seek.first.rec = FALSE,
  use.names = TRUE)

# generate phyloseq object
Ps_obj <- phyloseq(otu_table(abundance_mat, taxa_are_rows = FALSE),
                   sample_data(Metadata),
                   tax_table(Taxonomy),
                   refseq(ASV_seqs))
```
### Inspect data structure

```r
Ps_obj %>% 
  get_variable() %>% 
  vis_dat()
```

```
## Warning: attributes are not identical across measure variables;
## they will be dropped
```

![](02_Decontamination_files/figure-html/data-1.png)<!-- -->

```r
Ps_obj %>% 
  get_variable() %>% 
  select_if(is.numeric) %>% 
  vis_cor()
```

![](02_Decontamination_files/figure-html/data-2.png)<!-- -->

### Inspect Library sizes

```r
Ps_obj_df <-
  as.data.frame(sample_data(Ps_obj)) # Put sample_data into a ggplot-friendly data.frame
Ps_obj_df <- Ps_obj_df[order(Ps_obj_df$`Library size`), ]
Ps_obj_df$Index <- seq(nrow(Ps_obj_df))
ggplot(data = Ps_obj_df, 
       aes(x = Index, y = Library.size, color = Control)) + 
  geom_point() +
  scale_y_log10(breaks = c(
    min(Ps_obj_df$Library.size),
    10,
    100,
    1000,
    5000,
    10000,
    ceiling(max(Ps_obj_df$Library.size) / 10000) * 10000
    )) + 
  scale_color_brewer(type = 'qual', palette = 'Set1', direction = -1)
```

```
## Warning: Transformation introduced infinite values in continuous y-axis
```

![](02_Decontamination_files/figure-html/Library Sizes-1.png)<!-- -->

```r
summary(sample_sums(Ps_obj))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##       0   19959   45598   38869   53196   76544
```

```r
summary(taxa_sums(Ps_obj))
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##      1.0      4.0     14.0    895.5    120.8 579695.0
```


```r
Ps_obj %<>%
  prune_samples(names(which(sample_sums(Ps_obj) != 0)), .)
# Ps_obj <-
#   subset_samples(Ps_obj, Sample != "CTRL")
# summary(sample_sums(Ps_obj))
```

###  Identify contaminants - Frequency
Use the distribution of the frequency of each sequence feature as a function of the input DNA concentration to identify contaminants.


```r
contamdf.freq <-
  isContaminant(Ps_obj, method = "frequency", conc = "X16S.copies")
# print(contamdf.freq)
# How many contaminants are found?
table(contamdf.freq$contaminant)
```

```
## 
## FALSE  TRUE 
## 14954    20
```

```r
# Which ones
which(contamdf.freq$contaminant)
```

```
##  [1]  2752  3799  4049  4326  5490  5951  6616  6965  7201  7972  8601 10101
## [13] 10188 10500 11115 11252 11352 12399 12419 12517
```

Plot the frequency of sequnce 1 and 3 (non-contaminants) against the DNA concentration, as an example.

```r
plot_frequency(Ps_obj, taxa_names(Ps_obj)[c(1, 3)], conc = "X16S.copies")
```

```
## Warning: Transformation introduced infinite values in continuous y-axis
```

```
## Warning: Removed 1 row(s) containing missing values (geom_path).
```

![](02_Decontamination_files/figure-html/plot frequency 1-1.png)<!-- -->

Plot the frequency of the contaminant sequences against the DNA concentration.

```r
plot_frequency(Ps_obj, taxa_names(Ps_obj)[which(contamdf.freq$contaminant)[1:20]], conc = "X16S.copies")
```

```
## Warning: Transformation introduced infinite values in continuous y-axis
```

```
## Warning: Removed 1 row(s) containing missing values (geom_path).
```

![](02_Decontamination_files/figure-html/plot frequency 2-1.png)<!-- -->

The frequency analysis detected $20$ sequences as contaminants.

###  Identify contaminants - Prevalence
Use the prevalence of sequences found in the control samples (no-template controls) to identify contaminants.

```r
contamdf.prev <- isContaminant(Ps_obj, method = "prevalence", neg = "Control")
# How many contaminants are found?
table(contamdf.prev$contaminant)
```

```
## 
## FALSE  TRUE 
## 14934    40
```

```r
# Which ones
which(contamdf.prev$contaminant)
```

```
##  [1]    35  2427  2983  3367  3783  4294  4469  5105  6815  6881  6903  7226
## [13]  7354  7385  7705  7728  7986  8444  8460  8508  8810  8831  8835  8869
## [25]  9031  9078  9151  9183  9240  9603  9881 10033 10089 10100 10522 10653
## [37] 11275 11291 11972 12261
```

```r
# And using a more aggressive threshold
contamdf.prev05 <- isContaminant(Ps_obj, method = "prevalence", neg = "Control", threshold = 0.5)
table(contamdf.prev05$contaminant)
```

```
## 
## FALSE  TRUE 
## 14753   221
```

```r
# Make phyloseq object of presence-absence in negative controls
Ps_obj.pa <-
  transform_sample_counts(Ps_obj, function(abund)
    1 * (abund > 0))
Ps_obj.pa.neg <-
  prune_samples(sample_data(Ps_obj.pa)$Control == "TRUE", Ps_obj.pa)
Ps_obj.pa.pos <-
  prune_samples(sample_data(Ps_obj.pa)$Control == "FALSE", Ps_obj.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <-
  data.frame(
    pa.pos = taxa_sums(Ps_obj.pa.pos),
    pa.neg = taxa_sums(Ps_obj.pa.neg),
    contaminant = contamdf.prev$contaminant
  )
ggplot(data = df.pa, aes(x = pa.neg, y = pa.pos, color = contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
```

![](02_Decontamination_files/figure-html/prevalence-1.png)<!-- -->

The frequency analysis detected $40$ sequences as contaminants.
In total $60$ were detected as contaminants and will be removed.

### Save contaminant sequence names and decontaminated data

```r
c(taxa_names(Ps_obj)[which(contamdf.freq$contaminant)],
  taxa_names(Ps_obj)[which(contamdf.prev$contaminant)]) ->
  contaminant_seqs
  
write_csv(as_tibble(contaminant_seqs), 
            paste0(data_path, "decontam_contaminants.csv"), 
            col_names = FALSE)


good_seqs <- setdiff(taxa_names(Ps_obj), contaminant_seqs)
Ps_obj_clean <- prune_taxa(good_seqs, Ps_obj)

# save decontaminated seqtab
Ps_obj_clean %>% 
  t() %>% 
  get_taxa() %>%
  as_tibble(rownames = "ASV") %>%
  write_tsv(., 
            paste0(data_path, str_remove(Seq_table, ".tsv"), "_decontam.tsv"), 
            col_names = TRUE)

Ps_obj_clean %>% 
  t() %>% 
  tax_table() %>%
  as_tibble(rownames = "ASV") %>%
  write_tsv(., 
            paste0(data_path, str_remove(Tax_table, ".tsv"), "_decontam.tsv"), 
            col_names = TRUE)

# save decontaminated metadata (just in case some samples were dropped)
Ps_obj_clean %>% 
  t() %>% 
  get_variable() %>% 
  setNames(., colnames(Metadata)) %>% 
  # as_tibble(rownames = "ASV") %>%
  write_csv(., 
            paste0("./", str_remove(Metadata_table, ".csv"), "_decontam.csv"), 
            col_names = TRUE)

# save decontaminated seqs
Ps_obj_clean %>% 
   refseq() %>% 
  writeXStringSet(., filepath = paste0(data_path, str_remove(Seq_file, ".fa*"), "_decontam.fa"), format = "fasta", width = 1000)
 
# save R obj
saveRDS(Ps_obj_clean, file = paste0(data_path, "Ps_obj_decontam.Rds"))
```


```r
sessioninfo::session_info() %>%
  details::details(
    summary = 'Current session info',
    open    = TRUE
  )
```

<details open>
<summary> <span title='Click to Expand'> Current session info </span> </summary>

```r

─ Session info ───────────────────────────────────────────────────────────────
 setting  value                       
 version  R version 4.0.3 (2020-10-10)
 os       Ubuntu 18.04.5 LTS          
 system   x86_64, linux-gnu           
 ui       X11                         
 language (EN)                        
 collate  en_US.UTF-8                 
 ctype    en_US.UTF-8                 
 tz       Europe/Prague               
 date     2021-02-01                  

─ Packages ───────────────────────────────────────────────────────────────────
 package      * version    date       lib source                          
 ade4           1.7-16     2020-10-28 [1] CRAN (R 4.0.2)                  
 ape            5.4-1      2020-08-13 [1] CRAN (R 4.0.2)                  
 assertthat     0.2.1      2019-03-21 [1] CRAN (R 4.0.2)                  
 backports      1.2.1      2020-12-09 [1] CRAN (R 4.0.2)                  
 Biobase        2.48.0     2020-04-27 [1] Bioconductor                    
 BiocGenerics * 0.34.0     2020-04-27 [1] Bioconductor                    
 biomformat     1.16.0     2020-04-27 [1] Bioconductor                    
 Biostrings   * 2.56.0     2020-04-27 [1] Bioconductor                    
 broom          0.7.4      2021-01-29 [1] CRAN (R 4.0.3)                  
 cellranger     1.1.0      2016-07-27 [1] CRAN (R 4.0.2)                  
 cli            2.3.0      2021-01-31 [1] CRAN (R 4.0.3)                  
 clipr          0.7.1      2020-10-08 [1] CRAN (R 4.0.2)                  
 cluster        2.1.0      2019-06-19 [1] CRAN (R 4.0.2)                  
 codetools      0.2-18     2020-11-04 [1] CRAN (R 4.0.2)                  
 colorspace     2.0-0      2020-11-11 [1] CRAN (R 4.0.2)                  
 crayon         1.4.0      2021-01-30 [1] CRAN (R 4.0.3)                  
 data.table     1.13.6     2020-12-30 [1] CRAN (R 4.0.2)                  
 DBI            1.1.1      2021-01-15 [1] CRAN (R 4.0.3)                  
 dbplyr         2.0.0      2020-11-03 [1] CRAN (R 4.0.2)                  
 decontam     * 1.8.0      2020-04-27 [1] Bioconductor                    
 desc           1.2.0      2018-05-01 [1] CRAN (R 4.0.2)                  
 details        0.2.1      2020-01-12 [1] CRAN (R 4.0.2)                  
 digest         0.6.27     2020-10-24 [1] CRAN (R 4.0.2)                  
 dplyr        * 1.0.3      2021-01-15 [1] CRAN (R 4.0.3)                  
 ellipsis       0.3.1      2020-05-15 [1] CRAN (R 4.0.2)                  
 evaluate       0.14       2019-05-28 [1] CRAN (R 4.0.2)                  
 extrafont    * 0.17       2014-12-08 [1] CRAN (R 4.0.2)                  
 extrafontdb    1.0        2012-06-11 [1] CRAN (R 4.0.2)                  
 farver         2.0.3      2020-01-16 [1] CRAN (R 4.0.2)                  
 forcats      * 0.5.1      2021-01-27 [1] CRAN (R 4.0.3)                  
 foreach        1.5.1      2020-10-15 [1] CRAN (R 4.0.2)                  
 fs             1.5.0      2020-07-31 [1] CRAN (R 4.0.2)                  
 gdtools        0.2.3      2021-01-06 [1] CRAN (R 4.0.2)                  
 generics       0.1.0      2020-10-31 [1] CRAN (R 4.0.2)                  
 ggplot2      * 3.3.3      2020-12-30 [1] CRAN (R 4.0.2)                  
 glue           1.4.2      2020-08-27 [1] CRAN (R 4.0.2)                  
 gtable         0.3.0      2019-03-25 [1] CRAN (R 4.0.2)                  
 haven          2.3.1      2020-06-01 [1] CRAN (R 4.0.2)                  
 highr          0.8        2019-03-20 [1] CRAN (R 4.0.2)                  
 hms            1.0.0      2021-01-13 [1] CRAN (R 4.0.3)                  
 htmltools      0.5.1.1    2021-01-22 [1] CRAN (R 4.0.3)                  
 httr           1.4.2      2020-07-20 [1] CRAN (R 4.0.2)                  
 igraph         1.2.6      2020-10-06 [1] CRAN (R 4.0.2)                  
 IRanges      * 2.22.2     2020-05-21 [1] Bioconductor                    
 iterators      1.0.13     2020-10-15 [1] CRAN (R 4.0.2)                  
 jsonlite       1.7.2      2020-12-09 [1] CRAN (R 4.0.2)                  
 knitr          1.31       2021-01-27 [1] CRAN (R 4.0.3)                  
 labeling       0.4.2      2020-10-20 [1] CRAN (R 4.0.2)                  
 lattice        0.20-41    2020-04-02 [1] CRAN (R 4.0.2)                  
 lifecycle      0.2.0      2020-03-06 [1] CRAN (R 4.0.2)                  
 lubridate      1.7.9.2    2020-11-13 [1] CRAN (R 4.0.2)                  
 magrittr     * 2.0.1      2020-11-17 [1] CRAN (R 4.0.2)                  
 MASS           7.3-53     2020-09-09 [1] CRAN (R 4.0.2)                  
 Matrix         1.3-2      2021-01-06 [1] CRAN (R 4.0.2)                  
 mgcv           1.8-33     2020-08-27 [1] CRAN (R 4.0.2)                  
 modelr         0.1.8      2020-05-19 [1] CRAN (R 4.0.2)                  
 multtest       2.44.0     2020-04-27 [1] Bioconductor                    
 munsell        0.5.0      2018-06-12 [1] CRAN (R 4.0.2)                  
 nlme           3.1-151    2020-12-10 [1] CRAN (R 4.0.2)                  
 permute        0.9-5      2019-03-12 [1] CRAN (R 4.0.2)                  
 phyloseq     * 1.32.0     2020-04-27 [1] Bioconductor                    
 pillar         1.4.7      2020-11-20 [1] CRAN (R 4.0.2)                  
 pkgconfig      2.0.3      2019-09-22 [1] CRAN (R 4.0.2)                  
 plyr           1.8.6      2020-03-03 [1] CRAN (R 4.0.2)                  
 png            0.1-7      2013-12-03 [1] CRAN (R 4.0.2)                  
 prettyunits    1.1.1      2020-01-24 [1] CRAN (R 4.0.2)                  
 progress       1.2.2      2019-05-16 [1] CRAN (R 4.0.2)                  
 purrr        * 0.3.4      2020-04-17 [1] CRAN (R 4.0.2)                  
 R6             2.5.0      2020-10-28 [1] CRAN (R 4.0.2)                  
 RColorBrewer   1.1-2      2014-12-07 [1] CRAN (R 4.0.2)                  
 Rcpp           1.0.6      2021-01-15 [1] CRAN (R 4.0.3)                  
 readr        * 1.4.0      2020-10-05 [1] CRAN (R 4.0.2)                  
 readxl         1.3.1      2019-03-13 [1] CRAN (R 4.0.2)                  
 reprex         1.0.0      2021-01-27 [1] CRAN (R 4.0.3)                  
 reshape2       1.4.4      2020-04-09 [1] CRAN (R 4.0.2)                  
 rhdf5          2.32.4     2020-10-05 [1] Bioconductor                    
 Rhdf5lib       1.10.1     2020-07-09 [1] Bioconductor                    
 rlang          0.4.10     2020-12-30 [1] CRAN (R 4.0.2)                  
 rmarkdown      2.6        2020-12-14 [1] CRAN (R 4.0.2)                  
 rprojroot      2.0.2      2020-11-15 [1] CRAN (R 4.0.2)                  
 rstudioapi     0.13       2020-11-12 [1] CRAN (R 4.0.2)                  
 Rttf2pt1       1.3.8      2020-01-10 [1] CRAN (R 4.0.2)                  
 rvest          0.3.6      2020-07-25 [1] CRAN (R 4.0.2)                  
 S4Vectors    * 0.26.1     2020-05-16 [1] Bioconductor                    
 scales         1.1.1      2020-05-11 [1] CRAN (R 4.0.2)                  
 sessioninfo    1.1.1      2018-11-05 [1] CRAN (R 4.0.2)                  
 stringi        1.5.3      2020-09-09 [1] CRAN (R 4.0.2)                  
 stringr      * 1.4.0      2019-02-10 [1] CRAN (R 4.0.2)                  
 survival       3.2-7      2020-09-28 [1] CRAN (R 4.0.2)                  
 svglite      * 1.2.3.2    2020-07-07 [1] CRAN (R 4.0.2)                  
 systemfonts    0.3.2      2020-09-29 [1] CRAN (R 4.0.2)                  
 tibble       * 3.0.6      2021-01-29 [1] CRAN (R 4.0.3)                  
 tidyr        * 1.1.2      2020-08-27 [1] CRAN (R 4.0.2)                  
 tidyselect     1.1.0      2020-05-11 [1] CRAN (R 4.0.2)                  
 tidyverse    * 1.3.0      2019-11-21 [1] CRAN (R 4.0.2)                  
 vctrs          0.3.6      2020-12-17 [1] CRAN (R 4.0.2)                  
 vegan          2.5-7      2020-11-28 [1] CRAN (R 4.0.2)                  
 visdat       * 0.6.0.9000 2021-02-01 [1] Github (ropensci/visdat@8121dfe)
 withr          2.4.1      2021-01-26 [1] CRAN (R 4.0.3)                  
 xfun           0.20       2021-01-06 [1] CRAN (R 4.0.2)                  
 xml2           1.3.2      2020-04-23 [1] CRAN (R 4.0.2)                  
 XVector      * 0.28.0     2020-04-27 [1] Bioconductor                    
 yaml           2.2.1      2020-02-01 [1] CRAN (R 4.0.2)                  
 zlibbioc       1.34.0     2020-04-27 [1] Bioconductor                    

[1] /home/angel/R/library
[2] /usr/local/lib/R/site-library
[3] /usr/lib/R/site-library
[4] /usr/lib/R/library

```

</details>
<br>

## References
