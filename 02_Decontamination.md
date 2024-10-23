Anaerobic CUE
================
Roey Angel
2024-06-27

- <a href="#identify-and-remove-contaminant-asvs"
  id="toc-identify-and-remove-contaminant-asvs">Identify and remove
  contaminant ASVs</a>
  - <a href="#setting-general-parameters"
    id="toc-setting-general-parameters">Setting general parameters:</a>
  - <a href="#reading-in-raw-data-and-generate-phyloseq-object"
    id="toc-reading-in-raw-data-and-generate-phyloseq-object">Reading in raw
    data and generate phyloseq object</a>
  - <a href="#inspect-data-structure"
    id="toc-inspect-data-structure">Inspect data structure</a>
  - <a href="#inspect-library-sizes" id="toc-inspect-library-sizes">Inspect
    Library sizes</a>
  - <a href="#identify-contaminants---frequency"
    id="toc-identify-contaminants---frequency">Identify contaminants -
    Frequency</a>
  - <a href="#identify-contaminants---prevalence"
    id="toc-identify-contaminants---prevalence">Identify contaminants -
    Prevalence</a>
  - <a href="#save-contaminant-sequence-names-and-decontaminated-data"
    id="toc-save-contaminant-sequence-names-and-decontaminated-data">Save
    contaminant sequence names and decontaminated data</a>
- <a href="#references" id="toc-references">References</a>

## Identify and remove contaminant ASVs

Decontamination of sequence library based on [Introduction to
decontam](https://benjjneb.github.io/decontam/vignettes/decontam_intro.html)
and Davis and colleagues ([2018](#ref-davis_simple_2018)).
Decontamination is based on correlating sequence abundance frequencies
to initial DNA concentrations used for PCR and also on examining
sequence abundance prevalence patterns in the negative control samples.

### Setting general parameters:

``` r
set.seed(1000)
samples_prep_path <- "./"
data_path <- "./DADA2_pseudo/"
Metadata_table <- "./AnCUE_Metadata.csv"
Seq_table <- "DADA2.seqtab_nochim.tsv"
Tax_table <- "DADA2.taxa_silva.tsv"
Seq_file <- "DADA2.Seqs.fa"
```

### Reading in raw data and generate phyloseq object

``` r
# read OTU mat from data file
read_tsv(paste0(data_path, Seq_table), 
                        trim_ws = TRUE) %>% 
  column_to_rownames("ASV") %>% 
  t() %>% 
  as.data.frame() -> # not tibble because we need row names
  abundance_mat # convert to abundance matrix

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

# Order abundance_mat samples according to the metadata
sample_order <- match(str_remove(rownames(abundance_mat), "_L001"), rownames(Metadata))
abundance_mat %<>% arrange(sample_order)
rownames(abundance_mat) <- rownames(Metadata) # needed for pyhloseq

# read taxonomy from data file
Raw_tax_data <- read_tsv(paste0(data_path, Tax_table), 
                        trim_ws = TRUE, col_names = TRUE)
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

``` r
Ps_obj %>% 
  get_variable() %>% 
  vis_dat()
```

![](02_Decontamination_figures/data-1.png)<!-- -->

``` r
Ps_obj %>% 
  get_variable() %>% 
  select_if(is.numeric) %>% 
  vis_cor()
```

![](02_Decontamination_figures/data-2.png)<!-- -->

### Inspect Library sizes

``` r
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

![](02_Decontamination_figures/Library%20Sizes-1.png)<!-- -->

``` r
summary(sample_sums(Ps_obj))
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##       0   19959   45598   38869   53196   76544

``` r
summary(taxa_sums(Ps_obj))
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ##      1.0      4.0     14.0    895.5    120.8 579695.0

``` r
Ps_obj %<>%
  prune_samples(names(which(sample_sums(Ps_obj) != 0)), .)
# Ps_obj <-
#   subset_samples(Ps_obj, Sample != "CTRL")
# summary(sample_sums(Ps_obj))
```

### Identify contaminants - Frequency

Use the distribution of the frequency of each sequence feature as a
function of the input DNA concentration to identify contaminants.

``` r
contamdf.freq <-
  isContaminant(Ps_obj, method = "frequency", conc = "X16S.copies")
# print(contamdf.freq)
# How many contaminants are found?
table(contamdf.freq$contaminant)
```

    ## 
    ## FALSE  TRUE 
    ## 14954    20

``` r
# Which ones
which(contamdf.freq$contaminant)
```

    ##  [1]  2752  3799  4049  4326  5490  5951  6616  6965  7201  7972  8601 10101 10188 10500
    ## [15] 11115 11252 11352 12399 12419 12517

Plot the frequency of sequnce 1 and 3 (non-contaminants) against the DNA
concentration, as an example.

``` r
plot_frequency(Ps_obj, taxa_names(Ps_obj)[c(1, 3)], conc = "X16S.copies")
```

![](02_Decontamination_figures/plot%20frequency%201-1.png)<!-- -->

Plot the frequency of the contaminant sequences against the DNA
concentration.

``` r
plot_frequency(Ps_obj, taxa_names(Ps_obj)[which(contamdf.freq$contaminant)[1:20]], conc = "X16S.copies")
```

![](02_Decontamination_figures/plot%20frequency%202-1.png)<!-- -->

The frequency analysis detected $20$ sequences as contaminants.

### Identify contaminants - Prevalence

Use the prevalence of sequences found in the control samples
(no-template controls) to identify contaminants.

``` r
contamdf.prev <- isContaminant(Ps_obj, method = "prevalence", neg = "Control")
# How many contaminants are found?
table(contamdf.prev$contaminant)
```

    ## 
    ## FALSE  TRUE 
    ## 14934    40

``` r
# Which ones
which(contamdf.prev$contaminant)
```

    ##  [1]    35  2427  2983  3367  3783  4294  4469  5105  6815  6881  6903  7226  7354  7385
    ## [15]  7705  7728  7986  8444  8460  8508  8810  8831  8835  8869  9031  9078  9151  9183
    ## [29]  9240  9603  9881 10033 10089 10100 10522 10653 11275 11291 11972 12261

``` r
# And using a more aggressive threshold
contamdf.prev05 <- isContaminant(Ps_obj, method = "prevalence", neg = "Control", threshold = 0.5)
table(contamdf.prev05$contaminant)
```

    ## 
    ## FALSE  TRUE 
    ## 14753   221

``` r
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

![](02_Decontamination_figures/prevalence-1.png)<!-- -->

The frequency analysis detected $40$ sequences as contaminants. In total
$60$ were detected as contaminants and will be removed.

### Save contaminant sequence names and decontaminated data

``` r
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

``` r
sessioninfo::session_info() %>%
  details::details(
    summary = 'Current session info',
    open    = TRUE
  )
```

<details open>
<summary>
<span title="Click to Expand"> Current session info </span>
</summary>

``` r

─ Session info ─────────────────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.4.0 (2024-04-24)
 os       Ubuntu 22.04.4 LTS
 system   x86_64, linux-gnu
 ui       X11
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Europe/Prague
 date     2024-06-27
 pandoc   2.19.2 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/ (via rmarkdown)

─ Packages ─────────────────────────────────────────────────────────────────────────────
 package          * version    date (UTC) lib source
 ade4               1.7-22     2023-02-06 [1] CRAN (R 4.2.2)
 ape                5.8        2024-04-11 [1] CRAN (R 4.4.0)
 Biobase            2.60.0     2023-04-25 [1] Bioconductor
 BiocGenerics     * 0.46.0     2023-04-25 [1] Bioconductor
 biomformat         1.28.0     2023-04-25 [1] Bioconductor
 Biostrings       * 2.68.1     2023-05-16 [1] Bioconductor
 bit                4.0.5      2022-11-15 [1] CRAN (R 4.2.2)
 bit64              4.0.5      2020-08-30 [1] CRAN (R 4.0.2)
 bitops             1.0-7      2021-04-24 [1] CRAN (R 4.0.3)
 cli                3.6.2      2023-12-11 [1] CRAN (R 4.3.2)
 clipr              0.8.0      2022-02-22 [1] CRAN (R 4.1.2)
 cluster            2.1.6      2023-12-01 [1] CRAN (R 4.3.2)
 codetools          0.2-20     2024-03-31 [1] CRAN (R 4.3.3)
 colorspace         2.1-0      2023-01-23 [1] CRAN (R 4.2.2)
 crayon             1.5.2      2022-09-29 [1] CRAN (R 4.2.1)
 data.table         1.15.4     2024-03-30 [1] CRAN (R 4.3.3)
 decontam         * 1.20.0     2023-04-25 [1] Bioconductor
 desc               1.4.3      2023-12-10 [1] CRAN (R 4.3.2)
 details            0.3.0      2022-03-27 [1] CRAN (R 4.1.3)
 digest             0.6.35     2024-03-11 [1] CRAN (R 4.3.3)
 dplyr            * 1.1.4      2023-11-17 [1] CRAN (R 4.3.2)
 evaluate           0.24.0     2024-06-10 [1] CRAN (R 4.4.0)
 extrafont        * 0.19       2023-01-18 [1] CRAN (R 4.2.2)
 extrafontdb        1.0        2012-06-11 [1] CRAN (R 4.0.2)
 fansi              1.0.6      2023-12-08 [1] CRAN (R 4.3.2)
 farver             2.1.2      2024-05-13 [1] CRAN (R 4.4.0)
 fastmap            1.2.0      2024-05-15 [1] CRAN (R 4.4.0)
 forcats          * 1.0.0      2023-01-29 [1] CRAN (R 4.2.2)
 foreach            1.5.2      2022-02-02 [1] CRAN (R 4.1.2)
 generics           0.1.3      2022-07-05 [1] CRAN (R 4.2.0)
 GenomeInfoDb     * 1.36.2     2023-08-25 [1] Bioconductor
 GenomeInfoDbData   1.2.10     2023-05-31 [1] Bioconductor
 ggplot2          * 3.5.1      2024-04-23 [1] CRAN (R 4.4.0)
 glue               1.7.0      2024-01-09 [1] CRAN (R 4.3.2)
 gtable             0.3.5      2024-04-22 [1] CRAN (R 4.4.0)
 highr              0.11       2024-05-26 [1] CRAN (R 4.4.0)
 hms                1.1.3      2023-03-21 [1] CRAN (R 4.2.3)
 htmltools          0.5.8.1    2024-04-04 [1] CRAN (R 4.4.0)
 httr               1.4.7      2023-08-15 [1] CRAN (R 4.3.1)
 igraph             2.0.3      2024-03-13 [1] CRAN (R 4.3.3)
 IRanges          * 2.34.1     2023-06-22 [1] Bioconductor
 iterators          1.0.14     2022-02-05 [1] CRAN (R 4.1.2)
 jsonlite           1.8.8      2023-12-04 [1] CRAN (R 4.3.2)
 knitr              1.47       2024-05-29 [1] CRAN (R 4.4.0)
 labeling           0.4.3      2023-08-29 [1] CRAN (R 4.3.1)
 lattice            0.22-6     2024-03-20 [1] CRAN (R 4.3.3)
 lifecycle          1.0.4      2023-11-07 [1] CRAN (R 4.3.1)
 lubridate        * 1.9.3      2023-09-27 [1] CRAN (R 4.3.1)
 magrittr         * 2.0.3      2022-03-30 [1] CRAN (R 4.1.3)
 MASS               7.3-61     2024-06-13 [1] CRAN (R 4.4.0)
 Matrix             1.7-0      2024-04-26 [1] CRAN (R 4.4.0)
 mgcv               1.9-1      2023-12-21 [1] CRAN (R 4.3.2)
 multtest           2.56.0     2023-04-25 [1] Bioconductor
 munsell            0.5.1      2024-04-01 [1] CRAN (R 4.4.0)
 nlme               3.1-165    2024-06-06 [1] CRAN (R 4.4.0)
 permute            0.9-7      2022-01-27 [1] CRAN (R 4.1.2)
 phyloseq         * 1.44.0     2023-04-25 [1] Bioconductor
 pillar             1.9.0      2023-03-22 [1] CRAN (R 4.2.3)
 pkgconfig          2.0.3      2019-09-22 [1] CRAN (R 4.0.2)
 plyr               1.8.9      2023-10-02 [1] CRAN (R 4.3.1)
 png                0.1-8      2022-11-29 [1] CRAN (R 4.2.2)
 purrr            * 1.0.2      2023-08-10 [1] CRAN (R 4.3.1)
 R6                 2.5.1      2021-08-19 [1] CRAN (R 4.1.1)
 RColorBrewer       1.1-3      2022-04-03 [1] CRAN (R 4.1.3)
 Rcpp               1.0.12     2024-01-09 [1] CRAN (R 4.3.2)
 RCurl              1.98-1.14  2024-01-09 [1] CRAN (R 4.3.2)
 readr            * 2.1.5      2024-01-10 [1] CRAN (R 4.3.2)
 reshape2           1.4.4      2020-04-09 [1] CRAN (R 4.0.2)
 rhdf5              2.44.0     2023-04-25 [1] Bioconductor
 rhdf5filters       1.12.1     2023-04-30 [1] Bioconductor
 Rhdf5lib           1.22.0     2023-04-25 [1] Bioconductor
 rlang              1.1.4      2024-06-04 [1] CRAN (R 4.4.0)
 rmarkdown          2.27       2024-05-17 [1] CRAN (R 4.4.0)
 rstudioapi         0.16.0     2024-03-24 [1] CRAN (R 4.3.3)
 Rttf2pt1           1.3.12     2023-01-22 [1] CRAN (R 4.2.2)
 S4Vectors        * 0.38.1     2023-05-02 [1] Bioconductor
 scales             1.3.0      2023-11-28 [1] CRAN (R 4.3.2)
 sessioninfo        1.2.2      2021-12-06 [1] CRAN (R 4.1.2)
 stringi            1.8.4      2024-05-06 [1] CRAN (R 4.4.0)
 stringr          * 1.5.1      2023-11-14 [1] CRAN (R 4.3.2)
 survival           3.7-0      2024-06-05 [1] CRAN (R 4.4.0)
 svglite          * 2.1.3      2023-12-08 [1] CRAN (R 4.3.2)
 systemfonts        1.1.0      2024-05-15 [1] CRAN (R 4.4.0)
 tibble           * 3.2.1      2023-03-20 [1] CRAN (R 4.2.3)
 tidyr            * 1.3.1      2024-01-24 [1] CRAN (R 4.3.2)
 tidyselect         1.2.1      2024-03-11 [1] CRAN (R 4.3.3)
 tidyverse        * 2.0.0      2023-02-22 [1] CRAN (R 4.2.3)
 timechange         0.3.0      2024-01-18 [1] CRAN (R 4.3.2)
 tzdb               0.4.0      2023-05-12 [1] CRAN (R 4.2.3)
 utf8               1.2.4      2023-10-22 [1] CRAN (R 4.3.1)
 vctrs              0.6.5      2023-12-01 [1] CRAN (R 4.3.2)
 vegan              2.6-6.1    2024-05-21 [1] CRAN (R 4.4.0)
 visdat           * 0.6.0.9000 2023-05-31 [1] Github (ropensci/visdat@9199906)
 vroom              1.6.5      2023-12-05 [1] CRAN (R 4.3.2)
 withr              3.0.0      2024-01-16 [1] CRAN (R 4.3.2)
 xfun               0.45       2024-06-16 [1] CRAN (R 4.4.0)
 xml2               1.3.6      2023-12-04 [1] CRAN (R 4.3.2)
 XVector          * 0.40.0     2023-04-25 [1] Bioconductor
 yaml               2.3.8      2023-12-11 [1] CRAN (R 4.3.2)
 zlibbioc           1.46.0     2023-04-25 [1] Bioconductor

 [1] /home/angel/R/library
 [2] /usr/local/lib/R/site-library
 [3] /usr/lib/R/site-library
 [4] /usr/lib/R/library

────────────────────────────────────────────────────────────────────────────────────────
```

</details>

<br>

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-davis_simple_2018" class="csl-entry">

Davis NM, Proctor DM, Holmes SP *et al.* [Simple statistical
identification and removal of contaminant sequences in marker-gene and
metagenomics data](https://doi.org/gfxzx5). *Microbiome* 2018;**6**:226.

</div>

</div>
