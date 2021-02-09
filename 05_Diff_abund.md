---
title: "Anaerobic CUE"
subtitle: "05 Differential abundance modelling"
author: "Roey Angel"
email: "roey.angel@bc.cas.cz"
date: "2021-02-09"
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







## Differential abundance modelling of SIP gradients
Here we attempt to detect ASVs that were labelled with 13C in our soil incubations using differential abundance modelling.
Using DESeq2 \citep(love_moderated_2014) we compare the relative abundance of each ASV in the fractions where ^^13^^C label is supposed to be found ()

### Setting general parameters:

```r
set.seed(2021)
alpha_thresh <- 0.05
LFC_thresh <- 0.26
samples_prep_path <- "./"
data_path <- "./DADA2_pseudo/"
# Metadata_table <- "./AnCUE_Metadata_decontam.csv"
# Seq_table <- "DADA2.seqtab_nochim_decontam.tsv"
# Seq_file <- "DADA2.Seqs_decontam.fa"
Ps_file <- "Ps_obj_decontam_filt3.Rds"
Tree_file <- "./Tree/DADA2.Seqs_decontam_filtered.filtered.align.treefile"
```

### Read phyloseq object

```r
# Load phylogenetic tree
Tree <- read_tree(paste0(data_path, Tree_file))

# load and merge  phyloseq object
readRDS(paste0(data_path, Ps_file)) %>% 
  merge_phyloseq(.,
                 phy_tree(Tree)
  ) -> Ps_obj_SIP
```

### Beta diversity analysis
Let us look first at 

```r
Ord <- ordinate(Ps_obj_SIP, "CAP", "horn", formula =  ~ Oxygen + Density.zone)
explained <- as.numeric(format(round(eigenvals(Ord)/sum(eigenvals(Ord)) * 100, 1), nsmall = 1))
Ord_plt <- plot_ordination(Ps_obj_SIP, Ord, type = "samples", color = "Label..13C.", justDF = TRUE)

p_ord_joint <- ggplot(Ord_plt) +
  geom_point(aes(
               x = CAP1,
               y = CAP2,
               color = Label..13C.,
               size = Density..g.ml.1.,
               shape = Oxygen
             ), alpha = 2 / 3) +
  guides(colour = guide_legend(title = "Labelling"), 
         size = guide_legend(title = "Density (g ml<sup>-1</sup>)"),
         shape = guide_legend(title = "Oxygen")) +
  ggsci::scale_colour_locuszoom() +
  # scale_colour_manual(values = Gradient.colours) +
  # scale_fill_manual(values = Gradient.colours, guide = "none") +
  labs(x = sprintf("CAP1 (%s%%)", explained[1]),
  y = sprintf("CAP2 (%s%%)", explained[2])) +
  coord_fixed(ratio = sqrt(explained[2] / explained[1])) +
   theme(legend.justification = "top",
         legend.title = element_markdown(size = 11)
         ) +
  scale_size_continuous(breaks = c(seq(min(Ord_plt$Density..g.ml.1.), 
                                       max(Ord_plt$Density..g.ml.1.), 
                                       length.out = 5), 
                                   1),
                        range = c(0.1, 5)) +
  facet_grid(Site ~ Hours) +
  ggtitle("Joint analysis")
print(p_ord_joint)
```

![](05_Diff_abund_figures/beta div ord joint-1.png)<!-- -->


```r
test_expr_1 <- "(Site == '${Site}' & Oxygen == '${Oxygen}' & Label..13C. == 'Unlabelled') | (Site == '${Site}'  & Oxygen == '${Oxygen}' & Label..13C. == '${Label..13C.}')"
params_1 <- get_treatment_params(Ps_obj_SIP, c("Site",
                                   "Oxygen",
                                   "Glucose",
                                   "Label..13C."),
                     "Label..13C. != 'Unlabelled'")
Ps_obj_SIP_noTime_l <- phyloseq_subset(Ps_obj_SIP, params_1, test_expr_1) 
# names(Ps_obj_SIP_noTime_l) %<>% 
#   map(., ~str_remove_all(.x, "\\s\\|\\s.*")) %>% 
#   map(., ~str_remove_all(.x, "\\(|\\)|Site == |Hours == |Oxygen == |Label..13C. == |'")) %>% 
#   map(., ~str_replace_all(.x, "([0-9]+)", "\\1 h")) 

test_expr_2 <- "(Site == '${Site}' & Hours == '${Hours}' & Oxygen == '${Oxygen}' & Label..13C. == '${Label..13C.}') | (Site == '${Site}' & Hours == '${Hours}' & Oxygen == '${Oxygen}' & Label..13C. == '${Label..13C.}')"
params_2 <- get_treatment_params(Ps_obj_SIP, c("Site", 
                                   "Hours", 
                                   "Oxygen",
                                   "Glucose",
                                   "Label..13C."))

# Generate a list of subsetted phyloseq objects
Ps_obj_SIP_byTime_l <- phyloseq_subset(Ps_obj_SIP, params_2, test_expr_2) 
names(Ps_obj_SIP_byTime_l) %<>% 
  map(., ~str_remove_all(.x, "\\s\\|\\s.*")) %>% 
  map(., ~str_remove_all(.x, "\\(|\\)|Site == |Hours == |Oxygen == |Label..13C. == |'")) %>% 
  map(., ~str_replace_all(.x, "([0-9]+)", "\\1 h")) 

test_expr_3 <- "(Site == '${Site}' & Oxygen == '${Oxygen}' & Label..13C. == '${Label..13C.}') | (Site == '${Site}' & Oxygen == '${Oxygen}' & Label..13C. == '${Label..13C.}')"
params_3 <- get_treatment_params(Ps_obj_SIP, c("Site", 
                                   "Oxygen",
                                   "Glucose",
                                   "Label..13C."))

# Generate a list of subsetted phyloseq objects
Ps_obj_SIP_noTime_split_l <- phyloseq_subset(Ps_obj_SIP, params_3, test_expr_3) 
names(Ps_obj_SIP_noTime_split_l) %<>% 
  map(., ~str_remove_all(.x, "\\s\\|\\s.*")) %>% 
  map(., ~str_remove_all(.x, "\\(|\\)|Site == |Hours == |Oxygen == |Label..13C. == |'")) %>% 
  map(., ~str_replace_all(.x, "([0-9]+)", "\\1 h")) 
```


```r
Ord_l <- mclapply(Ps_obj_SIP_byTime_l, function(x) {ordinate(x, "CAP", "horn", formula = ~ Oxygen + Density.zone)}, mc.cores = nrow(params_2))

explained <- mclapply(Ord_l, function(x) {as.numeric(format(round(eigenvals(x)/sum( eigenvals(x)) * 100, 1), nsmall = 1))}, mc.cores = nrow(params_2))

bind_rows(map(seq(1, length(Ps_obj_SIP_byTime_l)), ~plot_ordination(Ps_obj_SIP_byTime_l[[.x]], Ord_l[[.x]], type = "samples", color = "Label..13C.", justDF = TRUE))) ->
  # mutate(MDS1 = coalesce(MDS1, MDS2)) ->
  Ord_plt_time_df 

p_ord_split_time <- ggplot(Ord_plt_time_df) +
  geom_point(aes(
               x = CAP1,
               y = MDS1,
               color = Label..13C.,
               size = Density..g.ml.1.,
               shape = Oxygen
             ), alpha = 2 / 3) +
  # guides(colour = guide_legend(title = "Labelling"), shape = guide_legend(title = "Oxygen")) +
  ggsci::scale_colour_locuszoom() +
  # labs(x = sprintf("CAP1 (%s%%)", explained[1]), 
       # y = sprintf("CAP2 (%s%%)", explained[2])) +
  # coord_fixed(ratio = sqrt(explained[2] / explained[1])) +
   theme(legend.justification = "top") +
  facet_grid(Site ~ Hours) +
  ggtitle("Split analysis (incl. time)")
print(p_ord_split_time)
```

![](05_Diff_abund_figures/beta div ord split by time-1.png)<!-- -->


```r
Ord_l <- mclapply(Ps_obj_SIP_noTime_l, function(x) {ordinate(x, "CAP", "horn", formula = ~ Oxygen + Density.zone)}, mc.cores = nrow(params_1))

explained <- mclapply(Ord_l, function(x) {as.numeric(format(round(eigenvals(x)/sum( eigenvals(x)) * 100, 1), nsmall = 1))}, mc.cores = nrow(params_1))

bind_rows(map(seq(1, length(Ps_obj_SIP_noTime_l)), ~plot_ordination(Ps_obj_SIP_noTime_l[[.x]], Ord_l[[.x]], type = "samples", color = "Label..13C.", justDF = TRUE))) ->
  # mutate(MDS1 = coalesce(MDS1, MDS2)) ->
  Ord_plt_df 

p_ord_split_time <- ggplot(Ord_plt_time_df) +
  geom_point(aes(
               x = CAP1,
               y = MDS1,
               color = Label..13C.,
               size = Density..g.ml.1.,
               shape = Oxygen
             ), alpha = 2 / 3) +
  # guides(colour = guide_legend(title = "Labelling"), shape = guide_legend(title = "Oxygen")) +
  ggsci::scale_colour_locuszoom() +
  # labs(x = sprintf("CAP1 (%s%%)", explained[1]), 
       # y = sprintf("CAP2 (%s%%)", explained[2])) +
  # coord_fixed(ratio = sqrt(explained[2] / explained[1])) +
   theme(legend.justification = "top") +
  facet_grid(Site ~ Oxygen) +
  ggtitle("Split analysis (excl. time)")
print(p_ord_split_time)
```

![](05_Diff_abund_figures/beta div ord split no time-1.png)<!-- -->


```r
# generate a deseq object
DESeq_obj_SIP_byTime_l <- mclapply(Ps_obj_SIP_byTime_l, 
                                   function(x) {phyloseq_to_deseq2_safe(x, 
                                                                        test_condition = "Density.zone", 
                                                                        ref_level = "Light")}, 
                                   mc.cores = nrow(params_2))


# run dds pipeline
DESeq_obj_SIP_byTime_l %<>% mclapply(., 
                                     function(x) {DESeq(x, 
                                                        test = "Wald", 
                                                        fitType = "local")}, 
                                     mc.cores = nrow(params_2)) # run dds pipeline

# extract results from a DESeq analysis
DESeq_res_SIP_byTime_l  <- mclapply(DESeq_obj_SIP_byTime_l, 
                             function(x) {
                               results(x, 
                                       altHypothesis = "greater",
                                       alpha = alpha_thresh, 
                                       contrast = c("Density.zone", "Heavy", "Light"))}, # redundant if phyloseq_to_deseq2_safe() was used but doesn't hurt
                             mc.cores = nrow(params_2)) 

DESeq_res_SIP_byTime_LFC_l <- mclapply(DESeq_obj_SIP_byTime_l, 
                                     function(x) {
                                       results(x,
                                               lfcThreshold = LFC_thresh,
                                               altHypothesis = "greater",
                                               alpha = alpha_thresh,
                                               contrast = c("Density.zone", "Heavy", "Light"))}, # redundant if phyloseq_to_deseq2_safe() was used but doesn't hurt
                                     mc.cores = nrow(params_2)) # Extract results from a DESeq analysis


DESeq_res_SIP_byTime_LFC_shrink_l <- map(seq(length(DESeq_obj_SIP_byTime_l)), 
                                             ~lfcShrink(DESeq_obj_SIP_byTime_l[[.x]],
                                                         res = DESeq_res_SIP_byTime_LFC_l[[.x]],
                                                         coef = "Density.zone_Heavy_vs_Light",
                                                         type = "ashr"))
names(DESeq_res_SIP_byTime_LFC_shrink_l) <- names(DESeq_res_SIP_byTime_LFC_l)

# Compare
plotMA(DESeq_res_SIP_byTime_l[[2]], ylim = c(-2,2))
```

![](05_Diff_abund_figures/DESeq2 models by time-1.png)<!-- -->

```r
plotMA(DESeq_res_SIP_byTime_LFC_l[[2]], ylim = c(-2,2))
```

![](05_Diff_abund_figures/DESeq2 models by time-2.png)<!-- -->

```r
plotMA(DESeq_res_SIP_byTime_LFC_shrink_l[[2]], ylim = c(-2,2))
```

![](05_Diff_abund_figures/DESeq2 models by time-3.png)<!-- -->

```r
# summarise results (lfcShrink doesn't change the values)
# map2(DESeq_res_SIP_byTime_l, print(names(DESeq_res_SIP_byTime_l)), ~summary(.x)) # summarise results
for (i in seq(1, length(DESeq_res_SIP_byTime_l))) { # didn't manage with map
  print(names(DESeq_res_SIP_byTime_l[i]))
  summary(DESeq_res_SIP_byTime_l[[i]])
}
```

```
## [1] "Certovo & 12 h & Anoxic & Labelled"
## 
## out of 3026 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 19, 0.63%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 3, 0.099%
## low counts [2]     : 2846, 94%
## (mean count < 41)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Certovo & 216 h & Anoxic & Labelled"
## 
## out of 2878 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 127, 4.4%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 1374, 48%
## (mean count < 2)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Certovo & 24 h & Anoxic & Labelled"
## 
## out of 2755 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 3, 0.11%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 4, 0.15%
## (mean count < 0)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Certovo & 48 h & Anoxic & Labelled"
## 
## out of 2787 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 31, 1.1%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 1541, 55%
## (mean count < 2)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Certovo & 216 h & Anoxic & Unlabelled"
## 
## out of 2825 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 0, 0%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 5, 0.18%
## (mean count < 0)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Certovo & 12 h & Oxic & Labelled"
## 
## out of 3053 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 14, 0.46%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 2408, 79%
## (mean count < 5)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Certovo & 24 h & Oxic & Labelled"
## 
## out of 2954 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 167, 5.7%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 1188, 40%
## (mean count < 1)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Certovo & 48 h & Oxic & Labelled"
## 
## out of 2882 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 160, 5.6%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 992, 34%
## (mean count < 1)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Certovo & 72 h & Oxic & Labelled"
## 
## out of 2898 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 198, 6.8%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 1, 0.035%
## low counts [2]     : 1274, 44%
## (mean count < 1)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Certovo & 72 h & Oxic & Unlabelled"
## 
## out of 2979 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 6, 0.2%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 6, 0.2%
## low counts [2]     : 0, 0%
## (mean count < 0)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Plesne & 12 h & Anoxic & Labelled"
## 
## out of 3031 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 104, 3.4%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 1562, 52%
## (mean count < 2)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Plesne & 216 h & Anoxic & Labelled"
## 
## out of 3210 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 33, 1%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 116, 3.6%
## low counts [2]     : 2445, 76%
## (mean count < 11)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Plesne & 24 h & Anoxic & Labelled"
## 
## out of 2848 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 69, 2.4%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 1362, 48%
## (mean count < 2)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Plesne & 48 h & Anoxic & Labelled"
## 
## out of 2758 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 13, 0.47%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 53, 1.9%
## (mean count < 0)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Plesne & 216 h & Anoxic & Unlabelled"
## 
## out of 2862 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 1, 0.035%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 2, 0.07%
## (mean count < 0)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Plesne & 12 h & Oxic & Labelled"
## 
## out of 2923 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 81, 2.8%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 1, 0.034%
## low counts [2]     : 1846, 63%
## (mean count < 2)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Plesne & 24 h & Oxic & Labelled"
## 
## out of 2548 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 137, 5.4%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 1206, 47%
## (mean count < 1)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Plesne & 48 h & Oxic & Labelled"
## 
## out of 2531 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 27, 1.1%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 1485, 59%
## (mean count < 1)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Plesne & 72 h & Oxic & Labelled"
## 
## out of 2556 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 158, 6.2%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 1307, 51%
## (mean count < 1)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Plesne & 72 h & Oxic & Unlabelled"
## 
## out of 2798 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 20, 0.71%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 7, 0.25%
## low counts [2]     : 1546, 55%
## (mean count < 2)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

```r
for (i in seq(1, length(DESeq_res_SIP_byTime_LFC_l))) { # didn't manage with map
  print(names(DESeq_res_SIP_byTime_LFC_l[i]))
  summary(DESeq_res_SIP_byTime_LFC_l[[i]])
}
```

```
## [1] "Certovo & 12 h & Anoxic & Labelled"
## 
## out of 3026 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0.26 (up)    : 0, 0%
## LFC < -0.26 (down) : 0, 0%
## outliers [1]       : 3, 0.099%
## low counts [2]     : 7, 0.23%
## (mean count < 0)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Certovo & 216 h & Anoxic & Labelled"
## 
## out of 2878 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0.26 (up)    : 37, 1.3%
## LFC < -0.26 (down) : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 1484, 52%
## (mean count < 2)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Certovo & 24 h & Anoxic & Labelled"
## 
## out of 2755 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0.26 (up)    : 1, 0.036%
## LFC < -0.26 (down) : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 4, 0.15%
## (mean count < 0)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Certovo & 48 h & Anoxic & Labelled"
## 
## out of 2787 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0.26 (up)    : 15, 0.54%
## LFC < -0.26 (down) : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 1222, 44%
## (mean count < 1)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Certovo & 216 h & Anoxic & Unlabelled"
## 
## out of 2825 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0.26 (up)    : 0, 0%
## LFC < -0.26 (down) : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 5, 0.18%
## (mean count < 0)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Certovo & 12 h & Oxic & Labelled"
## 
## out of 3053 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0.26 (up)    : 7, 0.23%
## LFC < -0.26 (down) : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 8, 0.26%
## (mean count < 0)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Certovo & 24 h & Oxic & Labelled"
## 
## out of 2954 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0.26 (up)    : 107, 3.6%
## LFC < -0.26 (down) : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 905, 31%
## (mean count < 1)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Certovo & 48 h & Oxic & Labelled"
## 
## out of 2882 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0.26 (up)    : 122, 4.2%
## LFC < -0.26 (down) : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 331, 11%
## (mean count < 0)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Certovo & 72 h & Oxic & Labelled"
## 
## out of 2898 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0.26 (up)    : 142, 4.9%
## LFC < -0.26 (down) : 0, 0%
## outliers [1]       : 1, 0.035%
## low counts [2]     : 1329, 46%
## (mean count < 1)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Certovo & 72 h & Oxic & Unlabelled"
## 
## out of 2979 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0.26 (up)    : 0, 0%
## LFC < -0.26 (down) : 0, 0%
## outliers [1]       : 6, 0.2%
## low counts [2]     : 0, 0%
## (mean count < 0)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Plesne & 12 h & Anoxic & Labelled"
## 
## out of 3031 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0.26 (up)    : 27, 0.89%
## LFC < -0.26 (down) : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 1678, 55%
## (mean count < 3)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Plesne & 216 h & Anoxic & Labelled"
## 
## out of 3210 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0.26 (up)    : 0, 0%
## LFC < -0.26 (down) : 0, 0%
## outliers [1]       : 116, 3.6%
## low counts [2]     : 3, 0.093%
## (mean count < 0)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Plesne & 24 h & Anoxic & Labelled"
## 
## out of 2848 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0.26 (up)    : 22, 0.77%
## LFC < -0.26 (down) : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 2340, 82%
## (mean count < 15)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Plesne & 48 h & Anoxic & Labelled"
## 
## out of 2758 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0.26 (up)    : 12, 0.44%
## LFC < -0.26 (down) : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 1527, 55%
## (mean count < 2)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Plesne & 216 h & Anoxic & Unlabelled"
## 
## out of 2862 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0.26 (up)    : 0, 0%
## LFC < -0.26 (down) : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 2, 0.07%
## (mean count < 0)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Plesne & 12 h & Oxic & Labelled"
## 
## out of 2923 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0.26 (up)    : 57, 2%
## LFC < -0.26 (down) : 0, 0%
## outliers [1]       : 1, 0.034%
## low counts [2]     : 1790, 61%
## (mean count < 1)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Plesne & 24 h & Oxic & Labelled"
## 
## out of 2548 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0.26 (up)    : 101, 4%
## LFC < -0.26 (down) : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 1061, 42%
## (mean count < 1)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Plesne & 48 h & Oxic & Labelled"
## 
## out of 2531 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0.26 (up)    : 13, 0.51%
## LFC < -0.26 (down) : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 1820, 72%
## (mean count < 3)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Plesne & 72 h & Oxic & Labelled"
## 
## out of 2556 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0.26 (up)    : 115, 4.5%
## LFC < -0.26 (down) : 0, 0%
## outliers [1]       : 0, 0%
## low counts [2]     : 1113, 44%
## (mean count < 1)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
## 
## [1] "Plesne & 72 h & Oxic & Unlabelled"
## 
## out of 2798 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0.26 (up)    : 3, 0.11%
## LFC < -0.26 (down) : 0, 0%
## outliers [1]       : 7, 0.25%
## low counts [2]     : 0, 0%
## (mean count < 0)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

```r
# Store labelled OTUs and save them to a file
# DESeq_res_SIP_byTime_l %>% 
#   map(., ~subset(.x, padj < alpha_thresh & log2FoldChange > LFC_thresh)) %>% 
#   map(., ~as.data.frame(.x)) %>% 
#   map(., ~rownames_to_column(.x, "ASV")) %>% 
#   bind_rows(., .id = "Comparison") %>% 
#   arrange(Comparison, desc(baseMean)) %T>% 
#   write_csv(., file = "DESeq2_byTime_a-0.05.txt") ->
#   DESeq_res_SIP_byTime_df

DESeq_res_SIP_byTime_LFC_l %>% 
  map(., ~subset(.x, padj < alpha_thresh & log2FoldChange > LFC_thresh)) %>% 
  map(., ~as.data.frame(.x)) %>% 
  map(., ~rownames_to_column(.x, "ASV")) %>% 
  bind_rows(., .id = "Comparison") %>% 
  arrange(Comparison, desc(baseMean)) %>% 
  separate(., "Comparison" ,c("Site","Hours", "Oxygen", "Label"), sep = " & ") %T>% 
  write_csv(., file = "DESeq2_byTime_a-0.05_LFC0-322.txt") ->
  DESeq_res_SIP_byTime_LFC_sig_df
```


```r
DESeq_res_SIP_byTime_LFC_sig_df %>% 
  get_variable() %>% 
  select_if(is.numeric) %>% 
  vis_value()
```

![](05_Diff_abund_figures/vis DES res-1.png)<!-- -->

```r
DESeq_res_SIP_byTime_LFC_sig_df %>% 
  get_variable() %>% 
  select_if(is.numeric) %>% 
  vis_cor()
```

![](05_Diff_abund_figures/vis DES res-2.png)<!-- -->


```r
# ps_obj <- Ps_obj_SIP
# DESeq_results <- DESeq_res_SIP_byTime_LFC0.322_l[9]
# plot_DESeq(DESeq_results, ps_obj, plot_title = names(DESeq_results))

DESeq_plots <- lapply(seq(length(DESeq_res_SIP_byTime_LFC_shrink_l)), 
                        function(x) {plot_DESeq(DESeq_res_SIP_byTime_LFC_shrink_l[x],  
                                                Ps_obj_SIP, plot_title = names(DESeq_res_SIP_byTime_LFC_shrink_l[x]))})


(DESeq_plots[[6]] + 
    theme(legend.position = "none") +
    theme(axis.text.x = element_blank())) +
  (DESeq_plots[[1]] + 
     theme(legend.position = "none", 
           axis.text.x = element_blank(), 
           axis.title.y = element_blank())) +
  (DESeq_plots[[7]] + 
     theme(legend.position = "none",
           axis.text.x = element_blank())) +
  (DESeq_plots[[3]] + 
     theme(legend.position = "none", 
           axis.text.x = element_blank(), 
           axis.title.y = element_blank())) +
  (DESeq_plots[[8]] + 
     theme(legend.position = "none") +
     theme(legend.position = "none",
           axis.text.x = element_blank())) +
  (DESeq_plots[[4]] + 
     theme(legend.position = "none") +
     theme(legend.position = "none", 
           axis.text.x = element_blank(), 
           axis.title.y = element_blank()) +
     ylim(NA, 5)) +
  (DESeq_plots[[9]] + 
     theme(legend.position = "none",
           axis.text.x = element_blank())) +
  (DESeq_plots[[2]] + 
     theme(legend.position = "none", 
           axis.text.x = element_blank(), 
           axis.title.y = element_blank())) +
  (DESeq_plots[[10]] + 
     theme(legend.position = "none")) +
  (DESeq_plots[[5]] + 
     theme(legend.position = "none", 
           axis.title.y = element_blank())) + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(legend.position = 'bottom')
```

![](05_Diff_abund_figures/plot DESeq2 models-1.png)<!-- -->

```r
(DESeq_plots[[16]] + 
    theme(legend.position = "none") +
    theme(axis.text.x = element_blank())) +
  (DESeq_plots[[11]] + 
     theme(legend.position = "none", 
           axis.text.x = element_blank(), 
           axis.title.y = element_blank())) +
  (DESeq_plots[[17]] + 
     theme(legend.position = "none",
           axis.text.x = element_blank())) +
  (DESeq_plots[[13]] + 
     theme(legend.position = "none", 
           axis.text.x = element_blank(), 
           axis.title.y = element_blank())) +
  (DESeq_plots[[18]] + 
     theme(legend.position = "none") +
     theme(legend.position = "none",
           axis.text.x = element_blank())) +
  (DESeq_plots[[14]] + 
     theme(legend.position = "none") +
     theme(legend.position = "none", 
           axis.text.x = element_blank(), 
           axis.title.y = element_blank())) +
  (DESeq_plots[[19]] + 
     theme(legend.position = "none",
           axis.text.x = element_blank())) +
  (DESeq_plots[[12]] + 
     theme(legend.position = "none", 
           axis.text.x = element_blank(), 
           axis.title.y = element_blank()) +
     ylim(-3, NA)) +
  (DESeq_plots[[20]] + 
     theme(legend.position = "none")) +
  (DESeq_plots[[15]] + 
     theme(legend.position = "none", 
           axis.title.y = element_blank())) + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(legend.position = 'bottom')
```

![](05_Diff_abund_figures/plot DESeq2 models-2.png)<!-- -->


```r
plot_otus_by_density(ps_obj = Ps_obj_SIP_noTime_l[[2]], 
                     ASV2plot = filter(DESeq_res_SIP_byTime_LFC_sig_df, Site == "Certovo", Oxygen == "Oxic"))
```

![](05_Diff_abund_figures/Plot incorporators-1.png)<!-- -->

```r
plot_otus_by_density(ps_obj = Ps_obj_SIP_noTime_l[[1]], 
                     ASV2plot = filter(DESeq_res_SIP_byTime_LFC_sig_df, Site == "Certovo", Oxygen == "Anoxic"))
```

![](05_Diff_abund_figures/Plot incorporators-2.png)<!-- -->

```r
plot_otus_by_density(ps_obj = Ps_obj_SIP_noTime_l[[4]], 
                     ASV2plot = filter(DESeq_res_SIP_byTime_LFC_sig_df, Site == "Plesne", Oxygen == "Oxic"))
```

![](05_Diff_abund_figures/Plot incorporators-3.png)<!-- -->

```r
plot_otus_by_density(ps_obj = Ps_obj_SIP_noTime_l[[3]],
                     ASV2plot = filter(DESeq_res_SIP_byTime_LFC_sig_df, Site == "Plesne", Oxygen == "Anoxic"))
```

![](05_Diff_abund_figures/Plot incorporators-4.png)<!-- -->


```r
c("Certovo Oxic 12 h", "Certovo Oxic 24 h", "Certovo Oxic 48 h", "Certovo Oxic 72 h", "Certovo Anoxic 12 h", "Certovo Anoxic 24 h", "Certovo Anoxic 48 h", "Certovo Anoxic 216 h", "Plesne Oxic 12 h", "Plesne Oxic 24 h", "Plesne Oxic 48 h", "Plesne Oxic 72 h", "Plesne Anoxic 12 h",  "Plesne Anoxic 24 h",  "Plesne Anoxic 48 h",  "Plesne Anoxic 216 h" )  ->
  col_order

DESeq_res_SIP_byTime_LFC_l %>% 
  map(., ~as.data.frame(.x)) %>% 
  map(., ~rownames_to_column(.x, "ASV")) %>% 
  bind_rows(., .id = "Comparison") %>% 
  filter(str_detect(Comparison, "Labelled")) %>% # remove unlabelled samples [c(-5, -10, -15, -20)]
  mutate(Labelled = ifelse(padj < alpha_thresh & log2FoldChange > LFC_thresh, "Labelled", "Unlabelled")) %>% 
  # arrange(Comparison, desc(baseMean)) %>% 
  separate(., "Comparison" ,c("Site","Hours", "Oxygen", "Label"), sep = " & ") %>% 
  mutate(Site_Oxygen_Hours = paste(Site, Oxygen, Hours)) %>% 
  mutate(across(Site_Oxygen_Hours, ~factor(., levels = col_order))) ->
  # mutate(Site_oxygen = paste(Site, Oxygen)) ->
  DESeq_res_SIP_byTime_all_df

# Summarise number of labelled and unlabelled ASVs
DESeq_res_SIP_byTime_all_df %>% 
  group_by(Labelled) %>% 
  summarise(n = n()) 
```

<div class="kable-table">

<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Labelled </th>
   <th style="text-align:right;"> n </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Labelled </td>
   <td style="text-align:right;"> 778 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Unlabelled </td>
   <td style="text-align:right;"> 37162 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> 21148 </td>
  </tr>
</tbody>
</table>

</div>

```r
# detect taxa with NA from DESeq analysis
DESeq_res_SIP_byTime_all_df %<>% 
  filter(!is.na(Labelled)) #%>% 
  # pull(Labelled) -> 
  # bad_seqs

# remove NA taxa from PS obj
Ps_obj_SIP %>% 
  prune_taxa(DESeq_res_SIP_byTime_all_df$ASV, .) -> 
  Ps_obj_SIP4tree_plot

p_actino_t <- plot_ggtree(Ps_obj_SIP4tree_plot, 
                          rank = "Class",
                          subrank = "Order",
                          Taxa2plot = "Actinobacteria")
p_actino_HM <- plot_ggtree_heatmap(p_actino_t, 
                                 Ps_obj_SIP4tree_plot,
                                 DESeq_res_SIP_byTime_all_df, 
                                 rank = "Class",
                                 Taxa2plot = "Actinobacteria",
                                 sample_names = "Site_Oxygen_Hours",
                                 sample_colours = rep(brewer.pal(n = 11, 
                                                                 "RdYlGn")[c(1, 3, 11, 9)],
                                                      each = 4))
p_actino_t + p_actino_HM
```

![](05_Diff_abund_figures/Plot trees-1.png)<!-- -->

```r
p_Aproteo_t <- plot_ggtree(Ps_obj_SIP4tree_plot, 
                          rank = "Class",
                          subrank = "Order",
                          Taxa2plot = "Alphaproteobacteria")
p_Aproteo_HM <- plot_ggtree_heatmap(p_Aproteo_t, 
                                 Ps_obj_SIP4tree_plot,
                                 DESeq_res_SIP_byTime_all_df, 
                                 rank = "Class",
                                 Taxa2plot = "Alphaproteobacteria",
                                 sample_names = "Site_Oxygen_Hours",
                                 sample_colours = rep(brewer.pal(n = 11, 
                                                                 "RdYlGn")[c(1, 3, 11, 9)],
                                                      each = 4))
p_Aproteo_t + p_Aproteo_HM
```

![](05_Diff_abund_figures/Plot trees-2.png)<!-- -->

```r
p_Gproteo_t <- plot_ggtree(Ps_obj_SIP4tree_plot, 
                          rank = "Class",
                          subrank = "Order",
                          Taxa2plot = "Gammaproteobacteria")
p_Gproteo_HM <- plot_ggtree_heatmap(p_Gproteo_t, 
                                 Ps_obj_SIP4tree_plot,
                                 DESeq_res_SIP_byTime_all_df, 
                                 rank = "Class",
                                 Taxa2plot = "Gammaproteobacteria",
                                 sample_names = "Site_Oxygen_Hours",
                                 sample_colours = rep(brewer.pal(n = 11, 
                                                                 "RdYlGn")[c(1, 3, 11, 9)],
                                                      each = 4))
p_Gproteo_t + p_Gproteo_HM
```

![](05_Diff_abund_figures/Plot trees-3.png)<!-- -->

```r
p_acido_t <- plot_ggtree(Ps_obj_SIP4tree_plot, 
                          rank = "Phylum",
                          subrank = "Class",
                          Taxa2plot = "Acidobacteriota")
p_acido_HM <- plot_ggtree_heatmap(p_acido_t, 
                                 Ps_obj_SIP4tree_plot,
                                 DESeq_res_SIP_byTime_all_df, 
                                 rank = "Phylum",
                                 Taxa2plot = "Acidobacteriota",
                                 sample_names = "Site_Oxygen_Hours",
                                 sample_colours = rep(brewer.pal(n = 11, 
                                                                 "RdYlGn")[c(1, 3, 11, 9)],
                                                      each = 4))
p_acido_t + p_acido_HM
```

![](05_Diff_abund_figures/Plot trees-4.png)<!-- -->

```r
p_verruc_t <- plot_ggtree(Ps_obj_SIP4tree_plot, 
                          rank = "Phylum",
                          subrank = "Class",
                          Taxa2plot = "Verrucomicrobiota")
p_verruc_HM <- plot_ggtree_heatmap(p_verruc_t, 
                                 Ps_obj_SIP4tree_plot,
                                 DESeq_res_SIP_byTime_all_df, 
                                 rank = "Phylum",
                                 Taxa2plot = "Verrucomicrobiota",
                                 sample_names = "Site_Oxygen_Hours",
                                 sample_colours = rep(brewer.pal(n = 11, 
                                                                 "RdYlGn")[c(1, 3, 11, 9)],
                                                      each = 4))
p_verruc_t + p_verruc_HM
```

![](05_Diff_abund_figures/Plot trees-5.png)<!-- -->

```r
p_bacteroi_t <- plot_ggtree(Ps_obj_SIP4tree_plot, 
                          rank = "Phylum",
                          subrank = "Class",
                          Taxa2plot = "Bacteroidota")
p_bacteroi_HM <- plot_ggtree_heatmap(p_bacteroi_t, 
                                 Ps_obj_SIP4tree_plot,
                                 DESeq_res_SIP_byTime_all_df, 
                                 rank = "Phylum",
                                 Taxa2plot = "Bacteroidota",
                                 sample_names = "Site_Oxygen_Hours",
                                 sample_colours = rep(brewer.pal(n = 11, 
                                                                 "RdYlGn")[c(1, 3, 11, 9)],
                                                      each = 4))
p_bacteroi_t + p_bacteroi_HM
```

![](05_Diff_abund_figures/Plot trees-6.png)<!-- -->

```r
p_firmi_t <- plot_ggtree(Ps_obj_SIP4tree_plot, 
                          rank = "Phylum",
                          subrank = "Class",
                          Taxa2plot = "Firmicutes")
p_firmi_HM <- plot_ggtree_heatmap(p_firmi_t, 
                                 Ps_obj_SIP4tree_plot,
                                 DESeq_res_SIP_byTime_all_df, 
                                 rank = "Phylum",
                                 Taxa2plot = "Firmicutes",
                                 sample_names = "Site_Oxygen_Hours",
                                 sample_colours = rep(brewer.pal(n = 11, 
                                                                 "RdYlGn")[c(1, 3, 11, 9)],
                                                      each = 4))
p_firmi_t + p_firmi_HM
```

![](05_Diff_abund_figures/Plot trees-7.png)<!-- -->

```r
c("Certovo Oxic 12 h", "Certovo Oxic 24 h", "Certovo Oxic 48 h", "Certovo Oxic 72 h", "Certovo Anoxic 12 h", "Certovo Anoxic 24 h", "Certovo Anoxic 48 h", "Certovo Anoxic 216 h", "Plesne Oxic 12 h", "Plesne Oxic 24 h", "Plesne Oxic 48 h", "Plesne Oxic 72 h", "Plesne Anoxic 12 h",  "Plesne Anoxic 24 h",  "Plesne Anoxic 48 h",  "Plesne Anoxic 216 h" ) %>% 
 map_dfr(~tibble(!!.x := double())) ->
  col_order

DESeq_res_SIP_byTime_all_df_actino %>% 
  filter(Labelled == "Labelled") %>% 
  select(ASV, Site, Oxygen, Hours, log2FoldChange) %>% 
  pivot_wider(names_from = c(Site, Oxygen, Hours), 
              values_from = log2FoldChange,
              names_sep = " ") %>% 
  bind_rows(col_order, .) %>% 
  # add_column(., !!!heatmap_blnk_df[setdiff(names(heatmap_blnk_df), names(.))])
  # set_names() %>% 
  # select(c("ASV", "Certovo Oxic_12 h", "Certovo Oxic_24 h", "Certovo Oxic_48 h", "Certovo Oxic_72 h", "Certovo Anoxic_12 h", "Certovo Anoxic_24 h", "Certovo Anoxic_48 h", "Certovo Anoxic_216 h", "Plesne Oxic_12 h", "Plesne Oxic_24 h", "Plesne Oxic_48 h", "Plesne Oxic_72 h", "Plesne Anoxic_12 h",  "Plesne Anoxic_24 h",  "Plesne Anoxic_48 h",  "Plesne Anoxic_216 h" )) %>%
  # modify_if(is.double, ~replace_na(.x, "")) %>%
  column_to_rownames("ASV") %>% 
  set_colnames(labels) ->
  DESeq_res_SIP_byTime_all_df_actino_hm 

gheatmap(
  p = p_actino_t, 
  data = DESeq_res_SIP_byTime_all_df_actino_hm,
  offset = 0, width = 1.6, font.size = 3, 
  low = "grey", 
  high = "darkred",
  color = "grey95",
  colnames_position = "bottom",
  colnames_offset_x = -0.05,
  colnames_offset_y = 1,
  colnames_angle = 25, 
  hjust = 0,
  legend_title = "L2FC") +
  theme_tree(legend.position = "bottom") +
  theme(
    axis.text.x = element_markdown(color = "red", size = 11)
  ) +
  guides(colour = guide_legend(nrow = 3, 
                               byrow = TRUE))


p <- ggtree(Actino_ps) + 
  # geom_treescale() + 
  geom_tiplab(size=2) +
  # geom_treescale(x=2008, y=1, offset=2) 
gheatmap(p, DESeq_res_SIP_byTime_all_df_actino_hm) +
  # scale_fill_brewer()
  scale_fill_manual(breaks=c("Unlabelled", "Labelled"), 
                    values=c("firebrick", "darkgreen"), name="DESeq_res_SIP_byTime_all_df_actino_hm")


gheatmap(p_actino_t, 
         DESeq_res_SIP_byTime_all_df_actino_hm, offset = 10, color = NULL, 
         colnames_position = "top", 
         colnames_angle = 90, colnames_offset_y = 1, 
         hjust = 0, font.size = 2)


+
  scale_fill_manual(values = heatmap.colours, breaks=0:15)




library(ggnewscale)

# p_actino_t %<+%  DESeq_res_SIP_byTime_all_df_actino# add data
# p_actino_t + 
#         new_scale_colour() + 
#   geom_fruit(
#           data = DESeq_res_SIP_byTime_all_df_actino,
#           geom = geom_point,
#           mapping = aes(y = ASV, x = Site_oxygen, fill = Labelled),
#           offset=0.08,   # The distance between layers, default is 0.03 of x range of tree.
#           pwidth=0.25, # width of the layer, default is 0.2 of x range of tree.
#           axis.params=list(                
#                         axis="x",  # add axis text of the layer.
#                         text.angle=-45, # the text size of axis.
#                         hjust=0  # adust the horizontal position of text of axis.
#                       )
#       ) 

beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
beast_tree <- read.beast(beast_file)

genotype_file <- system.file("examples/Genotype.txt", package="ggtree")
genotype <- read.table(genotype_file, sep="\t", stringsAsFactor=F)
colnames(genotype) <- sub("\\.$", "", colnames(genotype))
p <- ggtree(beast_tree, mrsd="2013-01-01") + 
  geom_treescale(x=2008, y=1, offset=2) + 
  geom_tiplab(size=2)
gheatmap(p, genotype, offset=5, width=0.5, font.size=3, 
         colnames_angle=-45, hjust=0) +
  scale_fill_manual(breaks=c("HuH3N2", "pdm", "trig"), 
                    values=c("steelblue", "firebrick", "darkgreen"), name="genotype")
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

─ Session info ─────────────────────────────────────────────────────────────────────────
 setting  value                       
 version  R version 4.0.3 (2020-10-10)
 os       Ubuntu 18.04.5 LTS          
 system   x86_64, linux-gnu           
 ui       X11                         
 language (EN)                        
 collate  en_US.UTF-8                 
 ctype    en_US.UTF-8                 
 tz       Europe/Prague               
 date     2021-02-09                  

─ Packages ─────────────────────────────────────────────────────────────────────────────
 package              * version    date       lib
 ade4                   1.7-16     2020-10-28 [1]
 annotate               1.66.0     2020-04-27 [1]
 AnnotationDbi          1.50.3     2020-07-25 [1]
 ape                    5.4-1      2020-08-13 [1]
 aplot                  0.0.6      2020-09-03 [1]
 ashr                   2.2-47     2020-02-20 [1]
 assertthat             0.2.1      2019-03-21 [1]
 backports              1.2.1      2020-12-09 [1]
 Biobase              * 2.48.0     2020-04-27 [1]
 BiocGenerics         * 0.34.0     2020-04-27 [1]
 BiocManager            1.30.10    2019-11-16 [1]
 BiocParallel           1.22.0     2020-04-27 [1]
 biomformat             1.16.0     2020-04-27 [1]
 Biostrings           * 2.56.0     2020-04-27 [1]
 bit                    4.0.4      2020-08-04 [1]
 bit64                  4.0.5      2020-08-30 [1]
 bitops                 1.0-6      2013-08-17 [1]
 blob                   1.2.1      2020-01-20 [1]
 broom                  0.7.4      2021-01-29 [1]
 cachem                 1.0.3      2021-02-04 [1]
 cellranger             1.1.0      2016-07-27 [1]
 cli                    2.3.0      2021-01-31 [1]
 clipr                  0.7.1      2020-10-08 [1]
 cluster                2.1.0      2019-06-19 [1]
 codetools              0.2-18     2020-11-04 [1]
 colorspace             2.0-0      2020-11-11 [1]
 crayon                 1.4.0      2021-01-30 [1]
 data.table             1.13.6     2020-12-30 [1]
 DBI                    1.1.1      2021-01-15 [1]
 dbplyr                 2.1.0      2021-02-03 [1]
 DelayedArray         * 0.14.1     2020-07-14 [1]
 desc                   1.2.0      2018-05-01 [1]
 DESeq2               * 1.28.1     2020-05-12 [1]
 details                0.2.1      2020-01-12 [1]
 digest                 0.6.27     2020-10-24 [1]
 dplyr                * 1.0.4      2021-02-02 [1]
 ellipsis               0.3.1      2020-05-15 [1]
 evaluate               0.14       2019-05-28 [1]
 extrafont            * 0.17       2014-12-08 [1]
 extrafontdb            1.0        2012-06-11 [1]
 farver                 2.0.3      2020-01-16 [1]
 fastmap                1.1.0      2021-01-25 [1]
 forcats              * 0.5.1      2021-01-27 [1]
 foreach                1.5.1      2020-10-15 [1]
 fs                     1.5.0      2020-07-31 [1]
 gdtools                0.2.3      2021-01-06 [1]
 genefilter             1.70.0     2020-04-27 [1]
 geneplotter            1.66.0     2020-04-27 [1]
 generics               0.1.0      2020-10-31 [1]
 GenomeInfoDb         * 1.24.2     2020-06-15 [1]
 GenomeInfoDbData       1.2.3      2020-08-13 [1]
 GenomicRanges        * 1.40.0     2020-04-27 [1]
 ggnewscale             0.4.5      2021-01-11 [1]
 ggplot2              * 3.3.3      2020-12-30 [1]
 ggpomological        * 0.1.2      2020-08-13 [1]
 ggrepel              * 0.9.1      2021-01-15 [1]
 ggsci                  2.9        2018-05-14 [1]
 ggtext               * 0.1.1      2020-12-17 [1]
 ggtree               * 2.2.4      2020-07-28 [1]
 ggtreeExtra          * 1.1.4      2021-02-01 [1]
 glue                 * 1.4.2      2020-08-27 [1]
 gridExtra              2.3        2017-09-09 [1]
 gridtext               0.1.4      2020-12-10 [1]
 gtable                 0.3.0      2019-03-25 [1]
 haven                  2.3.1      2020-06-01 [1]
 highr                  0.8        2019-03-20 [1]
 hms                    1.0.0      2021-01-13 [1]
 htmltools              0.5.1.1    2021-01-22 [1]
 HTSSIP               * 1.4.1      2021-01-15 [1]
 httr                   1.4.2      2020-07-20 [1]
 igraph                 1.2.6      2020-10-06 [1]
 invgamma               1.1        2017-05-07 [1]
 IRanges              * 2.22.2     2020-05-21 [1]
 irlba                  2.3.3      2019-02-05 [1]
 iterators              1.0.13     2020-10-15 [1]
 jsonlite               1.7.2      2020-12-09 [1]
 kableExtra           * 1.3.1      2020-10-22 [1]
 knitr                  1.31       2021-01-27 [1]
 labeling               0.4.2      2020-10-20 [1]
 lattice              * 0.20-41    2020-04-02 [1]
 lazyeval               0.2.2      2019-03-15 [1]
 lifecycle              0.2.0      2020-03-06 [1]
 locfit                 1.5-9.4    2020-03-25 [1]
 lubridate              1.7.9.2    2020-11-13 [1]
 magrittr             * 2.0.1      2020-11-17 [1]
 markdown               1.1        2019-08-07 [1]
 MASS                   7.3-53     2020-09-09 [1]
 Matrix                 1.3-2      2021-01-06 [1]
 matrixStats          * 0.58.0     2021-01-29 [1]
 memoise                2.0.0      2021-01-26 [1]
 mgcv                   1.8-33     2020-08-27 [1]
 mixsqp                 0.3-43     2020-05-14 [1]
 modelr                 0.1.8      2020-05-19 [1]
 multtest               2.44.0     2020-04-27 [1]
 munsell                0.5.0      2018-06-12 [1]
 nlme                   3.1-152    2021-02-04 [1]
 patchwork            * 1.1.1      2020-12-17 [1]
 permute              * 0.9-5      2019-03-12 [1]
 phyloseq             * 1.32.0     2020-04-27 [1]
 pillar                 1.4.7      2020-11-20 [1]
 pkgconfig              2.0.3      2019-09-22 [1]
 plyr                   1.8.6      2020-03-03 [1]
 png                    0.1-7      2013-12-03 [1]
 prettyunits            1.1.1      2020-01-24 [1]
 progress               1.2.2      2019-05-16 [1]
 purrr                * 0.3.4      2020-04-17 [1]
 R6                     2.5.0      2020-10-28 [1]
 RColorBrewer         * 1.1-2      2014-12-07 [1]
 Rcpp                   1.0.6      2021-01-15 [1]
 RCurl                  1.98-1.2   2020-04-18 [1]
 readr                * 1.4.0      2020-10-05 [1]
 readxl                 1.3.1      2019-03-13 [1]
 reprex                 1.0.0      2021-01-27 [1]
 reshape2               1.4.4      2020-04-09 [1]
 rhdf5                  2.32.4     2020-10-05 [1]
 Rhdf5lib               1.10.1     2020-07-09 [1]
 rlang                  0.4.10     2020-12-30 [1]
 rmarkdown              2.6        2020-12-14 [1]
 rprojroot              2.0.2      2020-11-15 [1]
 RSQLite                2.2.3      2021-01-24 [1]
 rstudioapi             0.13       2020-11-12 [1]
 Rttf2pt1               1.3.8      2020-01-10 [1]
 rvcheck                0.1.8      2020-03-01 [1]
 rvest                  0.3.6      2020-07-25 [1]
 S4Vectors            * 0.26.1     2020-05-16 [1]
 scales               * 1.1.1      2020-05-11 [1]
 sessioninfo            1.1.1      2018-11-05 [1]
 speedyseq            * 0.5.3.9001 2020-10-27 [1]
 SQUAREM                2021.1     2021-01-13 [1]
 stringi                1.5.3      2020-09-09 [1]
 stringr              * 1.4.0      2019-02-10 [1]
 SummarizedExperiment * 1.18.2     2020-07-09 [1]
 survival               3.2-7      2020-09-28 [1]
 svglite              * 1.2.3.2    2020-07-07 [1]
 systemfonts            1.0.0      2021-02-01 [1]
 tibble               * 3.0.6      2021-01-29 [1]
 tidyr                * 1.1.2      2020-08-27 [1]
 tidyselect             1.1.0      2020-05-11 [1]
 tidytree               0.3.3      2020-04-02 [1]
 tidyverse            * 1.3.0      2019-11-21 [1]
 treeio                 1.12.0     2020-04-27 [1]
 truncnorm              1.0-8      2018-02-27 [1]
 vctrs                  0.3.6      2020-12-17 [1]
 vegan                * 2.5-7      2020-11-28 [1]
 viridis              * 0.5.1      2018-03-29 [1]
 viridisLite          * 0.3.0      2018-02-01 [1]
 visdat               * 0.6.0.9000 2021-02-01 [1]
 webshot                0.5.2      2019-11-22 [1]
 withr                  2.4.1      2021-01-26 [1]
 xfun                   0.20       2021-01-06 [1]
 XML                    3.99-0.5   2020-07-23 [1]
 xml2                   1.3.2      2020-04-23 [1]
 xtable                 1.8-4      2019-04-21 [1]
 XVector              * 0.28.0     2020-04-27 [1]
 yaml                   2.2.1      2020-02-01 [1]
 zlibbioc               1.34.0     2020-04-27 [1]
 source                                  
 CRAN (R 4.0.2)                          
 Bioconductor                            
 Bioconductor                            
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 Bioconductor                            
 Bioconductor                            
 CRAN (R 4.0.2)                          
 Bioconductor                            
 Bioconductor                            
 Bioconductor                            
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.3)                          
 Bioconductor                            
 CRAN (R 4.0.2)                          
 Bioconductor                            
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 Bioconductor                            
 Bioconductor                            
 CRAN (R 4.0.2)                          
 Bioconductor                            
 Bioconductor                            
 Bioconductor                            
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.2)                          
 Github (gadenbuie/ggpomological@69f3815)
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 Bioconductor                            
 Github (xiangpin/ggtreeExtra@0c3e16c)   
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.3)                          
 Github (buckleylab/HTSSIP@29ec56b)      
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 Bioconductor                            
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 Bioconductor                            
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 Bioconductor                            
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.2)                          
 Bioconductor                            
 Bioconductor                            
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.2)                          
 Bioconductor                            
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 Github (mikemc/speedyseq@8daed32)       
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 Bioconductor                            
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.2)                          
 Bioconductor                            
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 Github (ropensci/visdat@8121dfe)        
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 Bioconductor                            
 CRAN (R 4.0.2)                          
 Bioconductor                            

[1] /home/angel/R/library
[2] /usr/local/lib/R/site-library
[3] /usr/lib/R/site-library
[4] /usr/lib/R/library

```

</details>
<br>

## References

