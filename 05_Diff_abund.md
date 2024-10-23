Anaerobic CUE
================
Roey Angel
`2021-02-20`

- <a href="#differential-abundance-modelling-of-sip-gradients"
  id="toc-differential-abundance-modelling-of-sip-gradients">Differential
  abundance modelling of SIP gradients</a>
  - <a href="#setting-general-parameters"
    id="toc-setting-general-parameters">Setting general parameters:</a>
  - <a href="#read-phyloseq-object" id="toc-read-phyloseq-object">Read
    phyloseq object</a>
  - <a href="#subset-the-dataset" id="toc-subset-the-dataset">Subset the
    dataset</a>
  - <a href="#beta-diversity-analysis" id="toc-beta-diversity-analysis">Beta
    diversity analysis</a>
  - <a href="#differential-abundance-models"
    id="toc-differential-abundance-models">Differential abundance models</a>
    - <a href="#inspect-results" id="toc-inspect-results">Inspect results</a>
    - <a href="#plot-differential-abundance-models"
      id="toc-plot-differential-abundance-models">Plot differential abundance
      models</a>
  - <a href="#plot-labelled-asvs" id="toc-plot-labelled-asvs">Plot labelled
    ASVs</a>
    - <a href="#plot-phylogenetic-trees-with-heatmaps"
      id="toc-plot-phylogenetic-trees-with-heatmaps">Plot phylogenetic trees
      with heatmaps</a>
    - <a href="#how-abundant-were-the-labelled-asvs"
      id="toc-how-abundant-were-the-labelled-asvs">How abundant were the
      labelled ASVs?</a>
    - <a href="#which-labelled-asvs-are-shared"
      id="toc-which-labelled-asvs-are-shared">Which labelled ASVs are
      shared?</a>
    - <a href="#calculate-nti" id="toc-calculate-nti">Calculate NTI</a>
- <a href="#references" id="toc-references">References</a>

## Differential abundance modelling of SIP gradients

Here we attempt to detect ASVs that were labelled with <sup>13</sup>C in
our soil incubations using differential abundance modelling. Using
DESeq2 ([Love, Huber and Anders 2014](#ref-love_moderated_2014)) we
compare the relative abundance of each ASV in the fractions where
<sup>13</sup>C-labelled RNA is expected to be found (\>1.795 g
ml<sup>-1</sup>; AKA ‘heavy’ fractions) to the fractions where
unlabelled RNA is expected to be found (\<1.795 g ml<sup>-1</sup>; AKA
‘light’ fractions). The method has been previously described in Angel et
al., ([2018](#ref-angel_application_2018)).

### Setting general parameters:

``` r
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

``` r
# Load phylogenetic tree
Tree <- read_tree(paste0(data_path, Tree_file))

# load and merge  phyloseq object
readRDS(paste0(data_path, Ps_file)) %>% 
  merge_phyloseq(.,
                 phy_tree(Tree)
  ) -> Ps_obj_SIP

sample_data(Ps_obj_SIP)$Site %<>% recode_factor(Plesne  = "Plešné",
                                               Certovo = "Čertovo",
                                               .ordered = TRUE)
```

### Subset the dataset

Because the DESeq2 models will be run on each gradient separately, we
need to subset This is easily done using `HTSSIP::phyloseq_subset`
([Youngblut, Barnett and Buckley 2018](#ref-youngblut_htssip_2018))

``` r
# split, ignore time points (for labelled ASV plots)
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

# split, include time points (for DESeq2 modelling)
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
```

### Beta diversity analysis

Let us look first at the dissimilarity in community composition between
the different fractions. If the labelling was strong enough we should
see a deivation of (some of) the heavy fractions from the light ones.
However, a lack of a significant deviation does not mean unsuccesful
labelling because if only a small minority of the community was labelled
we might not see it here (but we will, hopefully, see it using DESeq2
modelling).

``` r
(mod1 <- adonis2(vegdist(otu_table(Ps_obj_SIP), method = "horn") ~ Lib.size + Site * Oxygen * Hours ,
  data = as(sample_data(Ps_obj_SIP), "data.frame"),
  permutations = 999
))
```

<div class="kable-table">

<table>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
Df
</th>
<th style="text-align:right;">
SumOfSqs
</th>
<th style="text-align:right;">
R2
</th>
<th style="text-align:right;">
F
</th>
<th style="text-align:right;">
Pr(\>F)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Model
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
34.27409
</td>
<td style="text-align:right;">
0.6734282
</td>
<td style="text-align:right;">
75.0094
</td>
<td style="text-align:right;">
0.001
</td>
</tr>
<tr>
<td style="text-align:left;">
Residual
</td>
<td style="text-align:right;">
291
</td>
<td style="text-align:right;">
16.62085
</td>
<td style="text-align:right;">
0.3265718
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
Total
</td>
<td style="text-align:right;">
299
</td>
<td style="text-align:right;">
50.89494
</td>
<td style="text-align:right;">
1.0000000
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
</tbody>
</table>

</div>

``` r
plot_lib_dist(Ps_obj_SIP)
```

![](05_Diff_abund_figures/beta%20div%20joint-1.png)<!-- -->

``` r
Ps_obj_SIP %>%
  scale_libraries(round = "round") ->
  Ps_obj_SIP_scaled
  
plot_lib_dist(Ps_obj_SIP_scaled)
```

![](05_Diff_abund_figures/beta%20div%20joint-2.png)<!-- -->

``` r
(mod2 <- adonis2(vegdist(otu_table(Ps_obj_SIP_scaled), method = "horn") ~ Lib.size + Site * Oxygen * Hours,
  data = as(sample_data(Ps_obj_SIP_scaled), "data.frame"),
  permutations = 999
))
```

<div class="kable-table">

<table>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
Df
</th>
<th style="text-align:right;">
SumOfSqs
</th>
<th style="text-align:right;">
R2
</th>
<th style="text-align:right;">
F
</th>
<th style="text-align:right;">
Pr(\>F)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Model
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
34.61687
</td>
<td style="text-align:right;">
0.6748572
</td>
<td style="text-align:right;">
75.49891
</td>
<td style="text-align:right;">
0.001
</td>
</tr>
<tr>
<td style="text-align:left;">
Residual
</td>
<td style="text-align:right;">
291
</td>
<td style="text-align:right;">
16.67824
</td>
<td style="text-align:right;">
0.3251428
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
Total
</td>
<td style="text-align:right;">
299
</td>
<td style="text-align:right;">
51.29511
</td>
<td style="text-align:right;">
1.0000000
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
</tbody>
</table>

</div>

``` r
(mod3 <- adonis2(vegdist(otu_table(Ps_obj_SIP_scaled), method = "horn") ~ Lib.size + Site * Oxygen * Hours * Density.zone,
  data = as(sample_data(Ps_obj_SIP_scaled), "data.frame"),
  permutations = 999
))
```

<div class="kable-table">

<table>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
Df
</th>
<th style="text-align:right;">
SumOfSqs
</th>
<th style="text-align:right;">
R2
</th>
<th style="text-align:right;">
F
</th>
<th style="text-align:right;">
Pr(\>F)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Model
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
40.35622
</td>
<td style="text-align:right;">
0.786746
</td>
<td style="text-align:right;">
65.2535
</td>
<td style="text-align:right;">
0.001
</td>
</tr>
<tr>
<td style="text-align:left;">
Residual
</td>
<td style="text-align:right;">
283
</td>
<td style="text-align:right;">
10.93889
</td>
<td style="text-align:right;">
0.213254
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
Total
</td>
<td style="text-align:right;">
299
</td>
<td style="text-align:right;">
51.29511
</td>
<td style="text-align:right;">
1.000000
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
</tbody>
</table>

</div>

``` r
Site_disp <- betadisper(vegdist(otu_table(Ps_obj_SIP_scaled), method = "horn"), get_variable(Ps_obj_SIP_scaled, "Site"))
permutest(Site_disp)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##            Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
    ## Groups      1  0.2375 0.237546 5.5678    999  0.016 *
    ## Residuals 298 12.7139 0.042664                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(Site_disp)
```

![](05_Diff_abund_figures/beta%20div%20joint-3.png)<!-- -->

``` r
Oxygen_disp <- betadisper(vegdist(otu_table(Ps_obj_SIP_scaled), method = "horn"), get_variable(Ps_obj_SIP_scaled, "Oxygen"))
permutest(Oxygen_disp)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##            Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)  
    ## Groups      1 0.03023 0.0302278 3.1681    999  0.087 .
    ## Residuals 298 2.84331 0.0095413                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(Oxygen_disp)
```

![](05_Diff_abund_figures/beta%20div%20joint-4.png)<!-- -->

``` r
Hours_disp <- betadisper(vegdist(otu_table(Ps_obj_SIP_scaled), method = "horn"), get_variable(Ps_obj_SIP_scaled, "Hours"))
permutest(Hours_disp)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##            Df Sum Sq  Mean Sq      F N.Perm Pr(>F)   
    ## Groups      4 0.2595 0.064877 4.7232    999  0.002 **
    ## Residuals 295 4.0520 0.013736                        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(Hours_disp)
```

![](05_Diff_abund_figures/beta%20div%20joint-5.png)<!-- -->

``` r
Density_disp <- betadisper(vegdist(otu_table(Ps_obj_SIP_scaled), method = "horn"), get_variable(Ps_obj_SIP_scaled, "Density.zone"))
permutest(Density_disp)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##            Df  Sum Sq Mean Sq      F N.Perm Pr(>F)    
    ## Groups      1 0.34276 0.34276 33.288    999  0.001 ***
    ## Residuals 298 3.06845 0.01030                         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(Density_disp)
```

![](05_Diff_abund_figures/beta%20div%20joint-6.png)<!-- -->

``` r
Ord <- ordinate(Ps_obj_SIP_scaled, "CAP", "horn", 
                formula =  ~ Site * Oxygen * Hours * Density.zone)
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
  scale_colour_locuszoom() +
  # scale_colour_manual(values = Gradient.colours) +
  # scale_fill_manual(values = Gradient.colours, guide = "none") +
  labs(x = sprintf("CAP1 (%s%%)", explained[1]),
  y = sprintf("CAP2 (%s%%)", explained[2])) +
  coord_fixed(ratio = sqrt(explained[2] / explained[1])) +
   theme(legend.justification = "top",
         legend.title = element_markdown(size = 11)
         ) +
  scale_size_continuous(breaks = round(c(seq(min(Ord_plt$Density..g.ml.1.), 
                                       max(Ord_plt$Density..g.ml.1.), 
                                       length.out = 5), 
                                   1), 4),
                        range = c(0.1, 5)) +
  facet_grid(Site ~ Hours) +
  # ggtitle("Joint analysis") +
  NULL

save_figure(paste0(fig.path, "Oridnation"), 
            p_ord_joint, 
            pwidth = 10, 
            pheight = 8,
            dpi = 600)

knitr::include_graphics(paste0(fig.path, "Oridnation", ".png"))
```

<img src="05_Diff_abund_figures/Oridnation.png" width="1920" />

### Differential abundance models

Now run the differential abundance models using DESeq2. We then filter
the resutls to include only ASVs with Log_2\_ fold change \>`LFC_thresh`
and significant at P\<`alpha_thresh`. Lastly, we run ‘LFC-shrinking’
based on Stephens ([2016](#ref-stephens_fdr_2016)).

``` r
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

![](05_Diff_abund_figures/DESeq2%20models%20by%20time-1.png)<!-- -->

``` r
plotMA(DESeq_res_SIP_byTime_LFC_l[[2]], ylim = c(-2,2))
```

![](05_Diff_abund_figures/DESeq2%20models%20by%20time-2.png)<!-- -->

``` r
plotMA(DESeq_res_SIP_byTime_LFC_shrink_l[[2]], ylim = c(-2,2))
```

![](05_Diff_abund_figures/DESeq2%20models%20by%20time-3.png)<!-- -->

``` r
# summarise results (lfcShrink doesn't change the values)
# map2(DESeq_res_SIP_byTime_l, print(names(DESeq_res_SIP_byTime_l)), ~summary(.x)) # summarise results
# for (i in seq(1, length(DESeq_res_SIP_byTime_l))) { # didn't manage with map
#   print(names(DESeq_res_SIP_byTime_l[i]))
#   summary(DESeq_res_SIP_byTime_l[[i]])
# }

for (i in seq(1, length(DESeq_res_SIP_byTime_LFC_l))) { # didn't manage with map
  print(names(DESeq_res_SIP_byTime_LFC_l[i]))
  summary(DESeq_res_SIP_byTime_LFC_l[[i]])
}
```

    ## [1] "Čertovo & 12 h & Anoxic & Labelled"
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
    ## [1] "Čertovo & 216 h & Anoxic & Labelled"
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
    ## [1] "Čertovo & 24 h & Anoxic & Labelled"
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
    ## [1] "Čertovo & 48 h & Anoxic & Labelled"
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
    ## [1] "Čertovo & 216 h & Anoxic & Unlabelled"
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
    ## [1] "Čertovo & 12 h & Oxic & Labelled"
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
    ## [1] "Čertovo & 24 h & Oxic & Labelled"
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
    ## [1] "Čertovo & 48 h & Oxic & Labelled"
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
    ## [1] "Čertovo & 72 h & Oxic & Labelled"
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
    ## [1] "Čertovo & 72 h & Oxic & Unlabelled"
    ## 
    ## out of 2979 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0.26 (up)    : 11, 0.37%
    ## LFC < -0.26 (down) : 0, 0%
    ## outliers [1]       : 6, 0.2%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## [1] "Plešné & 12 h & Anoxic & Labelled"
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
    ## [1] "Plešné & 216 h & Anoxic & Labelled"
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
    ## [1] "Plešné & 24 h & Anoxic & Labelled"
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
    ## [1] "Plešné & 48 h & Anoxic & Labelled"
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
    ## [1] "Plešné & 216 h & Anoxic & Unlabelled"
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
    ## [1] "Plešné & 12 h & Oxic & Labelled"
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
    ## [1] "Plešné & 24 h & Oxic & Labelled"
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
    ## [1] "Plešné & 48 h & Oxic & Labelled"
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
    ## [1] "Plešné & 72 h & Oxic & Labelled"
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
    ## [1] "Plešné & 72 h & Oxic & Unlabelled"
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

``` r
# Store labelled ASVs and save them to a file
# DESeq_res_SIP_byTime_l %>% 
#   map(., ~subset(.x, padj < alpha_thresh & log2FoldChange > LFC_thresh)) %>% 
#   map(., ~as.data.frame(.x)) %>% 
#   map(., ~rownames_to_column(.x, "ASV")) %>% 
#   bind_rows(., .id = "Comparison") %>% 
#   arrange(Comparison, desc(baseMean)) %T>% 
#   write_csv(., file = "DESeq2_byTime_a-0.05.txt") ->
#   DESeq_res_SIP_byTime_df

# Store labelled ASVs and save them to a file
DESeq_res_SIP_byTime_LFC_shrink_l %>% 
  map(., ~subset(.x, padj < alpha_thresh & log2FoldChange > LFC_thresh)) %>% 
  map(., ~as.data.frame(.x)) %>% 
  map(., ~rownames_to_column(.x, "ASV")) %>% 
  bind_rows(., .id = "Comparison") %>% 
  arrange(Comparison, desc(baseMean)) %>% 
  separate(., "Comparison" ,c("Site","Hours", "Oxygen", "Label"), sep = " & ") ->
  DESeq_res_SIP_byTime_LFC_sig_df

# grab the taxonomy of the ASVs
prune_taxa(DESeq_res_SIP_byTime_LFC_sig_df$ASV, Ps_obj_SIP) %>% 
  tax_table() %>% 
  as("data.frame") %>% 
  rownames_to_column("ASV") %>% 
  merge(., DESeq_res_SIP_byTime_LFC_sig_df, by = "ASV") %>%  # watch out: this merge recycles values!
  arrange(Site, Oxygen, Hours, baseMean) %>% 
  write_csv(., file = "DESeq2_byTime_a-0.05_LFC0-322.txt")
```

#### Inspect results

``` r
DESeq_res_SIP_byTime_LFC_sig_df %>% 
  get_variable() %>% 
  select_if(is.numeric) %>% 
  vis_value()
```

![](05_Diff_abund_figures/vis%20DES%20res-1.png)<!-- -->

``` r
DESeq_res_SIP_byTime_LFC_sig_df %>% 
  get_variable() %>% 
  select_if(is.numeric) %>% 
  vis_cor()
```

![](05_Diff_abund_figures/vis%20DES%20res-2.png)<!-- -->

#### Plot differential abundance models

``` r
# DESeq_results <- DESeq_res_SIP_byTime_LFC_shrink_l[2]
# plot_DESeq(DESeq_results, Ps_obj_SIP, plot_title = names(DESeq_results))

DESeq_plots <- map(seq(length(DESeq_res_SIP_byTime_LFC_shrink_l)), 
                        ~plot_DESeq(DESeq_res_SIP_byTime_LFC_shrink_l[.x],  
                                                Ps_obj_SIP, plot_title = names(DESeq_res_SIP_byTime_LFC_shrink_l[.x])))

Certovo_DESeq <- ((DESeq_plots[[6]] + 
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
                       theme(legend.position = "none") +
                       ylim(NA, 5)) +
                    (DESeq_plots[[5]] + 
                       theme(legend.position = "none", 
                             axis.title.y = element_blank())) + 
                    plot_layout(ncol = 2, guides = "collect") & 
                    theme(legend.position = 'bottom'))

save_figure(paste0(fig.path, "Certovo_DESeq2"), 
            Certovo_DESeq, 
            pwidth = 8, 
            pheight = 10,
            dpi = 600)

knitr::include_graphics(paste0(fig.path, "Certovo_DESeq2", ".png"))
```

<img src="05_Diff_abund_figures/Certovo_DESeq2.png" width="1536" />

``` r
Plesne_DESeq <- ((DESeq_plots[[16]] + 
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
                   theme(legend.position = 'bottom'))

save_figure(paste0(fig.path, "Plesne_DESeq2"), 
            Plesne_DESeq, 
            pwidth = 8, 
            pheight = 10,
            dpi = 600)

knitr::include_graphics(paste0(fig.path, "Plesne_DESeq2", ".png"))
```

<img src="05_Diff_abund_figures/Plesne_DESeq2.png" width="1536" />

### Plot labelled ASVs

``` r
plot_combintions <- crossing(Site = c("Čertovo", "Plešné"), 
         Oxygen = c("Oxic", "Anoxic"))

Labelled_ASVs <- map(seq(length(Ps_obj_SIP_noTime_l)), 
                     ~plot_otus_by_density(Ps_obj_SIP_noTime_l[[.x]], 
                     ASV2plot = filter(DESeq_res_SIP_byTime_LFC_sig_df, 
                                       Site == plot_combintions$Site[.x], 
                                       Oxygen == plot_combintions$Oxygen[.x])))

map(seq(length(Ps_obj_SIP_noTime_l)), 
    ~save_figure(paste0(fig.path, "Labelled_ASVs_", paste(plot_combintions[.x, ], collapse = "_")), 
                 Labelled_ASVs[[.x]], 
                 pwidth = 16, 
                 pheight = 12,
                 dpi = 600))
```

    ## [[1]]
    ## [1] "05_Diff_abund_figures/Labelled_ASVs_Čertovo_Anoxic.svgz"
    ## 
    ## [[2]]
    ## [1] "05_Diff_abund_figures/Labelled_ASVs_Čertovo_Oxic.svgz"
    ## 
    ## [[3]]
    ## [1] "05_Diff_abund_figures/Labelled_ASVs_Plešné_Anoxic.svgz"
    ## 
    ## [[4]]
    ## [1] "05_Diff_abund_figures/Labelled_ASVs_Plešné_Oxic.svgz"

``` r
plots2display <- list.files(path = paste0(fig.path), 
                    pattern = "^Labelled_ASVs_(.*).png$",
                    full.names = TRUE)

knitr::include_graphics(plots2display)
```

<img src="05_Diff_abund_figures//Labelled_ASVs_Čertovo_Anoxic.png" width="3072" /><img src="05_Diff_abund_figures//Labelled_ASVs_Čertovo_Oxic.png" width="3072" /><img src="05_Diff_abund_figures//Labelled_ASVs_Plešné_Anoxic.png" width="3072" /><img src="05_Diff_abund_figures//Labelled_ASVs_Plešné_Oxic.png" width="3072" />

#### Plot phylogenetic trees with heatmaps

``` r
c("Čertovo Oxic 12 h", "Čertovo Oxic 24 h", "Čertovo Oxic 48 h", "Čertovo Oxic 72 h", "Čertovo Anoxic 12 h", "Čertovo Anoxic 24 h", "Čertovo Anoxic 48 h", "Čertovo Anoxic 216 h", "Plešné Oxic 12 h", "Plešné Oxic 24 h", "Plešné Oxic 48 h", "Plešné Oxic 72 h", "Plešné Anoxic 12 h",  "Plešné Anoxic 24 h",  "Plešné Anoxic 48 h",  "Plešné Anoxic 216 h" )  ->
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
  mutate(across(Site_Oxygen_Hours, ~factor(., levels = col_order))) %>% 
  mutate(Site_Oxygen = factor(paste0(Site, "-", Oxygen),
                              levels = c("Plešné-Oxic", "Plešné-Anoxic", "Čertovo-Oxic", "Čertovo-Anoxic"),
                              labels = c("Pl-Ox", "Pl-Anox", "Ct-Ox", "Ct-Anox"))) %>%
  mutate(across(c("Hours"), ~factor(., 
                                    levels = c("12 h", "24 h", "48 h", "72 h", "216 h"),
                                    labels = c("12", "24", "48", "72", "216")))) ->
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
<th style="text-align:left;">
Labelled
</th>
<th style="text-align:right;">
n
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Labelled
</td>
<td style="text-align:right;">
778
</td>
</tr>
<tr>
<td style="text-align:left;">
Unlabelled
</td>
<td style="text-align:right;">
37162
</td>
</tr>
<tr>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
21148
</td>
</tr>
</tbody>
</table>

</div>

``` r
# detect taxa with NA from DESeq analysis
DESeq_res_SIP_byTime_all_df %<>% 
  filter(!is.na(Labelled)) #%>% 
  # pull(Labelled) -> 
  # bad_seqs

# remove NA taxa from PS obj
Ps_obj_SIP %>% 
  prune_taxa(setdiff(taxa_names(Ps_obj_SIP), "Seq_2375"), .) %>% # outlier
  prune_taxa(DESeq_res_SIP_byTime_all_df$ASV, .) ->
  Ps_obj_SIP4tree_plot


# Remove long name
tax_table(Ps_obj_SIP4tree_plot)[, "Order"] %<>%  str_replace_all(., "Gammaproteobacteria Incertae Sedis", "Incertae Sedis")


taxa2plot <- tibble(rank = c(rep("Class", 3), rep("Phylum", 4)), 
                    subrank = c(rep("Order", 3), rep("Class", 4)), 
                    Taxa2plot = c("Actinobacteria", 
                                  "Alphaproteobacteria", 
                                  "Gammaproteobacteria", 
                                  "Acidobacteriota",
                                  "Verrucomicrobiota",
                                  "Bacteroidota",
                                  "Firmicutes"),
                    l_rows = c(4, 5, 6, 3, 3, 3, 3),
                    pwidth = c(5, 6, 8, 3, 3, 3, 3), 
                    pheight = c(rep(10, 7)),)

tree_p_l <- map(seq(nrow(taxa2plot)), 
                ~wrap_ggtree_heatmap(Ps_obj_SIP4tree_plot,
                                     DESeq_res_SIP_byTime_all_df,
                                     rank = taxa2plot$rank[.x],
                                     subrank = taxa2plot$subrank[.x],
                                     Taxa2plot = taxa2plot$Taxa2plot[.x],
                                     l_rows = 8,
                                     pwidth = 4,
                                     pheight = 10))

trees2display <- list.files(path = paste0(fig.path), 
                    pattern = "^Tree_HM_(.*).png$",
                    full.names = TRUE)

knitr::include_graphics(trees2display)
```

<img src="05_Diff_abund_figures//Tree_HM_Acidobacteriota.png" width="768" /><img src="05_Diff_abund_figures//Tree_HM_Actinobacteria.png" width="768" /><img src="05_Diff_abund_figures//Tree_HM_Alphaproteobacteria.png" width="768" /><img src="05_Diff_abund_figures//Tree_HM_Bacteroidota.png" width="768" /><img src="05_Diff_abund_figures//Tree_HM_Firmicutes.png" width="768" /><img src="05_Diff_abund_figures//Tree_HM_Gammaproteobacteria.png" width="768" /><img src="05_Diff_abund_figures//Tree_HM_Verrucomicrobiota.png" width="768" />

``` r
all_trees <- ((tree_p_l[[1]] | tree_p_l[[2]] + guides(fill = FALSE) | tree_p_l[[3]] + guides(fill = FALSE) | tree_p_l[[4]] + guides(fill = FALSE)) / (tree_p_l[[5]] + guides(fill = FALSE) | tree_p_l[[6]] + guides(fill = FALSE) | tree_p_l[[7]] + guides(fill = FALSE) | plot_spacer())) + plot_layout(heights = c(2, 1))

save_figure(paste0(fig.path, "all_trees"), 
            all_trees, 
            pwidth = 16, 
            pheight = 18,
            dpi = 900)
```

#### How abundant were the labelled ASVs?

``` r
# Average relative abundance by time point of all labelled ASVs
DESeq_res_SIP_byTime_LFC_sig_df %>% 
  get_variable(.,c("Site", "Oxygen", "Hours")) %>% 
  distinct %>% 
  arrange(Site, Oxygen, Hours) %>% 
  remove_rownames() -> 
  variables_table 

results_table <- tibble(variables_table[1,], Total = 0)
for(i in seq(nrow(variables_table))) {
  results_table[i, ] <- sum_labelled_ASVs(Ps_obj_SIP, DESeq_res_SIP_byTime_LFC_sig_df, site = variables_table$Site[i], oxygen = variables_table$Oxygen[i], Time = variables_table$Hours[i])[1, ]
}
results_table %>%
  mutate_at(., "Total", ~ . * 100) %>%  
  rename(., Total = "Total (%)") %>% 
  kable(., digits = c(3)) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), 
                full_width = F)
```

<table class="table table-striped table-hover table-condensed table-responsive" style="color: black; width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Site
</th>
<th style="text-align:left;">
Oxygen
</th>
<th style="text-align:left;">
Hours
</th>
<th style="text-align:right;">
Total (%)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Plesne
</td>
<td style="text-align:left;">
Anoxic
</td>
<td style="text-align:left;">
12
</td>
<td style="text-align:right;">
4.041
</td>
</tr>
<tr>
<td style="text-align:left;">
Plesne
</td>
<td style="text-align:left;">
Anoxic
</td>
<td style="text-align:left;">
24
</td>
<td style="text-align:right;">
4.958
</td>
</tr>
<tr>
<td style="text-align:left;">
Plesne
</td>
<td style="text-align:left;">
Anoxic
</td>
<td style="text-align:left;">
48
</td>
<td style="text-align:right;">
3.427
</td>
</tr>
<tr>
<td style="text-align:left;">
Plesne
</td>
<td style="text-align:left;">
Oxic
</td>
<td style="text-align:left;">
12
</td>
<td style="text-align:right;">
20.867
</td>
</tr>
<tr>
<td style="text-align:left;">
Plesne
</td>
<td style="text-align:left;">
Oxic
</td>
<td style="text-align:left;">
24
</td>
<td style="text-align:right;">
26.763
</td>
</tr>
<tr>
<td style="text-align:left;">
Plesne
</td>
<td style="text-align:left;">
Oxic
</td>
<td style="text-align:left;">
48
</td>
<td style="text-align:right;">
7.469
</td>
</tr>
<tr>
<td style="text-align:left;">
Plesne
</td>
<td style="text-align:left;">
Oxic
</td>
<td style="text-align:left;">
72
</td>
<td style="text-align:right;">
31.964
</td>
</tr>
<tr>
<td style="text-align:left;">
Certovo
</td>
<td style="text-align:left;">
Anoxic
</td>
<td style="text-align:left;">
216
</td>
<td style="text-align:right;">
2.173
</td>
</tr>
<tr>
<td style="text-align:left;">
Certovo
</td>
<td style="text-align:left;">
Anoxic
</td>
<td style="text-align:left;">
24
</td>
<td style="text-align:right;">
0.127
</td>
</tr>
<tr>
<td style="text-align:left;">
Certovo
</td>
<td style="text-align:left;">
Anoxic
</td>
<td style="text-align:left;">
48
</td>
<td style="text-align:right;">
0.929
</td>
</tr>
<tr>
<td style="text-align:left;">
Certovo
</td>
<td style="text-align:left;">
Oxic
</td>
<td style="text-align:left;">
12
</td>
<td style="text-align:right;">
4.167
</td>
</tr>
<tr>
<td style="text-align:left;">
Certovo
</td>
<td style="text-align:left;">
Oxic
</td>
<td style="text-align:left;">
24
</td>
<td style="text-align:right;">
24.781
</td>
</tr>
<tr>
<td style="text-align:left;">
Certovo
</td>
<td style="text-align:left;">
Oxic
</td>
<td style="text-align:left;">
48
</td>
<td style="text-align:right;">
28.204
</td>
</tr>
<tr>
<td style="text-align:left;">
Certovo
</td>
<td style="text-align:left;">
Oxic
</td>
<td style="text-align:left;">
72
</td>
<td style="text-align:right;">
24.938
</td>
</tr>
</tbody>
</table>

#### Which labelled ASVs are shared?

``` r
Ps_obj_SIP %>% 
  subset_samples(Site == "Plešné" & Oxygen == "Oxic") %>% 
  sum_base_means() %>% 
  mutate(Site = "Plešné", Oxygen = "Oxic") -> PO

Ps_obj_SIP %>% 
  subset_samples(Site == "Plešné" & Oxygen == "Anoxic") %>% 
  sum_base_means() %>% 
  mutate(Site = "Plešné", Oxygen = "Anoxic") -> PA

Ps_obj_SIP %>% 
  subset_samples(Site == "Čertovo" & Oxygen == "Oxic") %>% 
  sum_base_means() %>% 
  mutate(Site = "Čertovo", Oxygen = "Oxic") -> CO

Ps_obj_SIP %>% 
  subset_samples(Site == "Čertovo" & Oxygen == "Anoxic") %>% 
  sum_base_means() %>% 
  mutate(Site = "Čertovo", Oxygen = "Anoxic") -> CA

rel_abunds <- rbind(PO, PA, CO, CA)

DESeq_res_SIP_byTime_LFC_sig_df <- left_join(DESeq_res_SIP_byTime_LFC_sig_df, rel_abunds, by = c("Oxygen", "Site", "ASV")) 

DESeq_res_SIP_byTime_LFC_sig_df %>% 
  mutate(ASV = str_remove(ASV, "Seq_")) %>%
  mutate(site_condition = factor(paste0(Site, " ", Oxygen))) %>% 
  mutate(site_condition = fct_relevel(site_condition, c("Čertovo Anoxic", "Plešné Anoxic", "Čertovo Oxic", "Plešné Oxic"))) %>% 
  group_by(site_condition, Oxygen, Site, ASV) %>% 
  summarise(Abundance = mean(`Mean abundance (%)`), LFC = mean(log2FoldChange)) %>% 
  group_by(Oxygen, ASV) %>% mutate(Sn = n()) %>% # shared between sites
  group_by(Site, ASV) %>% mutate(Oxn = n()) %>% # shared between oxygen
  group_by(ASV) %>% mutate(Bn = n()) %>% 
  mutate(Shared = as.factor(case_when(Bn >= 3 ~ "Both sites and conditions",
                            Oxn == 1 & Sn == 1 ~ "Unique",
                            Oxn == 1 & Sn == 2 ~ "Both sites",
                            Oxn == 2 & Sn == 1 ~ "Both conditions"))) %>% 
  mutate(Shared = fct_relevel(Shared, c("Unique", "Both sites", "Both conditions", "Both sites and conditions"))) %>% 
  ungroup() %>% 
  arrange(desc(Abundance)) %>% 
  arrange(match(ASV, unique(ASV))) %>% 
  mutate(group = rep(seq(1:5), length.out = length(Abundance), each = ceiling((length(Abundance))/5))) %>% 
  ggplot(.) + 
    geom_point(aes(x = reorder(ASV, Abundance, decreasing = TRUE), y = site_condition, size = Abundance, colour = Shared)) +
  facet_wrap(~ group, nrow = 5, scales = "free_x") + theme(
  strip.background = element_blank(),
  strip.text.x = element_blank(),
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
)
```

![](05_Diff_abund_figures/shared%20labelled%20ASV-1.png)<!-- -->

#### Calculate NTI

``` r
# grab names of labelled ASVs
DESeq_res_SIP_byTime_all_df %>% 
  as.data.frame %>% 
  # rownames_to_column("ASV") %>% 
  filter(padj < alpha_thresh) %>% 
  filter(log2FoldChange > LFC_thresh) %>% 
  pull(ASV) %>% 
  unique() ->
  Labelled_ASVs_char

# keep only labelled gradients and heavy fractions  
Ps_obj_SIP4tree_plot %>% 
  subset_samples(Label..13C. == "Labelled") %>% 
  subset_samples(Density.zone == "Heavy") %>% 
  prune_taxa(Labelled_ASVs_char, .) ->
  Ps_obj_SIP4tree_plot_labelled 

Ps_obj_SIP4tree_plot_labelled %>% 
  otu_table() %>% 
  as(., "matrix") %>% 
  {if(taxa_are_rows(Ps_obj_SIP4tree_plot_labelled)) t(.) else .}  %>% 
  {ifelse(. > 0, 1, 0)} ->
  presence_matrix

# Ps_tree <- phy_tree(Ps_obj_top)
all.equal(phy_tree(Ps_obj_SIP4tree_plot_labelled)$tip.label, colnames(otu_table(Ps_obj_SIP4tree_plot_labelled))) # just a test (should be TRUE)
```

    ## [1] TRUE

``` r
weights <- colSums(otu_table(Ps_obj_SIP4tree_plot_labelled)) / sum(colSums(otu_table(Ps_obj_SIP4tree_plot_labelled)))

# # Calculate MNTD using PhyloMeasures.
# MNTD_PM <- tibble(Identifier = rownames(otu_table(Ps_obj_SIP4tree_plot_labelled)), MNTD = mntd.query(tree = phy_tree(Ps_obj_SIP4tree_plot_labelled), matrix = presence_matrix, standardize = F, null.model = "frequency.by.richness", abundance.weights = weights)) # faster but results differ from picante 

# Using picante
MNTD_Pic <- tibble(Identifier = rownames(otu_table(Ps_obj_SIP4tree_plot_labelled)), MNTD = mntd(samp = as(otu_table(Ps_obj_SIP4tree_plot_labelled), "matrix"), dis = cophenetic(phy_tree(Ps_obj_SIP4tree_plot_labelled)), abundance.weighted = TRUE))

# # are the values equal?
# (MNTD_PM$MNTD == MNTD_Pic$MNTD)
# # Unfortunately not!

# I'll stick to picante's version
# add MNTD data to metadata 
Ps_obj_SIP4tree_plot_labelled %>%
  get_variable() %>% 
  left_join(., MNTD_Pic, "Identifier") %>% 
  mutate(Site_Oxygen = paste(Site, Oxygen)) %>% 
  mutate(Site_Oxygen_Hours = paste(Site, Oxygen, Hours)) ->
  MNTD_DF
```

``` r
# Calculate NTI using PhyloMeasures.
# NTI_PM <- tibble(Sample = rownames(otu_table(Ps_obj_top)), MNTD = -1 * mntd.query(tree = phy_tree(Ps_obj_top), matrix = otu_table(Ps_obj_top), standardize = T, null.model = "frequency.by.richness", abundance.weights = colSums(otu_table(Ps_obj_top)))) # only works on rooted trees

# using picante
NTI_Pic <- bind_cols(Identifier = rownames(otu_table(Ps_obj_SIP4tree_plot_labelled)), 
                     ses.mntd(samp = as(otu_table(Ps_obj_SIP4tree_plot_labelled), "matrix"), 
                              dis = cophenetic(phy_tree(Ps_obj_SIP4tree_plot_labelled)), 
                              null.model = "taxa.labels",
                              runs = 999, 
                              iterations = 999, 
                              abundance.weighted = TRUE))

# transform z-scores in NTI by multiplying by -1
NTI_Pic %<>% mutate(NTI = mntd.obs.z * -1)

save(NTI_Pic, file = "NTI.RData")
# load("NTI.RData")

# add MNTD data to metadata
MNTD_DF %<>%
  left_join(., NTI_Pic, "Identifier")
```

``` r
NTI_fig <- ggplot(data = MNTD_DF, 
                  aes(x = as.factor(Hours), 
                      y = NTI,  
           colour = Density..g.ml.1.)) +
  see::geom_point2(position = position_jitterdodge(), 
             size = 4, 
             alpha = 4/5)  +
  # geom_violinhalf(colour = "grey", scale = "area", trim = FALSE) +
  stat_summary(aes(group = as.factor(Hours)),
               fun.data = mean_cl_normal,
               fun.args = list(mult = 1),
               geom = "pointrange",
               colour = "black",
               alpha = 2/5,
               position = position_dodge(width = 0.5)) +
  facet_grid(Site ~ Oxygen) +
  scale_colour_viridis_c(name = "Density g ml^-1") +
  labs(x = "Incubation time (h)", 
       y = "NTI") +
  theme_bw(base_size = f_size, base_family = f_name) +
  theme(legend.title = element_markdown(), 
        text = element_text(size = f_size + 2))

save_figure(paste0(fig.path, "NTI_fig"), 
            NTI_fig, 
            # pwidth = 16, 
            # pheight = 18,
            dpi = 300)
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
 version  R version 4.4.1 (2024-06-14)
 os       Ubuntu 22.04.5 LTS
 system   x86_64, linux-gnu
 ui       X11
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Europe/Prague
 date     2024-10-22
 pandoc   2.19.2 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/ (via rmarkdown)

─ Packages ─────────────────────────────────────────────────────────────────────────────
 package              * version    date (UTC) lib source
 abind                  1.4-8      2024-09-12 [1] CRAN (R 4.4.1)
 ade4                   1.7-22     2023-02-06 [1] CRAN (R 4.2.2)
 ape                  * 5.8        2024-04-11 [1] CRAN (R 4.4.0)
 aplot                  0.2.3      2024-06-17 [1] CRAN (R 4.4.0)
 ashr                   2.2-63     2023-08-21 [1] CRAN (R 4.3.1)
 backports              1.5.0      2024-05-23 [1] CRAN (R 4.4.0)
 base64enc              0.1-3      2015-07-28 [1] CRAN (R 4.0.2)
 Biobase              * 2.60.0     2023-04-25 [1] Bioconductor
 BiocGenerics         * 0.46.0     2023-04-25 [1] Bioconductor
 BiocParallel           1.34.2     2023-05-22 [1] Bioconductor
 biomformat             1.28.0     2023-04-25 [1] Bioconductor
 Biostrings           * 2.68.1     2023-05-16 [1] Bioconductor
 bit                    4.5.0      2024-09-20 [1] CRAN (R 4.4.1)
 bit64                  4.5.2      2024-09-22 [1] CRAN (R 4.4.1)
 bitops                 1.0-8      2024-07-29 [1] CRAN (R 4.4.1)
 checkmate              2.3.2      2024-07-29 [1] CRAN (R 4.4.1)
 cli                    3.6.3      2024-06-21 [1] CRAN (R 4.4.1)
 clipr                  0.8.0      2022-02-22 [1] CRAN (R 4.1.2)
 cluster                2.1.6      2023-12-01 [1] CRAN (R 4.3.2)
 codetools              0.2-20     2024-03-31 [1] CRAN (R 4.3.3)
 colorspace             2.1-1      2024-07-26 [1] CRAN (R 4.4.1)
 commonmark             1.9.1      2024-01-30 [1] CRAN (R 4.3.2)
 crayon                 1.5.3      2024-06-20 [1] CRAN (R 4.4.1)
 data.table             1.16.0     2024-08-27 [1] CRAN (R 4.4.1)
 DelayedArray           0.26.7     2023-07-28 [1] Bioconductor
 desc                   1.4.3      2023-12-10 [1] CRAN (R 4.3.2)
 DESeq2               * 1.40.2     2023-06-23 [1] Bioconductor
 details                0.3.0      2022-03-27 [1] CRAN (R 4.1.3)
 digest                 0.6.37     2024-08-19 [1] CRAN (R 4.4.1)
 dplyr                * 1.1.4      2023-11-17 [1] CRAN (R 4.3.2)
 evaluate               1.0.0      2024-09-17 [1] CRAN (R 4.4.1)
 extrafont            * 0.19       2023-01-18 [1] CRAN (R 4.2.2)
 extrafontdb            1.0        2012-06-11 [1] CRAN (R 4.0.2)
 fansi                  1.0.6      2023-12-08 [1] CRAN (R 4.3.2)
 farver                 2.1.2      2024-05-13 [1] CRAN (R 4.4.0)
 fastmap                1.2.0      2024-05-15 [1] CRAN (R 4.4.0)
 forcats              * 1.0.0      2023-01-29 [1] CRAN (R 4.2.2)
 foreach                1.5.2      2022-02-02 [1] CRAN (R 4.1.2)
 foreign                0.8-87     2024-06-26 [1] CRAN (R 4.4.1)
 Formula                1.2-5      2023-02-24 [1] CRAN (R 4.2.3)
 fs                     1.6.4      2024-04-25 [1] CRAN (R 4.4.0)
 generics               0.1.3      2022-07-05 [1] CRAN (R 4.2.0)
 GenomeInfoDb         * 1.36.2     2023-08-25 [1] Bioconductor
 GenomeInfoDbData       1.2.10     2023-05-31 [1] Bioconductor
 GenomicRanges        * 1.52.0     2023-04-25 [1] Bioconductor
 ggfun                  0.1.6      2024-08-28 [1] CRAN (R 4.4.1)
 ggplot2              * 3.5.1      2024-04-23 [1] CRAN (R 4.4.0)
 ggplotify              0.1.2      2023-08-09 [1] CRAN (R 4.3.1)
 ggpomological        * 0.1.2      2020-08-13 [1] Github (gadenbuie/ggpomological@69f3815)
 ggrepel              * 0.9.6      2024-09-07 [1] CRAN (R 4.4.1)
 ggsci                * 3.2.0      2024-06-18 [1] CRAN (R 4.4.1)
 ggtext               * 0.1.2      2022-09-16 [1] CRAN (R 4.2.1)
 ggtree               * 3.8.2      2023-07-24 [1] Bioconductor
 glue                 * 1.7.0      2024-01-09 [1] CRAN (R 4.3.2)
 gridExtra              2.3        2017-09-09 [1] CRAN (R 4.0.2)
 gridGraphics           0.5-1      2020-12-13 [1] CRAN (R 4.0.2)
 gridtext               0.1.5      2022-09-16 [1] CRAN (R 4.2.1)
 gtable                 0.3.5      2024-04-22 [1] CRAN (R 4.4.0)
 highr                  0.11       2024-05-26 [1] CRAN (R 4.4.0)
 Hmisc                  5.1-3      2024-05-28 [1] CRAN (R 4.4.0)
 hms                    1.1.3      2023-03-21 [1] CRAN (R 4.2.3)
 htmlTable              2.4.3      2024-07-21 [1] CRAN (R 4.4.1)
 htmltools              0.5.8.1    2024-04-04 [1] CRAN (R 4.4.0)
 htmlwidgets            1.6.4      2023-12-06 [1] CRAN (R 4.3.2)
 HTSSIP               * 1.4.1      2021-01-15 [1] Github (buckleylab/HTSSIP@29ec56b)
 httr                   1.4.7      2023-08-15 [1] CRAN (R 4.3.1)
 igraph                 2.0.3      2024-03-13 [1] CRAN (R 4.3.3)
 invgamma               1.1        2017-05-07 [1] CRAN (R 4.0.2)
 IRanges              * 2.34.1     2023-06-22 [1] Bioconductor
 irlba                  2.3.5.1    2022-10-03 [1] CRAN (R 4.2.1)
 iterators              1.0.14     2022-02-05 [1] CRAN (R 4.1.2)
 jsonlite               1.8.9      2024-09-20 [1] CRAN (R 4.4.1)
 kableExtra           * 1.4.0      2024-01-24 [1] CRAN (R 4.3.2)
 knitr                  1.48       2024-07-07 [1] CRAN (R 4.4.1)
 labeling               0.4.3      2023-08-29 [1] CRAN (R 4.3.1)
 lattice              * 0.22-6     2024-03-20 [1] CRAN (R 4.3.3)
 lazyeval               0.2.2      2019-03-15 [1] CRAN (R 4.0.2)
 lifecycle              1.0.4      2023-11-07 [1] CRAN (R 4.3.1)
 locfit                 1.5-9.10   2024-06-24 [1] CRAN (R 4.4.1)
 lubridate            * 1.9.3      2023-09-27 [1] CRAN (R 4.3.1)
 magrittr             * 2.0.3      2022-03-30 [1] CRAN (R 4.1.3)
 markdown               1.13       2024-06-04 [1] CRAN (R 4.4.0)
 MASS                   7.3-61     2024-06-13 [1] CRAN (R 4.4.0)
 Matrix                 1.7-0      2024-04-26 [1] CRAN (R 4.4.0)
 MatrixGenerics       * 1.12.3     2023-07-30 [1] Bioconductor
 matrixStats          * 1.4.1      2024-09-08 [1] CRAN (R 4.4.1)
 mgcv                   1.9-1      2023-12-21 [1] CRAN (R 4.3.2)
 mixsqp                 0.3-54     2023-12-20 [1] CRAN (R 4.3.2)
 multtest               2.56.0     2023-04-25 [1] Bioconductor
 munsell                0.5.1      2024-04-01 [1] CRAN (R 4.4.0)
 nlme                 * 3.1-166    2024-08-14 [1] CRAN (R 4.4.1)
 nnet                   7.3-19     2023-05-03 [1] CRAN (R 4.2.3)
 patchwork            * 1.3.0      2024-09-16 [1] CRAN (R 4.4.1)
 permute              * 0.9-7      2022-01-27 [1] CRAN (R 4.1.2)
 phyloseq             * 1.44.0     2023-04-25 [1] Bioconductor
 picante              * 1.8.2      2020-06-10 [1] CRAN (R 4.0.2)
 pillar                 1.9.0      2023-03-22 [1] CRAN (R 4.2.3)
 pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.0.2)
 plyr                   1.8.9      2023-10-02 [1] CRAN (R 4.3.1)
 png                    0.1-8      2022-11-29 [1] CRAN (R 4.2.2)
 purrr                * 1.0.2      2023-08-10 [1] CRAN (R 4.3.1)
 R6                     2.5.1      2021-08-19 [1] CRAN (R 4.1.1)
 ragg                 * 1.3.3      2024-09-11 [1] CRAN (R 4.4.1)
 RColorBrewer         * 1.1-3      2022-04-03 [1] CRAN (R 4.1.3)
 Rcpp                   1.0.13     2024-07-17 [1] CRAN (R 4.4.1)
 RCurl                  1.98-1.16  2024-07-11 [1] CRAN (R 4.4.1)
 readr                * 2.1.5      2024-01-10 [1] CRAN (R 4.3.2)
 reshape2               1.4.4      2020-04-09 [1] CRAN (R 4.0.2)
 rhdf5                  2.44.0     2023-04-25 [1] Bioconductor
 rhdf5filters           1.12.1     2023-04-30 [1] Bioconductor
 Rhdf5lib               1.22.0     2023-04-25 [1] Bioconductor
 rlang                  1.1.4      2024-06-04 [1] CRAN (R 4.4.0)
 rmarkdown              2.28       2024-08-17 [1] CRAN (R 4.4.1)
 rpart                  4.1.23     2023-12-05 [1] CRAN (R 4.3.2)
 rstudioapi             0.16.0     2024-03-24 [1] CRAN (R 4.3.3)
 Rttf2pt1               1.3.12     2023-01-22 [1] CRAN (R 4.2.2)
 S4Arrays               1.0.6      2023-08-30 [1] Bioconductor
 S4Vectors            * 0.38.1     2023-05-02 [1] Bioconductor
 scales               * 1.3.0      2023-11-28 [1] CRAN (R 4.3.2)
 see                    0.9.0      2024-09-06 [1] CRAN (R 4.4.1)
 sessioninfo            1.2.2      2021-12-06 [1] CRAN (R 4.1.2)
 speedyseq            * 0.5.3.9021 2024-06-17 [1] Github (mikemc/speedyseq@0057652)
 SQUAREM                2021.1     2021-01-13 [1] CRAN (R 4.0.2)
 stringi                1.8.4      2024-05-06 [1] CRAN (R 4.4.0)
 stringr              * 1.5.1      2023-11-14 [1] CRAN (R 4.3.2)
 SummarizedExperiment * 1.30.2     2023-06-06 [1] Bioconductor
 survival               3.7-0      2024-06-05 [1] CRAN (R 4.4.0)
 svglite              * 2.1.3      2023-12-08 [1] CRAN (R 4.3.2)
 systemfonts            1.1.0      2024-05-15 [1] CRAN (R 4.4.0)
 textshaping            0.4.0      2024-05-24 [1] CRAN (R 4.4.0)
 tibble               * 3.2.1      2023-03-20 [1] CRAN (R 4.2.3)
 tidyr                * 1.3.1      2024-01-24 [1] CRAN (R 4.3.2)
 tidyselect             1.2.1      2024-03-11 [1] CRAN (R 4.3.3)
 tidytree               0.4.6      2023-12-12 [1] CRAN (R 4.3.2)
 tidyverse            * 2.0.0      2023-02-22 [1] CRAN (R 4.2.3)
 timechange             0.3.0      2024-01-18 [1] CRAN (R 4.3.2)
 treeio                 1.24.3     2023-07-24 [1] Bioconductor
 truncnorm              1.0-9      2023-03-20 [1] CRAN (R 4.2.3)
 tzdb                   0.4.0      2023-05-12 [1] CRAN (R 4.2.3)
 utf8                   1.2.4      2023-10-22 [1] CRAN (R 4.3.1)
 vctrs                  0.6.5      2023-12-01 [1] CRAN (R 4.3.2)
 vegan                * 2.6-8      2024-08-28 [1] CRAN (R 4.4.1)
 viridis              * 0.6.5      2024-01-29 [1] CRAN (R 4.3.2)
 viridisLite          * 0.4.2      2023-05-02 [1] CRAN (R 4.2.3)
 visdat               * 0.6.0.9000 2023-05-31 [1] Github (ropensci/visdat@9199906)
 vroom                  1.6.5      2023-12-05 [1] CRAN (R 4.3.2)
 withr                  3.0.1      2024-07-31 [1] CRAN (R 4.4.1)
 xfun                   0.47       2024-08-17 [1] CRAN (R 4.4.1)
 xml2                   1.3.6      2023-12-04 [1] CRAN (R 4.3.2)
 XVector              * 0.40.0     2023-04-25 [1] Bioconductor
 yaml                   2.3.10     2024-07-26 [1] CRAN (R 4.4.1)
 yulab.utils            0.1.7      2024-08-26 [1] CRAN (R 4.4.1)
 zlibbioc               1.46.0     2023-04-25 [1] Bioconductor

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

<div id="ref-angel_application_2018" class="csl-entry">

Angel R, Panhölzl C, Gabriel R *et al.* [Application of stable-isotope
labelling techniques for the detection of active
diazotrophs](https://doi.org/10.1111/1462-2920.13954). *Environ
Microbiol* 2018;**20**:44–61.

</div>

<div id="ref-love_moderated_2014" class="csl-entry">

Love MI, Huber W, Anders S. [Moderated estimation of fold change and
dispersion for RNA-seq data with
DESeq2](https://doi.org/10.1186/s13059-014-0550-8). *Genome Biol*
2014;**15**:550.

</div>

<div id="ref-stephens_fdr_2016" class="csl-entry">

Stephens M. [<span class="nocase">False discovery rates: a new
deal</span>](https://doi.org/10.1093/biostatistics/kxw041).
*Biostatistics* 2016;**18**:275–94.

</div>

<div id="ref-youngblut_htssip_2018" class="csl-entry">

Youngblut ND, Barnett SE, Buckley DH. [HTSSIP: An R package for analysis
of high throughput sequencing data from nucleic acid stable isotope
probing (SIP)
experiments](https://doi.org/10.1371/journal.pone.0189616). *PLOS ONE*
2018;**13**:e0189616.

</div>

</div>
