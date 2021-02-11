---
title: "Anaerobic CUE"
subtitle: "03 Filtering out rare taxa"
author: "Roey Angel"
email: "roey.angel@bc.cas.cz"
date: "2021-02-11"
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







## Filter out rare taxa and those not classified as bacteria or archaea

### Setting general parameters:

```r
set.seed(2021)
prev_thresh <- 0.1 # remove ASVs with prevalence less than that
samples_prep_path <- "./"
data_path <- "./DADA2_pseudo/"
# Metadata_table <- "./AnCUE_Metadata_decontam.csv"
# Seq_table <- "DADA2.seqtab_nochim_decontam.tsv"
# Seq_file <- "DADA2.Seqs_decontam.fa"
Ps_file <- "Ps_obj_decontam.Rds"
Seq_file <- "DADA2.Seqs_decontam.fa"
```

### Read phyloseq object
Also remove controls

```r
readRDS(paste0(data_path, Ps_file)) ->
  Ps_obj

Ps_obj %>% 
  subset_samples(., Hours == 0) ->
  Ps_obj_t0

Ps_obj  %>% 
  subset_samples(., Hours != 0 & !is.na(Site)) -> 
  Ps_obj_SIP 
```

Combine repeated runs

```r
sample_data(Ps_obj_SIP) %<>% 
  as("data.frame") %>% 
  rownames_to_column() %>% 
  mutate(Identifier = paste(Site, Oxygen, Glucose, Label..13C., Hours, Fraction.no., 
                            sep = "_")) %>% 
  column_to_rownames("rowname") 
  
phyloseq_merge_samples(Ps_obj_SIP, "Identifier") %>% 
  filter_taxa(., function(x) sum(x) > 0, TRUE) ->
  Ps_obj_SIP_merged 

# Compute lib sizes
sample_data(Ps_obj_SIP_merged)$Lib.size <- rowSums(otu_table(Ps_obj_SIP_merged))
```

### Exploring dataset features
First let's look at the count data distribution

```r
plot_lib_size(Ps_obj_SIP_merged, x = "Fraction.no.", fill = "Oxygen", facet1 = "Site", facet2 = "Hours")
```

![](03_Filter_taxa_figures/plot abundance-1.png)<!-- -->

I will test now the effect of library size and all other experimental factors on the community composition and also plot 

```r
(mod1 <- adonis(vegdist(otu_table(Ps_obj_SIP_merged), method = "bray") ~ Site * Oxygen * Hours + Lib.size,
  data = as(sample_data(Ps_obj_SIP_merged), "data.frame"),
  permutations = 999
))
```

```
## 
## Call:
## adonis(formula = vegdist(otu_table(Ps_obj_SIP_merged), method = "bray") ~      Site * Oxygen * Hours + Lib.size, data = as(sample_data(Ps_obj_SIP_merged),      "data.frame"), permutations = 999) 
## 
## Permutation: free
## Number of permutations: 999
## 
## Terms added sequentially (first to last)
## 
##                    Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## Site                1    22.555 22.5553 314.924 0.39126  0.001 ***
## Oxygen              1     6.377  6.3767  89.033 0.11061  0.001 ***
## Hours               1     0.333  0.3332   4.652 0.00578  0.002 ** 
## Lib.size            1     4.943  4.9430  69.016 0.08575  0.001 ***
## Site:Oxygen         1     1.483  1.4830  20.706 0.02573  0.001 ***
## Site:Hours          1     0.471  0.4707   6.572 0.00817  0.002 ** 
## Oxygen:Hours        1     0.281  0.2809   3.922 0.00487  0.005 ** 
## Site:Oxygen:Hours   1     0.363  0.3631   5.070 0.00630  0.001 ***
## Residuals         291    20.842  0.0716         0.36154           
## Total             299    57.648                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
plot_lib_dist(Ps_obj_SIP_merged)
```

![](03_Filter_taxa_figures/mod abundance 1-1.png)<!-- -->

```r
plot_read_dist(Ps_obj_SIP_merged)
```

![](03_Filter_taxa_figures/mod abundance 1-2.png)<!-- -->

```r
plot_mean_SD(Ps_obj_SIP_merged)
```

![](03_Filter_taxa_figures/mod abundance 1-3.png)<!-- -->

Modelling library size shows a significant effect of read depth on the community structure, but explaining only 9% of the variance.
The reads histogram shows as expected a highly sparse and skewed sequence matrix.
The mean vs SD also shows as expected large dependency of SD on the mean reads of a sequence across all samples.

#### Taxa-based filtering 
Frist let's look at the taxonomic distribution

```r
table(tax_table(Ps_obj_SIP_merged)[, "Kingdom"], exclude = NULL)
```

```
## 
##  Archaea Bacteria 
##      185    14280
```

```r
table(tax_table(Ps_obj_SIP_merged)[, "Class"], exclude = NULL)
```

```
## 
##          Abditibacteria          Acidimicrobiia          Acidobacteriae 
##                      21                     679                    1199 
##          Actinobacteria                     AD3     Alphaproteobacteria 
##                     784                      27                    2031 
##            Anaerolineae           Armatimonadia                Babeliae 
##                       3                      84                     680 
##                 Bacilli             Bacteroidia                  BD7-11 
##                     459                     503                      16 
##         Bdellovibrionia          Blastocatellia         Campylobacteria 
##                     110                      20                       5 
##              Chlamydiae            Chloroflexia        Chthonomonadetes 
##                     583                       4                      63 
##              Clostridia          Coriobacteriia          Cyanobacteriia 
##                     276                       4                     178 
##         Dehalococcoidia              Deinococci      Desulfitobacteriia 
##                       1                      10                      12 
##           Desulfobulbia        Desulfovibrionia           Elusimicrobia 
##                       1                       2                      14 
##            Endomicrobia           Fibrobacteria          Fimbriimonadia 
##                       2                       5                      58 
##           Fusobacteriia     Gammaproteobacteria        Gemmatimonadetes 
##                      21                    2250                      25 
##             Gitt-GS-136         Gracilibacteria            Halobacteria 
##                       1                       1                       1 
##              Holophagae         Hydrogenedentia            JG30-KF-CM66 
##                       2                       1                      25 
##            Kapabacteria                  KD4-96         Ktedonobacteria 
##                      10                       1                      82 
##           Lentisphaeria             Lineage IIa           Longimicrobia 
##                       2                       5                       2 
##         Methanobacteria            Micrarchaeia              Myxococcia 
##                       1                       4                     111 
##           Negativicutes         Nitrososphaeria             Oligoflexia 
##                      13                     129                     106 
##                   OM190             Omnitrophia           Parcubacteria 
##                       1                      11                      11 
##           Phycisphaerae          Planctomycetes               Polyangia 
##                     154                    1703                     146 
##            Rhodothermia S0134 terrestrial group         Saccharimonadia 
##                       2                       1                      49 
##       Sericytochromatia                  SJA-28            Spirochaetia 
##                       1                       1                       1 
##             Subgroup 22              Subgroup 5       Syntrophobacteria 
##                       1                      12                       1 
##     Thermoanaerobaculia         Thermoleophilia          Thermoplasmata 
##                       2                     223                      50 
##                    TK10            Unclassified               vadinHA49 
##                      10                     615                      11 
##        Vampirivibrionia        Verrucomicrobiae        Vicinamibacteria 
##                     206                     599                      27
```

```r
# table(tax_table(Ps_obj)[, "Family"], exclude = NULL)
```
Now let's remove some taxa which are obvious artefacts or those which aren't bacteria or archaea

```r
kingdoms2remove <- c("", "Eukaryota", "Unclassified")
orders2remove <- c("Chloroplast")
families2remove <- c("Mitochondria")

Ps_obj_kingdoms <- tax_glom(Ps_obj_SIP_merged, "Kingdom")
Ps_obj_orders <- tax_glom(Ps_obj_SIP_merged, "Order")
Ps_obj_families <- tax_glom(Ps_obj_SIP_merged, "Family")

Ps_obj_SIP_merged_filt <- subset_taxa(Ps_obj_SIP_merged, !is.na(Phylum) &
                        !Kingdom %in% kingdoms2remove &
                      !Order %in% orders2remove &
                      !Family %in% families2remove)
```


```r
Summary_pruned <- tibble(
  Level = c("Kingdom", "Order", "Family"),
  ASVs.removed = c(
    table(tax_table(Ps_obj)[, "Kingdom"], exclude = NULL) %>% as.data.frame() %>% .[.$Var1 == "" | .$Var1 == "Eukaryota" | .$Var1 == "Unclassified", 2] %>% sum(),
    table(tax_table(Ps_obj)[, "Order"], exclude = NULL) %>% as.data.frame() %>% .[.$Var1 == "Chloroplast", 2] %>% sum(),
    table(tax_table(Ps_obj)[, "Family"], exclude = NULL) %>% as.data.frame() %>% .[.$Var1 == "Mitochondria", 2] %>% sum()
                     ),
  Seqs.removed = c(
    psmelt(Ps_obj_kingdoms) %>%
      group_by(Kingdom) %>%
      filter(Kingdom == "" |
               Kingdom == "Eukaryota" | Kingdom == "Unclassified") %>%
      summarise(sum = sum(Abundance)) %>% .$sum %>% sum(),
    psmelt(Ps_obj_orders) %>%
      group_by(Order) %>%
      filter(Order == orders2remove) %>%
      summarise(sum = sum(Abundance)) %>% .$sum %>% sum(),
    psmelt(Ps_obj_families) %>%
      group_by(Family) %>%
      filter(Family == families2remove) %>%
      summarise(sum = sum(Abundance)) %>% .$sum %>% sum()
    )
  )

Summary_pruned %>% 
  kable(., digits = c(0, 1, 0)) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)
```

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Level </th>
   <th style="text-align:right;"> ASVs.removed </th>
   <th style="text-align:right;"> Seqs.removed </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Kingdom </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Order </td>
   <td style="text-align:right;"> 136 </td>
   <td style="text-align:right;"> 4294 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Family </td>
   <td style="text-align:right;"> 138 </td>
   <td style="text-align:right;"> 15725 </td>
  </tr>
</tbody>
</table>

Removed 0.0335% of the sequences.

Now let's explore the prevalence of different taxa in the database.
Prevalence is the number of samples in which a taxa appears at least once. So "Mean prevalence" refers to in how many samples does a sequence belonging to the phylum appears on average, and "Sum prevalence" is the sum of all samples where any sequence from the taxon appears.

```r
prevdf <- apply(X = otu_table(Ps_obj_SIP_merged_filt),
                 MARGIN = ifelse(taxa_are_rows(Ps_obj_SIP_merged_filt), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(Ps_obj_SIP_merged_filt),
                      tax_table(Ps_obj_SIP_merged_filt))

prevdf %>%
  group_by(Phylum) %>%
  summarise(`Mean prevalence` = mean(Prevalence),
            `Sum prevalence` = sum(Prevalence)) ->
  Prevalence_phylum_summary

Prevalence_phylum_summary %>% 
  kable(., digits = c(0, 1, 0)) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)
```

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Phylum </th>
   <th style="text-align:right;"> Mean prevalence </th>
   <th style="text-align:right;"> Sum prevalence </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Abditibacteriota </td>
   <td style="text-align:right;"> 21.3 </td>
   <td style="text-align:right;"> 448 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Acidobacteriota </td>
   <td style="text-align:right;"> 58.2 </td>
   <td style="text-align:right;"> 73619 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Actinobacteriota </td>
   <td style="text-align:right;"> 51.3 </td>
   <td style="text-align:right;"> 86754 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Armatimonadota </td>
   <td style="text-align:right;"> 43.5 </td>
   <td style="text-align:right;"> 9178 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Bacteroidota </td>
   <td style="text-align:right;"> 20.9 </td>
   <td style="text-align:right;"> 10783 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Bdellovibrionota </td>
   <td style="text-align:right;"> 7.5 </td>
   <td style="text-align:right;"> 1625 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Campilobacterota </td>
   <td style="text-align:right;"> 3.6 </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Chloroflexi </td>
   <td style="text-align:right;"> 29.5 </td>
   <td style="text-align:right;"> 4632 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Crenarchaeota </td>
   <td style="text-align:right;"> 38.9 </td>
   <td style="text-align:right;"> 5012 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cyanobacteria </td>
   <td style="text-align:right;"> 17.5 </td>
   <td style="text-align:right;"> 4458 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Deinococcota </td>
   <td style="text-align:right;"> 1.5 </td>
   <td style="text-align:right;"> 15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Dependentiae </td>
   <td style="text-align:right;"> 11.7 </td>
   <td style="text-align:right;"> 7981 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Desulfobacterota </td>
   <td style="text-align:right;"> 6.7 </td>
   <td style="text-align:right;"> 40 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Elusimicrobiota </td>
   <td style="text-align:right;"> 13.0 </td>
   <td style="text-align:right;"> 273 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Euryarchaeota </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FCPU426 </td>
   <td style="text-align:right;"> 44.6 </td>
   <td style="text-align:right;"> 1606 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fibrobacterota </td>
   <td style="text-align:right;"> 9.0 </td>
   <td style="text-align:right;"> 45 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Firmicutes </td>
   <td style="text-align:right;"> 22.1 </td>
   <td style="text-align:right;"> 18329 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fusobacteriota </td>
   <td style="text-align:right;"> 4.4 </td>
   <td style="text-align:right;"> 93 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GAL15 </td>
   <td style="text-align:right;"> 40.1 </td>
   <td style="text-align:right;"> 281 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Gemmatimonadota </td>
   <td style="text-align:right;"> 16.7 </td>
   <td style="text-align:right;"> 468 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Halobacterota </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hydrogenedentes </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Micrarchaeota </td>
   <td style="text-align:right;"> 1.5 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Myxococcota </td>
   <td style="text-align:right;"> 46.7 </td>
   <td style="text-align:right;"> 11996 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Patescibacteria </td>
   <td style="text-align:right;"> 9.0 </td>
   <td style="text-align:right;"> 588 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Planctomycetota </td>
   <td style="text-align:right;"> 36.0 </td>
   <td style="text-align:right;"> 67790 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Proteobacteria </td>
   <td style="text-align:right;"> 28.2 </td>
   <td style="text-align:right;"> 119347 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RCP2-54 </td>
   <td style="text-align:right;"> 105.0 </td>
   <td style="text-align:right;"> 8924 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAR324 clade(Marine group B) </td>
   <td style="text-align:right;"> 33.0 </td>
   <td style="text-align:right;"> 33 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Spirochaetota </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Thermoplasmatota </td>
   <td style="text-align:right;"> 49.3 </td>
   <td style="text-align:right;"> 2464 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Unclassified </td>
   <td style="text-align:right;"> 10.3 </td>
   <td style="text-align:right;"> 1825 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Verrucomicrobiota </td>
   <td style="text-align:right;"> 29.7 </td>
   <td style="text-align:right;"> 35502 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WPS-2 </td>
   <td style="text-align:right;"> 63.1 </td>
   <td style="text-align:right;"> 8456 </td>
  </tr>
</tbody>
</table>

```r
prevdf %>%
  group_by(Order) %>%
  summarise(`Mean prevalence` = mean(Prevalence),
            `Sum prevalence` = sum(Prevalence)) ->
  Prevalence_order_summary

Prevalence_order_summary %>% 
  kable(., digits = c(0, 1, 0)) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)
```

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Order </th>
   <th style="text-align:right;"> Mean prevalence </th>
   <th style="text-align:right;"> Sum prevalence </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 0319-6G20 </td>
   <td style="text-align:right;"> 8.2 </td>
   <td style="text-align:right;"> 655 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 11-24 </td>
   <td style="text-align:right;"> 7.5 </td>
   <td style="text-align:right;"> 15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Abditibacteriales </td>
   <td style="text-align:right;"> 21.3 </td>
   <td style="text-align:right;"> 448 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Absconditabacteriales (SR1) </td>
   <td style="text-align:right;"> 2.0 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Acetobacterales </td>
   <td style="text-align:right;"> 69.1 </td>
   <td style="text-align:right;"> 16435 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Acholeplasmatales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Acidaminococcales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Acidimicrobiales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Acidobacteriales </td>
   <td style="text-align:right;"> 51.6 </td>
   <td style="text-align:right;"> 30887 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Actinomycetales </td>
   <td style="text-align:right;"> 3.5 </td>
   <td style="text-align:right;"> 38 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Aeromonadales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 540 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Alicyclobacillales </td>
   <td style="text-align:right;"> 26.6 </td>
   <td style="text-align:right;"> 186 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Alteromonadales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ardenticatenales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Armatimonadales </td>
   <td style="text-align:right;"> 43.8 </td>
   <td style="text-align:right;"> 3676 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Azospirillales </td>
   <td style="text-align:right;"> 38.2 </td>
   <td style="text-align:right;"> 992 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> B12-WMSP1 </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Babeliales </td>
   <td style="text-align:right;"> 11.7 </td>
   <td style="text-align:right;"> 7981 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Bacillales </td>
   <td style="text-align:right;"> 35.0 </td>
   <td style="text-align:right;"> 3151 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Bacteriovoracales </td>
   <td style="text-align:right;"> 1.2 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Bacteroidales </td>
   <td style="text-align:right;"> 2.0 </td>
   <td style="text-align:right;"> 117 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Bacteroidetes VC2.1 Bac22 </td>
   <td style="text-align:right;"> 1.5 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Bdellovibrionales </td>
   <td style="text-align:right;"> 6.5 </td>
   <td style="text-align:right;"> 687 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Bifidobacteriales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Blastocatellales </td>
   <td style="text-align:right;"> 1.3 </td>
   <td style="text-align:right;"> 22 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Blfdi19 </td>
   <td style="text-align:right;"> 22.5 </td>
   <td style="text-align:right;"> 45 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Bryobacterales </td>
   <td style="text-align:right;"> 76.7 </td>
   <td style="text-align:right;"> 11734 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Burkholderiales </td>
   <td style="text-align:right;"> 23.8 </td>
   <td style="text-align:right;"> 6741 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Caedibacterales </td>
   <td style="text-align:right;"> 11.7 </td>
   <td style="text-align:right;"> 82 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Caldalkalibacillales </td>
   <td style="text-align:right;"> 3.5 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Campylobacterales </td>
   <td style="text-align:right;"> 3.6 </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Candidatus Jorgensenbacteria </td>
   <td style="text-align:right;"> 3.8 </td>
   <td style="text-align:right;"> 23 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Candidatus Kaiserbacteria </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cardiobacteriales </td>
   <td style="text-align:right;"> 4.3 </td>
   <td style="text-align:right;"> 13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Catenulisporales </td>
   <td style="text-align:right;"> 97.8 </td>
   <td style="text-align:right;"> 2151 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Caulobacterales </td>
   <td style="text-align:right;"> 57.1 </td>
   <td style="text-align:right;"> 9421 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cellvibrionales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Chitinophagales </td>
   <td style="text-align:right;"> 34.7 </td>
   <td style="text-align:right;"> 5136 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Chlamydiales </td>
   <td style="text-align:right;"> 17.1 </td>
   <td style="text-align:right;"> 9989 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Chloroflexales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Christensenellales </td>
   <td style="text-align:right;"> 1.4 </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Chthoniobacterales </td>
   <td style="text-align:right;"> 26.0 </td>
   <td style="text-align:right;"> 4752 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Chthonomonadales </td>
   <td style="text-align:right;"> 19.0 </td>
   <td style="text-align:right;"> 1199 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Clostridia UCG-014 </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Clostridiales </td>
   <td style="text-align:right;"> 16.8 </td>
   <td style="text-align:right;"> 606 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Coriobacteriales </td>
   <td style="text-align:right;"> 1.2 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Corynebacteriales </td>
   <td style="text-align:right;"> 41.8 </td>
   <td style="text-align:right;"> 5059 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Coxiellales </td>
   <td style="text-align:right;"> 13.5 </td>
   <td style="text-align:right;"> 3232 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cyanobacteriales </td>
   <td style="text-align:right;"> 3.9 </td>
   <td style="text-align:right;"> 126 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cytophagales </td>
   <td style="text-align:right;"> 5.1 </td>
   <td style="text-align:right;"> 471 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Deinococcales </td>
   <td style="text-align:right;"> 1.5 </td>
   <td style="text-align:right;"> 15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Desulfitobacteriales </td>
   <td style="text-align:right;"> 23.2 </td>
   <td style="text-align:right;"> 279 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Desulfobulbales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Desulfovibrionales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Diplorickettsiales </td>
   <td style="text-align:right;"> 10.1 </td>
   <td style="text-align:right;"> 2318 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Dongiales </td>
   <td style="text-align:right;"> 2.0 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Elev-1554 </td>
   <td style="text-align:right;"> 14.5 </td>
   <td style="text-align:right;"> 29 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Elev-16S-1166 </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Elsterales </td>
   <td style="text-align:right;"> 57.4 </td>
   <td style="text-align:right;"> 19463 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Endomicrobiales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Enterobacterales </td>
   <td style="text-align:right;"> 4.5 </td>
   <td style="text-align:right;"> 150 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Entomoplasmatales </td>
   <td style="text-align:right;"> 8.2 </td>
   <td style="text-align:right;"> 41 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Erysipelotrichales </td>
   <td style="text-align:right;"> 6.7 </td>
   <td style="text-align:right;"> 100 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Exiguobacterales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FCPU453 </td>
   <td style="text-align:right;"> 32.5 </td>
   <td style="text-align:right;"> 130 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fibrobacterales </td>
   <td style="text-align:right;"> 9.0 </td>
   <td style="text-align:right;"> 45 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fimbriimonadales </td>
   <td style="text-align:right;"> 70.2 </td>
   <td style="text-align:right;"> 4072 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Flavobacteriales </td>
   <td style="text-align:right;"> 2.6 </td>
   <td style="text-align:right;"> 124 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Frankiales </td>
   <td style="text-align:right;"> 83.3 </td>
   <td style="text-align:right;"> 28396 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fusobacteriales </td>
   <td style="text-align:right;"> 4.4 </td>
   <td style="text-align:right;"> 93 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Gaiellales </td>
   <td style="text-align:right;"> 71.8 </td>
   <td style="text-align:right;"> 933 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Gammaproteobacteria Incertae Sedis </td>
   <td style="text-align:right;"> 32.0 </td>
   <td style="text-align:right;"> 6088 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Gastranaerophilales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Gemmatales </td>
   <td style="text-align:right;"> 34.4 </td>
   <td style="text-align:right;"> 49765 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Gemmatimonadales </td>
   <td style="text-align:right;"> 18.6 </td>
   <td style="text-align:right;"> 465 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Group 1.1c </td>
   <td style="text-align:right;"> 39.2 </td>
   <td style="text-align:right;"> 3335 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Haliangiales </td>
   <td style="text-align:right;"> 37.5 </td>
   <td style="text-align:right;"> 488 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Halobacterales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Holosporales </td>
   <td style="text-align:right;"> 22.9 </td>
   <td style="text-align:right;"> 1280 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hungateiclostridiaceae </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hydrogenedentiales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> I3A </td>
   <td style="text-align:right;"> 10.6 </td>
   <td style="text-align:right;"> 53 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IMCC26256 </td>
   <td style="text-align:right;"> 37.0 </td>
   <td style="text-align:right;"> 7954 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Isosphaerales </td>
   <td style="text-align:right;"> 55.7 </td>
   <td style="text-align:right;"> 5624 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> JG36-TzT-191 </td>
   <td style="text-align:right;"> 82.3 </td>
   <td style="text-align:right;"> 3375 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Kapabacteriales </td>
   <td style="text-align:right;"> 38.3 </td>
   <td style="text-align:right;"> 383 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KF-JG30-C25 </td>
   <td style="text-align:right;"> 20.0 </td>
   <td style="text-align:right;"> 60 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Kineosporiales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ktedonobacterales </td>
   <td style="text-align:right;"> 27.6 </td>
   <td style="text-align:right;"> 2232 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Lachnospirales </td>
   <td style="text-align:right;"> 5.6 </td>
   <td style="text-align:right;"> 767 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Lactobacillales </td>
   <td style="text-align:right;"> 14.4 </td>
   <td style="text-align:right;"> 1209 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Legionellales </td>
   <td style="text-align:right;"> 19.7 </td>
   <td style="text-align:right;"> 3175 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Leptolyngbyales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Lineage IV </td>
   <td style="text-align:right;"> 10.9 </td>
   <td style="text-align:right;"> 98 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Longimicrobiales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Methanobacteriales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Methanomassiliicoccales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Methylacidiphilales </td>
   <td style="text-align:right;"> 33.8 </td>
   <td style="text-align:right;"> 2163 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Methylococcales </td>
   <td style="text-align:right;"> 9.0 </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Micavibrionales </td>
   <td style="text-align:right;"> 57.0 </td>
   <td style="text-align:right;"> 57 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Micrarchaeales </td>
   <td style="text-align:right;"> 1.5 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Micrococcales </td>
   <td style="text-align:right;"> 7.4 </td>
   <td style="text-align:right;"> 891 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Micromonosporales </td>
   <td style="text-align:right;"> 1.5 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Micropepsales </td>
   <td style="text-align:right;"> 56.0 </td>
   <td style="text-align:right;"> 4645 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Microtrichales </td>
   <td style="text-align:right;"> 2.1 </td>
   <td style="text-align:right;"> 27 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mle1-27 </td>
   <td style="text-align:right;"> 11.0 </td>
   <td style="text-align:right;"> 44 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Monoglobales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Myxococcales </td>
   <td style="text-align:right;"> 22.2 </td>
   <td style="text-align:right;"> 2464 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Nitrosococcales </td>
   <td style="text-align:right;"> 37.0 </td>
   <td style="text-align:right;"> 74 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Nitrososphaerales </td>
   <td style="text-align:right;"> 60.5 </td>
   <td style="text-align:right;"> 1150 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Nitrosotaleales </td>
   <td style="text-align:right;"> 25.6 </td>
   <td style="text-align:right;"> 436 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Obscuribacterales </td>
   <td style="text-align:right;"> 21.2 </td>
   <td style="text-align:right;"> 4153 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Oceanospirillales </td>
   <td style="text-align:right;"> 2.7 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Oligoflexales </td>
   <td style="text-align:right;"> 16.6 </td>
   <td style="text-align:right;"> 216 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Omnitrophales </td>
   <td style="text-align:right;"> 25.5 </td>
   <td style="text-align:right;"> 281 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Opitutales </td>
   <td style="text-align:right;"> 38.4 </td>
   <td style="text-align:right;"> 3957 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Oscillospirales </td>
   <td style="text-align:right;"> 4.6 </td>
   <td style="text-align:right;"> 146 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Oxyphotobacteria Incertae Sedis </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Paenibacillales </td>
   <td style="text-align:right;"> 53.1 </td>
   <td style="text-align:right;"> 8012 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Paracaedibacterales </td>
   <td style="text-align:right;"> 11.6 </td>
   <td style="text-align:right;"> 792 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Pasteurellales </td>
   <td style="text-align:right;"> 14.1 </td>
   <td style="text-align:right;"> 127 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Pedosphaerales </td>
   <td style="text-align:right;"> 70.4 </td>
   <td style="text-align:right;"> 13805 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Peptococcales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Peptostreptococcales-Tissierellales </td>
   <td style="text-align:right;"> 6.4 </td>
   <td style="text-align:right;"> 256 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Phormidesmiales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Phycisphaerales </td>
   <td style="text-align:right;"> 21.0 </td>
   <td style="text-align:right;"> 524 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Pirellulales </td>
   <td style="text-align:right;"> 34.5 </td>
   <td style="text-align:right;"> 4455 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Planctomycetales </td>
   <td style="text-align:right;"> 85.8 </td>
   <td style="text-align:right;"> 2230 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Polyangiales </td>
   <td style="text-align:right;"> 70.8 </td>
   <td style="text-align:right;"> 8783 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Propionibacteriales </td>
   <td style="text-align:right;"> 10.5 </td>
   <td style="text-align:right;"> 221 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Pseudanabaenales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Pseudomonadales </td>
   <td style="text-align:right;"> 15.1 </td>
   <td style="text-align:right;"> 1782 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Pseudonocardiales </td>
   <td style="text-align:right;"> 36.8 </td>
   <td style="text-align:right;"> 699 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Pyrinomonadales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Reyranellales </td>
   <td style="text-align:right;"> 29.3 </td>
   <td style="text-align:right;"> 352 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rhizobiales </td>
   <td style="text-align:right;"> 62.3 </td>
   <td style="text-align:right;"> 23562 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rhodobacterales </td>
   <td style="text-align:right;"> 2.6 </td>
   <td style="text-align:right;"> 84 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rhodospirillales </td>
   <td style="text-align:right;"> 13.8 </td>
   <td style="text-align:right;"> 1658 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rhodothermales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rickettsiales </td>
   <td style="text-align:right;"> 15.6 </td>
   <td style="text-align:right;"> 1485 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S-BQ2-57 soil group </td>
   <td style="text-align:right;"> 11.8 </td>
   <td style="text-align:right;"> 177 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S085 </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Saccharimonadales </td>
   <td style="text-align:right;"> 11.0 </td>
   <td style="text-align:right;"> 537 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Salinisphaerales </td>
   <td style="text-align:right;"> 48.2 </td>
   <td style="text-align:right;"> 1013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAR11 clade </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SBR1031 </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Silvanigrellales </td>
   <td style="text-align:right;"> 4.8 </td>
   <td style="text-align:right;"> 62 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Solibacterales </td>
   <td style="text-align:right;"> 89.8 </td>
   <td style="text-align:right;"> 10771 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Solirubrobacterales </td>
   <td style="text-align:right;"> 75.9 </td>
   <td style="text-align:right;"> 15946 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sphingobacteriales </td>
   <td style="text-align:right;"> 29.4 </td>
   <td style="text-align:right;"> 4468 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sphingomonadales </td>
   <td style="text-align:right;"> 4.6 </td>
   <td style="text-align:right;"> 705 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Spirochaetales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Staphylococcales </td>
   <td style="text-align:right;"> 21.6 </td>
   <td style="text-align:right;"> 410 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Steroidobacterales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Streptomycetales </td>
   <td style="text-align:right;"> 63.1 </td>
   <td style="text-align:right;"> 883 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Streptosporangiales </td>
   <td style="text-align:right;"> 2.9 </td>
   <td style="text-align:right;"> 105 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Subgroup 12 </td>
   <td style="text-align:right;"> 10.8 </td>
   <td style="text-align:right;"> 43 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Subgroup 13 </td>
   <td style="text-align:right;"> 61.6 </td>
   <td style="text-align:right;"> 1725 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Subgroup 15 </td>
   <td style="text-align:right;"> 5.1 </td>
   <td style="text-align:right;"> 72 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Subgroup 2 </td>
   <td style="text-align:right;"> 73.8 </td>
   <td style="text-align:right;"> 16985 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Subgroup 7 </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Synechococcales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Syntrophobacterales </td>
   <td style="text-align:right;"> 35.0 </td>
   <td style="text-align:right;"> 35 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Tepidisphaerales </td>
   <td style="text-align:right;"> 35.1 </td>
   <td style="text-align:right;"> 4533 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Thalassobaculales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Thermicanales </td>
   <td style="text-align:right;"> 5.0 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Thermoactinomycetales </td>
   <td style="text-align:right;"> 30.0 </td>
   <td style="text-align:right;"> 1951 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Thermoanaerobaculales </td>
   <td style="text-align:right;"> 15.5 </td>
   <td style="text-align:right;"> 31 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Thermomicrobiales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TSBb06 </td>
   <td style="text-align:right;"> 10.2 </td>
   <td style="text-align:right;"> 61 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Unclassified </td>
   <td style="text-align:right;"> 31.5 </td>
   <td style="text-align:right;"> 56145 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Vampirovibrionales </td>
   <td style="text-align:right;"> 17.3 </td>
   <td style="text-align:right;"> 104 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Veillonellales-Selenomonadales </td>
   <td style="text-align:right;"> 6.7 </td>
   <td style="text-align:right;"> 80 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Verrucomicrobiales </td>
   <td style="text-align:right;"> 10.1 </td>
   <td style="text-align:right;"> 375 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Vibrionales </td>
   <td style="text-align:right;"> 1.7 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Vicinamibacterales </td>
   <td style="text-align:right;"> 31.3 </td>
   <td style="text-align:right;"> 846 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Victivallales </td>
   <td style="text-align:right;"> 1.0 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WD260 </td>
   <td style="text-align:right;"> 158.8 </td>
   <td style="text-align:right;"> 3335 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Xanthomonadales </td>
   <td style="text-align:right;"> 29.6 </td>
   <td style="text-align:right;"> 2336 </td>
  </tr>
</tbody>
</table>

Based on that I'll remove all orders with a sum prevalence of under 5% (15) of all samples

```r
Prevalence_order_summary %>% 
  filter(`Sum prevalence` < (0.05 * nsamples(Ps_obj_SIP_merged_filt))) %>% 
  dplyr::select(Order) %>% 
  map(as.character) %>% 
  unlist() ->
  filterOrder

Ps_obj_SIP_merged_filt2 <- subset_taxa(Ps_obj_SIP_merged_filt, !Order %in% filterOrder)

sample_data(Ps_obj_SIP_merged_filt2)$Lib.size <- rowSums(otu_table(Ps_obj_SIP_merged_filt2))
print(Ps_obj_SIP_merged_filt)
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:          [ 14199 taxa and 300 samples ]:
## sample_data() Sample Data:        [ 300 samples by 29 sample variables ]:
## tax_table()   Taxonomy Table:     [ 14199 taxa by 6 taxonomic ranks ]:
## taxa are columns
```

```r
print(Ps_obj_SIP_merged_filt2)
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:          [ 14099 taxa and 300 samples ]:
## sample_data() Sample Data:        [ 300 samples by 29 sample variables ]:
## tax_table()   Taxonomy Table:     [ 14099 taxa by 6 taxonomic ranks ]:
## taxa are columns
```

This removed 100 or 1% of the ESVs, and 0.046% of the reads.

Plot general prevalence features of the phyla

```r
# Subset to the remaining phyla
prevdf_phylum_filt <- subset(prevdf, 
                             Phylum %in% get_taxa_unique(Ps_obj_SIP_merged_filt2, 
                                                         "Phylum"))
ggplot(prevdf_phylum_filt,
       aes(TotalAbundance, 
           Prevalence / nsamples(Ps_obj_SIP_merged_filt2), 
           color = Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05,
             alpha = 0.5,
             linetype = 2) + 
  geom_point2(size = 2, alpha = 0.7) +
  scale_x_log10() +  
  xlab("Total Abundance") + 
  ylab("Prevalence [Frac. Samples]") +
  facet_wrap( ~ Phylum) + 
  theme(legend.position = "none")
```

![](03_Filter_taxa_figures/prevalence phylum-1.png)<!-- -->

Plot general prevalence features of the top 20 orders

```r
# Subset to the remaining phyla
prevdf_order_filt <- subset(prevdf, 
                            Order %in% get_taxa_unique(Ps_obj_SIP_merged_filt2, "Order"))

# grab the top 30 most abundant orders
prevdf_order_filt %>% 
  group_by(Order) %>%
  summarise(Combined.abundance = sum(TotalAbundance)) %>% 
  arrange(desc(Combined.abundance)) %>% 
  .[1:30, "Order"]  ->
  Orders2plot

prevdf_order_filt2 <- subset(prevdf,
                             Order %in% Orders2plot$Order)

ggplot(prevdf_order_filt2,
       aes(TotalAbundance,
           Prevalence / nsamples(Ps_obj_SIP_merged_filt2), 
           color = Order)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05,
             alpha = 0.5,
             linetype = 2) + 
  geom_point2(size = 2, alpha = 0.7) +
  scale_x_log10() +  
  xlab("Total Abundance") +
  ylab("Prevalence [Frac. Samples]") +
  facet_wrap( ~ Order) + 
  theme(legend.position = "none")
```

![](03_Filter_taxa_figures/prevalence order-1.png)<!-- -->

#### Unsupervised filtering by prevalence
I'll remove all sequences which appear in less than 5% of the samples

```r
# Define prevalence threshold as 5% of total samples
prevalenceThreshold <- prev_thresh * nsamples(Ps_obj_SIP_merged_filt2)
prevalenceThreshold
```

```
## [1] 30
```

```r
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa <-
  row.names(prevdf_phylum_filt)[(prevdf_phylum_filt$Prevalence >= prevalenceThreshold)]
Ps_obj_SIP_merged_filt3 <- prune_taxa(keepTaxa, Ps_obj_SIP_merged_filt2)

sample_data(Ps_obj_SIP_merged_filt3)$Lib.size <- rowSums(otu_table(Ps_obj_SIP_merged_filt3))
print(Ps_obj_SIP_merged_filt2)
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:          [ 14099 taxa and 300 samples ]:
## sample_data() Sample Data:        [ 300 samples by 29 sample variables ]:
## tax_table()   Taxonomy Table:     [ 14099 taxa by 6 taxonomic ranks ]:
## taxa are columns
```

```r
print(Ps_obj_SIP_merged_filt3)
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:          [ 3693 taxa and 300 samples ]:
## sample_data() Sample Data:        [ 300 samples by 29 sample variables ]:
## tax_table()   Taxonomy Table:     [ 3693 taxa by 6 taxonomic ranks ]:
## taxa are columns
```
This removed 10406 or 74% of the ESVs! 

However all these removed ESVs accounted for only: 

```r
prevdf_phylum_filt %>% 
  arrange(., Prevalence) %>% 
  group_by(Prevalence > prevalenceThreshold) %>% 
  summarise(Abundance = sum(TotalAbundance)) %>%
  mutate(`Rel. Ab.` = percent(Abundance / sum(Abundance), accuracy = 0.01)) %>% 
  kable(., digits = c(0, 1, 0)) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)
```

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Prevalence &gt; prevalenceThreshold </th>
   <th style="text-align:right;"> Abundance </th>
   <th style="text-align:left;"> Rel. Ab. </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:right;"> 303010 </td>
   <td style="text-align:left;"> 2.37% </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:right;"> 12498417 </td>
   <td style="text-align:left;"> 97.63% </td>
  </tr>
</tbody>
</table>
So it's fine to remove them.

Test again the effect of library size and all other experimental factors on the community composition after filtering

```r
(mod1 <- adonis(vegdist(otu_table(Ps_obj_SIP_merged_filt3), method = "bray") ~ Site * Oxygen * Hours + Lib.size,
  data = as(sample_data(Ps_obj_SIP_merged_filt3), "data.frame"),
  permutations = 999
))
```

```
## 
## Call:
## adonis(formula = vegdist(otu_table(Ps_obj_SIP_merged_filt3),      method = "bray") ~ Site * Oxygen * Hours + Lib.size, data = as(sample_data(Ps_obj_SIP_merged_filt3),      "data.frame"), permutations = 999) 
## 
## Permutation: free
## Number of permutations: 999
## 
## Terms added sequentially (first to last)
## 
##                    Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## Site                1    22.638 22.6379  342.15 0.40294  0.001 ***
## Oxygen              1     6.358  6.3575   96.09 0.11316  0.001 ***
## Hours               1     0.322  0.3223    4.87 0.00574  0.005 ** 
## Lib.size            1     5.069  5.0686   76.61 0.09022  0.001 ***
## Site:Oxygen         1     1.475  1.4747   22.29 0.02625  0.001 ***
## Site:Hours          1     0.458  0.4576    6.92 0.00814  0.001 ***
## Oxygen:Hours        1     0.264  0.2643    3.99 0.00470  0.007 ** 
## Site:Oxygen:Hours   1     0.345  0.3448    5.21 0.00614  0.002 ** 
## Residuals         291    19.254  0.0662         0.34270           
## Total             299    56.181                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
plot_lib_dist(Ps_obj_SIP_merged_filt3)
```

![](03_Filter_taxa_figures/mod abundance 2-1.png)<!-- -->

```r
plot_read_dist(Ps_obj_SIP_merged_filt3)
```

![](03_Filter_taxa_figures/mod abundance 2-2.png)<!-- -->

```r
plot_mean_SD(Ps_obj_SIP_merged_filt3)
```

![](03_Filter_taxa_figures/mod abundance 2-3.png)<!-- -->

#### Save filtered phyloseq object

```r
saveRDS(Ps_obj_SIP_merged_filt3, file = paste0(data_path, str_remove(Ps_file, ".Rds"), "_filt3.Rds"))
readDNAStringSet(paste0(data_path, Seq_file)) %>% 
  .[taxa_names(Ps_obj_SIP_merged_filt3)] %>%  
  writeXStringSet(., filepath = paste0(data_path, str_remove(Seq_file, ".fa*"), "_filtered.fa"), format = "fasta", width = 1000)
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
 date     2021-02-11                  

─ Packages ─────────────────────────────────────────────────────────────────────────────
 package        * version    date       lib source                           
 ade4             1.7-16     2020-10-28 [1] CRAN (R 4.0.2)                   
 affy             1.66.0     2020-04-27 [1] Bioconductor                     
 affyio           1.58.0     2020-04-27 [1] Bioconductor                     
 ape              5.4-1      2020-08-13 [1] CRAN (R 4.0.2)                   
 assertthat       0.2.1      2019-03-21 [1] CRAN (R 4.0.2)                   
 backports        1.2.1      2020-12-09 [1] CRAN (R 4.0.2)                   
 bayestestR       0.8.2      2021-01-26 [1] CRAN (R 4.0.3)                   
 Biobase        * 2.48.0     2020-04-27 [1] Bioconductor                     
 BiocGenerics   * 0.34.0     2020-04-27 [1] Bioconductor                     
 BiocManager      1.30.10    2019-11-16 [1] CRAN (R 4.0.2)                   
 biomformat       1.16.0     2020-04-27 [1] Bioconductor                     
 Biostrings     * 2.56.0     2020-04-27 [1] Bioconductor                     
 broom            0.7.4      2021-01-29 [1] CRAN (R 4.0.3)                   
 cellranger       1.1.0      2016-07-27 [1] CRAN (R 4.0.2)                   
 cli              2.3.0      2021-01-31 [1] CRAN (R 4.0.3)                   
 clipr            0.7.1      2020-10-08 [1] CRAN (R 4.0.2)                   
 cluster          2.1.0      2019-06-19 [1] CRAN (R 4.0.2)                   
 codetools        0.2-18     2020-11-04 [1] CRAN (R 4.0.2)                   
 colorspace       2.0-0      2020-11-11 [1] CRAN (R 4.0.2)                   
 crayon           1.4.1      2021-02-08 [1] CRAN (R 4.0.3)                   
 data.table       1.13.6     2020-12-30 [1] CRAN (R 4.0.2)                   
 DBI              1.1.1      2021-01-15 [1] CRAN (R 4.0.3)                   
 dbplyr           2.1.0      2021-02-03 [1] CRAN (R 4.0.3)                   
 desc             1.2.0      2018-05-01 [1] CRAN (R 4.0.2)                   
 details          0.2.1      2020-01-12 [1] CRAN (R 4.0.2)                   
 digest           0.6.27     2020-10-24 [1] CRAN (R 4.0.2)                   
 dplyr          * 1.0.4      2021-02-02 [1] CRAN (R 4.0.3)                   
 effectsize       0.4.3      2021-01-18 [1] CRAN (R 4.0.3)                   
 ellipsis         0.3.1      2020-05-15 [1] CRAN (R 4.0.2)                   
 evaluate         0.14       2019-05-28 [1] CRAN (R 4.0.2)                   
 extrafont      * 0.17       2014-12-08 [1] CRAN (R 4.0.2)                   
 extrafontdb      1.0        2012-06-11 [1] CRAN (R 4.0.2)                   
 farver           2.0.3      2020-01-16 [1] CRAN (R 4.0.2)                   
 forcats        * 0.5.1      2021-01-27 [1] CRAN (R 4.0.3)                   
 foreach          1.5.1      2020-10-15 [1] CRAN (R 4.0.2)                   
 fs               1.5.0      2020-07-31 [1] CRAN (R 4.0.2)                   
 gdtools          0.2.3      2021-01-06 [1] CRAN (R 4.0.2)                   
 generics         0.1.0      2020-10-31 [1] CRAN (R 4.0.2)                   
 ggplot2        * 3.3.3      2020-12-30 [1] CRAN (R 4.0.2)                   
 ggridges         0.5.3      2021-01-08 [1] CRAN (R 4.0.2)                   
 glue             1.4.2      2020-08-27 [1] CRAN (R 4.0.2)                   
 gtable           0.3.0      2019-03-25 [1] CRAN (R 4.0.2)                   
 haven            2.3.1      2020-06-01 [1] CRAN (R 4.0.2)                   
 hexbin           1.28.2     2021-01-08 [1] CRAN (R 4.0.2)                   
 highr            0.8        2019-03-20 [1] CRAN (R 4.0.2)                   
 hms              1.0.0      2021-01-13 [1] CRAN (R 4.0.3)                   
 htmltools        0.5.1.1    2021-01-22 [1] CRAN (R 4.0.3)                   
 httr             1.4.2      2020-07-20 [1] CRAN (R 4.0.2)                   
 igraph           1.2.6      2020-10-06 [1] CRAN (R 4.0.2)                   
 insight          0.12.0     2021-01-14 [1] CRAN (R 4.0.3)                   
 IRanges        * 2.22.2     2020-05-21 [1] Bioconductor                     
 iterators        1.0.13     2020-10-15 [1] CRAN (R 4.0.2)                   
 jsonlite         1.7.2      2020-12-09 [1] CRAN (R 4.0.2)                   
 kableExtra     * 1.3.1      2020-10-22 [1] CRAN (R 4.0.2)                   
 knitr            1.31       2021-01-27 [1] CRAN (R 4.0.3)                   
 labeling         0.4.2      2020-10-20 [1] CRAN (R 4.0.2)                   
 lattice        * 0.20-41    2020-04-02 [1] CRAN (R 4.0.2)                   
 lifecycle        0.2.0      2020-03-06 [1] CRAN (R 4.0.2)                   
 limma            3.44.3     2020-06-12 [1] Bioconductor                     
 lubridate        1.7.9.2    2020-11-13 [1] CRAN (R 4.0.2)                   
 magrittr       * 2.0.1      2020-11-17 [1] CRAN (R 4.0.2)                   
 MASS             7.3-53     2020-09-09 [1] CRAN (R 4.0.2)                   
 Matrix           1.3-2      2021-01-06 [1] CRAN (R 4.0.2)                   
 mgcv             1.8-33     2020-08-27 [1] CRAN (R 4.0.2)                   
 modelr           0.1.8      2020-05-19 [1] CRAN (R 4.0.2)                   
 multtest         2.44.0     2020-04-27 [1] Bioconductor                     
 munsell          0.5.0      2018-06-12 [1] CRAN (R 4.0.2)                   
 nlme             3.1-152    2021-02-04 [1] CRAN (R 4.0.3)                   
 parameters       0.11.0     2021-01-15 [1] CRAN (R 4.0.3)                   
 permute        * 0.9-5      2019-03-12 [1] CRAN (R 4.0.2)                   
 phyloseq       * 1.32.0     2020-04-27 [1] Bioconductor                     
 pillar           1.4.7      2020-11-20 [1] CRAN (R 4.0.2)                   
 pkgconfig        2.0.3      2019-09-22 [1] CRAN (R 4.0.2)                   
 plyr             1.8.6      2020-03-03 [1] CRAN (R 4.0.2)                   
 png              0.1-7      2013-12-03 [1] CRAN (R 4.0.2)                   
 preprocessCore   1.50.0     2020-04-27 [1] Bioconductor                     
 prettyunits      1.1.1      2020-01-24 [1] CRAN (R 4.0.2)                   
 progress         1.2.2      2019-05-16 [1] CRAN (R 4.0.2)                   
 purrr          * 0.3.4      2020-04-17 [1] CRAN (R 4.0.2)                   
 R6               2.5.0      2020-10-28 [1] CRAN (R 4.0.2)                   
 RColorBrewer     1.1-2      2014-12-07 [1] CRAN (R 4.0.2)                   
 Rcpp             1.0.6      2021-01-15 [1] CRAN (R 4.0.3)                   
 readr          * 1.4.0      2020-10-05 [1] CRAN (R 4.0.2)                   
 readxl           1.3.1      2019-03-13 [1] CRAN (R 4.0.2)                   
 reprex           1.0.0      2021-01-27 [1] CRAN (R 4.0.3)                   
 reshape2         1.4.4      2020-04-09 [1] CRAN (R 4.0.2)                   
 rhdf5            2.32.4     2020-10-05 [1] Bioconductor                     
 Rhdf5lib         1.10.1     2020-07-09 [1] Bioconductor                     
 rlang            0.4.10     2020-12-30 [1] CRAN (R 4.0.2)                   
 rmarkdown        2.6        2020-12-14 [1] CRAN (R 4.0.2)                   
 rprojroot        2.0.2      2020-11-15 [1] CRAN (R 4.0.2)                   
 rstudioapi       0.13       2020-11-12 [1] CRAN (R 4.0.2)                   
 Rttf2pt1         1.3.8      2020-01-10 [1] CRAN (R 4.0.2)                   
 rvest            0.3.6      2020-07-25 [1] CRAN (R 4.0.2)                   
 S4Vectors      * 0.26.1     2020-05-16 [1] Bioconductor                     
 scales         * 1.1.1      2020-05-11 [1] CRAN (R 4.0.2)                   
 see            * 0.6.2      2021-02-04 [1] CRAN (R 4.0.3)                   
 sessioninfo      1.1.1      2018-11-05 [1] CRAN (R 4.0.2)                   
 speedyseq      * 0.5.3.9001 2020-10-27 [1] Github (mikemc/speedyseq@8daed32)
 stringi          1.5.3      2020-09-09 [1] CRAN (R 4.0.2)                   
 stringr        * 1.4.0      2019-02-10 [1] CRAN (R 4.0.2)                   
 survival         3.2-7      2020-09-28 [1] CRAN (R 4.0.2)                   
 svglite        * 1.2.3.2    2020-07-07 [1] CRAN (R 4.0.2)                   
 systemfonts      1.0.1      2021-02-09 [1] CRAN (R 4.0.3)                   
 tibble         * 3.0.6      2021-01-29 [1] CRAN (R 4.0.3)                   
 tidyr          * 1.1.2      2020-08-27 [1] CRAN (R 4.0.2)                   
 tidyselect       1.1.0      2020-05-11 [1] CRAN (R 4.0.2)                   
 tidyverse      * 1.3.0      2019-11-21 [1] CRAN (R 4.0.2)                   
 vctrs            0.3.6      2020-12-17 [1] CRAN (R 4.0.2)                   
 vegan          * 2.5-7      2020-11-28 [1] CRAN (R 4.0.2)                   
 viridisLite      0.3.0      2018-02-01 [1] CRAN (R 4.0.2)                   
 vsn            * 3.56.0     2020-04-27 [1] Bioconductor                     
 webshot          0.5.2      2019-11-22 [1] CRAN (R 4.0.2)                   
 withr            2.4.1      2021-01-26 [1] CRAN (R 4.0.3)                   
 xfun             0.20       2021-01-06 [1] CRAN (R 4.0.2)                   
 xml2             1.3.2      2020-04-23 [1] CRAN (R 4.0.2)                   
 XVector        * 0.28.0     2020-04-27 [1] Bioconductor                     
 yaml             2.2.1      2020-02-01 [1] CRAN (R 4.0.2)                   
 zlibbioc         1.34.0     2020-04-27 [1] Bioconductor                     

[1] /home/angel/R/library
[2] /usr/local/lib/R/site-library
[3] /usr/lib/R/site-library
[4] /usr/lib/R/library

```

</details>
<br>

## References
