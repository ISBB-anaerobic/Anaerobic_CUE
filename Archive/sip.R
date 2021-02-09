.libPaths( c( "/home/lara/librerias" , .libPaths() ) )
setwd("./SIP")
# 1 corre los dos primeras lineas cada vez antes de hacer cuaquier cosa

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("dplyr")
BiocManager::install("HTSSIP")
BiocManager::install("phyloseq")

# 2 carga las librerias
library(dplyr)
library(HTSSIP)
library(phyloseq)
library(ggplot2)
# Crea el phylo object
otu_file<-read.delim(file = 'DADA2.seqtab_nochim.tsv', header=T)
tax_tab<-read.delim(file='DADA2.taxa_silva.tsv',header=T)
meta<-read.csv(file = 'Metadata.csv', header=T )
otu_file$OTU
row.names(otu_file)<-otu_file$OTU
View(otu_file)
otu_file <- otu_file %>% select (-OTU) 
row.names(tax_tab)<-tax_tab$OTU
View(tax_tab)
View(otu_file)
tax_tab <- tax_tab %>% select (-OTU) 
View(tax_tab)
meta$sample_name
row.names(meta)<-meta$sample_name
meta <- meta %>% select(-sample_name)
otu_file <- as.matrix(otu_file)
tax_tab <- as.matrix(tax_tab)
otu_file <- as.matrix(otu_file)
tax_tab <- as.matrix(tax_tab)
OTU = otu_table(otu_file, taxa_are_rows = TRUE)
TAX = tax_table(tax_tab)
samples = sample_data(meta)
cue <- phyloseq(OTU, TAX, samples)
cue #phylo-object
cue_rel <- transform_sample_counts(cue, function(x) x / sum(x))# convert to relative abundance
as.data.frame(otu_table(cue_rel))[1,1]
cue_rel_filt = filter_taxa(cue_rel, function(x) mean(x) > .1,TRUE)# filtra por abundancias relativas mayores a 0.1
cue_rel_filt = filter_taxa(cue_rel, function(x) sum(x) > .005, TRUE)#filtra solo las abundancias relativas mayores a 0.005
cue_rel_filt
cue_rel
params = get_treatment_params(cue_rel, c('Sample', 'Time.point', 'fraction_weight'))
params
params = get_treatment_params(cue_rel, c('Sample', 'Time.point', 'fraction_weight'))
params = dplyr::filter(params, fraction_weight!='light')
params
ex = "(fraction_weight=='light' & Sample=='${Sample}') | (fraction_weight=='${fraction_weight}' & Sample == '${Sample}')"
cue_parsed = phyloseq_subset(cue_rel, params, ex)
cue_parsed
padj_cutoff = 0.1
ncores = 2
physeq = cue_parsed[[1]]
physeq
physeq %>% sample_data %>% .$fraction_weight %>% table


# esta dando un error relacionado con la exitencia de ceros en el objeto phy. Falta arbol? falta seq-ref? 
df_l2fc = HRSIP(physeq,
design = ~fraction_weight,
padj_cutoff = padj_cutoff,
sparsity_threshold = c(0,0.15,0.3))

## grafica de *** para las fracciones pesadas y ligeras de la muestra ***
m = sample_data(cue_rel)
sample_variables(cue_rel)
cue_subset1 <- subset_samples(cue_rel, Sample=='AP9')
cue_subset1 <- merge_samples(cue_subset1, "fraction_weight")
p <- plot_bar(cue_subset1, fill = "Family")+
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  theme(legend.position="bottom")+
  guides(fill=guide_legend(nrow=5, byrow=TRUE))
print(p)
