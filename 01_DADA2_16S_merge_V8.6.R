################################################################################
# DADA2 pipeline - 16S                                                         #
# This implementation of the DADA2 pipeline merges 2 read files per sample:    #
# *_R1_* and *_R2_*.                                                           #
# Optimised for both MiSeq and MiniSeq runs                                    #
# Usage: preferebly run on a server using run_DADA2.sh                         #
#                                                                              #
# V8.6                                                                         #
# Roey Angel - August 2020                                                     #
# roey.angel@bc.cas.cz                                                         #
# Based on http://benjjneb.github.io/dada2/tutorial.html                       #
################################################################################


############################################################
# TO DO                                                    #
# 1. Which mock seqs match                                 #
############################################################

# Load libraries ----------------------------
Sys.setenv(R_LIBS_USER = "~/R/library") #
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())) # Uncomment if you have no write access to R path
# library("tidyverse") # "No used functions found"
# library("broom") # "No used functions found"
library("magrittr") # %>% %<>%
# library("scales") # "No used functions found"
# library("svglite") # "No used functions found"
library("seqinr") # write.fasta 
library("Biostrings") # strsplit sort as.data.frame %in% head reverseComplement nchar duplicated as.list
library("dada2") # plotQualityProfile filterAndTrim learnErrors plotErrors derepFastq dada mergePairs makeSequenceTable getSequences removeBimeraDenovo getUniques assignTaxonomy addSpecies
library("ShortRead") # %in% reverseComplement readFasta sapply sread
library("microseq") # reverseComplement readFasta
# library("DECIPHER") 

# Set general parameters ----------------------------
# set default arguments
pooling <- FALSE # defualt - single sample analysis
data_path <- "../" # directory containing the zipped fastq files after untaring
db_path <- "~/Resources/"
mock_id <- ""
# file_suffix <- "/*(.*)_L001_R1_001_noPrimers.fastq.gz$"
file_suffix <- "/*(.*)_R1_(.*).fastq.gz$"
# file_prefix <- "^Angel[0-9]+-(.*)"

# change them if arguments were supplied 
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  writeLines("No args supplied; using default values")
  writeLines(paste0("pooling=", pooling))
  writeLines(paste0("data_path=", data_path))
  writeLines(paste0("db_path=", db_path))
  writeLines(paste0("mock_id=", mock_id))
} else if (length(args) == 1) {
  writeLines("Pooling type supplied")
  pooling <- args[1]
  writeLines(paste0("pooling=", pooling))
  writeLines(paste0("data_path=", data_path))
  writeLines(paste0("db_path=", db_path))
  writeLines(paste0("mock_id=", mock_id))
} else if (length(args) == 2) {
  writeLines("Pooling type and data path supplied")
  pooling <- args[1]
  data_path <- args[2]
  writeLines(paste0("pooling=", pooling))
  writeLines(paste0("data_path=", data_path))
  writeLines(paste0("db_path=", db_path))
  writeLines(paste0("mock_id=", mock_id))
} else if (length(args) == 3) {
  writeLines("Pooling type, data path and db path supplied")
  pooling <- args[1]
  data_path <- args[2]
  db_path <- args[3]
  writeLines(paste0("pooling=", pooling))
  writeLines(paste0("data_path=", data_path))
  writeLines(paste0("db_path=", db_path))
  writeLines(paste0("mock_id=", mock_id))
} else if (length(args) == 4) {
  writeLines("Pooling type, data path, db path and mock id args supplied")
  pooling <- args[1]
  data_path <- args[2]
  db_path <- args[3]
  mock_id <- args[4]
  writeLines(paste0("pooling=", pooling))
  writeLines(paste0("data_path=", data_path))
  writeLines(paste0("db_path=", db_path))
  writeLines(paste0("mock_id=", mock_id))
} else {
  stop("Only up to 3 arguments are allowed.", call. = FALSE)
}

if (pooling != "pseudo") pooling <- as.logical(pooling)

# list.files(db_path)
# list.files(db_path)
silva_db <- "SILVA/silva_nr99_v138_wSpecies_train_set.fa.gz"
writeLines(paste("Using SILVA file:", unlist(strsplit(silva_db, "/"))[2]))
#silva_sp_db <- "SILVA/silva_species_assignment_v138.fa.gz"
gtdb_db <- "GTDB/GTDB_bac120_arc122_ssu_r95_fullTaxo.fa.gz"
writeLines(paste("Using GTDB file:", unlist(strsplit(gtdb_db, "/"))[2]))
# silva_decipher_db <- "SILVA_SSU_r132_March2018.RData"
# GTDB_decipher_db <- "GTDB_r86-mod_September2018.RData"
mock_ref_seqs <- "mock/ZymoBIOMICS.STD.refseq.v2.ssrRNA_16S_515F-806R.pcr.fasta"
writeLines(paste("Using mock community file:", unlist(strsplit(mock_ref_seqs, "/"))[2]))
# theme_set(theme_bw(base_size = 14, base_family = "DejaVu Sans"))
amp_size <- 253
writeLines(paste("Assuming amplicon length:", amp_size))
# truncate <- c(130, 130) # set trimming points for read 1 and read 2 NOT NEEDED FOR MINISEQ

# Filtering and trimming ----------------------------------------------------------------------
writeLines("\nFilter and trim raw sequences:")
fastqs <- list.files(data_path, pattern = "*.fastq.gz$", full.names = TRUE) # list gz fastq files only
if (isEmpty(fastqs)) { # make sure files are there
  stop("could not find any fastq.gz files in ", data_path)
}

fastqs <-
  sort(fastqs) # Sort ensures forward/reverse reads are in same order
fastq1s <-
  fastqs[str_detect(fastqs, file_suffix)] # Just the read 1 files
fastq2s <-
  fastqs[str_detect(fastqs, str_replace(file_suffix, "_R1_", "_R2_"))] # Just the read 2 files

# Get sample names, assuming files named as so: SAMPLENAME_XXX.fastq
fastq1s %>% 
  # str_replace(., file.path(data_path, "/*(.*)\\.[0-9]\\.fastq\\.gz$"), "\\1") %>% # or
  str_replace(., file.path(data_path, file_suffix), "\\1") %>% # and maybe 
  str_replace(., paste0(data_path, "/*(.*)\\.[0-9]\\.fastq\\.gz$"), "\\1") -> 
  # str_replace(.,  file_prefix, "\\1") ->
  sample_names
# writeLines(sample_names)
if (isEmpty(sample_names)) {
  stop("could not match sample name pattern in fastq.gz file names in ", 
       data_path
  )
}

pdf("QualProf.pdf")
suppressMessages(plotQualityProfile(fastq1s[1:2]))
suppressMessages(plotQualityProfile(fastq2s[1:2]))
graphics.off()

# Make directory and filenames for the filtered fastqs
filt_path <- file.path(data_path, "filtered")
if (!file_test("-d", filt_path))
  dir.create(filt_path)
if (length(fastq1s) != length(fastq2s))
  stop("Forward and reverse files do not match.")
filt1s <-
  file.path(filt_path, paste0(sample_names, "_R1_filt.fastq.gz"))
filt2s <-
  file.path(filt_path, paste0(sample_names, "_R2_filt.fastq.gz"))

# Filter
writeLines("filter_output <- filterAndTrim(
  fastq1s,
  filt1s,
  fastq2s,
  filt2s,
  truncLen = truncate, 
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)\n")

filter_output <- filterAndTrim(
  fastq1s,
  filt1s,
  fastq2s,
  filt2s,
  truncLen = truncate,
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)
filter_output <- bind_cols(sample_name = sample_names, file_name = rownames(filter_output), as.data.frame(filter_output)) 
print(filter_output)

# Scan filter path for files that passed filtering
final_filt1s <- list.files(file.path(data_path, "filtered"),  pattern = "*_R1_filt.fastq.gz", full.names = TRUE)
final_filt2s <- list.files(file.path(data_path, "filtered"), pattern = "*_R2_filt.fastq.gz", full.names = TRUE)
final_filt1s %>% 
  sub("_R1_filt.fastq.gz$", "", .) %>% 
  sub((file.path(data_path, "filtered/")), "", .) ->
  filt_sample_names 

filtered_out <- sample_names[!(sample_names %in% filt_sample_names)]
if (length(filtered_out) == 0) {
  writeLines("No sample was filtered out")
} else if (length(filtered_out) == 1) {
  print(paste0("Sample: ", filtered_out, " was filtered out"))
  filter_output %<>% filter(!sample_name %in% filtered_out)
} else if (length(filtered_out) > 1) {
  print(paste0("Samples: ", paste(filtered_out, collapse = ", "), " were filtered out"))
  filter_output %<>% filter(!sample_name %in% filtered_out)
}

# Learn the error rates -----------------------------------------------------------------------
writeLines("\nLearn the error rates:\n")
err1 <- learnErrors(final_filt1s, randomize = TRUE, multithread = TRUE)
err2 <- learnErrors(final_filt2s, randomize = TRUE, multithread = TRUE)

suppressMessages(ggsave("Read1Errors.pdf", plotErrors(err1, nominalQ = TRUE)))
suppressMessages(ggsave("Read2Errors.pdf", plotErrors(err2, nominalQ = TRUE)))

# Dereplication -------------------------------------------------------------------------------
writeLines("\nDereplicate:\n")
derep1s <- derepFastq(final_filt1s, verbose = TRUE)
derep2s <- derepFastq(final_filt2s, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derep1s) <- filt_sample_names
names(derep2s) <- filt_sample_names

# Sample inference ----------------------------------------------------------------------------
writeLines("\nSample inference with DADA2:\n")
dada1s <- dada(derep1s, err = err1, pool = pooling, multithread = TRUE)
dada2s <- dada(derep2s, err = err2, pool = pooling, multithread = TRUE)

writeLines("\nSummarise DADA2 results:\n")
dada1s[[1]]

# Merge paired reads --------------------------------------------------------------------------
writeLines("\nMerge paired reads:\n")
mergers <-
  mergePairs(dada1s, derep1s, dada2s, derep2s, minOverlap = 6, verbose = TRUE)
# Inspect the merger data.frame from the first sample
writeLines("First 10 merged sequences:")
head(mergers[[1]], n = 10)

# reverseComplement rev-primer reads (if analysing a 4 files per sample MiSeq run (i.e. DOME method))
for (rev_sample in grep("\\_R", names(mergers))) {
 mergers[[rev_sample]]$sequence <-
   microseq::reverseComplement(mergers[[rev_sample]]$sequence)
}

# Construct sequence table and remove chimeras ------------------------------------------------
seqtab <- makeSequenceTable(mergers)
#seqtab <- makeSequenceTable(mergers[grep("^Zymo", names(mergers), orderBy = TRUE)])

writeLines("\nSequence table dimensions:")
dim(seqtab)
writeLines("\nInspect distribution of sequence lengths:")
table(nchar(getSequences(seqtab)))

# Now merge fwd-primer and rev-primer samples (if analysing a 4-files-per-sample-MiSeq run (i.e. DOME method)) 
samples2merge <- gsub('\\_[FR]$', '', rownames(seqtab)) # remove FR to make identical names

seqtab %>% 
  rowsum(., as_vector(samples2merge)) ->
  seqtab_merged

duplicated_samples <-
  duplicated(samples2merge) |
  duplicated(samples2merge, fromLast = TRUE)

if (!all(duplicated_samples)) {
  cat(
    "The following sample(s) were missing either fwd or reverse primer reads:",
    rownames(seqtab_merged)[!duplicated_samples]
  )
}


writeLines("\nRemove sequences differing from excected length by more than 5 bp:")
seqtab_clean <- seqtab_merged[, nchar(colnames(seqtab_merged)) %in% seq(amp_size - 5, amp_size + 5)]

## Chimeras
writeLines("\nRemove chimeras:")
seqtab_nochim <-
  removeBimeraDenovo(
    seqtab_clean,
    method = "consensus", 
    allowOneOff = TRUE,
    verbose = TRUE,
    multithread = TRUE
  )
writeLines("\nSequence table dimensions after chimera removal:")
dim(seqtab_nochim)
writeLines("\nWhat proportions passed through chimera check?")
sum(seqtab_nochim) / sum(seqtab_clean)

# Save sequences and seq_table --------------------------------------------
# save ASVs to a fasta file
Seq_names <- paste0("Seq_", seq(1, ncol(seqtab_nochim)))
Seqs_nochim <- data.frame(
  Seq.name = Seq_names,
  Sequence = colnames(seqtab_nochim),
  stringsAsFactors = FALSE
)
write.fasta(
  as.list(Seqs_nochim$Sequence),
  Seqs_nochim$Seq.name,
  "DADA2.Seqs.fa",
  as.string = TRUE,
  nbchar = 1000
)

saveRDS(seqtab_nochim, "DADA2.seqtab_nochim.RDS") # save also an RDS object for easy merging

# colnames(seqtab_nochim) %<>% 
# gsub(".*([mM]ock)", "\\1", .) # this is needed for uncross

# save seq table without sequences
seqtab_nochim %>% 
  t() %>% 
  as_tibble(.name_repair = "check_unique")  %>% 
  #  bind_cols("Seq#" = Seq_names, .) 
  bind_cols("ASV" = Seq_names, .) %>% 
  write_delim(.,
              path = "DADA2.seqtab_nochim.tsv",
              delim = "\t",
              col_names = TRUE
  )

# Track reads through the pipeline ------------------------------------------------------------
filter_output %>% 
  select_if(is.numeric) %>% 
  bind_cols(merged_sample_name = as_vector(samples2merge), .) %>% 
  group_by(merged_sample_name) %>% 
  summarise_all(., "sum", na.rm = TRUE) %>% 
  select(-merged_sample_name) -> 
  numeric_cols

filter_output %>% 
  select_if(negate(is.numeric)) %>% 
  bind_cols(merged_sample_name = as_vector(samples2merge), .) %>% 
  group_by(merged_sample_name) %>% 
  summarise_all(list(~paste(unique(.), collapse = ",")))  -> 
  other_cols

getN <- function(x) sum(getUniques(x))
track <- cbind(other_cols, 
               numeric_cols,
               rowsum(map_dbl(dada1s, getN), as_vector(samples2merge)), 
               rowsum(map_dbl(mergers, getN), as_vector(samples2merge)), 
               rowSums(seqtab_clean), 
               rowSums(seqtab_nochim))
colnames(track) <- c("merged_sample_name", "sample_name", "Read1_file", "input", "filtered", "denoised", "merged", "tabled", "nonchim")
# track <- cbind(filt_sample_names, track)
writeLines("\nTrack reads through the pipeline:")
head(track)

write_delim(track, 
            path = "DADA2.track_each.tsv",
            delim = "\t",
            col_names = TRUE)

# Assign taxonomy -----------------------------------------------------------------------------
writeLines("\nAssign taxonomy:")
writeLines("Assigning SILVA taxonomy:")
# Using assignTaxonomy - SILVA
taxa_silva <- assignTaxonomy(
  seqs = seqtab_nochim,
  refFasta = file.path(db_path, silva_db),
  tryRC = TRUE,
  outputBootstraps = TRUE, 
  multithread = TRUE
)

# Add species
# taxa_silva[[1]] <- addSpecies(taxa_silva[[1]], file.path(db_path, silva_sp_db))

colnames(taxa_silva$boot) %<>% 
  paste(., "(BS)")

# save taxa table without sequences
taxa_silva %>% 
  do.call(cbind, .) %>% # only needed if bootstraps are included
  as_tibble(.name_repair = "minimal")  %>% 
  bind_cols("ASV" = Seq_names, .) %>% 
  write_delim(., 
              path = "DADA2.taxa_silva.tsv",
              delim = "\t",
              col_names = TRUE)
writeLines("Done")

writeLines("Assigning GTDB taxonomy:")
# Using assignTaxonomy - GTDB
taxa_gtdb <- assignTaxonomy(
  seqs = seqtab_nochim,
  refFasta = file.path(db_path, gtdb_db),
  tryRC = TRUE,
  outputBootstraps = TRUE, 
  multithread = TRUE
)

# Add species
# taxa_gtdb[[1]] <- addSpecies(taxa_gtdb[[1]], file.path(db_path, gtdb_sp_db))

colnames(taxa_gtdb$boot) %<>% 
  paste(., "(BS)")

# save taxa table without sequences
taxa_gtdb %>% 
  do.call(cbind, .) %>% # only needed if bootstraps are included
  as_tibble(.name_repair = "minimal")  %>%  
  bind_cols("ASV" = Seq_names, .) %>% 
  write_delim(., 
              path = "DADA2.taxa_gtdb.tsv",
              delim = "\t",
              col_names = TRUE)
writeLines("Done")
# Using Decipher
# doesn't work as well as assignTaxonomy 07/08/2019
# dna <-
#   DNAStringSet(getSequences(seqtab_nochim)) # Create a DNAStringSet from the ASVs
# load(file.path(db_path, silva_decipher_db)) # load training set
# ids_silva <-
#   IdTaxa(
#     dna,
#     trainingSet,
#     strand = "top",
#     threshold = 50,
#     processors = NULL,
#     verbose = FALSE
#   ) # use all processors
# ranks <-
#   c("domain",
#     "phylum",
#     "class",
#     "order",
#     "family",
#     "genus",
#     "species") # ranks of interest
# # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
# taxid_silva <- t(sapply(ids_silva, function(x) {
#   m <- match(ranks, x$rank)
#   taxa <- x$taxon[m]
#   taxa[startsWith(taxa, "unclassified_")] <- NA
#   taxa
# }))
# rownames(taxid_silva) <- getSequences(seqtab_nochim)
# taxid_silva <- addSpecies(taxid_silva, file.path(db_path, silva_sp_db))
# colnames(taxid_silva) <- c("Domain",
#                            "Phylum",
#                            "Class",
#                            "Order",
#                            "Family",
#                            "Genus",
#                            "Species",
#                            "Exact species")
# 
# taxid_silva %<>% 
#   as_tibble() %>% 
#   bind_cols("Seq#" = Seq_names, .) 
# 
# write_delim(taxid_silva, 
#             path = "DADA2.taxid_silva.tsv",
#             delim = "\t",
#             col_names = TRUE)

# load(file.path(db_path, GTDB_decipher_db)) # load training set
# ids_GTDB <-
#   IdTaxa(
#     dna,
#     trainingSet,
#     strand = "top",
#     threshold = 50,
#     processors = NULL,
#     verbose = FALSE
#   ) # use all processors
# ranks <-
#   c("domain",
#     "phylum",
#     "class",
#     "order",
#     "family",
#     "genus",
#     "species") # ranks of interest
# # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
# taxid_GTDB <- t(sapply(ids_GTDB, function(x) {
#   m <- match(ranks, x$rank)
#   taxa <- x$taxon[m]
#   taxa[startsWith(taxa, "unclassified_")] <- NA
#   taxa
# }))
# colnames(taxid_GTDB) <- ranks
# rownames(taxid_GTDB) <- getSequences(seqtab_nochim)
# 
# taxid_GTDB %<>% 
#   as_tibble() %>% 
#   bind_cols("Seq#" = Seq_names, .) 
# 
# write_delim(taxid_GTDB, 
#             path = "DADA2.taxid_GTDB.tsv",
#             delim = "\t",
#             col_names = TRUE)


# Construct phylogenetic trees ----------------------
## run seperately through 05_calc_FastTree.sh

# seqs2align <- Seqs_nochim$Sequence
# names(seqs2align) <- seqs2align # This propagates to the tip labels of the tree
# alignment <- AlignSeqs(DNAStringSet(seqs2align), anchor = NA)
# 
# phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
# dm <- dist.ml(phang.align)
# treeNJ <- NJ(dm) # Note, tip order != sequence order
# fit = pml(treeNJ, data = phang.align)
# 
# ## negative edges length changed to 0!
# 
# fitGTR <- update(fit, k = 4, inv = 0.2)
# fitGTR <- optim.pml(
#   fitGTR,
#   model = "GTR",
#   optInv = TRUE,
#   optGamma = TRUE,
#   rearrangement = "stochastic",
#   control = pml.control(trace = 0)
# )

# Evaluate accuracy ---------------------------------------------------------------------------
writeLines("\nEvaluate mock community:")
if (nchar(mock_id) != 0) {# if mock id was supplied
  mock_names <-
    seqtab_nochim[grep(mock_id, rownames(seqtab_nochim), ignore.case = TRUE)[1], ] # first mock community
  mock_names <-
    sort(mock_names[mock_names > 0], decreasing = TRUE) # Drop ASVs that are absent in the Mock
  mock_seqs <-
    Seqs_nochim$Sequence[Seqs_nochim$Sequence %in% names(mock_names)]
  print(paste(
    "DADA2 inferred",
    length(mock_seqs),
    "sequences present in the Mock community sample:",
    rownames(seqtab_nochim)[grep(mock_id, rownames(seqtab_nochim), ignore.case = TRUE)[1]]
  ))
  
  mock_ref <- ShortRead::readFasta(file.path(db_path, mock_ref_seqs))
  match_ref <-
    sum(sapply(mock_seqs, function(x)
      any(grepl(
        x, as.character(sread(mock_ref))
      ))))
  print(paste("Of those,",
              sum(match_ref),
              "were exact matches to the expected reference sequences."))
} else {writeLines("No mock community identifier supplied")}



