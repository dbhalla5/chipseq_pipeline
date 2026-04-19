
# conda install -c bioconda bioconductor-biomaRt

library(DiffBind)
library(tidyverse)
library(BiocParallel)

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene) 

library(biomaRt)

# Use parallel processing for the analysis
library(BiocParallel)
param <- MulticoreParam(workers = 6)  # Set the number of workers to use (adjust based on your system)


path_to_main <- "/path/to/chip_seq_analysis"
path_to_dir <- "/path/to/chip_seq_analysis/files"
path_to_csv <- paste0(path_to_dir, "/path/to/diff_bd_input.csv")


samples <- read.csv(path_to_csv)
dbObj <- dba(sampleSheet=samples)

dbObj <- dba.blacklist(dbObj, NULL)
dbObj <- dba.greylist(dbObj, NULL)

dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)
dbObj <- dba.contrast(dbObj, categories=DBA_CONDITION, minMembers = 2)

# dbObj <- dba.blacklist(dbObj, NULL)
# dbObj <- dba.greylist(dbObj, NULL)

# Now, run the analysis
dbObj <- dba.analyze(dbObj, BPPARAM = param)

########

report <- dba.report(dbObj)
annotated_peaks <- annotatePeak(report, TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene)
peak_anno_df <- as.data.frame(annotated_peaks)
ensembl <- useEnsembl(biomart="genes", dataset="mmusculus_gene_ensembl")  # For mouse

################
# Split the cleaned list into chunks of 1000 IDs
split_list <- function(lst, chunk_size) {
  split(lst, ceiling(seq_along(lst) / chunk_size))
}

# Split the transcript IDs into smaller chunks (e.g., 1000 per batch)

# Reduce chunk size to 500 or 200
chunk_size <- 500  # or 200

####

# Check if the chunk looks correct
#print(transcript_id_chunks_8207[[1]])  # This should print the first chunk of unique transcript IDs

# Retry logic with delay
retry_limit <- 5  # Set a maximum retry limit
retry_delay <- 5  # Delay in seconds between retries

# Function to handle retries
get_with_retry <- function(chunk, mart) {
  retries <- 0
  success <- FALSE
  while (!success && retries < retry_limit) {
    try({
      gene_info_chunk <- getBM(attributes = c('ensembl_transcript_id', 'ensembl_gene_id', 'external_gene_name'),
                               filters = 'ensembl_transcript_id',
                               values = chunk,
                               mart = mart)
      success <- TRUE  # If no error occurs, set success to TRUE
    }, silent = TRUE)
    
    if (!success) {
      retries <- retries + 1
      message(paste("Retrying", retries, "of", retry_limit, "for chunk"))
      Sys.sleep(retry_delay)  # Pause before retrying
    }
  }
  return(gene_info_chunk)
}

################

# transcript_ids_df <- dplyr::select(peak_anno_df, transcriptId)
# transcript_ids <- unlist(transcript_ids_df)
# # Remove any leading/trailing spaces if present
# transcript_ids <- trimws(transcript_ids)
# # Remove version numbers (if any)
# transcript_ids_clean <- sub("\\.\\d+$", "", transcript_ids)
# transcript_id_chunks <- split_list(transcript_ids_clean, chunk_size)
# 
# #################
# 
# # Use the function in the loop
# gene_info_list <- list()
# for (chunk in transcript_id_chunks) {
#   gene_info_chunk <- get_with_retry(chunk, ensembl)  # ENSEMBL 102 is GRCm38 = mm10 
#   gene_info_list[[length(gene_info_list) + 1]] <- gene_info_chunk
# }
# 
# # Combine the results
# gene_info <- do.call(rbind, gene_info_list)
# 
# 
# # Step 1: Remove the version number from the transcriptId column (i.e., remove everything after the '.')
# peak_anno_df$transcriptId_cleaned <- sub("\\.\\d+$", "", peak_anno_df$transcriptId)
# 
# # Step 2: Merge the two dataframes on the cleaned transcriptId column and ensembl_transcript_id
# 
# merged_df_0 <- merge(peak_anno_df, gene_info, by.x = "transcriptId_cleaned", by.y = "ensembl_transcript_id", all.x = TRUE)
# merged_df <-  merged_df_0  %>% relocate(transcriptId_cleaned, .after=transcriptId) # move "transcriptId_cleaned" column to a different location 


#########

# mapped_file <- "/DiffBind_contrast_gene_name_mapped.txt"
# output_file_mapped <- paste0(path_to_dir, mapped_file)
# write.table(merged_df, file = output_file_mapped, sep = "\t", row.names = FALSE, quote = FALSE)

#########


#DESeq2

# https://hbctraining.github.io/Intro-to-ChIPseq/lessons/08_diffbind_differential_peaks.html
# Additionally, we will want to create BED files for each set of significant regions identified by DESeq2,
# separating them based on the gain or loss of enrichment. We will write these regions to file and 
# use as input for downstream visualization.


# Convert res_deseq to a data frame
res_deseq_df <- as.data.frame(res_deseq)

# Convert res_deseq to GRanges if it is not already in that format
res_deseq_gr <- GRanges(seqnames = res_deseq_df$seqnames, 
                        ranges = IRanges(start = res_deseq_df$start, end = res_deseq_df$end))

# Annotate the peaks using ChIPseeker
annotated_peaks <- annotatePeak(res_deseq_gr, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)


# Convert to data frame
annotated_peaks_df <- as.data.frame(annotated_peaks)

# Merge with gene_info to get gene name and gene ID
annotated_peaks_df$transcriptId_cleaned <- sub("\\.\\d+$", "", annotated_peaks_df$transcriptId)


merged_peaks_with_deseq <- merge(annotated_peaks_df, res_deseq_df, 
                                 by = c("seqnames", "start", "end"), 
                                 all.x = TRUE)

genotype_2_enrich <- merged_peaks_with_deseq %>%
  filter(FDR < 0.05 & Fold > 0)


genotype_1_enrich <- merged_peaks_with_deseq  %>%
  filter(FDR < 0.05 & Fold < 0) #%>%
#  select(seqnames, start, end)


# Clean transcriptId columns by removing version numbers
genotype_1_enrich$transcriptId_cleaned <- sub("\\.\\d+$", "", genotype_1_enrich$transcriptId_cleaned)
genotype_2_enrich$transcriptId_cleaned <- sub("\\.\\d+$", "", genotype_2_enrich$transcriptId_cleaned)

# Merge with gene_info for genotype_1_enrich
genotype_1_enrich_merged_0 <- merge(genotype_1_enrich, gene_info, 
                                  by.x = "transcriptId_cleaned", 
                                  by.y = "ensembl_transcript_id", 
                                  all.x = TRUE)

genotype_1_enrich_merged <-  genotype_1_enrich_merged_0  %>% relocate(transcriptId_cleaned, .after=transcriptId)

# Merge with gene_info for genotype_2_enrich
genotype_2_enrich_merged_0 <- merge(genotype_2_enrich, gene_info, 
                                  by.x = "transcriptId_cleaned", 
                                  by.y = "ensembl_transcript_id", 
                                  all.x = TRUE)
genotype_2_enrich_merged <-  genotype_2_enrich_merged_0  %>% relocate(transcriptId_cleaned, .after=transcriptId)

# path_to_g_type_2_file <-  paste0(path_to_dir, "/Diff_bind_plots/genotype_2_enriched.bed")
# write.table(genotype_2_enrich, file=path_to_g_type_2_file, sep="\t", quote=F, row.names=F, col.names=F)
# 
# 
# path_to_g_type_1_file <-  paste0(path_to_dir, "/Diff_bind_plots/genotype_1_enriched.bed")
# write.table(genotype_1_enrich, file=path_to_g_type_1_file, sep="\t", quote=F, row.names=F, col.names=F)


# Save to file
path_to_mapped_res_genotype_1 <- paste0(path_to_dir, "/DiffBind_DESeq2_gene_mapped_genotype_1_enriched.txt")
write.table(genotype_1_enrich_merged, file = path_to_mapped_res_genotype_1, sep = "\t", row.names = FALSE, quote = FALSE)


path_to_mapped_res_genotype_2 <- paste0(path_to_dir, "/DiffBind_DESeq2_gene_mapped_genotype_2_enriched.txt")
write.table(genotype_2_enrich_merged, file = path_to_mapped_res_genotype_2, sep = "\t", row.names = FALSE, quote = FALSE)






