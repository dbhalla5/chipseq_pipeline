library(DiffBind)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene) 
library(biomaRt)
library(tidyverse)
library(dplyr)


# Path to your data
path_to_main <- "/path/to/chip_seq_analysis"
path_to_csv <- paste0(path_to_main, "/diff_bd_input.csv")

# Load your sample sheet and create the dbObj (make sure this is before the count step)
samples <- read.csv(path_to_csv)
samples$Condition <- fct_relevel(samples$Condition, "S")
dbObj <- dba(sampleSheet = samples)

counts <- dba.count(dbObj,fragmentSize=0,bUseSummarizeOverlaps=T,summits=0 , bParallel=FALSE)
# counts <- dba.count(dbObj, fragmentSize=0, bUseSummarizeOverlaps=F, summits=0 , bParallel=FALSE)


# Create contrasts based on conditions
#dbObj_contra <- dba.contrast(counts, categories = DBA_CONDITION, minMembers = 2)

dbObj_contra <- dba.contrast(counts, contrast = c("Condition","N","S"), minMembers = 2)

# Perform the analysis
dbObj <- dba.analyze(dbObj_contra,bGreylist=FALSE , bBlacklist=FALSE)
########

report <- dba.report(dbObj)
annotated_peaks_t <- annotatePeak(report, TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene)
peak_anno_df <- as.data.frame(annotated_peaks_t)
ensembl <- useEnsembl(biomart="genes", dataset="mmusculus_gene_ensembl")  # For mouse

split_list <- function(lst, chunk_size) {
  split(lst, ceiling(seq_along(lst) / chunk_size))
}

# Split the transcript IDs into smaller chunks (e.g., 1000 per batch)

# Reduce chunk size to 500 or 200
chunk_size <- 500  # or 200

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

res_deseq_df <- as.data.frame(report)

# Convert res_deseq to GRanges if it is not already in that format
res_deseq_gr <- GRanges(seqnames = res_deseq_df$seqnames, 
                        ranges = IRanges(start = res_deseq_df$start, end = res_deseq_df$end))

# Annotate the peaks using ChIPseeker
annotated_peaks_c <- annotatePeak(res_deseq_gr, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)

# Convert to data frame
annotated_peaks_df_n <- as.data.frame(annotated_peaks_c)

##########
# Ensure the transcript IDs are extracted correctly
transcript_ids_df <- dplyr::select(annotated_peaks_df_n, transcriptId)
transcript_ids <- unlist(transcript_ids_df)

# Clean up transcript IDs
transcript_ids <- trimws(transcript_ids)  # Remove any leading/trailing spaces
transcript_ids_clean <- sub("\\.\\d+$", "", transcript_ids)  # Remove version numbers

# Split the transcript IDs into chunks
transcript_id_chunks <- split_list(transcript_ids_clean, chunk_size)

# Initialize a list to hold gene_info chunks
gene_info_list <- list()

# Process each chunk and populate gene_info_list
for (chunk in transcript_id_chunks) {
  gene_info_chunk <- get_with_retry(chunk, ensembl)
  if (!is.null(gene_info_chunk)) {
    gene_info_list[[length(gene_info_list) + 1]] <- gene_info_chunk
  }
}

# Combine all chunks into the final `gene_info` data frame
if (length(gene_info_list) > 0) {
  gene_info <- do.call(rbind, gene_info_list)
} else {
  stop("No gene information was retrieved. Check your transcript IDs or connection to Ensembl.")
}

##########

# Merge with gene_info to get gene name and gene ID
annotated_peaks_df_n$transcriptId_cleaned <- sub("\\.\\d+$", "", annotated_peaks_df_n$transcriptId)


merged_peaks_with_deseq_0 <- merge(annotated_peaks_df_n, res_deseq_df, 
                                 by = c("seqnames", "start", "end"), 
                                 all.x = TRUE)

merged_peaks_with_deseq_1 <- merged_peaks_with_deseq_0 %>% relocate(transcriptId_cleaned, .after=transcriptId)



merged_peaks_with_deseq_2 <- merge(merged_peaks_with_deseq_1, gene_info, 
                                 by.x = "transcriptId_cleaned", 
                                 by.y = "ensembl_transcript_id", 
                                 all.x = TRUE)

merged_peaks_with_deseq_3 <- merged_peaks_with_deseq_2 %>% relocate(transcriptId_cleaned, .after=transcriptId)

merged_peaks_with_deseq_4 <- merged_peaks_with_deseq_3 %>% relocate(transcriptId_cleaned, .before=ensembl_gene_id)
merged_peaks_with_deseq_5 <- merged_peaks_with_deseq_4 %>% relocate(transcriptId, .before=transcriptId_cleaned)
merged_peaks_with_deseq <- merged_peaks_with_deseq_5

path_to_mapped_res_new <- paste0(path_to_main, "/DiffBind_DESeq2_S_vs_N_gene_mapped_results.txt")
write.table(merged_peaks_with_deseq, file = path_to_mapped_res_new, sep = "\t")



genotype_N_enrich <- merged_peaks_with_deseq %>%
  filter(FDR < 0.05 & Fold > 0)


genotype_S_enrich <- merged_peaks_with_deseq  %>%
  filter(FDR < 0.05 & Fold < 0) #%>%
#  select(seqnames, start, end)


# Clean transcriptId columns by removing version numbers
genotype_S_enrich$transcriptId_cleaned <- sub("\\.\\d+$", "", genotype_S_enrich$transcriptId_cleaned)
genotype_N_enrich$transcriptId_cleaned <- sub("\\.\\d+$", "", genotype_N_enrich$transcriptId_cleaned)

# Merge with gene_info for genotype_1_enrich
genotype_S_enrich_merged_0 <- merge(genotype_S_enrich, gene_info, 
                                    by.x = "transcriptId_cleaned", 
                                    by.y = "ensembl_transcript_id", 
                                    all.x = TRUE)

genotype_S_enrich_merged <-  genotype_S_enrich_merged_0  %>% relocate(transcriptId_cleaned, .after=transcriptId)



# Merge with gene_info for genotype_2_enrich
genotype_N_enrich_merged_0 <- merge(genotype_N_enrich, gene_info, 
                                    by.x = "transcriptId_cleaned", 
                                    by.y = "ensembl_transcript_id", 
                                    all.x = TRUE)
genotype_N_enrich_merged <-  genotype_N_enrich_merged_0  %>% relocate(transcriptId_cleaned, .after=transcriptId)

# Save to file
path_to_mapped_res_genotype_S <- paste0(path_to_main, "/DiffBind_DESeq2_Decreased_binding_in_N.txt")
write.table(genotype_S_enrich_merged, file = path_to_mapped_res_genotype_S, sep = "\t", row.names = FALSE, quote = FALSE)


path_to_mapped_res_genotype_N <- paste0(path_to_main, "/DiffBind_DESeq2_Increased_binding_in_N.txt")
write.table(genotype_N_enrich_merged, file = path_to_mapped_res_genotype_N, sep = "\t", row.names = FALSE, quote = FALSE)




