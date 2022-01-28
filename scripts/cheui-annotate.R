#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=4) {
  stop("\nUsage: Rscript liftoverCustomTranscriptome.R /path/to/bar.gtf /path/to/foo.bed /path/to/output.bed cheui-model", call.=FALSE)
}

library(GenomicFeatures)
library(rtracklayer)
library(tidyverse)

# interactive args for testing p-val
# args <- c("~/localGadiData/2021-12-24_test-custom-metaplot/Mus_musculus.GRCm39.104.chr.gtf", "/Users/asethi/localGadiData/2021-12-25_develop-cheu-annotate/WT_E18_m6a.tempbed", "/Users/asethi/localGadiData/2021-12-25_develop-cheu-annotate/WT_E18_m6a.bed","pval")

################################################################################
################################################################################
################################################################################

# start 

# tstart <- print(paste("start time is", as.character(Sys.time())))
# tend <- print(paste("end time is", as.character(Sys.time())))
# tstart
# tend

##################################################

# process the annotation 

# read in reference transcripts
gtf <- makeTxDbFromGFF(file=args[1], format = "gtf")

# make an exon database from the reference transcripts 
exons <- exonsBy(gtf, "tx", use.names=TRUE)

# convert the transcripts to a tibble
exons_tib <- as_tibble(as(exons, "data.frame"))

# fetch the lengthe of each transcript segment from the gtf 
txlen <- transcriptLengths(gtf, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE) %>% 
  as_tibble() %>% 
  mutate(diff = tx_len - cds_len - utr5_len - utr3_len) %>% 
  dplyr::rename(transcript_id = tx_name) %>% 
  dplyr::select(-tx_id, -nexon, -gene_id, -diff)

# the last command doesn't store biotype, so we read in the gtf again using another package 
tx_biotype <- rtracklayer::import(args[1]) %>% 
  as_tibble() %>% 
  dplyr::select(transcript_id, transcript_biotype, gene_name, gene_id) %>% 
  na.omit() %>% 
  distinct()

# merge the biotypes with the transcript segment lengths 
merged_metadata <- inner_join(tx_biotype, txlen, by = "transcript_id")

##################################################

# lift-over the bed-like CHEUI output 

# import bed file of transcriptome alignments
mappedLocus <- read_tsv(file = args[2], col_names = F) %>%
  rename(transcript = 1, start = 2, end = 3, name = 4, score = 5, strand = 6) %>%
  dplyr::select(-strand) %>% 
  separate(transcript, into=c("transcript_id", "transcript_version"), sep = "([.])", extra = "merge") %>% 
  dplyr::select(-transcript_version)

# make lookup table for strand
print("preparing strand lookup table")
strand_lookup <- exons_tib %>%
  rename(transcript_id = group_name) %>%
  dplyr::select(transcript_id, strand) %>% 
  distinct()

# attach the correct strand to the bed sites
print("reparing strand")
mappedLocus_fixedStrand <- inner_join(mappedLocus, strand_lookup, by = "transcript_id") 

# write out the strand-repaired file as a temporary file
print("writing strand bedfile")
write_tsv(mappedLocus_fixedStrand, args[3], col_names = F)

# read in the bed sites with corrected strand
print("importing strand bedfile")
mappedLocus <- import.bed(args[3])

# map transcript coordinates to genome
print("mapping transcript coordinates to genome")
genomeLocus <- mapFromTranscripts(x=mappedLocus, transcripts=exons)

# bind score to output
# the score column contains cheui-specific output (e.g. stoich, prob, coverage, which we aggreagate into the score column and delimit using semicolons)
print("binding output")
mcols(genomeLocus)<-cbind(mcols(genomeLocus),DataFrame(mappedLocus[genomeLocus$xHits]))

# convert output to tibble
genome_coordinates = as_tibble(as(genomeLocus, "data.frame"))

# prepare the output by selecting bed-like coordinates from
print("filtering output")
output <- genome_coordinates %>% dplyr::select(seqnames, start, end, X.name, X.seqnames, strand) %>%
  unique() %>% 
  dplyr::rename(chr = seqnames, data = X.name, transcript = X.seqnames) %>% 
  mutate(score = ".") %>% 
  dplyr::select(chr, start, end, transcript, score, strand, data)

# based on the CHEUI runmode, separate the data column into it's actual values 
if (args[4] == "pval"){
  output <- output %>% 
    separate(data, into=c("tx_coord", "motif", "coverage", "stoich", "prob", "p.raw"), sep = ";") %>% 
    type_convert(col_types = "fiiffficiddd")
  
  output$p.adj <- p.adjust(output$p.raw, method = "fdr")

  
}else{print("no")}


##################################################

# attach the transcript lengths and biotype to output using the transcriptome coordinate 
merge_out <- inner_join(output, merged_metadata %>% dplyr::rename(transcript = transcript_id), by = "transcript")

meta <- merge_out %>% 
  mutate(cds_start = utr5_len,
         cds_end = utr5_len + cds_len, 
         tx_end = cds_end + utr3_len) %>% 
  mutate(rel_pos = ifelse(tx_coord < cds_start, # if the site is in the 5' utr 
         ((tx_coord)/(cds_start)), # the relative position is simply the position of the site in the UTR 
          ifelse(tx_coord < cds_end, # else, if the site is in the CDS
          (1 + (tx_coord - utr5_len)/(cds_len)), # the relative position is betwee 1 and 2 and represents the fraction of the cds the site is in 
          (2 + (tx_coord - utr5_len - cds_len) / utr3_len))),  # no final condition, the site must be in the 3' utr, similar as before but the rel_pos is between 2 and 3 
         abs_cds_start = tx_coord - cds_start, # absolute pos relative to CDS start 
         abs_cds_end = tx_coord - cds_end) # absolute pos relative to CDS end 
         
##################################################

# write the output
print("writing final output")
write_tsv(output, args[3], col_names = T, append = FALSE)

quit()

################################################################################
################################################################################
################################################################################
