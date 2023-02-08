#!/usr/bin/env Rscript
# Gennaro Calendo
# 2023-02-02
#
# Create a tx2gene file for tximport and a GRangeslist object from the indexed
# transcripts and TE sequences
#
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if ( (length(args) != 3) | any(c("-h", '--help') %in% args)) {
  msg <- "Usage:
  
  ./05_rmsk_annotation.R rmsk-duplicateInfo.json transcripts.fa out_dir

Help:
  duplicateInfo.json : RepeatMasker duplicate information produced by the extract_fasta script
  transcripts.fa     : GTF file containing transcript information for all transcripts. Can be gzipped.
  out_dir            : Directory to save the annotation files to.
  
  "
  cat(msg)
  stop()
}

info_json <- args[1]
gtf <- args[2]
out_dir <- args[3]

# Processing --------------------------------------------------------------

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("data.table library is necessary")
}
suppressPackageStartupMessages(library(data.table))

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("jsonlite library is necessary")
}
if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
  stop("GenomicRanges library is necessary")
}
if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  stop("rtracklayer library is necessary")
}

message("Reading in the duplicate information from ", basename(info_json), "...")
info <- jsonlite::fromJSON(info_json)

message("Creating data.table from json information...")
dt <- data.table(Hash = names(info), Instance = info)

message("Unnesting json information...")
dt <- dt[, .(Instance = as.character(unlist(Instance))), by = Hash]

message("Extracting repetitive element names from location strings...")
dt[, c("ID", "Location") := tstrsplit(Instance, "::", fixed=TRUE)]
dt[, `:=`(RepID = gsub("\\..*", "", ID, perl=TRUE),
          RepName = gsub("^[0-9]+\\.", "", ID, perl=TRUE))]

message("Extracting position information...")
dt[, c("seqnames", "position") := tstrsplit(Location, ":",  fixed = TRUE)]
dt[, strand := regmatches(position, regexec("\\(.\\)", position, perl=TRUE))][, 
     strand := gsub("\\(", "", strand, perl=TRUE)][, 
     strand := gsub("\\)", "", strand, perl=TRUE)]
dt[, position := gsub("\\(.\\)", "", position, perl=TRUE)][,
     c("start", "end") := tstrsplit(position, "-", fixed = TRUE)]
dt[, `:=`(Location = NULL, position = NULL, Instance = NULL, ID = NULL)]

message("Writing out RepeatMasker Hash ID to Location mapping to ", file.path(out_dir, "rmsk-annotation.tsv"))
fwrite(dt, file.path(out_dir, "rmsk-annotation.tsv"), sep="\t")

message("Creating tx2gene.tsv file for tximport...")
tx2gene <- unique(dt[, .(Hash, RepName)])
setkey(tx2gene, "Hash")

# Find all hashes that have more than one element contained
message("Finding hash values that map to multiple elements...")
dups <- tx2gene[, .N, by = Hash][N > 1]
message("Found ", nrow(dups), "such elements.")
setkey(dups, "Hash")

# For all Hashes with multiple entries, collapse RepNames into single strings
message("Collapsing multi-element hashes into single comma-delimited strings...")
dups2 <- tx2gene[dups]
dups2[,  c("Class", "Family", "Subfamily") := tstrsplit(RepName, ".", fixed=TRUE)]
dedup <- dups2[, .(RepName = paste0(paste(unique(Class), collapse = ","), ".", 
                                    paste(unique(Family), collapse = ","), ".",
                                    paste(unique(Subfamily), collapse = ","))), 
               by = Hash]

# remove duplicated hashes from tx2gene before adding cleaned hashes
tx2gene2 <- tx2gene[!Hash %chin% dups[, unique(Hash)]]

# Now add the de-dupped data back
tx2gene2 <- rbind(tx2gene2, dedup)

# Check that hashes with multiple different subfamilies have been properly collapsed
stopifnot("The number of hashes does not equal the number of features" = length(info) == nrow(tx2gene2))

message("Reading in transcript information from ", basename(gtf), "...")
tx_gtf <- rtracklayer::import(gtf)
tx_dt <- setDT(data.frame(tx_gtf))
tx <- unique(tx_dt[, .(transcript_id, gene_id)])
tx <- tx[!is.na(transcript_id)]

setnames(tx2gene2, c("TXNAME", "GENEID"))
setnames(tx, c("TXNAME", "GENEID"))
message("Combining transcript and repetitive element information to create tx2gene...")
tx2gene3 <- rbind(tx, tx2gene2)

message("Writing out tx2gene file to ", file.path(out_dir, "tx2gene.tsv"), "...")
fwrite(tx2gene3, file.path(out_dir, "tx2gene.tsv"), sep = "\t")

# GRangesList of rmsk elements --------------------------------------------

message("Creating GRangesList object for all RepeatMasker Hashes...")
rmsk_gr <- GenomicRanges::makeGRangesFromDataFrame(dt, keep.extra.columns = TRUE)
rmsk_grl <- S4Vectors::splitAsList(rmsk_gr, rmsk_gr$Hash)

message("Saving GRangesList object to ", file.path(out_dir, "rmsk_grl.rds"))
saveRDS(rmsk_grl, file.path(out_dir, "rmsk_grl.rds"))
message("Done.")
