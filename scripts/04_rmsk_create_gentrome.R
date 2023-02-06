#!/usr/bin/env Rscript
# Gennaro Calendo
# 2023-02-02
#
# Create the gentrome (transcripts fasta + RepeatMasker fasta + genome fasta)
# required for Salmon indexing. Also return the decoys.txt file.
#
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if ( (length(args) != 4) | any(c("-h", '--help') %in% args)) {
  msg <- "Usage:
  
  ./04_rmsk_create_gentrome.R genome.fa tx.fa rmsk.fa out_dir

Help:
  genome.fa : The genome fasta file for the given organism. e.g. GRCh38.primary_assembly.fa.gz
  tx.fa     : Transcripts fasta file. e.g. gencode.v24.transcripts.fa.gz
  rmsk.fa   : Fasta file containing the unique TE sequences. e.g. rmsk-unique.fa
  out_dir   : Directory to save the output files.
  "
  cat(msg)
  stop()
}

genome_fa <- args[1]
tx_fa <- args[2]
rmsk_fa <- args[3]
out_dir <- args[4]

# Gentrome generation ---------------------------------------------------------
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Biostrings library is necessary for script execution!")
}

message("Reading in genome fasta...")
genome_seqs <- Biostrings::readDNAStringSet(genome_fa, format = "fasta")

# Subset for only the primary chromosomes
genome_seqs <- genome_seqs[grepl("chr[0-9]+|chr[XY]", names(genome_seqs))]

# Fix the names of the DNAStringSet (they import as "chr1 1", "chr2 2", etc.)
names(genome_seqs) <- regmatches(names(genome_seqs), regexpr("chr[0-9]+|chr[XY]", names(genome_seqs)))

message("Reading in transcripts fasta...")
tx_seqs <- Biostrings::readDNAStringSet(tx_fa, format = "fasta")

message("Reading in unique RepeatMasker instances fasta...")
rmsk_seqs <- Biostrings::readDNAStringSet(rmsk_fa, format = "fasta")

# Create the gentrome from combined seqs and write out
gentrome <- c(tx_seqs, rmsk_seqs, genome_seqs)
gentrome_fa <- file.path(out_dir, "rmsk-gentrome.fa.gz")

message("Writing out gentrome to", gentrome_fa)
Biostrings::writeXStringSet(gentrome, filepath = gentrome_fa, compress = TRUE)

# Decoy generation -------------------------------------------------------------
# Get the names of the genome fasta headers for the decoys file
decoys <- names(genome_seqs)
decoy_file <- file.path(out_dir, "decoys.txt")
message("Writing out decoys to", decoy_file)
write.table(
  decoys,
  file = decoy_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

message("Done.\nCreate a Salmon index with:\n\n")
message("salmon index -t", gentrome_fa, "-d", decoy_file, "-p 12 -i GRCh38.rmsk.salmon_index --gencode\n")
