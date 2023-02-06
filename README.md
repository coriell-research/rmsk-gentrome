# RepeatMasker Salmnon Index and Analysis

Create a new gentrome using the RepeatMasker (rmsk) input locations from all 
unique rmsk instances in the genome from the 5 main TE families only. 
Additionally track their locations so that multi-location subfamilies can be 
tracked and then later re-estimated.

## Structure

Currently, only a series of scripts has been implemented for the generation of
the Salmon index and annotation files. Each script is to be run in succession:

1. `01_rmsk_download.py` Retrieves the resources needed for 'gentrome' construction. 
It downloads the RepeatMasker.out file from [RepeatMasker.org](http://repeatmasker.org/species/hg.html), 
transcript sequences from [GENCODE](https://www.gencodegenes.org/human/), 
primary assembly genome sequences from GENCODE, and the gene annotation GTF 
from GENCODE. 

Only human sequences are supported at this time. The GENCODE version is v42.

2. `02_rmsk_to_bed.py` Filters out Simple_repeat", "Low_complexity", "Satellite", 
"RNA", "rRNA", "snRNA", "scRNA", "srpRNA", "tRNA", "Unknown" elements from the
RepeatMasker.out file as well as removes any sequences that are shorter than the
standard Salmon kmer length of 31. the exclusion list and kmer can be modified.

3. `03_rmsk_extract_fasta.py` Extracts the unique RepeatMasker TE sequences using 
the rmsk BED file and a genomic fasta file. The output of thus script is a fasta
file with each header corresponding to the SHA1 hash of the sequence. This script
also outputs a json file mapping rmsk hash values (keys) to the RepeatMasker 
elements and their locations (values). 

4. `04_rmsk_create_gentrome.R` Creates the gentrome (Transcripts + 
RepeatMasker elements + decoy sequences) fasta file used to create the 
[Salmon](https://salmon.readthedocs.io/en/latest/) index. A decoys.txt file is 
also produced for use in the indexing step. Indexing with `Salmon` can be run 
with:

`salmon index -t rmsk-gentrome.fa.gz -d decoys.txt -p 12 -i GRCh38.rmsk.salmon_index --gencode`

5. `05_rmsk_annotation.R` Produces data files that can be used in downstream 
analysis with R. The script produces an "rmsk-annotation.tsv" file that maps 
every hash value to every RepeatMasker element and it's genomic location. This
file is essentially a long version of the "rmsk-duplicateInfo.json" file. It also
produces a "tx2gene.tsv" which maps transcript IDs / RepeatMasker sequence hashes
to the Gene Id / RepeatMasker element. This file is used by the `tximport` library
in order to import transcript / hash counts and summarize to the gene / sub-family 
level. In the cases where one hash maps to multiple RepeatMasker elements, the 
elements are concatenated at their lowest common level by a ",". The script also 
produces "rmsk-grl.rds", a `GRangesList` object with keys as RepeatMasker hashes 
and values as `GRanges` objects with all of the ranges for sequences with the 
given hash.

## Requirements

- Python requirements : Python 3.10+ and `pybedtools`
- R requirements : `Biostrings`, `data.table`, `jsonlite`, `GenomicRanges`, and
`rtracklayer`
