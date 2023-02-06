#!/opt/miniconda3/envs/rmsk-gentrome/bin/python
"""
Gennaro Calendo
2023-01-04

This module includes functions to download the neccesary files from the remote 
sources. The program requires:

- RepeatMasker.out file
- Transcript sequences from GENCODE
- Primary assembly genome sequences
- Gene annotation GTF

"""
import argparse
import sys
import urllib.request
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Download the file necessary for downstream index creation")
    parser.add_argument('-o', '--out_dir', help='Where to save the downloaded files.', default='.')
    parser.add_argument('-s', '--species', help='What species to download files for.', choices=['Hs', 'Mm'], default='Hs')
    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    species = args.species
    
    if species != "Hs":
        sys.exit("Error: Only Homo sapiens (Hs) is available at this time.")

    urls = [
        "http://repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz",
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.transcripts.fa.gz",
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/GRCh38.primary_assembly.genome.fa.gz",
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz"
        ]
    fnames = [
        "hg38.fa.out.gz", 
        "gencode.v42.transcripts.fa.gz", 
        "GRCh38.primary_assembly.genome.fa.gz", 
        "gencode.v42.annotation.gtf.gz"
        ]
    resources = list(zip(urls, fnames))

    for resource in resources:
        url = resource[0]
        fpath = Path(out_dir, resource[1])
        if fpath.exists():
            print(f"{fpath} already exists. Skipping...")
            continue
        print(f"Downloading {url} to {fpath}")
        urllib.request.urlretrieve(url, fpath)


if __name__ == '__main__':
    main()


