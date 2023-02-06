#!/opt/miniconda3/envs/rmsk-gentrome/bin/python
"""
Convert RepeatMasker out file to BED format.
"""
from pathlib import Path
import argparse
import gzip


def main():
    parser = argparse.ArgumentParser(description="Convert RepeatMasker out to BED file")
    parser.add_argument("rmsk", help="RepeatMasker out file")
    parser.add_argument(
        "-e",
        "--exclude",
        help="Features to exclude from the BED output",
        default=["Simple_repeat", "Low_complexity", "Satellite", "RNA", "rRNA",
                "snRNA", "scRNA", "srpRNA", "tRNA", "Unknown"],
    )
    parser.add_argument("-l", "--min_length", help="Minimum length of the sequence", default=31)
    parser.add_argument("--gzipped", help="The input file is gzip compressed", action="store_true")

    # Setup parameters ------------------------------------------------------------
    args = parser.parse_args()
    rmsk = Path(args.rmsk)
    exclude = args.exclude
    gzipped = args.gzipped
    min_len = args.min_length
    outfile = rmsk.with_suffix(".bed")
    out = open(outfile, "w")

    # Define canonical chromosomes to keep in output bedfile
    keep_chroms = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", 
                "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
                "chrX", "chrY"}

    if gzipped:
        f = gzip.open(rmsk, "rt")
    else:
        f = open(rmsk, "r")

    # Skip header
    for _ in range(3):
        next(f)

    # Processing ------------------------------------------------------------------
    line_counter = 1
    processed_lines = 0
    bad = 0
    too_short = 0
    excluded = 0
    non_canonical = 0
    for line in f:
        line_counter += 1
        l = line.strip().split()

        # Skip bad lines
        if len(l) != 15:
            bad += 1
            print(f"Bad line found at line number: {line_counter}")
            print("The offending line is:")
            print(line)
            print("Skipping...")
            continue

        score = l[0]
        chrom = l[4]
        start = int(l[5]) - 1  # Convert to 0-based BED style
        end = l[6]
        strand = "-" if l[8] == "C" else "+"
        rep_elem = l[9]
        class_fam = l[10].split("/")
        rep_class = class_fam[0]
        rep_id = l[14]

        # Skip non-canconical chromosomes
        if chrom not in keep_chroms:
            non_canonical += 1
            continue

        # Fill repeat family name with class name if absent
        if len(class_fam) == 1:
            rep_fam = rep_class
        else:
            rep_fam = class_fam[1]

        # Skip lines that have element in exclusion list
        if (rep_elem in exclude) or (rep_class in exclude) or (rep_fam in exclude):
            excluded += 1
            continue

        # Skip lines where the sequence isnt long enough to hash
        if abs(int(start) - int(end)) < min_len:
            too_short += 1
            continue

        # Create a BED formatted record with a score column
        name = ".".join([rep_id, rep_class, rep_fam, rep_elem])
        out.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
        processed_lines += 1

    f.close()
    out.close()

    # Summary ---------------------------------------------------------------------
    print(f"\nMalformed RepeatMasker input lines (skipped): {bad}")
    print(f"Excluded records: {excluded}")
    print(f"Records shorter than {min_len} bp: {too_short}")
    print(f"Records in non-canonical chromosomes (skipped): {non_canonical}")
    print(f"Records written to file: {processed_lines}")


if __name__ == '__main__':
    main()


