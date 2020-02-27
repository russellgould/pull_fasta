from argparse import ArgumentParser, FileType
from numpy import where
from os import system
from .utils import *

parser = ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-gff", action="store_true", help="indicate the input file is GFF")
group.add_argument("-bed", action="store_true", help="indicate the input file is BED")
group.add_argument("-peak", action="store_true", help="indicate the input file is PEAK")
parser.add_argument("input", type=FileType("r"), help="file containing regions")
parser.add_argument(
    "-nu",
    "--nucs_up",
    required=True,
    dest="nucs_up",
    help="number of nucleotides upstream",
)
parser.add_argument(
    "-nd",
    "--nucs_down",
    required=True,
    dest="nucs_down",
    help="number of nucleotides downstream",
)
args = parser.parse_args()

if __name__ == "__main__":
    bedtools_cols = ["chrom", "start", "end", "id", "score", "strand"]

    if args.peak:
        regions = read_peak(args.input)
        regions["start"] = where(
            regions["strand"] == "-",
            regions["loc"] - args.nucs_down - 1,
            regions["loc"] - args.nucs_up - 1,
        )
        regions["end"] = where(
            regions["strand"] == "-",
            regions["loc"] + args.nucs_up,
            regions["loc"] + args.nucs_down,
        )
        regions["score"] = 0
        regions = regions[bedtools_cols]
    else:
        if args.gff:
            regions = read_gff(args.input)
        elif args.bed:
            regions = read_bed(args.input)

    regions = regions[bedtools_cols]
