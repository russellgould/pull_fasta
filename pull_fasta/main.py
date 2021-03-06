#!/usr/bin/env python

import sys

assert sys.version_info >= (3, 7), "Wrong Python version! Requires at least 3.7.0"

try:
    from argparse import ArgumentParser, FileType
    from pathlib import Path
    from subprocess import run
    from itertools import zip_longest
    from numpy import where
    from pandas import read_csv, Series, DataFrame
except:
    print(
        "Problem with imports! Ensure you are not running this from the base Conda environment."
    )
    sys.exit(1)


def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks
    Taken from https://docs.python.org/3/library/itertools.html#itertools-recipes
    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def read_peak(file_path: str):
    """Load TSSseq data file and return a pandas DataFrame."""
    return (
        read_csv(file_path)
        .rename(
            columns={
                "Chromosome": "chrom",
                "Strand": "strand",
                "Start": "start",
                "End": "end",
                "ReadCount": "total_count",
                "ModeLocation": "loc",
                "ModeReadCount": "mode_count",
                "Shape": "shape",
                "TranscriptLocation": "script_loc",
                "TranscriptID": "id",
                "GeneName": "name",
                "GeneType": "type",
                "%-Capped": "capped",
            }
        )
        .drop(
            columns=[
                "start",
                "end",
                "total_count",
                "mode_count",
                "capped",
                "shape",
                "script_loc",
            ]
        )
    )


def read_gff(file_path: str) -> DataFrame:
    """Load a GFF file and return a pandas DataFrame."""

    def format_attributes(df):
        for item in df["attributes"]:
            split_on_semi = item.split(";")
            if "=" in split_on_semi[0]:
                yield {
                    x.split("=")[0]: x.split("=")[1] for x in split_on_semi if "=" in x
                }
            else:
                yield "."

    def get_name(row):
        if type(row["attributes"]) == dict:
            if "ID" in row["attributes"].keys():
                return row["attributes"]["ID"]
            elif "Name" in row["attributes"].keys():
                return row["attributes"]["Name"]
            else:
                return "."
        else:
            return "."

    cols = [
        "chrom",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]
    df = read_csv(file_path, sep="\t", header=None)
    if df.shape[1] > 9:
        df = df.drop(df.columns[[i for i in range(9, df.shape[1])]], axis=1)

    df = df.rename(columns={i: cols[i] for i in range(df.shape[1])})

    if df["attributes"].all() != ".":
        df["attributes"] = Series(format_attributes(df))
        df["name"] = df.apply(get_name, axis=1)
    return df


def read_bed(file_path: str):
    """Load a BED file and return a pandas DataFrame."""
    cols = [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "thick_start",
        "thick_end",
        "item_rgb",
        "block_count",
        "block_sizes",
        "block_starts",
    ]
    df = read_csv(file_path, sep="\t", header=None)
    return df.rename(columns={i: cols[i] for i in range(df.shape[1])})


parser = ArgumentParser(
    usage="Pull FASTA sequences",
    description="This tool pulls FASTA sequences from an input file containing regions in a number of file formats.",
)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-gff", action="store_true", help="indicate the input file is GFF")
group.add_argument("-bed", action="store_true", help="indicate the input file is BED")
group.add_argument("-peak", action="store_true", help="indicate the input file is PEAK")
parser.add_argument("regions", type=FileType("r"), help="file containing regions")
parser.add_argument("output", type=FileType("w"), help="FASTA file to store output")
parser.add_argument(
    "-ref",
    "--reference",
    required=True,
    type=Path,
    help="reference fasta file",
    dest="ref",
)
parser.add_argument(
    "-nu", type=int, default=0, dest="nucs_up", help="number of nucleotides upstream",
)
parser.add_argument(
    "-nd",
    type=int,
    default=0,
    dest="nucs_down",
    help="number of nucleotides downstream",
)
args = parser.parse_args()

if args.peak:
    regions = read_peak(args.regions)
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
else:
    if args.gff:
        regions = read_gff(args.regions)
        regions["start"] -= 1
        if "name" not in regions.columns:
            regions["name"] = "."
    else:
        regions = read_bed(args.regions)

bedtools_cols = ["chrom", "start", "end", "name", "score", "strand"]
regions = regions[bedtools_cols]

tmp_path = Path("/tmp/tmp.bed").resolve()
regions.to_csv(tmp_path, sep="\t", index=False, header=False)

cmd = run(
    ["bedtools", "getfasta", "-fi", args.ref, "-bed", tmp_path, "-s"],
    capture_output=True,
)

output = cmd.stdout.decode().split()

with args.output as f:
    for idx, group in enumerate(grouper(output, 2)):
        chrom = regions.iloc[idx, 0]
        start = regions.iloc[idx, 1] + 1
        end = regions.iloc[idx, 2]
        name = regions.iloc[idx, 3]
        strand = regions.iloc[idx, 5]
        seq = group[1]
        f.write(f">{name} {chrom}:{start}-{end}({strand})\n")
        f.write(f"{seq}\n")
