from pandas import read_csv, Series


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
                "GeneName": "gene_name",
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


def read_gff(file_path: str):
    """Load a GFF file and return a pandas DataFrame."""

    def format_attributes(df):
        for item in df["attributes"]:
            yield {
                x.split("=")[0]: x.split("=")[1]
                for x in [sub_item for sub_item in item.split(";")]
            }

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
        df["name"] = Series(x["ID"] for x in df["attributes"])
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
