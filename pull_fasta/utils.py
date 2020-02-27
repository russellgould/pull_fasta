from pandas import read_csv


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
    return df.rename(columns={i: cols[i] for i in range(df.shape[1])})


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
