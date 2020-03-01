import pytest
import pull_fasta as app
import subprocess
from pathlib import Path

nu = "2"
nd = "4"
ref_path = Path("tests/data/reference.fa").resolve()


def test_version():
    assert app.__version__ == "0.1.0"


def test_exit_status():
    cmd = subprocess.run(["python", "-m", "pull_fasta.main", "-h"])
    assert cmd.check_returncode


def test_peak():
    infile = Path("tests/data/example.peaks").resolve()
    cmd = subprocess.run(
        [
            "python",
            "-m",
            "pull_fasta.main",
            "-ref",
            ref_path,
            "-peak",
            infile,
            "-nu",
            nu,
            "-nd",
            nd,
        ],
        capture_output=True,
    )
    test = ">Chr1:4-11(+)AAACCCT>Chr1:2-9(-)GGTTTAG"
    assert "".join(cmd.stdout.decode().split()) == test


def test_gff_regions():
    infile = Path("tests/data/example_regions.gff").resolve()
    cmd = subprocess.run(
        ["python", "-m", "pull_fasta.main", "-gff", infile, "-nu", nu, "-nd", nd],
        capture_output=True,
    )
    test = ">Chr1:4-11(+)AAACCCT>Chr1:2-9(-)GGTTTAG"
    assert "".join(cmd.stdout.decode().split()) == test


def test_bed_regions():
    infile = Path("tests/data/example_regions.bed").resolve()
    cmd = subprocess.run(
        ["python", "-m", "pull_fasta.main", "-bed", infile, "-nu", nu, "-nd", nd],
        capture_output=True,
    )
    test = ">Chr1:4-11(+)AAACCCT>Chr1:2-9(-)GGTTTAG"
    assert "".join(cmd.stdout.decode().split()) == test
