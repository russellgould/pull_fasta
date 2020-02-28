import pytest
import pull_fasta as app
import subprocess
from pathlib import Path

nu = "2"
nd = "4"


def test_version():
    assert app.__version__ == "0.1.0"


def test_exit_status():
    cmd = subprocess.run(["python", "-m", "pull_fasta.main", "-h"])
    assert cmd.check_returncode


def test_peak():
    infile = Path("tests/data/regions.peaks").resolve()
    cmd = subprocess.run(
        ["python", "-m", "pull_fasta.main", "-peak", infile, "-nu", nu, "-nd", nd],
        capture_output=True,
    )
    test_str = "chromstartendnamescorestrand0Chr1411AT3G092600+1Chr129AT3G092600-"
    assert "".join(cmd.stdout.decode().split()) == test_str


def test_gff():
    infile = Path("tests/data/regions.gff").resolve()
    cmd = subprocess.run(
        ["python", "-m", "pull_fasta.main", "-gff", infile, "-nu", nu, "-nd", nd],
        capture_output=True,
    )
    test_str = "chromstartendnamescorestrand0Chr167..+1Chr167..-"
    assert "".join(cmd.stdout.decode().split()) == test_str


def test_bed():
    infile = Path("tests/data/regions.bed").resolve()
    cmd = subprocess.run(
        ["python", "-m", "pull_fasta.main", "-bed", infile, "-nu", nu, "-nd", nd],
        capture_output=True,
    )
    test_str = "chromstartendnamescorestrand0Chr167.0+1Chr167.0-"
    assert "".join(cmd.stdout.decode().split()) == test_str
