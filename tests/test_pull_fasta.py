import pytest
import pull_fasta as app
import subprocess
from pathlib import Path


def test_version():
    assert app.__version__ == "0.1.0"


def test_exit_status():
    cmd = subprocess.run(["python", "-m", "pull_fasta.main", "-h"])
    assert cmd.check_returncode


def test_peak():
    infile = "/Users/user/Documents/Lab/pull_fasta/tests/data/regions.peaks"
    cmd = subprocess.run(
        ["python", "-m", "pull_fasta.main", "-peak", infile, "-nu", "2", "-nd", "4"],
        capture_output=True,
    )
    test_str = "chromstartendnamescorestrand0Chr1411AT3G092600+1Chr129AT3G092600-"
    assert "".join(cmd.stdout.decode().split()) == test_str
