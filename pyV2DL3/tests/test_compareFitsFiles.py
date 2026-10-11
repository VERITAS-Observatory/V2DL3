import numpy as np
from astropy.io import fits
from click.testing import CliRunner

from pyV2DL3.script.compareFitsFiles import cli


def write_fits(path, value):
    fits.PrimaryHDU(data=np.array([value], dtype=np.int16)).writeto(path)


def test_compare_fits_returns_zero_for_identical_files(tmp_path):
    file1 = tmp_path / "one.fits"
    file2 = tmp_path / "two.fits"
    diff_file = tmp_path / "diff.txt"
    write_fits(file1, 1)
    write_fits(file2, 1)

    result = CliRunner().invoke(
        cli, ["--file_pair", str(file1), str(file2), "--diff_file", str(diff_file)]
    )

    assert result.exit_code == 0
    assert "FITS files are identical" in result.output


def test_compare_fits_returns_nonzero_for_different_files(tmp_path):
    file1 = tmp_path / "one.fits"
    file2 = tmp_path / "two.fits"
    diff_file = tmp_path / "diff.txt"
    write_fits(file1, 1)
    write_fits(file2, 2)

    result = CliRunner().invoke(
        cli, ["--file_pair", str(file1), str(file2), "--diff_file", str(diff_file)]
    )

    assert result.exit_code == 1
    assert "FITS files differ" in result.output
    assert str(diff_file) in result.output
    assert "Data contains differences" in diff_file.read_text()
