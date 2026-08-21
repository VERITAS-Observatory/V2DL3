from unittest.mock import Mock

import pytest
from astropy.io import fits
from click.testing import CliRunner

from pyV2DL3.script import v2dl3_for_Eventdisplay as eventdisplay_cli


@pytest.mark.parametrize(
    ("response_option", "expected_irfs"),
    [
        ([], {"full-enclosure": False, "point-like": True}),
        (["--full-enclosure"], {"full-enclosure": True, "point-like": False}),
    ],
)
def test_cli_forwards_conversion_options(
    tmp_path, monkeypatch, response_option, expected_irfs
):
    anasum_file = tmp_path / "run.anasum.root"
    effective_area_file = tmp_path / "effective_area.root"
    filter_file = tmp_path / "selection.yaml"
    db_fits_file = tmp_path / "run.db.fits"
    for path in (anasum_file, effective_area_file, filter_file, db_fits_file):
        path.touch()

    datasource = Mock()
    hdu_list = fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU()])
    hdu_list[1].header["OBS_ID"] = 12345
    load_root_files = Mock(return_value=datasource)
    gen_hdu_list = Mock(return_value=hdu_list)
    monkeypatch.setattr(eventdisplay_cli, "loadROOTFiles", load_root_files)
    monkeypatch.setattr(eventdisplay_cli, "genHDUlist", gen_hdu_list)

    result = CliRunner().invoke(
        eventdisplay_cli.cli,
        [
            "--file_pair", str(anasum_file), str(effective_area_file),
            "--evt_filter", str(filter_file),
            "--force_extrapolation",
            "--fuzzy_boundary", "zenith", "0.1",
            "--db_fits_file", str(db_fits_file),
            "--interpolator_name", "RegularGridInterpolator",
            "--instrument_epoch", "V6",
            "--save_multiplicity",
            *response_option,
            str(tmp_path / "output.fits.gz"),
        ],
    )

    assert result.exit_code == 0, result.output
    load_root_files.assert_called_once_with(
        str(anasum_file), str(effective_area_file), "Eventdisplay"
    )
    datasource.set_irfs_to_store.assert_called_once_with(expected_irfs)
    datasource.fill_data.assert_called_once_with(
        evt_filter=str(filter_file),
        db_fits_file=str(db_fits_file),
        force_extrapolation=True,
        fuzzy_boundary=(("zenith", 0.1),),
        interpolator_name="RegularGridInterpolator",
    )
    gen_hdu_list.assert_called_once_with(
        datasource, save_multiplicity=True, instrument_epoch="V6"
    )
    assert (tmp_path / "output.fits.gz").is_file()
