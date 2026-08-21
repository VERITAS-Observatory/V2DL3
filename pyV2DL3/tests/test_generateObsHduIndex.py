
import os
from unittest.mock import Mock

import pytest
from astropy.io import fits
from click.testing import CliRunner

from pyV2DL3.generateObsHduIndex import (
    _add_table_units,
    _check_unit_consistency,
    _default_null_value,
    _fill_table_data,
    gen_hdu_index,
    get_unit_string_from_comment,
)
from pyV2DL3.script import generate_index_file


class TestGetUnitStringFromComment:

    # Returns the unit string when it is enclosed in square brackets.
    def test_returns_unit_string_when_enclosed_in_square_brackets(self):
        comment_string = 'average pointing azimuth [deg]'
        assert get_unit_string_from_comment(comment_string) == 'deg'

    # Returns None when the unit string is not enclosed in square brackets.
    def test_returns_none_when_unit_string_not_enclosed_in_square_brackets(self):
        comment_string = 'average pointing azimuth'
        assert get_unit_string_from_comment(comment_string) is None

    # Returns 'uA' when the unit string is 'muA'.
    def test_returns_uA_when_unit_string_is_muA(self):
        comment_string = 'current [muA]'
        assert get_unit_string_from_comment(comment_string) == 'uA'


class Test_FillTableData:

    # Appends value to the list for a given key in the header
    def test_appends_value_to_list(self):
        key = "RA_PNT"
        value = []
        data_type = ">f4"
        header = fits.Header()
        header["RA_PNT"] = 180.0

        result = _fill_table_data(key, value, data_type, header, dqm_header=False)

        assert result is True
        assert value == [180.0]

    # Handles non-null values in the header
    def test_handles_non_null_values(self):
        key = "RA_PNT"
        value = []
        data_type = ">f4"
        header = fits.Header()
        header["RA_PNT"] = 180.0

        result = _fill_table_data(key, value, data_type, header, dqm_header=False)

        assert result is True
        assert value == [180.0]

    # Handles null values in the header
    def test_handles_null_values(self):
        key = "RA_PNT"
        value = []
        data_type = ">f4"
        header = fits.Header()
        header["RA_PNT"] = 'NULL'

        result = _fill_table_data(key, value, data_type, header, dqm_header=False)

        assert result is True
        assert value == [-9999.0]

    # Handles missing key in the header when dqm_header is False
    def test_handles_missing_key_when_dqm_header_is_false(self):
        key = "RA_PNT"
        value = []
        data_type = ">f4"
        header = fits.Header()

        result = _fill_table_data(key, value, data_type, header, dqm_header=False)

        assert result is False
        assert value == []

    # Handles missing key in the header when dqm_header is True
    def test_handles_missing_key_when_dqm_header_is_true(self):
        key = "RA_PNT"
        value = []
        data_type = ">f4"
        header = fits.Header()

        result = _fill_table_data(key, value, data_type, header, dqm_header=True)

        assert result is True
        assert value == [-9999.0]

    # Handles key 'ZEN_PNT' by calculating its value from 'ALT_PNT' in the header
    def test_handles_key_zen_pnt(self):
        key = "ZEN_PNT"
        value = []
        data_type = ">f4"
        header = fits.Header()
        header["ALT_PNT"] = 40.0

        result = _fill_table_data(key, value, data_type, header, dqm_header=False)

        assert result is True
        assert value == [50.0]


class Test_AddTableUnits:

    # Unit string is found in header comments and added to _table_units dictionary
    def test_unit_string_found(self):
        header = fits.Header()
        header['ALT_PNT'] = 'altitude pointing'
        header.comments['ALT_PNT'] = 'altitude pointing comment [deg]'

        _table_units = {}
        key = 'ALT_PNT'

        _add_table_units(key, _table_units, header)

        assert _table_units == {'ALT_PNT': 'deg'}

    # Key is "ZEN_PNT" and corresponding unit string is found in
    # header comments and added to _table_units dictionary
    def test_zen_pnt_unit_string_found(self):
        header = fits.Header()
        header['ALT_PNT'] = 'altitude pointing'
        header.comments['ALT_PNT'] = 'altitude pointing comment [deg]'
        header['ZEN_PNT'] = 'zenith pointing'
        header.comments['ZEN_PNT'] = 'zenith pointing comment [deg]'

        _table_units = {}
        key = 'ZEN_PNT'

        _add_table_units(key, _table_units, header)

        assert _table_units == {'ZEN_PNT': 'deg'}

    # Key is not in header comments and not in _table_units dictionary,
    # but no error is raised and _table_units remains unchanged
    def test_key_not_in_comments_or_table_units(self):
        header = fits.Header()
        header['ALT_PNT'] = 'altitude pointing'
        header.comments['ALT_PNT'] = 'altitude pointing comment [deg]'

        _table_units = {}
        key = 'RA_PNT'

        _add_table_units(key, _table_units, header)

        assert _table_units == {'RA_PNT': None}

    # Key is not in header comments and not in _table_units dictionary,
    # and None is added to _table_units dictionary
    def test_key_not_in_comments_or_table_units_none(self):
        header = fits.Header()
        header['ALT_PNT'] = 'altitude pointing'
        header.comments['ALT_PNT'] = 'altitude pointing comment [deg]'

        _table_units = {}
        key = 'RA_PNT'

        _add_table_units(key, _table_units, header)

        assert _table_units == {'RA_PNT': None}


class Test_CheckUnitConsistency:

    # Function returns the same value and units when no conversion is needed
    def test_same_value_and_units(self):
        key = "RA_PNT"
        value = [1.0, 2.0, 3.0]
        units = "deg"

        result = _check_unit_consistency(key, value, units)

        assert result == (value, units)

    # Function converts mph to km/h when key contains "WIND"
    def test_convert_mph_to_kmh(self):
        key = "WINDSPE"
        value = [10, 20, 30]
        units = "mph"

        result = _check_unit_consistency(key, value, units)

        expected_value = [v * 1.60934 for v in value]
        expected_units = "km/h"

        assert result[0] == pytest.approx(expected_value, abs=1.e-3)
        assert result[1] == expected_units

    # Function returns None when value is None
    def test_value_is_none(self):
        key = "RA_PNT"
        value = None
        units = "deg"

        result = _check_unit_consistency(key, value, units)

        assert result == (None, units)

    # Function returns None when units is None
    def test_units_is_none(self):
        key = "RA_PNT"
        value = [1.0, 2.0, 3.0]
        units = None

        result = _check_unit_consistency(key, value, units)

        assert result == (value, None)

    # Function returns empty list when value is empty
    def test_value_is_empty(self):
        key = "RA_PNT"
        value = []
        units = "deg"

        result = _check_unit_consistency(key, value, units)

        assert result == ([], units)

    # Function returns empty list when units is empty
    def test_units_is_empty(self):
        key = "RA_PNT"
        value = [1.0, 2.0, 3.0]
        units = ""

        result = _check_unit_consistency(key, value, units)

        assert result == (value, "")


class Test_DefaultNullValue:

    # Returns -9999.0 for data_type '>f4'
    def test_returns_minus_9999_0_for_data_type_f4(self):
        result = _default_null_value('>f4')
        assert result == -9999.0

    # Returns -9999 for data_type '>i8'
    def test_returns_minus_9999_for_data_type_i8(self):
        result = _default_null_value('>i8')
        assert result == -9999

    # Returns empty string for any other data_type
    def test_returns_empty_string_for_any_other_data_type(self):
        result = _default_null_value('any_other_data_type')
        assert result == ''

    # Returns empty string for invalid data_type
    def test_returns_empty_string_for_invalid_data_type(self):
        result = _default_null_value('invalid_data_type')
        assert result == ''

    # Returns empty string for None data_type
    def test_returns_empty_string_for_none_data_type(self):
        result = _default_null_value(None)
        assert result == ''

    # Covers all possible data_type values
    def test_covers_all_possible_data_type_values(self):
        data_types = ['>f4', '>i8', 'any_other_data_type', 'invalid_data_type', None]
        for data_type in data_types:
            result = _default_null_value(data_type)
            assert isinstance(result, (float, int, str))


def test_hdu_index_preserves_long_relative_paths(tmp_path):
    source_dir = tmp_path / ("nested_" + "x" * 45)
    source_dir.mkdir()
    source = source_dir / ("observation_" + "y" * 50 + ".fits")
    events = fits.BinTableHDU.from_columns([])
    events.name = "EVENTS"
    events.header["OBS_ID"] = 42
    events.header["HDUCLAS1"] = "EVENTS"
    fits.HDUList([fits.PrimaryHDU(), events]).writeto(source)

    index = gen_hdu_index([str(source)], str(tmp_path))
    row = index.data[0]
    expected_dir = os.path.dirname(os.path.relpath(source, tmp_path))
    expected_name = os.path.basename(source)
    actual_dir = row["FILE_DIR"].decode() if isinstance(row["FILE_DIR"], bytes) else row["FILE_DIR"]
    actual_name = row["FILE_NAME"].decode() if isinstance(row["FILE_NAME"], bytes) else row["FILE_NAME"]

    assert actual_dir == expected_dir
    assert actual_name == expected_name


def test_index_cli_checks_requested_output_names(tmp_path, monkeypatch):
    (tmp_path / "obs-index.fits.gz").touch()
    monkeypatch.setattr(generate_index_file.glob, "glob", lambda pattern: ["input.fits"])
    create_mock = Mock()
    monkeypatch.setattr(generate_index_file, "create_obs_hdu_index_file", create_mock)
    result = CliRunner().invoke(
        generate_index_file.cli,
        [
            "--folder_location", str(tmp_path),
            "--index_file_dir", str(tmp_path),
            "--obs_index_file", "custom-obs.fits.gz",
            "--hdu_index_file", "custom-hdu.fits.gz",
        ],
    )

    assert result.exit_code == 0, result.output
    create_mock.assert_called_once()
