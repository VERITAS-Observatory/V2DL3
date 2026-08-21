from unittest.mock import Mock

import numpy as np
import pytest
from astropy.table import MaskedColumn, Table

import pyV2DL3.eventdisplay.fillEVENTS as fill_events
from pyV2DL3.eventdisplay import DBFitsFile, EventDisplayDataSource, IrfExtractor
from pyV2DL3.eventdisplay.fillRESPONSE import (
    __fill_response__,
    fill_direction_migration,
    fill_effective_area,
    fill_energy_migration,
    find_camera_offsets,
)
from pyV2DL3.eventdisplay.IrfInterpolator import IrfInterpolator
from pyV2DL3.eventdisplay.util import (
    ZeroLengthEventList,
    bin_centers_to_edges,
    duplicate_dimensions,
    get_root_log_lines,
    getGTI,
    produce_tel_list,
)


class Branch:
    def __init__(self, values):
        self.values = np.asarray(values)

    def array(self, library=None):
        return self.values


class RootFile(dict):
    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


class RootLog:
    def __init__(self, lines):
        self.lines = lines

    def member(self, name):
        if name != "fLines":
            raise KeyError(name)
        return self.lines


def effective_area_tree():
    values = {
        "azMin": [-180.0, 0.0, 180.0], "azMax": [-180.0, 0.0, 180.0],
        "az": [0], "pedvar": [4.0], "ze": [20.0], "Woff": [0.5],
        "e0": [[0.0, 1.0]], "eff": [[10.0, 20.0]],
        "hEsysMCRelative2D_binsx": [2], "hEsysMCRelative2D_minx": [-1.0],
        "hEsysMCRelative2D_maxx": [1.0], "hEsysMCRelative2D_binsy": [2],
        "hEsysMCRelative2D_miny": [0.0], "hEsysMCRelative2D_maxy": [2.0],
        "hEsysMCRelative2D_value": [[1.0, 2.0, 3.0, 4.0]],
    }
    return {name: Branch(value) for name, value in values.items()}


def test_eventdisplay_utilities_cover_telescope_bins_dimensions_and_empty_gti():
    assert produce_tel_list({"TelType": [1, 3, 4]}) == "T1,T3,T4"
    edges, low, high = bin_centers_to_edges(np.array([1.0, 3.0]), logaxis=False)
    assert np.array_equal(edges, [0.0, 2.0, 4.0])
    assert np.array_equal(low, [0.0, 2.0])
    assert np.array_equal(high, [2.0, 4.0])
    assert duplicate_dimensions(np.array([[7.0]])).shape == (2, 2)

    start, stop, ontime = getGTI(np.array([0, 0], dtype=np.uint8), 100.0)
    assert ontime == 0
    assert start.size == stop.size == 0


def test_root_log_adapter_uses_public_member_and_decodes_lines():
    assert get_root_log_lines(RootLog([b"first", "second"])) == ["first", "second"]


def test_eventdisplay_version_parser_reports_missing_version_line(monkeypatch):
    root_file = RootFile({"anasumLog;1": RootLog(["unrelated line"])})
    monkeypatch.setattr(EventDisplayDataSource.uproot, "open", Mock(return_value=root_file))
    source = EventDisplayDataSource.EventDisplayDataSource("events.root", "irf.root")

    with pytest.raises(ValueError, match="VERITAS Analysis Summary"):
        source.get_version()


def test_missing_time_mask_falls_back_to_full_run_interval():
    file = {"run_42": {"stereo": {}}}

    start, stop, ontime = fill_events.__get_ontime(file, 42, 10.0, 20.0)

    assert start == [10.0]
    assert stop == [20.0]
    assert ontime == 10.0


def test_event_list_helpers_apply_selection_and_reject_an_empty_result():
    tree = {
        "eventNumber": np.array([1, 2]),
        "timeOfDay": np.array([2.0, 3.0]),
        "RA": np.array([10.0, 20.0]),
        "DEC": np.array([30.0, 40.0]),
        "El": np.array([50.0, 60.0]),
        "Az": np.array([70.0, 80.0]),
        "Energy": np.array([1.0, 2.0]),
        "NImages": np.array([2, 3]),
        "Xoff": np.array([0.1, 0.2]),
        "Yoff": np.array([0.3, 0.4]),
        "ImgSel": np.array([3, 7]),
        "MeanPedvar": np.array([4.0, 6.0]),
    }

    file = {"run_42/stereo/DL3EventTree": Mock()}
    file["run_42/stereo/DL3EventTree"].arrays.return_value = tree
    events, max_img_sel, pedvar = fill_events.__fill_event_list(
        file, 42, {"Energy": [1.5, 2.5]}, 100.0
    )
    assert events["EVENT_ID"].tolist() == [2]
    assert events["TIME"].tolist() == [103.0]
    assert max_img_sel == 7
    assert pedvar == 6.0
    metadata = fill_events.__get_run_event_metadata(file, 42)
    assert metadata["altitude"] == 55.0
    assert metadata["azimuth"] == pytest.approx(75.0)
    assert metadata["max_img_sel"] == 7
    assert metadata["pedvar"] == 5.0

    with pytest.raises(ZeroLengthEventList):
        fill_events.__fill_event_list(file, 42, {"Energy": [3.0, 4.0]}, 100.0)


def test_db_fits_reader_converts_masked_values(monkeypatch):
    table = Table()
    table["runNumber"] = [41, 42]
    table["weather"] = ["B", "A"]
    table.add_column(
        MaskedColumn([1.0, 2.0], mask=[False, True], name="l3_rate_mean")
    )
    table["unsupported"] = [1, 2]
    monkeypatch.setattr(DBFitsFile.Table, "read", Mock(return_value=table))

    assert DBFitsFile.read_db_fits_file("run.db.fits", 42) == {
        "weather": "A", "l3_rate_mean": None
    }
    assert DBFitsFile.read_db_fits_file(None) == {}


def test_db_fits_reader_rejects_missing_or_ambiguous_run_metadata(monkeypatch):
    table = Table()
    table["runNumber"] = [41, 43]
    table["weather"] = ["B", "A"]
    monkeypatch.setattr(DBFitsFile.Table, "read", Mock(return_value=table))

    with pytest.raises(ValueError, match="expected exactly one"):
        DBFitsFile.read_db_fits_file("run.db.fits", 42)

    table["runNumber"] = [42, 42]
    with pytest.raises(ValueError, match="expected exactly one"):
        DBFitsFile.read_db_fits_file("run.db.fits", 42)


def test_db_fits_reader_rejects_core_metadata_collisions(monkeypatch):
    table = Table()
    table["runNumber"] = [42]
    table["OBS_ID"] = [999]
    monkeypatch.setattr(DBFitsFile.Table, "read", Mock(return_value=table))

    with pytest.raises(ValueError, match="overwrite core"):
        DBFitsFile.read_db_fits_file(
            "run.db.fits", 42, protected_keys={"OBS_ID"}
        )


def test_db_fits_reader_preserves_read_errors(monkeypatch):
    monkeypatch.setattr(
        DBFitsFile.Table, "read", Mock(side_effect=FileNotFoundError)
    )
    with pytest.raises(FileNotFoundError):
        DBFitsFile.read_db_fits_file("missing.db.fits", 42)


def test_data_source_passes_event_and_response_options(monkeypatch, tmp_path):
    event_result = (
        {"goodTimeStart": [1], "goodTimeStop": [2]},
        {"azimuth": 0.0, "zenith": 20.0, "pedvar": 4.0},
        {"OBS_ID": 42},
    )
    fill_events_mock = Mock(return_value=event_result)
    fill_response_mock = Mock(return_value={"EA": "response"})
    monkeypatch.setattr(EventDisplayDataSource, "__fillEVENTS__", fill_events_mock)
    monkeypatch.setattr(EventDisplayDataSource, "__fill_response__", fill_response_mock)

    source = EventDisplayDataSource.EventDisplayDataSource("events.root", "irf.root")
    source.fill_data(
        evt_filter=None,
        db_fits_file="run.db.fits",
        interpolator_name="RegularGridInterpolator",
        force_extrapolation=True,
    )

    fill_events_mock.assert_called_once_with("events.root", {}, "run.db.fits")
    fill_response_mock.assert_called_once_with(
        "events.root", "irf.root", 0.0, 20.0, 4.0,
        {"point-like": True, "full-enclosure": False},
        evt_filter=None,
        db_fits_file="run.db.fits",
        interpolator_name="RegularGridInterpolator",
        force_extrapolation=True,
    )
    assert source.get_evt_data() == {"OBS_ID": 42}
    assert source.get_response_data() == {"EA": "response"}


def test_event_helpers_keep_times_directions_and_filter_errors_explicit():
    start, stop, offset = fill_events.__get_times_since_reference_time(
        fill_events.Time(51544.0, format="mjd", scale="utc"),
        fill_events.Time(51544.5, format="mjd", scale="utc"),
    )
    assert stop - start == 43200.0
    assert offset == start
    altitude, azimuth = fill_events.__get_average_event_direction(
        np.array([40.0, 60.0]), np.array([359.0, 1.0])
    )
    assert altitude == 50.0
    assert azimuth == pytest.approx(360.0)
    assert fill_events.__get_time_vector(np.array([1.0, 2.0]), 10.0).tolist() == [11.0, 12.0]
    with pytest.raises(ValueError):
        fill_events.__get_time_vector(np.array([86401.0]), 0.0)
    with pytest.raises(TypeError):
        fill_events.__get_mask({"RA": np.array([1.0])}, {"RA": "bad"})


def test_event_builder_combines_run_metadata_events_gti_and_db_metadata(monkeypatch):
    run_summary = Mock()
    run_summary.arrays.return_value = {
        "runOn": np.array([42]), "DeadTimeFracOn": np.array([0.1]),
        "TargetName": np.array(["Crab"]), "TargetRAJ2000": np.array([83.0]),
        "TargetDecJ2000": np.array([22.0]),
    }
    tel_config = Mock()
    tel_config.arrays.return_value = {"TelType": np.array([1, 2])}
    root_file = RootFile({"total_1/stereo/tRunSummary": run_summary, "run_42/stereo/telconfig": tel_config})
    monkeypatch.setattr(fill_events.uproot, "open", Mock(return_value=root_file))
    start = fill_events.Time(60000.0, format="mjd", scale="utc")
    stop = fill_events.Time(60000.1, format="mjd", scale="utc")
    monkeypatch.setattr(fill_events, "__get_start_stop_times", Mock(return_value=(start, stop, start)))
    monkeypatch.setattr(fill_events, "__get_times_since_reference_time", Mock(return_value=(10.0, 20.0, 0.0)))
    monkeypatch.setattr(
        fill_events, "__fill_event_list",
        Mock(return_value=({"ALT": np.array([50.0]), "AZ": np.array([180.0])}, 3, 4.0)),
    )
    monkeypatch.setattr(
        fill_events,
        "__get_run_event_metadata",
        Mock(return_value={"altitude": 55.0, "azimuth": 190.0, "max_img_sel": 7, "pedvar": 5.0}),
    )
    monkeypatch.setattr(fill_events, "__get_average_pointing", Mock(return_value=(83.0, 22.0)))
    monkeypatch.setattr(fill_events, "__get_average_event_direction", Mock(return_value=(50.0, 180.0)))
    monkeypatch.setattr(fill_events, "__read_quality_flag_from_log", Mock(return_value=0))
    monkeypatch.setattr(fill_events, "__get_ontime", Mock(return_value=([10.0], [18.0], 8.0)))
    db_reader = Mock(return_value={"weather": "A"})
    monkeypatch.setattr(fill_events, "read_db_fits_file", db_reader)

    gti, irf_query, events = fill_events.__fillEVENTS__("events.root", db_fits_file="run.db.fits")

    assert gti == {"goodTimeStart": [10.0], "goodTimeStop": [18.0], "TSTART": 10.0, "TSTOP": 20.0}
    assert irf_query == {"azimuth": 190.0, "zenith": 35.0, "pedvar": 5.0}
    assert events["OBS_ID"] == 42
    assert events["LIVETIME"] == pytest.approx(7.2)
    assert events["TELLIST"] == "T1,T2"
    assert events["weather"] == "A"
    db_reader.assert_called_once_with(
        "run.db.fits", 42, protected_keys=events.keys()
    )


def test_irf_extraction_and_interpolation_accept_north_azimuth(tmp_path, monkeypatch):
    extract_1d = Mock(return_value=("data", "axes"))
    monkeypatch.setattr(IrfExtractor, "extract_irf_1d", extract_1d)
    assert IrfExtractor.extract_irf("irf.root", "eff", azimuth=0.0, irf1d=True) == (
        "data", "axes"
    )
    with pytest.raises(ValueError):
        IrfExtractor.extract_irf("irf.root", "eff", azimuth=None, irf1d=True)

    filename = tmp_path / "irf.root"
    filename.touch()
    interpolator = IrfInterpolator(str(filename), 0.0, "KNeighborsRegressor")
    interpolator.irf_name = "eff"
    interpolator.irf_axes = [np.array([1.0, 2.0])]
    interpolator.interpolator = Mock()
    interpolator.interpolator.predict.return_value = np.array([4.0, 5.0])

    values, axis = interpolator.interpolate([3.0, 20.0, 0.5])
    assert values.tolist() == [4.0, 5.0]
    assert np.array_equal(axis[0], [1.0, 2.0])
    assert interpolator.interpolator.predict.call_args.args[0].shape == (2, 4)


def test_irf_extractor_builds_regular_and_knn_input_arrays(monkeypatch):
    tree = effective_area_tree()
    monkeypatch.setattr(IrfExtractor.uproot, "open", Mock(return_value={"fEffAreaH2F": tree}))

    one_dimensional, one_axes = IrfExtractor.extract_irf_1d("irf.root", "eff", azimuth=0.0)
    two_dimensional, two_axes = IrfExtractor.extract_irf_2d(
        "irf.root", "hEsysMCRelative2D", azimuth=0.0
    )
    coordinates, values = IrfExtractor.extract_irf_for_knn(
        "irf.root", "eff", irf1d=True, azimuth=0.0
    )

    assert one_dimensional.shape == (2, 1, 1, 1)
    assert np.array_equal(one_dimensional[:, 0, 0, 0], [10.0, 20.0])
    assert one_axes["energies"].tolist() == [0.0, 1.0]
    assert two_dimensional.shape == (2, 2, 1, 1, 1)
    assert two_dimensional[:, :, 0, 0, 0].tolist() == [[1.0, 3.0], [2.0, 4.0]]
    assert two_axes["irf_dimension_1"].tolist() == [-0.5, 0.5]
    assert coordinates.shape == (2, 4)
    assert values.tolist() == [10.0, 20.0]


def test_irf_interpolator_loads_knn_coordinates(tmp_path, monkeypatch):
    filename = tmp_path / "irf.root"
    filename.touch()
    coordinates = np.array([
        [1.0, 1.0, 0.5, 0.0], [1.0, 1.0, 0.5, 1.0],
        [2.0, 1.0, 0.5, 0.0], [2.0, 1.0, 0.5, 1.0],
        [3.0, 1.0, 0.5, 0.0],
    ])
    extract = Mock(return_value=(coordinates, np.arange(5.0)))
    monkeypatch.setattr("pyV2DL3.eventdisplay.IrfInterpolator.extract_irf_for_knn", extract)

    interpolator = IrfInterpolator(str(filename), 10.0, "KNeighborsRegressor")
    interpolator.set_irf("eff")

    extract.assert_called_once_with(str(filename), "eff", irf1d=True, azimuth=10.0)
    assert np.array_equal(interpolator.irf_axes[0], [0.0, 1.0])


def test_regular_grid_interpolator_loads_and_interpolates(tmp_path, monkeypatch):
    filename = tmp_path / "irf.root"
    filename.touch()
    axes = {
        "energies": np.array([0.0, 1.0]), "pedvars": np.array([4.0]),
        "zeniths": np.array([20.0]), "woffs": np.array([0.5]),
    }
    data = np.array([[[[10.0]]], [[[20.0]]]])
    monkeypatch.setattr(
        "pyV2DL3.eventdisplay.IrfInterpolator.extract_irf", Mock(return_value=(data, axes))
    )

    interpolator = IrfInterpolator(str(filename), 0.0, "RegularGridInterpolator")
    interpolator.set_irf("eff", use_click=False)
    values, axis = interpolator.interpolate([4.0, 20.0, 0.5])

    assert np.allclose(values, [10.0, 20.0])
    assert np.array_equal(axis[0], [0.0, 1.0])


def test_response_builder_selects_the_requested_response_mode(monkeypatch):
    irf = Mock()
    monkeypatch.setattr("pyV2DL3.eventdisplay.fillRESPONSE.IrfInterpolator", Mock(return_value=irf))
    fast_eff_area = {
        "Woff": Mock(), "ze": Mock(), "pedvar": Mock(),
    }
    fast_eff_area["Woff"].array.return_value = np.array([0.5, 1.0])
    fast_eff_area["ze"].array.return_value = np.array([20.0])
    fast_eff_area["pedvar"].array.return_value = np.array([4.0])
    run_summary = Mock()
    run_summary.arrays.return_value = {"Theta2Max": np.array([0.04])}
    monkeypatch.setattr(
        "pyV2DL3.eventdisplay.fillRESPONSE.uproot.open",
        Mock(return_value={"fEffAreaH2F": fast_eff_area, "total_1/stereo/tRunSummary": run_summary}),
    )
    point_area = Mock(return_value=("EA", 0.1, 1.0))
    point_migration = Mock(return_value="MIGRATION")
    monkeypatch.setattr("pyV2DL3.eventdisplay.fillRESPONSE.fill_effective_area", point_area)
    monkeypatch.setattr("pyV2DL3.eventdisplay.fillRESPONSE.fill_energy_migration", point_migration)

    response = __fill_response__(
        "events.root", "irf.root", 10.0, 20.0, 4.0,
        {"point-like": True, "full-enclosure": False}, use_click=False,
    )

    assert response["EA"] == "EA"
    assert response["MIGRATION"] == "MIGRATION"
    assert point_area.call_args.args[0] == "eff"


def test_response_builder_creates_all_full_enclosure_products(monkeypatch):
    fast_eff_area = {"Woff": Mock(), "ze": Mock(), "pedvar": Mock()}
    fast_eff_area["Woff"].array.return_value = np.array([0.5, 1.0])
    fast_eff_area["ze"].array.return_value = np.array([20.0])
    fast_eff_area["pedvar"].array.return_value = np.array([4.0])
    monkeypatch.setattr(
        "pyV2DL3.eventdisplay.fillRESPONSE.uproot.open",
        Mock(return_value={"fEffAreaH2F": fast_eff_area}),
    )
    area = Mock(return_value=("EA", 0.1, 1.0))
    migration = Mock(return_value="MIGRATION")
    direction = Mock(return_value="PSF")
    monkeypatch.setattr("pyV2DL3.eventdisplay.fillRESPONSE.IrfInterpolator", Mock())
    monkeypatch.setattr("pyV2DL3.eventdisplay.fillRESPONSE.fill_effective_area", area)
    monkeypatch.setattr("pyV2DL3.eventdisplay.fillRESPONSE.fill_energy_migration", migration)
    monkeypatch.setattr("pyV2DL3.eventdisplay.fillRESPONSE.fill_direction_migration", direction)

    response = __fill_response__(
        "events.root", "irf.root", 10.0, 20.0, 4.0,
        {"point-like": False, "full-enclosure": True}, use_click=False,
    )

    assert response == {"FULL_EA": "EA", "LO_THRES": 0.1, "HI_THRES": 1.0,
                        "FULL_MIGRATION": "MIGRATION", "PSF": "PSF"}
    assert area.call_args.args[0] == "effNoTh2"
    assert migration.call_args.args[0] == "hEsysMCRelative2DNoDirectionCut"


def test_response_offset_coordinates_are_preserved_in_serialized_irfs():
    theta_low, theta_high = find_camera_offsets(np.array([0.5, 1.0]))
    assert np.array_equal(theta_low, [0.5, 1.0])
    assert np.array_equal(theta_high, [0.5, 1.0])

    interpolator = Mock()
    interpolator.interpolate.return_value = (np.array([10.0, 20.0]), [np.array([0.0, 1.0])])
    table, low_energy, high_energy = fill_effective_area(
        "eff", interpolator, np.array([0.5, 1.0]), 4.0, 20.0, theta_low, theta_high
    )
    assert interpolator.interpolate.call_count == 2
    assert low_energy < high_energy
    assert np.array_equal(table["THETA_LO"][0], theta_low)
    assert np.array_equal(table["THETA_HI"][0], theta_high)


def test_response_array_builders_normalize_and_preserve_axis_shapes():
    interpolator = Mock()
    interpolator.interpolate.return_value = (
        np.array([[1.0, 1.0], [1.0, 1.0]]),
        [np.array([0.0, 1.0]), np.array([0.0, 1.0])],
    )
    migration = fill_energy_migration(
        "hEsysMCRelative2D", interpolator, [0.5, 1.0], 4.0, 20.0,
        np.array([0.25, 0.75]), np.array([0.75, 1.25]),
    )
    assert migration["MATRIX"].shape == (1, 2, 2, 2)
    assert np.all(migration["MIGRA_HI"][0] > migration["MIGRA_LO"][0])

    interpolator.interpolate.return_value = (
        np.ones((2, 6)),
        [np.array([-2.0, -1.0, 0.0, 1.0, 2.0, 3.0]), np.array([-2.0, -1.0])],
    )
    psf = fill_direction_migration(
        interpolator, [0.5, 1.0], 4.0, 20.0,
        np.array([0.25, 0.75]), np.array([0.75, 1.25]),
    )
    assert psf["RPSF"].shape[1] == 2
    assert np.all(psf["RAD_HI"][0] > psf["RAD_LO"][0])
