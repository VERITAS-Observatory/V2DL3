from unittest.mock import Mock

import numpy as np
import pytest
from astropy.table import MaskedColumn, Table

import pyV2DL3.eventdisplay.fillEVENTS as fill_events
from pyV2DL3.eventdisplay import DBFitsFile, EventDisplayDataSource, IrfExtractor
from pyV2DL3.eventdisplay.fillRESPONSE import (
    __fill_response__,
    fill_effective_area,
    find_camera_offsets,
)
from pyV2DL3.eventdisplay.IrfInterpolator import IrfInterpolator
from pyV2DL3.eventdisplay.util import (
    ZeroLengthEventList,
    bin_centers_to_edges,
    duplicate_dimensions,
    getGTI,
    produce_tel_list,
)


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

    with pytest.raises(ZeroLengthEventList):
        fill_events.__fill_event_list(file, 42, {"Energy": [3.0, 4.0]}, 100.0)


def test_db_fits_reader_converts_masked_values(monkeypatch):
    table = Table()
    table["weather"] = ["A"]
    table.add_column(MaskedColumn([1.0], mask=[True], name="l3_rate_mean"))
    monkeypatch.setattr(DBFitsFile.Table, "read", Mock(return_value=table))

    assert DBFitsFile.read_db_fits_file("run.db.fits") == {
        "weather": "A", "l3_rate_mean": None
    }
    assert DBFitsFile.read_db_fits_file(None) == {}


def test_db_fits_reader_preserves_read_errors(monkeypatch):
    monkeypatch.setattr(DBFitsFile.Table, "read", Mock(side_effect=FileNotFoundError))
    with pytest.raises(FileNotFoundError):
        DBFitsFile.read_db_fits_file("missing.db.fits")


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


def test_response_offset_edges_are_non_empty_and_used_separately_from_centers():
    theta_low, theta_high = find_camera_offsets(np.array([0.5, 1.0]))
    assert np.all(theta_high > theta_low)
    assert np.allclose(theta_low, [0.25, 0.75])
    assert np.allclose(theta_high, [0.75, 1.25])

    interpolator = Mock()
    interpolator.interpolate.return_value = (np.array([10.0, 20.0]), [np.array([0.0, 1.0])])
    table, low_energy, high_energy = fill_effective_area(
        "eff", interpolator, np.array([0.5, 1.0]), 4.0, 20.0, theta_low, theta_high
    )
    assert interpolator.interpolate.call_count == 2
    assert low_energy < high_energy
    assert np.all(table["THETA_HI"][0] > table["THETA_LO"][0])
