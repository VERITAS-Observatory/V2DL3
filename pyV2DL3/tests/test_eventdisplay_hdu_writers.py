import numpy as np

import pyV2DL3.genHDUList as hdu_list_module
from pyV2DL3.fillEVENTS import fillEVENTS
from pyV2DL3.fillGTI import fillGTI
from pyV2DL3.fillRESPONSE import fillRESPONSE
from pyV2DL3.genHDUList import genHDUlist


class DataSource:
    __irf_to_store__ = {"point-like": True, "full-enclosure": False}

    def __init__(self):
        self.events = {
            "EVENT_ID": np.array([1, 2]), "TIME": np.array([10.0, 11.0]),
            "RA": np.array([20.0, 21.0]), "DEC": np.array([30.0, 31.0]),
            "ENERGY": np.array([1.0, 2.0]), "EVENT_TYPE": np.array([2, 3]),
            "OBS_ID": 42, "DATE-OBS": "2020-01-01T00:00:00",
            "DATE-AVG": "2020-01-01T00:05:00", "DATE-END": "2020-01-01T00:10:00",
            "TSTART": 10.0, "TSTOP": 20.0, "ONTIME": 8.0, "LIVETIME": 7.0,
            "DEADC": 0.875, "OBJECT": "Crab", "RA_OBJ": 83.0, "DEC_OBJ": 22.0,
            "RA_PNT": 83.0, "DEC_PNT": 22.0, "ALT_PNT": 60.0, "AZ_PNT": 180.0,
            "TELLIST": "T1,T2", "N_TELS": 2,
        }
        self.gti = {"goodTimeStart": [10.0], "goodTimeStop": [18.0], "TSTART": 10.0, "TSTOP": 20.0}
        response_table = np.array([([0.1], [1.0])], dtype=[("ENERG_LO", "f4", (1,)), ("ENERG_HI", "f4", (1,))])
        self.response = {
            "EA": response_table, "MIGRATION": response_table,
            "LO_THRES": 0.1, "HI_THRES": 1.0, "RAD_MAX": 0.2,
        }

    def get_evt_data(self):
        return self.events

    def get_gti_data(self):
        return self.gti

    def get_response_data(self):
        return self.response

    def get_source_name(self):
        return "Eventdisplay"

    def get_version(self):
        return "v1"


def test_shared_writers_create_consistent_event_gti_and_response_hdus():
    datasource = DataSource()

    events = fillEVENTS(datasource, save_multiplicity=True, instrument_epoch="V6")
    gti = fillGTI(datasource)
    response = fillRESPONSE(datasource, instrument_epoch="V6")

    assert events.name == "EVENTS"
    assert events.header["OBS_ID"] == 42
    assert events.header["INSTRUME"] == "Epoch V6"
    assert "EVENT_TYPE" in events.columns.names
    assert gti.name == "GTI"
    assert gti.data["START"].tolist() == [10.0]
    assert [hdu.name for hdu in response] == ["EFFECTIVE AREA", "ENERGY DISPERSION"]
    assert {hdu.header["OBS_ID"] for hdu in response} == {events.header["OBS_ID"]}


def test_hdu_list_joins_the_shared_products_in_gadf_order():
    hdus = genHDUlist(DataSource(), save_multiplicity=True, instrument_epoch="V6")
    assert [hdu.name for hdu in hdus] == [
        "PRIMARY", "EVENTS", "GTI", "EFFECTIVE AREA", "ENERGY DISPERSION"
    ]


def test_response_writer_builds_all_full_enclosure_hdus():
    datasource = DataSource()
    datasource.__irf_to_store__ = {"point-like": False, "full-enclosure": True}
    datasource.response["FULL_EA"] = datasource.response["EA"]
    datasource.response["FULL_MIGRATION"] = datasource.response["MIGRATION"]
    datasource.response["PSF"] = datasource.response["EA"]

    response = fillRESPONSE(datasource)

    assert [hdu.name for hdu in response] == ["EFFECTIVE AREA", "ENERGY DISPERSION", "PSF"]
    assert all(hdu.header["HDUCLAS3"] == "FULL-ENCLOSURE" for hdu in response)


def test_root_file_loader_selects_eventdisplay_and_rejects_missing_effective_area(monkeypatch):
    eventdisplay_source = object()
    monkeypatch.setattr(
        "pyV2DL3.eventdisplay.EventDisplayDataSource.EventDisplayDataSource",
        lambda events, irfs: eventdisplay_source,
    )

    assert hdu_list_module.loadROOTFiles("events.root", "irf.root", "Eventdisplay") is eventdisplay_source
    try:
        hdu_list_module.loadROOTFiles("events.root", None, "Eventdisplay")
    except Exception as error:
        assert "effective area" in str(error)
    else:
        raise AssertionError("missing effective-area file must fail")
