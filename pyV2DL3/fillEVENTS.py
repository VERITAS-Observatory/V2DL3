import logging

from astropy.io import fits

import pyV2DL3.constant as constant
from pyV2DL3.addHDUClassKeyword import addHDUClassKeyword


def add_existing_column(columns, evt_dict, name, format, unit=None):
    """
    Test if key `name` exists in `evt_dict`. If `true` add a new column to columns inplace.

    Parameters
    ----------
    columns: list
    evt_dict: dict
    name: str
        Used as name for the new column and key in evt_dict.
    format: str
    unit: str or None

    Returns
    -------
    """

    if name in evt_dict:
        columns.append(
            fits.Column(name=name, format=format, array=evt_dict[name], unit=unit)
        )


def fillEVENTS(datasource, save_multiplicity=False, instrument_epoch=None, event_class_idx=None):
    logging.debug("Create EVENT HDU")
    evt_dict = datasource.get_evt_data()

    if event_class_idx is not None:
        # Add 'EV_CLASS' header key if there is more than one event group
        add_evclass = True if len(evt_dict) > 1 else False
        evt_dict = evt_dict[event_class_idx]

    # Columns to be saved
    columns = [
        fits.Column(name="EVENT_ID", format="1K", array=evt_dict["EVENT_ID"]),
        fits.Column(name="TIME", format="1D", array=evt_dict["TIME"], unit="s"),
        fits.Column(name="RA", format="1E", array=evt_dict["RA"], unit="deg"),
        fits.Column(name="DEC", format="1E", array=evt_dict["DEC"], unit="deg"),
        fits.Column(name="ENERGY", format="1E", array=evt_dict["ENERGY"], unit="TeV"),
    ]

    # Test if key exists in evt_dict. If yes, add these columns.
    add_existing_column(columns, evt_dict, name="ALT", format="1E", unit="deg")
    add_existing_column(columns, evt_dict, name="AZ", format="1E", unit="deg")
    add_existing_column(columns, evt_dict, name="Xoff", format="1E")
    add_existing_column(columns, evt_dict, name="Yoff", format="1E")
    add_existing_column(columns, evt_dict, name="MSW", format="1D")
    add_existing_column(columns, evt_dict, name="MSL", format="1D")
    add_existing_column(columns, evt_dict, name="IS_GAMMA", format="1L")
    add_existing_column(columns, evt_dict, name="GAMMANESS", format="1E")

    # Number of triggered telescope if necessary
    if save_multiplicity:
        columns.append(
            fits.Column(name="EVENT_TYPE", format="1J", array=evt_dict["EVENT_TYPE"])
        )

    # Create HDU
    hdu1 = fits.BinTableHDU.from_columns(columns)
    hdu1.name = "EVENTS"

    # Fill Standard HDUCLASS headers
    hdu1 = addHDUClassKeyword(hdu1, class1="EVENTS")

    # Fill Header
    hdu1.header.set("RADECSYS", constant.RADECSYS, "equatorial system type")
    hdu1.header.set("EQUINOX", constant.EQUINOX, "base equinox")
    hdu1.header.set(
        "CREATOR",
        "pyV2DL3 v{}::{}".format(constant.VERSION, datasource.get_source_name()),
    )
    hdu1.header.set("ORIGIN", "VERITAS Collaboration", "Data from VERITAS")
    hdu1.header.set("TELESCOP", "VERITAS")
    if instrument_epoch:
        hdu1.header.set("INSTRUME", "Epoch " + instrument_epoch)
    else:
        hdu1.header.set("INSTRUME", "VERITAS")

    hdu1.header.set("OBS_ID  ", evt_dict["OBS_ID"], "Run Number")

    hdu1.header.set(
        "DATE-OBS", evt_dict["DATE-OBS"], "start date (UTC) of obs yy-mm-dd hh:mm:ss"
    )
    hdu1.header.set(
        "DATE-AVG", evt_dict["DATE-AVG"], "average date (UTC) of obs"
    )
    hdu1.header.set(
        "DATE-END", evt_dict["DATE-END"], "end date (UTC) of obs yy-mm-dd hh:mm:ss"
    )

    hdu1.header.set("TSTART  ", evt_dict["TSTART"], "mission time of start of obs [s]")
    hdu1.header.set("TSTOP   ", evt_dict["TSTOP"], "mission time of end of obs [s]")
    hdu1.header.set(
        "MJDREFI ",
        constant.VTS_REFERENCE_MJD,
        "int part of reference MJD [days]",
    )
    hdu1.header.set("MJDREFF ", 0.0, "fractional part of reference MJD [days]")

    hdu1.header.set("TIMEUNIT", "s", "time unit is seconds since MET start")
    hdu1.header.set("TIMESYS ", "utc", "time scale is UTC")
    hdu1.header.set("TIMEREF ", "topocenter", "location from where the observation was made")

    hdu1.header.set("ONTIME  ", evt_dict["ONTIME"], "sum of good time intervals [s]")

    # Correct live time for time cuts
    hdu1.header.set(
        "LIVETIME", evt_dict["LIVETIME"], "(ontime * deadtime time correction) [s] "
    )

    hdu1.header.set(
        "DEADC   ", evt_dict["DEADC"], "Average dead time correction (LIVETIME/ONTIME)"
    )

    hdu1.header.set("OBJECT  ", evt_dict["OBJECT"], "observed object")
    hdu1.header.set("RA_OBJ  ", evt_dict["RA_OBJ"], "right ascension of object [deg]")
    hdu1.header.set("DEC_OBJ ", evt_dict["DEC_OBJ"], "declination of object [deg]")

    hdu1.header.set("RA_PNT  ", evt_dict["RA_PNT"], "pointing right ascension [deg]")
    hdu1.header.set("DEC_PNT ", evt_dict["DEC_PNT"], "pointing declination [deg]")
    hdu1.header.set(
        "ALT_PNT ", evt_dict["ALT_PNT"], "average pointing altitude [deg]"
    )
    hdu1.header.set("AZ_PNT  ", evt_dict["AZ_PNT"], "average pointing azimuth [deg]")

    # get the list of telescopes that participate in the event
    hdu1.header.set("TELLIST", evt_dict["TELLIST"], "comma-separated list of tel IDs")
    hdu1.header.set("N_TELS", evt_dict["N_TELS"], "number of telescopes in event list")

    hdu1.header.set("EUNIT   ", "TeV", "energy unit")
    hdu1.header.set(
        "GEOLON  ",
        constant.VTS_REFERENCE_LON,
        "longitude of array center [deg]",
    )
    hdu1.header.set(
        "GEOLAT  ", constant.VTS_REFERENCE_LAT, "latitude of array center [deg]"
    )
    hdu1.header.set(
        "ALTITUDE",
        constant.VTS_REFERENCE_HEIGHT,
        "altitude of array center [m]",
    )
    if event_class_idx is not None and add_evclass:
        hdu1.header.set("EV_CLASS", event_class_idx, "Event class number")

    if hasattr(datasource, "__pedvar__"):
        hdu1.header.set(
            "PED_VAR",
            datasource.__pedvar__,
            "average pedestal variance",
        )

    try:
        hdu1.header.set(
            "QUALITY",
            evt_dict["QUALITY"],
            "run quality flag based on vpm data used or not",
        )
    except KeyError:
        logging.debug("Keyword QUALITY not set in the EVENTS header")

    try:
        hdu1.header.set(
            "NSBLEVEL",
            evt_dict["NSBLEVEL"],
            "nsb level (mean of pedestal variations)",
        )
    except KeyError:
        logging.debug("Keyword NSBLEVEL not set in the EVENTS header")

    return hdu1
