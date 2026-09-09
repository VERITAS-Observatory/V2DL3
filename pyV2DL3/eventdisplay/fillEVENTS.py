import logging

import astropy.units as u
import numpy as np
import uproot
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from scipy.stats import circmean

from pyV2DL3.constant import (
    VTS_REFERENCE_HEIGHT,
    VTS_REFERENCE_LAT,
    VTS_REFERENCE_LON,
    VTS_REFERENCE_MJD,
)
from pyV2DL3.eventdisplay.DBFitsFile import read_db_fits_file
from pyV2DL3.eventdisplay.util import (
    ZeroLengthEventList,
    getGTI,
    getRunQuality,
    produce_tel_list,
)

logger = logging.getLogger(__name__)


def __fillEVENTS__(edFileIO, select=None, db_fits_file=None):
    """
    Fill event list and event header from anasum file

    """

    with uproot.open(edFileIO) as file:
        runSummary = file["total_1/stereo/tRunSummary"].arrays(library="np")
        runNumber = runSummary["runOn"][0]
        logger.info("Run number: %d", runNumber)

        t_start, t_stop, t_avg = __get_start_stop_times(file)
        t_start_from_reference, t_stop_from_reference, seconds_from_reference = \
            __get_times_since_reference_time(t_start, t_stop)

        event_tree = file[f"run_{runNumber}/stereo/DL3EventTree"].arrays(library="np")
        run_metadata = __get_run_event_metadata(file, runNumber, event_tree=event_tree)
        evt_dict, _, _ =  \
            __fill_event_list(file, runNumber, select, seconds_from_reference, event_tree=event_tree)
        pointing_ra, pointing_dec = __get_average_pointing(file, runNumber)
        pointing_altitude, pointing_azimuth = __get_pointing_altaz(
            pointing_ra, pointing_dec, t_avg
        )

        # Header info
        evt_dict["OBS_ID"] = runNumber
        evt_dict["DATE-OBS"] = t_start.to_value("fits")
        evt_dict["DATE-AVG"] = t_avg.to_value("fits")
        evt_dict["DATE-END"] = t_stop.to_value("fits")
        evt_dict["TSTART"] = t_start_from_reference
        evt_dict["TSTOP"] = t_stop_from_reference
        evt_dict["MJDREFI"] = int(VTS_REFERENCE_MJD)
        evt_dict["DEADC"] = 1 - runSummary["DeadTimeFracOn"][0]
        evt_dict["OBJECT"] = runSummary["TargetName"][0]
        evt_dict["RA_PNT"] = pointing_ra
        evt_dict["DEC_PNT"] = pointing_dec
        evt_dict["ALT_PNT"] = pointing_altitude
        evt_dict["AZ_PNT"] = pointing_azimuth
        evt_dict["RA_OBJ"] = runSummary["TargetRAJ2000"][0]
        evt_dict["DEC_OBJ"] = runSummary["TargetDecJ2000"][0]
        evt_dict["TELLIST"] = produce_tel_list(
            file[f"run_{runNumber}/stereo/telconfig"].arrays(library="np"))
        evt_dict["N_TELS"] = np.binary_repr(run_metadata["max_img_sel"]).count("1")
        logger.info("Number of Telescopes: %d", evt_dict["N_TELS"])
        evt_dict["GEOLON"] = VTS_REFERENCE_LON
        evt_dict["GEOLAT"] = VTS_REFERENCE_LAT
        evt_dict["ALTITUDE"] = VTS_REFERENCE_HEIGHT
        evt_dict["NSBLEVEL"] = run_metadata["pedvar"]
        evt_dict["QUALITY"] = __read_quality_flag_from_log(file, runNumber)
        gti_tstart_from_reference, gti_tstop_from_reference, evt_dict["ONTIME"] = \
            __get_ontime(file, runNumber, t_start_from_reference, t_stop_from_reference)
        evt_dict["LIVETIME"] = evt_dict["ONTIME"] * evt_dict["DEADC"]

    evt_dict.update(
        read_db_fits_file(
            db_fits_file, runNumber, protected_keys=evt_dict.keys()
        )
    )

    return (
        {
            "goodTimeStart": gti_tstart_from_reference,
            "goodTimeStop": gti_tstop_from_reference,
            "TSTART": t_start_from_reference,
            "TSTOP": t_stop_from_reference,
        },
        {
            "azimuth": evt_dict["AZ_PNT"],
            "zenith": (90.0 - evt_dict["ALT_PNT"]),
            "pedvar": evt_dict["NSBLEVEL"],
        },
        evt_dict,
    )


def __fill_event_list(file, runNumber, select, seconds_from_reference, event_tree=None):
    """
    Fill event list from DL3EventTree

    """

    if event_tree is None:
        event_tree = file[f"run_{runNumber}/stereo/DL3EventTree"].arrays(library="np")
    if len(event_tree["eventNumber"]) == 0:
        logger.error("Empty event list")
        raise ZeroLengthEventList

    mask = __get_mask(event_tree, select)
    if not np.any(mask):
        logger.error("Empty event list after selection")
        raise ZeroLengthEventList

    if np.sum(mask) == 0:
        logger.error("Empty event list after applying selection filter")
        raise ZeroLengthEventList

    evt_dict = {}
    evt_dict["EVENT_ID"] = event_tree["eventNumber"][mask]
    evt_dict["TIME"] = __get_time_vector(event_tree["timeOfDay"][mask], seconds_from_reference)
    evt_dict["RA"] = event_tree["RA"][mask]
    evt_dict["DEC"] = event_tree["DEC"][mask]
    evt_dict["ALT"] = event_tree["El"][mask]
    evt_dict["AZ"] = event_tree["Az"][mask]
    evt_dict["ENERGY"] = event_tree["Energy"][mask]
    evt_dict["EVENT_TYPE"] = event_tree["NImages"][mask]
    evt_dict["Xoff"] = event_tree["Xoff"][mask]
    evt_dict["Yoff"] = event_tree["Yoff"][mask]
    try:
        # Test if anasum file was created using the all events option.
        # In this case write out the additional output.
        evt_dict["GAMMANESS"] = event_tree["MVA"][mask]
        evt_dict["IS_GAMMA"] = event_tree["IsGamma"][mask]
    except KeyError:
        pass

    logger.info("Number of events: %d", len(evt_dict["EVENT_ID"]))

    return (
        evt_dict,
        np.max(event_tree["ImgSel"][mask]),
        np.mean(event_tree["MeanPedvar"][mask]),
    )


def __get_run_event_metadata(file, runNumber, event_tree=None):
    """Return run-level event metadata without applying an event selection."""

    if event_tree is None:
        event_tree = file[f"run_{runNumber}/stereo/DL3EventTree"].arrays(library="np")
    if len(event_tree["eventNumber"]) == 0:
        logger.error("Empty event list")
        raise ZeroLengthEventList

    altitude, azimuth = __get_average_event_direction(
        event_tree["El"], event_tree["Az"]
    )
    return {
        "altitude": altitude,
        "azimuth": azimuth,
        "max_img_sel": np.max(event_tree["ImgSel"]),
        "pedvar": np.mean(event_tree["MeanPedvar"]),
    }


def __get_start_stop_times(file):
    """
    Return run start and stop time read from tRunSummary tree

    """

    start_mjd = file["total_1/stereo/tRunSummary/MJDrunstart"].array(library="np")[0]
    stop_mjd = file["total_1/stereo/tRunSummary/MJDrunstop"].array(library="np")[0]

    # convert mjd to fits format
    t_start = Time(start_mjd, format="mjd", scale="utc")
    t_stop = Time(stop_mjd, format="mjd", scale="utc")
    t_avg = t_start + (t_stop - t_start) / 2.0

    return t_start, t_stop, t_avg


def __get_times_since_reference_time(t_start, t_stop):
    """
    Return time since reference time in seconds

    Returns
    -------
    t_start_from_reference
        Time since reference time in seconds at start of run
    t_stop_from_reference
        Time since reference time in seconds at end of run
    seconds_from_reference
        Seconds between reference time and run MJD at 00:00:00:

    """

    t_ref = Time(VTS_REFERENCE_MJD, format="mjd", scale="utc")
    seconds_from_reference = (Time(np.trunc(t_start.mjd), format="mjd", scale="utc") - t_ref).sec

    return (
        (t_start - t_ref).sec,
        (t_stop - t_ref).sec,
        seconds_from_reference
    )


def __get_average_event_direction(altArr, azArr):
    """
    Return average azimuth and elevation events

    """

    avAlt = np.mean(altArr)
    # Calculate average azimuth angle from average vector on a circle
    # https://en.wikipedia.org/wiki/Mean_of_circular_quantities
    avAz_rad = np.deg2rad(azArr)
    avAz = np.rad2deg(
        np.arctan2(np.sum(np.sin(avAz_rad)), np.sum(np.cos(avAz_rad)))
    )
    avAz = avAz if avAz > 0 else avAz + 360

    return avAlt, avAz


def __get_average_pointing(file, runNumber):
    """
    Return circular mean RA and average DEC of telescope pointing

    """
    pointingDataReduced = file[
        f"run_{runNumber}/stereo/pointingDataReduced"].arrays(library="np")
    avRA = np.rad2deg(circmean(pointingDataReduced["TelRAJ2000"]))
    avDec = np.mean(np.rad2deg(pointingDataReduced["TelDecJ2000"]))

    return avRA, avDec


def __get_pointing_altaz(pointing_ra, pointing_dec, obstime):
    """Return the telescope pointing altitude and azimuth at ``obstime``.

    Reconstructed event directions are affected by the event selection.  The
    telescope pointing stored in ``pointingDataReduced`` is independent of
    those cuts, so it is used for the zenith coordinate of the IRF query.
    The midpoint of a run is the appropriate representative time for the
    run-averaged pointing position.
    """

    location = EarthLocation.from_geodetic(
        lon=VTS_REFERENCE_LON * u.deg,
        lat=VTS_REFERENCE_LAT * u.deg,
        height=VTS_REFERENCE_HEIGHT * u.m,
    )
    pointing = SkyCoord(
        ra=pointing_ra * u.deg,
        dec=pointing_dec * u.deg,
        frame="fk5",
        equinox="J2000",
    )
    altaz = pointing.transform_to(AltAz(obstime=obstime, location=location))

    return altaz.alt.to_value(u.deg), altaz.az.to_value(u.deg)


def __read_quality_flag_from_log(file, runNumber):
    """
    Return quality flag read from evndispLog

    """
    try:
        return getRunQuality(file["run_{}/stereo/evndispLog".format(runNumber)].member("fLines"))
    except KeyError:
        logger.info("Eventdisplay logfile not found in anasum root file. Quality flag set to 0")
    return 0


def __get_ontime(file, runNumber, t_start_from_reference, t_stop_from_reference):
    """
    time on target in seconds, taking into account time masks

    """

    try:
        BitArray = file[f"run_{runNumber}"]["stereo"]["timeMask"]["maskBits"].member("fAllBits")
        gti_tstart_from_reference, gti_tstop_from_reference, ontime_s = getGTI(
            BitArray, t_start_from_reference
        )
    except KeyError:
        try:
            time_mask = file[f"run_{runNumber}"]["stereo"]["timeMask"]
        except KeyError:
            logger.info("Eventdisplay time mask not found; using the full run interval")
        else:
            for key in time_mask.keys():
                logger.info("maskBits not found, available key: %s", key)
        gti_tstart_from_reference = [t_start_from_reference]
        gti_tstop_from_reference = [t_stop_from_reference]
        ontime_s = t_stop_from_reference - t_start_from_reference

    return gti_tstart_from_reference, gti_tstop_from_reference, ontime_s


def __get_time_vector(time_of_day, seconds_from_reference):
    """
    Time vector in seconds since reference time

    This should already have microsecond resolution if stored with
    double precision. Max 24*60*60 seconds

    """

    if time_of_day.max() > 24 * 60 * 60:
        logger.error("Max value in time_of_day array exceeds length of a day")
        raise ValueError
    return seconds_from_reference + time_of_day


def __get_mask(DL3EventTree, select):
    """
    Apply a selection to the event list

    """

    mask = np.ones(len(DL3EventTree["RA"]), bool)
    if select is not None and len(select) > 0:
        logger.info("Applying event selection filter: %s", select)
        for key, value in select.items():
            if isinstance(value, (list, tuple)):
                mask = (
                    mask
                    & (DL3EventTree[key] >= value[0])
                    & (DL3EventTree[key] <= value[1])
                )
            elif isinstance(value, (int, float)):
                mask = mask & (DL3EventTree[key] == value)
            else:
                logger.error("select condition required a list or tuple of ranges")
                raise TypeError
        logger.info("%d of %d events after selection.", np.sum(mask), len(mask))

    return mask
