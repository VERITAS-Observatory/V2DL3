import logging

import numpy as np
import uproot
from astropy.time import Time

from pyV2DL3.constant import (
    VTS_REFERENCE_HEIGHT,
    VTS_REFERENCE_LAT,
    VTS_REFERENCE_LON,
    VTS_REFERENCE_MJD,
)
from pyV2DL3.eventdisplay.util import (
    ZeroLengthEventList,
    getGTI,
    getRunQuality,
    produce_tel_list,
)

logger = logging.getLogger(__name__)


def __fillEVENTS__(edFileIO, select=None):
    """
    Fill event list and event header from anasum file

    """

    evt_dict = {}

    with uproot.open(edFileIO) as file:
        runSummary = file["total_1/stereo/tRunSummary"].arrays(library="np")
        runNumber = runSummary["runOn"][0]
        logging.info("Run number: %d", runNumber)
        telConfig = file[f"run_{runNumber}/stereo/telconfig"].arrays(library="np")

        t_start, t_stop, t_avg = __get_start_stop_times(file)
        t_start_from_reference, t_stop_from_reference, seconds_from_reference = \
            __get_times_since_reference_time(t_start, t_stop)

        DL3EventTree = file["run_{}/stereo/DL3EventTree".format(runNumber)].arrays(
            library="np"
        )
        if len(DL3EventTree["eventNumber"]) == 0:
            logging.error("Empty event list")
            raise ZeroLengthEventList

        mask = __get_mask(DL3EventTree, select)

        try:
            # Test if anasum file was created using the all events option.
            # In this case write out the additional output.
            IsGamma = DL3EventTree["IsGamma"][mask]
            bdtScore = DL3EventTree["MVA"][mask]
            all_events = True
        except KeyError:
            all_events = False

        # Event List
        evt_dict["EVENT_ID"] = DL3EventTree["eventNumber"]
        evt_dict["TIME"] = __get_time_vector(DL3EventTree, seconds_from_reference)
        evt_dict["RA"] = DL3EventTree["RA"][mask]
        evt_dict["DEC"] = DL3EventTree["DEC"][mask]
        evt_dict["ALT"] = DL3EventTree["El"][mask]
        evt_dict["AZ"] = DL3EventTree["Az"][mask]
        evt_dict["ENERGY"] = DL3EventTree["Energy"][mask]
        evt_dict["EVENT_TYPE"] = DL3EventTree["NImages"][mask]
        if all_events:
            evt_dict["GAMMANESS"] = bdtScore
            evt_dict["IS_GAMMA"] = IsGamma

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
        evt_dict["RA_PNT"], evt_dict["DEC_PNT"] = __get_average_pointing(file, runNumber)
        evt_dict["ALT_PNT"], evt_dict["AZ_PNT"] = __get_average_event_direction(
            evt_dict["ALT"], evt_dict["AZ"])
        evt_dict["RA_OBJ"] = runSummary["TargetRAJ2000"][0]
        evt_dict["DEC_OBJ"] = runSummary["TargetDecJ2000"][0]
        evt_dict["TELLIST"] = produce_tel_list(telConfig)
        MaxImgSel = np.max(DL3EventTree["ImgSel"])
        evt_dict["N_TELS"] = np.binary_repr(MaxImgSel).count("1")
        logging.info("Number of Telescopes: {}".format(evt_dict["N_TELS"]))
        evt_dict["GEOLON"] = VTS_REFERENCE_LON
        evt_dict["GEOLAT"] = VTS_REFERENCE_LAT
        evt_dict["ALTITUDE"] = VTS_REFERENCE_HEIGHT
        evt_dict["Xoff"] = DL3EventTree["Xoff"][mask]
        evt_dict["Yoff"] = DL3EventTree["Yoff"][mask]
        evt_dict["NSBLEVEL"] = np.mean(DL3EventTree["MeanPedvar"][mask])

        evt_dict["QUALITY"] = __read_quality_flag_from_log(file, runNumber)
        gti_tstart_from_reference, gti_tstop_from_reference, evt_dict["ONTIME"] = \
            __get_ontime(file, runNumber, t_start_from_reference, t_stop_from_reference)
        evt_dict["LIVETIME"] = evt_dict["ONTIME"] * (1 - evt_dict["DEADC"])

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


def __get_start_stop_times(file):
    """
    Return run start and stop time read from tRunSummary tree

    """

    start_mjd = file["total_1/stereo/tRunSummary/MJDrunstart"].array(library="np")[0]
    stop_mjd = file["total_1/stereo/tRunSummary/MJDrunstop"].array(library="np")[0]

    # convert mjd to fits format
    t_start = Time(start_mjd, format="mjd", scale="utc")
    t_stop = Time(stop_mjd, format="mjd", scale="utc")
    t_avg = Time(start_mjd, format="mjd", scale="utc") + \
        (
            Time(stop_mjd, format="mjd", scale="utc") - Time(start_mjd, format="mjd", scale="utc")
        ) / 2.

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
    avRA = np.rad2deg(
        np.arctan2(
            np.sum(np.sin(pointingDataReduced["TelRAJ2000"])),
            np.sum(np.cos(pointingDataReduced["TelRAJ2000"])),
        )
    )
    avDec = np.mean(np.rad2deg(pointingDataReduced["TelDecJ2000"]))

    return avRA, avDec


def __read_quality_flag_from_log(file, runNumber):
    """
    Return quality flag read from evndispLog

    """
    try:
        return getRunQuality(file["run_{}/stereo/evndispLog".format(runNumber)].member(
                "fLines"
            )
        )
    except KeyError:
        logging.error("Eventdisplay logfile not found in anasum root file")
        logging.error("Please make sure to use Eventdisplay >= 486")
        raise


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
        for k in file["run_{}".format(runNumber)]["stereo"]["timeMask"].keys():
            logging.info("maskBits not found, Available keys: {0}".format(k))
        gti_tstart_from_reference = [t_start_from_reference]
        gti_tstop_from_reference = [t_stop_from_reference]
        ontime_s = t_stop_from_reference - t_start_from_reference

    return gti_tstart_from_reference, gti_tstop_from_reference, ontime_s


def __get_time_vector(DL3EventTree, seconds_from_reference):
    """
    Time vector in seconds since reference time

    This should already have microsecond resolution if stored with
    double precision. Max 24*60*60 seconds

    """

    time_of_day = DL3EventTree["timeOfDay"]
    if time_of_day.max() > 24 * 60 * 60:
        logging.error(
            "Max value in time_of_day  \
                            array exceeds length of a day"
        )
        raise ValueError
    return seconds_from_reference + time_of_day


def __get_mask(DL3EventTree, select):
    """

    """

    mask = np.ones(len(DL3EventTree["RA"]), bool)
    if select is not None:
        logging.info(select)
        for key, value in select.items():
            if isinstance(value, (list, tuple)):
                mask = (
                    mask
                    & (DL3EventTree[key] >= value[0])
                    & (DL3EventTree[key] <= value[1])
                )
            else:
                logging.error(
                    "select condition required a list or tuple of ranges"
                )
                raise TypeError
        logging.info(f"{np.sum(mask)} of {len(mask)} events after selection.")

    return mask
