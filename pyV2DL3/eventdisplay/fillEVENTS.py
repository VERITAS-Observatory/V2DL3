import logging

from astropy.time import Time
import numpy as np
import uproot

from pyV2DL3.constant import VTS_REFERENCE_HEIGHT
from pyV2DL3.constant import VTS_REFERENCE_LAT
from pyV2DL3.constant import VTS_REFERENCE_LON
from pyV2DL3.constant import VTS_REFERENCE_MJD

from pyV2DL3.eventdisplay.util import getGTI
from pyV2DL3.eventdisplay.util import getRunQuality
from pyV2DL3.eventdisplay.util import produce_tel_list

logger = logging.getLogger(__name__)


def __fillEVENTS__(edFileIO, select=None):
    if select is None:
        select = {}

    evt_dict = {}

    # reading variables with uproot
    with uproot.open(edFileIO) as file:
        runSummary = file["total_1/stereo/tRunSummary"].arrays(library="np")
        runNumber = runSummary["runOn"][0]
        logging.info("Run number: {}".format(runNumber))
        telConfig = file["run_{}/stereo/telconfig".format(runNumber)].arrays(
            library="np"
        )

        # Get start and stop time within the run.
        start_mjd = file["total_1/stereo/tRunSummary/MJDrunstart"].array(library="np")[
            0
        ]
        stop_mjd = file["total_1/stereo/tRunSummary/MJDrunstop"].array(library="np")[0]

        # convert mjd to fits format
        t_start_fits = Time(start_mjd, format="mjd", scale="utc").to_value("fits")
        t_stop_fits = Time(stop_mjd, format="mjd", scale="utc").to_value("fits")
        deadtime = file["total_1/stereo/tRunSummary/DeadTimeFracOn"].array(
            library="np"
        )[0]

        # Number of seconds between reference time and run MJD at 00:00:00:
        t_ref = Time(VTS_REFERENCE_MJD, format="mjd", scale="utc")
        seconds_from_reference = (
            Time(np.trunc(start_mjd), format="mjd", scale="utc") - t_ref
        ).sec

        tstart_from_reference = (Time(start_mjd, format="mjd", scale="utc") - t_ref).sec
        tstop_from_reference = (Time(stop_mjd, format="mjd", scale="utc") - t_ref).sec

        DL3EventTree = file["run_{}/stereo/DL3EventTree".format(runNumber)].arrays(
            library="np"
        )
        evNumArr = DL3EventTree["eventNumber"]

        # This should already have microsecond resolution if stored with
        # double precision. Max 24*60*60 seconds
        time_of_day = DL3EventTree["timeOfDay"]
        if time_of_day.max() > 24 * 60 * 60:
            raise ValueError(
                "Max value in time_of_day  \
                              array exceeds length of a day"
            )
        # mjd time in full days since reference time + event time
        timeArr = seconds_from_reference + time_of_day

        mask = np.ones(len(DL3EventTree["RA"]), bool)
        if select:
            logging.info(select)
            for key, value in select.items():
                if isinstance(value, (list, tuple)):
                    mask = (
                        mask
                        & (DL3EventTree[key] >= value[0])
                        & (DL3EventTree[key] <= value[1])
                    )
                else:
                    raise TypeError(
                        "select condition required \
                                     a list or tuple of ranges"
                    )
            logging.info(f"{np.sum(mask)} of {len(mask)} events after selection.")

        raArr = DL3EventTree["RA"][mask]
        decArr = DL3EventTree["DEC"][mask]
        azArr = DL3EventTree["Az"][mask]
        altArr = DL3EventTree["El"][mask]
        energyArr = DL3EventTree["Energy"][mask]
        nTelArr = DL3EventTree["NImages"][mask]
        try:
            # Test if anasum file was created using the all events option.
            # In this case write out the additional output.
            IsGamma = DL3EventTree["IsGamma"][mask]
            bdtScore = DL3EventTree["MVA"][mask]
            all_events = True
        except KeyError:
            all_events = False

        avAlt = np.mean(altArr)
        # Calculate average azimuth angle from average vector on a circle
        # https://en.wikipedia.org/wiki/Mean_of_circular_quantities
        avAz_rad = np.deg2rad(azArr)
        avAz = np.rad2deg(
            np.arctan2(np.sum(np.sin(avAz_rad)), np.sum(np.cos(avAz_rad)))
        )
        avAz = avAz if avAz > 0 else avAz + 360

        # Calculate circular mean RA and average DEC
        pointingDataReduced = file[
            "run_{}/stereo/pointingDataReduced".format(runNumber)
        ].arrays(library="np")
        avRA = np.rad2deg(
            np.arctan2(
                np.sum(np.sin(pointingDataReduced["TelRAJ2000"])),
                np.sum(np.cos(pointingDataReduced["TelRAJ2000"])),
            )
        )
        avDec = np.mean(np.rad2deg(pointingDataReduced["TelDecJ2000"]))

        # Filling Event List
        evt_dict["EVENT_ID"] = evNumArr
        evt_dict["TIME"] = timeArr
        evt_dict["RA"] = raArr
        evt_dict["DEC"] = decArr
        evt_dict["ALT"] = altArr
        evt_dict["AZ"] = azArr
        evt_dict["ENERGY"] = energyArr
        evt_dict["EVENT_TYPE"] = nTelArr
        if all_events:
            evt_dict["BDT_SCORE"] = bdtScore
            evt_dict["IS_GAMMA"] = IsGamma

        # Filling Header info
        evt_dict["OBS_ID"] = runNumber
        evt_dict["DATE-OBS"] = t_start_fits
        evt_dict["DATE-END"] = t_stop_fits
        evt_dict["TSTART"] = tstart_from_reference
        evt_dict["TSTOP"] = tstop_from_reference
        evt_dict["MJDREFI"] = int(VTS_REFERENCE_MJD)
        evt_dict["DEADC"] = 1 - deadtime
        evt_dict["OBJECT"] = runSummary["TargetName"][0]
        evt_dict["RA_PNT"] = avRA
        evt_dict["DEC_PNT"] = avDec
        evt_dict["ALT_PNT"] = avAlt
        evt_dict["AZ_PNT"] = avAz
        evt_dict["RA_OBJ"] = runSummary["TargetRAJ2000"][0]
        evt_dict["DEC_OBJ"] = runSummary["TargetDecJ2000"][0]
        evt_dict["TELLIST"] = produce_tel_list(telConfig)
        MaxImgSel = np.max(DL3EventTree["ImgSel"])
        evt_dict["N_TELS"] = np.binary_repr(MaxImgSel).count("1")
        logging.info("Number of Telescopes: {}".format(evt_dict["N_TELS"]))
        evt_dict["GEOLON"] = VTS_REFERENCE_LON
        evt_dict["GEOLAT"] = VTS_REFERENCE_LAT
        evt_dict["ALTITUDE"] = VTS_REFERENCE_HEIGHT

        # Read evndispLog which is stored as TMacro in anasum root file (ED >= 486)
        try:
            evndisplog_data = file["run_{}/stereo/evndispLog".format(runNumber)].member(
                "fLines"
            )
            evt_dict["QUALITY"] = getRunQuality(evndisplog_data)
        except KeyError:
            logging.exception("Eventdisplay logfile not found in anasum root file")
            logging.exception("Please make sure to use ED >= 486")

        avPedvar = runSummary["pedvarsOn"][0]

        try:
            BitArray = file["run_{}".format(runNumber)]["stereo"]["timeMask"][
                "maskBits"
            ].member("fAllBits")
            gti_tstart_from_reference, gti_tstop_from_reference, ontime_s = getGTI(
                BitArray, tstart_from_reference
            )
        except KeyError:
            for k in file["run_{}".format(runNumber)]["stereo"]["timeMask"].keys():
                logging.info("maskBits not found, Available keys: {0}".format(k))
            gti_tstart_from_reference = [tstart_from_reference]
            gti_tstop_from_reference = [tstop_from_reference]
            ontime_s = tstop_from_reference - tstart_from_reference

        evt_dict["ONTIME"] = ontime_s
        evt_dict["LIVETIME"] = ontime_s * (1 - deadtime)

    return (
        {
            "goodTimeStart": gti_tstart_from_reference,
            "goodTimeStop": gti_tstop_from_reference,
            "TSTART": tstart_from_reference,
            "TSTOP": tstop_from_reference,
        },
        {"azimuth": avAz, "zenith": (90.0 - avAlt), "pedvar": avPedvar},
        evt_dict,
    )
