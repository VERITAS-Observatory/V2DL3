import logging
from copy import deepcopy

from astropy.time import Time
import numpy as np

from pyV2DL3.constant import VTS_REFERENCE_MJD
from pyV2DL3.vegas.util import decodeConfigMask
from pyV2DL3.vegas.util import getGTArray
from pyV2DL3.vegas.util import getTimeCut
from pyV2DL3.vegas.util import mergeTimeCut
from pyV2DL3.vegas.util import produceTelList

logger = logging.getLogger(__name__)

windowSizeForNoise = 7


def __fillEVENTS_not_safe__(vegasFileIO, event_classes=None, save_msw_msl=False):
    # Load header ,array info and selected event tree ( vegas > v2.5.7)
    runHeader = vegasFileIO.loadTheRunHeader()
    selectedEventsTree = vegasFileIO.loadTheCutEventTree()
    qStatsData = vegasFileIO.loadTheQStatsData()
    pixelData = vegasFileIO.loadThePixelStatusData()
    arrayInfo = vegasFileIO.loadTheArrayInfo(0)
    cuts = vegasFileIO.loadTheCutsInfo()

    # Calculate time references
    startTime = runHeader.getStartTime()
    endTime = runHeader.getEndTime()

    mjd = startTime.getMJDInt()
    startTime_s = float(startTime.getDayNS()) / 1e9
    endTime_s = float(endTime.getDayNS()) / 1e9

    t_ref = Time(VTS_REFERENCE_MJD, format="mjd", scale="utc")
    seconds_from_reference_t0 = (Time(mjd, format="mjd", scale="utc") - t_ref).sec
    startTime_s = startTime_s + seconds_from_reference_t0
    endTime_s = endTime_s + seconds_from_reference_t0

    if event_classes is None:
        # We use a single "event class" when not in event class mode
        num_event_classes = 1
    else:
        # Set num_event_classes so we dont call len(event_classes)
        # thousands of times when filling events.
        num_event_classes = len(event_classes)
        if num_event_classes < 1:
            # Dev exception
            raise Exception("event_classes was passed in as an empty List")

    # These arrays are the same for every event class
    avAlt = []
    avAz = []
    avRA = []
    avDec = []

    # These arrays are unique to each event class
    event_arrays = {
        "evNumArr": [],
        "timeArr": [],
        "raArr": [],
        "decArr": [],
        "azArr": [],
        "altArr": [],
        "energyArr": [],
        "nTelArr": [],
    }

    # Arrays exclusive to event class mode
    if event_classes is not None or save_msw_msl:
        event_arrays["mswArr"] = []
        event_arrays["mslArr"] = []
        if event_classes is not None:
            event_arrays["eClassArr"] = []

    # Deep copy the dictionary for each event class
    event_class_dicts = [event_arrays]
    for i in range(num_event_classes - 1):
        event_class_dicts.append(deepcopy(event_arrays))

    logger.debug("Start filling events ...")

    for ev in selectedEventsTree:
        # When not using event class mode, this will stay 0
        event_class_idx = 0

        if event_classes is not None:
            fMSW = ev.S.fMSW
            """Determine which event class (if any) the event falls into.

            Simply loop through the event classes and break if this event meets all
            of an event class' parameters

            For now, we only do it based on the MSW intervals
            """
            event_class_idx = 0
            # Loop through event classes to check if this event satisfies the parameters
            for ec in event_classes:
                if ec.msw_lower <= fMSW < ec.msw_upper:
                    break
                event_class_idx += 1

            # If this event falls into an event classes
            if event_class_idx < num_event_classes:
                event_class_dicts[event_class_idx]["eClassArr"].append(
                    8 * 2**event_class_idx)
                event_class_dicts[event_class_idx]["mswArr"].append(fMSW)
                event_class_dicts[event_class_idx]["mslArr"].append(ev.S.fMSL)
            # Else skip to next event
            else:
                logger.debug("Event excluded: " + str(ev.S.fArrayEventNum)
                             + " MSW: " + str(fMSW))
                continue

        elif save_msw_msl:
            event_class_dicts[event_class_idx]["mswArr"].append(ev.S.fMSW)
            event_class_dicts[event_class_idx]["mslArr"].append(ev.S.fMSL)

        # seconds since first light
        time_relative_to_reference = (
            float(ev.S.fTime.getDayNS()) / 1e9 + seconds_from_reference_t0
        )
        this_event_class = event_class_dicts[event_class_idx]
        this_event_class["evNumArr"].append(ev.S.fArrayEventNum)
        this_event_class["timeArr"].append(time_relative_to_reference)
        this_event_class["raArr"].append(np.rad2deg(ev.S.fDirectionRA_J2000_Rad))
        this_event_class["decArr"].append(np.rad2deg(ev.S.fDirectionDec_J2000_Rad))
        this_event_class["azArr"].append(np.rad2deg(ev.S.fDirectionAzimuth_Rad))
        this_event_class["altArr"].append(np.rad2deg(ev.S.fDirectionElevation_Rad))
        this_event_class["energyArr"].append(ev.S.fEnergy_GeV / 1000.0)
        this_event_class["nTelArr"].append(ev.S.fImages)

        avAlt.append(ev.S.fArrayTrackingElevation_Deg)
        avAz.append(ev.S.fArrayTrackingAzimuth_Deg)
        avRA.append(ev.S.fArrayTrackingRA_J2000_Rad)
        avDec.append(ev.S.fArrayTrackingDec_J2000_Rad)

    avAlt = np.mean(avAlt)
    # Calculate average azimuth angle from average vector on a circle
    # https://en.wikipedia.org/wiki/Mean_of_circular_quantities
    avAz_rad = np.deg2rad(avAz)
    avAz = np.rad2deg(np.arctan2(np.sum(np.sin(avAz_rad)), np.sum(np.cos(avAz_rad))))
    avAz = avAz if avAz > 0 else avAz + 360

    avRA = np.rad2deg(np.mean(avRA))
    avDec = np.rad2deg(np.mean(avDec))

    # Get Time Cuts and build GTI start and stop time array
    # and calculate live time
    for k in cuts:
        tc = getTimeCut(k.fCutsFileText)

    goodTimeStart, goodTimeStop = getGTArray(startTime_s, endTime_s, mergeTimeCut(tc))
    real_live_time = np.sum(np.array(goodTimeStop) - np.array(goodTimeStart))
    
    # Construct an array to hold the event dict(s) to be returned:
    returned_dicts = []
    for i in range(num_event_classes):
        returned_dicts.append({})

    # Filling Event List(s)
    for index in range(num_event_classes):
        # Fill each event class's FITS format from its corresponding event arrays
        evt_dict = returned_dicts[index]
        arr_dict = event_class_dicts[index]
        evt_dict["EVENT_ID"] = arr_dict["evNumArr"]
        evt_dict["TIME"] = arr_dict["timeArr"]
        evt_dict["RA"] = arr_dict["raArr"]
        evt_dict["DEC"] = arr_dict["decArr"]
        evt_dict["ALT"] = arr_dict["altArr"]
        evt_dict["AZ"] = arr_dict["azArr"]
        evt_dict["ENERGY"] = arr_dict["energyArr"]
        evt_dict["EVENT_TYPE"] = arr_dict["nTelArr"]
        if "mswArr" in arr_dict:
            evt_dict["MSW"] = arr_dict["mswArr"]
        if "mslArr" in arr_dict:
            evt_dict["MSL"] = arr_dict["mslArr"]
        if "eclassArr" in arr_dict:
            evt_dict["EVENT_CLASS"] = arr_dict["eClassArr"]
        # Filling Header info
        evt_dict["OBS_ID"] = runHeader.getRunNumber()
        evt_dict["DATE-OBS"] = startTime.getString().split()[0]
        evt_dict["TIME-OBS"] = startTime.getString().split()[1]
        evt_dict["DATE-END"] = endTime.getString().split()[0]
        evt_dict["TIME-END"] = endTime.getString().split()[1]
        evt_dict["TSTART"] = startTime_s
        evt_dict["TSTOP"] = endTime_s
        evt_dict["ONTIME"] = endTime_s - startTime_s
        evt_dict["LIVETIME"] = (
            runHeader.getLiveTimeFrac(True) * real_live_time
        )  # True to suppress error warnings
        evt_dict["DEADC"] = evt_dict["LIVETIME"] / evt_dict["ONTIME"]
        evt_dict["OBJECT"] = runHeader.getSourceId()
        evt_dict["RA_PNT"] = avRA
        evt_dict["DEC_PNT"] = avDec
        evt_dict["ALT_PNT"] = avAlt
        evt_dict["AZ_PNT"] = avAz
        evt_dict["RA_OBJ"] = np.rad2deg(runHeader.getSourceRA())
        evt_dict["DEC_OBJ"] = np.rad2deg(runHeader.getSourceDec())
        evt_dict["TELLIST"] = produceTelList(runHeader.fRunInfo.fConfigMask)
        evt_dict["N_TELS"] = runHeader.pfRunDetails.fTels

    avNoise = 0
    nTels = 0
    for telID in decodeConfigMask(runHeader.fRunInfo.fConfigMask):
        avNoise += qStatsData.getCameraAverageTraceVarTimeIndpt(
            telID - 1, windowSizeForNoise, pixelData, arrayInfo
        )
        nTels += 1

    avNoise /= nTels
    return (
        {
            "goodTimeStart": goodTimeStart,
            "goodTimeStop": goodTimeStop,
            "TSTART": startTime_s,
            "TSTOP": endTime_s,
        },
        {"azimuth": avAz, "zenith": (90.0 - avAlt), "noise": avNoise},
        returned_dicts,
    )
