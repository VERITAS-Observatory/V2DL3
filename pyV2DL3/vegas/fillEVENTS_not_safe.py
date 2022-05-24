import logging
from copy import deepcopy

from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as units
import numpy as np

from pyV2DL3.constant import VTS_REFERENCE_MJD
from pyV2DL3.vegas.util import decodeConfigMask
from pyV2DL3.vegas.util import getGTArray
from pyV2DL3.vegas.util import getTimeCut
from pyV2DL3.vegas.util import mergeTimeCut
from pyV2DL3.vegas.util import produceTelList

logger = logging.getLogger(__name__)

windowSizeForNoise = 7


def __fillEVENTS_not_safe__(vegasFileIO, reco_type=1, user_cuts_dict=None, event_classes=None, store_msw_msl=False):
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

    # Default value when event_classes is None
    num_event_classes = 1
    if event_classes is not None:
        # Set num_event_classes so we dont call len(event_classes)
        # thousands of times when filling events.
        num_event_classes = len(event_classes)
        if num_event_classes < 0:
            # Dev exception
            raise Exception("event_classes was passed in as an empty List")

    # Arrays common to all event classes
    avAlt = []
    avAz = []
    avRA = []
    avDec = []

    arrays_dict = {
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
    if event_classes is not None:
        arrays_dict["eclassArr"] = []
        arrays_dict["mswArr"] = []
        arrays_dict["mslArr"] = []

    # Load an array with deep copies of our dictionary for each event class
    array_dicts = [arrays_dict]
    for i in range(num_event_classes - 1):
        array_dicts.append(deepcopy(arrays_dict))

    spatial_cuts = False
    # Load user cuts if provided
    if user_cuts_dict is not None:
        if "spatial_exclusion" in user_cuts_dict:
            (srcra, srcdec, srcsep) = user_cuts_dict["spatial_exclusion"]
            spatial_cuts = True

    logger.debug("Start filling events ...")

    for ev in selectedEventsTree:
        # For a single event class (or when not using EventClass), this will stay 0
        # so that all qualifying events will be added to our single dict
        event_class_idx = 0

        if reco_type == 1:
            reco = ev.S
        elif reco_type == 2:
            reco = ev.M3D
        else:
            raise Exception("Invalid reconstruction type!"
                            "\nSee --help for supported arguments")

        if spatial_cuts:
            for xra, xdec, xsep in zip(srcra, srcdec, srcsep):
                # Source position
                c1 = SkyCoord(xra, xdec, frame='icrs',
                              unit=(units.deg, units.deg))
                # Actual reconstructed shower direction
                c2 = SkyCoord(np.rad2deg(reco.fDirectionRA_J2000_Rad), np.rad2deg(
                    reco.fDirectionDec_J2000_Rad), frame='icrs', unit=(units.deg, units.deg))
                sep = c1.separation(c2)
                sep = sep.degree
                xsep = float(xsep)
                # If one event falls in an exclusion region, skip the other exclusion regions
                if (sep < xsep):
                    break
            # Remove event
            if (sep < xsep):
                continue

        if event_classes is not None:
            fMSW = reco.fMSW
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

            # If this event falls into one of our event classes
            if event_class_idx < num_event_classes:
                array_dicts[event_class_idx]["eclassArr"].append(
                    8 * 2**event_class_idx)
                array_dicts[event_class_idx]["mswArr"].append(fMSW)
                array_dicts[event_class_idx]["mslArr"].append(reco.fMSL)
            else:
                logger.debug("Event excluded: " + str(reco.fArrayEventNum)
                             + "MSW: " + str(fMSW))
                # Skip to next event
                continue

        elif store_msw_msl:
            array_dicts[event_class_idx]["mswArr"].append(reco.fMSW)
            array_dicts[event_class_idx]["mslArr"].append(reco.fMSL)

        # Seconds since first light
        time_relative_to_reference = (
            float(reco.fTime.getDayNS()) / 1e9 + seconds_from_reference_t0
        )
        arr_dict = array_dicts[event_class_idx]
        arr_dict["evNumArr"].append(reco.fArrayEventNum)
        arr_dict["timeArr"].append(time_relative_to_reference)
        arr_dict["raArr"].append(np.rad2deg(reco.fDirectionRA_J2000_Rad))
        arr_dict["decArr"].append(np.rad2deg(reco.fDirectionDec_J2000_Rad))
        arr_dict["azArr"].append(np.rad2deg(reco.fDirectionAzimuth_Rad))
        arr_dict["altArr"].append(np.rad2deg(reco.fDirectionElevation_Rad))
        arr_dict["energyArr"].append(reco.fEnergy_GeV / 1000.0)
        arr_dict["nTelArr"].append(reco.fImages)
        avAlt.append(reco.fArrayTrackingElevation_Deg)
        avAz.append(reco.fArrayTrackingAzimuth_Deg)
        avRA.append(reco.fArrayTrackingRA_J2000_Rad)
        avDec.append(reco.fArrayTrackingDec_J2000_Rad)

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
    event_dicts = []
    for i in range(num_event_classes):
        event_dicts.append({})

    # Filling Event List
    for i in range(num_event_classes):
        # Fill each 'event_dict' from its corresponding 'array_dict'
        evt_dict = event_dicts[i]
        arr_dict = array_dicts[i]
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
            evt_dict["EVENT_CLASS"] = arr_dict["eclassArr"]
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
        event_dicts,
    )
