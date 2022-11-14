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


def __fillEVENTS_not_safe__(vegasFileIO, effective_area_files,
                            fov_cut=True, event_class_mode=False, reco_type=1, save_msw_msl=False):
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
    seconds_from_reference_t0 = (
        Time(mjd, format="mjd", scale="utc") - t_ref).sec
    startTime_s = startTime_s + seconds_from_reference_t0
    endTime_s = endTime_s + seconds_from_reference_t0

    time_avg = (Time(startTime.getString(), format="iso", scale="utc")) + \
        (Time(endTime.getString(), format="iso", scale="utc") - Time(startTime.getString(), format="iso", scale="utc")
         ) / 2.

    # Set num_event_groups so we dont call len(effective_area_files)
    # thousands of times when filling events.
    num_event_groups = len(effective_area_files)
    if num_event_groups < 1:
        # Dev exception
        raise Exception("effective_area_files was passed in as an empty List")

    # Average arrays are independent of event group
    avAlt = []
    avAz = []
    avRA = []
    avDec = []

    # These arrays are unique to each event group
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
    if event_class_mode or save_msw_msl:
        event_arrays["mswArr"] = []
        event_arrays["mslArr"] = []

    # Deep copy the dictionary for each event class
    event_groups = [event_arrays]
    for i in range(num_event_groups - 1):
        event_groups.append(deepcopy(event_arrays))

    logger.debug("Start filling events ...")

    for ev in selectedEventsTree:
        # Event reconstruction method
        if reco_type == 1:
            reco = ev.S
        elif reco_type == 2:
            reco = ev.M3D
        else:
            raise Exception("Invalid reconstruction type!"
                            "\nSee --help for supported arguments")

        # Include event in averages regardless of whether it gets cut
        avAlt.append(reco.fArrayTrackingElevation_Deg)
        avAz.append(reco.fArrayTrackingAzimuth_Deg)
        avRA.append(reco.fArrayTrackingRA_J2000_Rad)
        avDec.append(reco.fArrayTrackingDec_J2000_Rad)

        event_class_idx = 0
        if event_class_mode:
            fMSW = reco.fMSW
            """Determine which event class (if any) the event falls into.

            Simply loop through the effective area files and break if this event meets all
            of one's cuts parameters

            For now, we only do it based on the MSW intervals
            """
            for ec in effective_area_files:
                if ec.msw_lower <= fMSW < ec.msw_upper:
                    break
                event_class_idx += 1

            # If this event falls outside of any event class
            if event_class_idx >= num_event_groups:
                # Skip to next event
                logger.debug("Event excluded: " + str(reco.fArrayEventNum)
                             + " MSW: " + str(fMSW) + " not within an event class' MSW range")
                continue

        # Check FoV if appropriate
        if fov_cut:
            fov_cut_upper = effective_area_files[event_class_idx].fov_cut_upper
            fov_cut_lower = effective_area_files[event_class_idx].fov_cut_lower
            if fov_cut_lower > 0 or fov_cut_upper <= 180:
                excluded, tel_sep = check_FoV_exclusion(
                    reco, fov_cut_upper, fov_cut_lower)
                if excluded:
                    logger.debug("Event excluded: " + str(reco.fArrayEventNum)
                                 + " separation: " +
                                 str(tel_sep) + " not within FoVCut range: "
                                 + str(fov_cut_lower) + "-" + str(fov_cut_upper))
                    continue

        if event_class_mode or save_msw_msl:
            event_groups[event_class_idx]["mswArr"].append(reco.fMSW)
            event_groups[event_class_idx]["mslArr"].append(reco.fMSL)

        # seconds since first light
        time_relative_to_reference = (
            float(reco.fTime.getDayNS()) / 1e9 + seconds_from_reference_t0
        )
        this_event_group = event_groups[event_class_idx]
        this_event_group["evNumArr"].append(reco.fArrayEventNum)
        this_event_group["timeArr"].append(time_relative_to_reference)
        this_event_group["raArr"].append(
            np.rad2deg(reco.fDirectionRA_J2000_Rad))
        this_event_group["decArr"].append(
            np.rad2deg(reco.fDirectionDec_J2000_Rad))
        this_event_group["azArr"].append(
            np.rad2deg(reco.fDirectionAzimuth_Rad))
        this_event_group["altArr"].append(
            np.rad2deg(reco.fDirectionElevation_Rad))
        this_event_group["energyArr"].append(reco.fEnergy_GeV / 1000.0)
        this_event_group["nTelArr"].append(reco.fImages)

    avAlt = np.mean(avAlt)
    # Calculate average azimuth angle from average vector on a circle
    # https://en.wikipedia.org/wiki/Mean_of_circular_quantities
    avAz_rad = np.deg2rad(avAz)
    avAz = np.rad2deg(np.arctan2(
        np.sum(np.sin(avAz_rad)), np.sum(np.cos(avAz_rad))))
    avAz = avAz if avAz > 0 else avAz + 360

    avRA = np.rad2deg(np.mean(avRA))
    avDec = np.rad2deg(np.mean(avDec))

    # Get Time Cuts and build GTI start and stop time array
    # and calculate live time
    for k in cuts:
        tc = getTimeCut(k.fCutsFileText)

    goodTimeStart, goodTimeStop = getGTArray(
        startTime_s, endTime_s, mergeTimeCut(tc))
    real_live_time = np.sum(np.array(goodTimeStop) - np.array(goodTimeStart))

    # Construct an array to hold the event dict(s) to be returned:
    returned_dicts = []
    for i in range(num_event_groups):
        returned_dicts.append({})

    # Filling Event List(s)
    for index in range(num_event_groups):
        # Fill each event class's FITS format from its corresponding event arrays
        evt_dict = returned_dicts[index]
        arr_dict = event_groups[index]
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
        # Filling Header info
        evt_dict["OBS_ID"] = runHeader.getRunNumber()
        evt_dict["DATE-OBS"] = Time(startTime.getString(),
                                    format="iso", scale="utc").to_value("fits")
        evt_dict["TIME-OBS"] = startTime.getString().split()[1]
        evt_dict["DATE-END"] = Time(endTime.getString(),
                                    format="iso", scale="utc").to_value("fits")
        evt_dict["TIME-END"] = endTime.getString().split()[1]
        evt_dict["DATE-AVG"] = time_avg.to_value("fits")
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


"""
Check event against the provided FoV upper and lower bounds

Arguments:
    event_skycoord  --  Reconstructed shower direction as an astropy.SkyCoord
    reco            --  The event reconstruction (e.g ev.S or ev.ITM)
    fov_cut_upper   --  The FoV upper limit in degrees

Returns:
    Bool  --  True if event falls outside of the FoV (event needs excluded)
    float --  The computed separation from the FoV center
"""


def check_FoV_exclusion(reco, fov_cut_upper, fov_cut_lower):
    event_skycoord = SkyCoord(
        np.rad2deg(reco.fDirectionRA_J2000_Rad),
        np.rad2deg(reco.fDirectionDec_J2000_Rad),
        frame='icrs', unit=(units.deg, units.deg)
    )

    pointing_position = SkyCoord(
        np.rad2deg(reco.fArrayTrackingRA_J2000_Rad),
        np.rad2deg(reco.fArrayTrackingDec_J2000_Rad),
        frame='icrs', unit=(units.deg, units.deg)
    )

    tel_sep = pointing_position.separation(event_skycoord).degree

    if tel_sep > fov_cut_upper or tel_sep < fov_cut_lower:
        return True, tel_sep

    return False, tel_sep
