import ctypes
import logging
from copy import deepcopy

import astropy.units as units
import numpy as np
import ROOT
from astropy.coordinates import SkyCoord
from astropy.time import Time

from pyV2DL3.constant import VTS_REFERENCE_MJD
from pyV2DL3.vegas.util import (
    decodeConfigMask,
    getGTArray,
    getTimeCut,
    mergeTimeCut,
    produceTelList,
)

logger = logging.getLogger(__name__)

windowSizeForNoise = 7


def __fillEVENTS_not_safe__(
    vegasFileIO,
    effective_area_files,
    irf_to_store,
    fov_cut=True,
    event_class_mode=False,
    reco_type=1,
    save_msw_msl=False,
    corr_EB=False,
    psf_king_params=None,
    st6_configs=None,
):
    if corr_EB and event_class_mode:
        raise Exception("Currently Energy Bias and multiple EAs not supported")
    # Load header ,array info and selected event tree ( vegas > v2.5.7)
    runHeader = vegasFileIO.loadTheRunHeader()
    selectedEventsTree = vegasFileIO.loadTheCutEventTree()
    cuts = vegasFileIO.loadTheCutsInfo()
    arrayInfo = vegasFileIO.loadTheArrayInfo()
    qStatsData = vegasFileIO.loadTheQStatsData()
    pixelData = vegasFileIO.loadThePixelStatusData()

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

    time_avg = (Time(startTime.getString(), format="iso", scale="utc")) + (
        Time(endTime.getString(), format="iso", scale="utc")
        - Time(startTime.getString(), format="iso", scale="utc")
    ) / 2.0

    # Threshold for total pixels suppressed across all telescopes to cause warning.
    n_suppresed_pixel_thresh = 200
    # Threshold of how many standard deviations from the mean to consider a noise value
    # in a run with more than n_suppresed_pixel_thresh suppressed pixels to be considered
    # artificially low/high and replaced with mean for that run.
    n_noise_stddev_thresh = 3

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
        "fNoise": [],
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
            raise Exception(
                "Invalid reconstruction type!" "\nSee --help for supported arguments"
            )

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
                logger.debug(
                    "Event excluded: "
                    + str(reco.fArrayEventNum)
                    + " MSW: "
                    + str(fMSW)
                    + " not within an event class' MSW range"
                )
                continue

        # Check FoV if appropriate
        if fov_cut:
            fov_cut_upper = effective_area_files[event_class_idx].fov_cut_upper
            fov_cut_lower = effective_area_files[event_class_idx].fov_cut_lower
            if fov_cut_lower > 0 or fov_cut_upper <= 180:
                excluded, tel_sep = check_FoV_exclusion(
                    reco, fov_cut_upper, fov_cut_lower
                )
                if excluded:
                    logger.debug(
                        "Event excluded: "
                        + str(reco.fArrayEventNum)
                        + " separation: "
                        + str(tel_sep)
                        + " not within FoVCut range: "
                        + str(fov_cut_lower)
                        + "-"
                        + str(fov_cut_upper)
                    )
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
        this_event_group["raArr"].append(np.rad2deg(reco.fDirectionRA_J2000_Rad))
        this_event_group["decArr"].append(np.rad2deg(reco.fDirectionDec_J2000_Rad))
        this_event_group["azArr"].append(np.rad2deg(reco.fDirectionAzimuth_Rad))
        this_event_group["altArr"].append(np.rad2deg(reco.fDirectionElevation_Rad))
        this_event_group["energyArr"].append(reco.fEnergy_GeV / 1000.0)
        this_event_group["nTelArr"].append(reco.fImages)

        avNoise = 0
        nTels = 0
        for telID in decodeConfigMask(runHeader.fRunInfo.fConfigMask):
            Noise = qStatsData.getCameraAverageTraceVar(
                telID - 1, windowSizeForNoise, reco.fTime, pixelData, arrayInfo
            )
            if Noise > 0:
                # Using mean not median to match VEGAS behaviour
                avNoise += Noise
                nTels += 1
        if nTels == 0:
            avNoise = -100
            logger.warning(
                f"Warning! {runHeader.getRunNumber()} has no (time dependent) noise found for this event. "
                "Will use average of other events in this run."
            )
        else:
            if 0 < nTels < 4:
                logger.warning(
                    f"Warning {runHeader.getRunNumber()} is missing (time dependent) "
                    f"noise for {4 - nTels} "
                    f"telescopes for this event. Using average of other {nTels} telescopes."
                )
            avNoise /= nTels
        this_event_group["fNoise"].append(avNoise)

    # Since avNoise set to -100 if no telescopes, need to replace the -100 values with the average for the run.
    nNonNegativeNoises = 0
    avNonNegativeNoises = 0
    for fNoise in this_event_group["fNoise"]:
        if fNoise > 0:
            nNonNegativeNoises += 1
            avNonNegativeNoises += fNoise

    if nNonNegativeNoises == 0:
        raise ValueError(f"No valid noises found for run {runHeader.getRunNumber()}.")

    else:
        avNonNegativeNoises /= nNonNegativeNoises
        for i, fNoise in enumerate(this_event_group["fNoise"]):
            if fNoise <= 0:  # replace negative noises with average
                this_event_group["fNoise"][i] = avNonNegativeNoises

    # If the high voltage was turned off for part of a time slice, the noise in that timeslice will be artificially low.
    # It should therefore be excluded via cutting the entire time slice in Stage 5.
    # This checks for the number of pixels excluded.
    ntels = runHeader.pfRunDetails.fTels
    isSuppressed = ctypes.c_bool(False)
    n_suppressed_all_tels = 0
    for tel_i in range(ntels):
        n_PMT = 0
        n_suppressed = 0
        tel_info = arrayInfo.telescope(tel_i)
        for chan_i in range(499):
            if tel_info.channel(chan_i).hasPMT():
                n_PMT += 1
                pixelData.getSuppressedTimeIndpt(tel_i, chan_i, isSuppressed)
                if isSuppressed.value:
                    n_suppressed += 1
                    isSuppressed.value = False
        n_suppressed_all_tels += n_suppressed
    if n_suppressed_all_tels > n_suppresed_pixel_thresh:
        MeanPerEventNoise = np.mean(list(this_event_group["fNoise"]))
        StDevPerEventNoise = np.std(list(this_event_group["fNoise"]))
        unique_noises = []  # Use this rather than set() so they are in order
        for x in this_event_group["fNoise"]:
            if x not in unique_noises:
                unique_noises.append(x)

        logger.warning(
            f"Warning! {n_suppressed_all_tels} Pixels Suppressed for Run {runHeader.getRunNumber()}: \n"
            "    This will make noise artificially low in that timeslice. \n"
            f"    Replacing any Time Dependent Noise values more than {n_noise_stddev_thresh} sigma "
            f"below the run mean with the mean of those within {n_noise_stddev_thresh} sigma \n"
            "Alternatively, consider cutting time slice in Stage 5. \n"
            f"All timeslice noises in run: \n {[format(x, '.4f') for x in unique_noises]}"
        )
        # Replace any noise values that are more than n_noise_stddev_thresh * sigma
        # below the mean with the average for the run.
        # Only looking for low noises as these are the result of suppressed pixels.
        nWithinXSigmaOfMeanNoises = 0
        avWithinXSigmaOfMeanNoises = 0
        ValuesReplaced = []
        for fNoise in this_event_group["fNoise"]:
            if fNoise > MeanPerEventNoise - n_noise_stddev_thresh * StDevPerEventNoise:
                nWithinXSigmaOfMeanNoises += 1
                avWithinXSigmaOfMeanNoises += fNoise

        if nWithinXSigmaOfMeanNoises == 0:
            logger.error(
                "Error! No valid noises found for this run. Setting TimeDependentNoise: -100. Do not use."
            )
            avNonNegativeNonSuppressedNoises = -100
        else:
            avWithinXSigmaOfMeanNoises /= nWithinXSigmaOfMeanNoises
            for i, fNoise in enumerate(this_event_group["fNoise"]):
                if (
                    fNoise
                    < MeanPerEventNoise - n_noise_stddev_thresh * StDevPerEventNoise
                ):
                    this_event_group["fNoise"][i] = avWithinXSigmaOfMeanNoises
                    ValuesReplaced.append(fNoise)
            avNonNegativeNonSuppressedNoises = avWithinXSigmaOfMeanNoises
        if len(ValuesReplaced) > 0:
            logger.warning(
                f"Warning! The following Time Dependent Noise values were more than {n_noise_stddev_thresh} "
                f"sigma below the mean (as well as having suppressed pixels): {[f'{n:.4f}' for n in set(ValuesReplaced)]} \n"
                f"These values have been replaced with the average of the other noise values within {n_noise_stddev_thresh} "
                f"sigma of the mean ({avWithinXSigmaOfMeanNoises:.4f}) \n"
            )
        else:  # len(ValuesReplaced) == 0
            logger.info(
                f"Info: No Time Dependent Noise values were more than {n_noise_stddev_thresh} sigma below the mean "
                f"despite {n_suppressed_all_tels} suppressed pixels. Continuing normally."
            )
    else:  # n_suppressed_all_tels < n_suppresed_pixel_thresh:
        avNonNegativeNonSuppressedNoises = avNonNegativeNoises

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

    # Load the L3 scalar tree and use it to calculate the livetime
    L3ScalarTree = vegasFileIO.getTFilePtr().Get("/Diagnostics/LiveTime/L3ScalarTree")
    L3ScalarTree_df = ROOT.RDataFrame(L3ScalarTree)

    lastElapsedTime = (
        L3ScalarTree_df
        .Max("ElapsedTimeNow")
        .GetValue()
    )

    goodTimeStart, goodTimeStop = getGTArray(startTime_s, endTime_s, mergeTimeCut(tc), lastElapsedTime)
    ontime_after_timecuts = np.sum(np.array(goodTimeStop) - np.array(goodTimeStart))

    L3ScalarTreeElapsedTime = 0.0
    L3ScalarTreeLiveTime = 0.0

    for start, stop in zip(goodTimeStart - startTime_s, goodTimeStop - startTime_s):
        df_filt = (
            L3ScalarTree_df.Define("entry", "rdfentry_")
            .Filter(f"ElapsedTimeNow >= {start} and ElapsedTimeNow <= {stop}")
        )
        n = df_filt.Count().GetValue()
        if n == 0:
            logger.warning(f"Warning: no events between seconds {start} and {stop}")
            continue

        # Get start values
        df_sel = df_filt.Range(0, 0)
        start_entry = df_sel.Take["ULong64_t"]("entry").GetValue()[0]
        L3ScalarTree.GetEntry(start_entry)
        elapsedtime_start = L3ScalarTree.ElapsedTimeNow
        livetime_start = L3ScalarTree.LiveTimeNow

        # Get end values
        df_sel = df_filt.Range(n - 1, n)
        end_entry = df_sel.Take["ULong64_t"]("entry").GetValue()[0]
        L3ScalarTree.GetEntry(end_entry)
        elapsedtime_end = L3ScalarTree.ElapsedTimeNow
        livetime_end = L3ScalarTree.LiveTimeNow

        # Add elapsed and live times
        L3ScalarTreeElapsedTime += (elapsedtime_end - elapsedtime_start)
        L3ScalarTreeLiveTime    += (livetime_end - livetime_start)
    if L3ScalarTreeElapsedTime <= 0:
        raise ValueError("Total elapsed time from L3 scaler tree is 0. Cannot comput livetime fraction")
    LiveTimeFraction = L3ScalarTreeLiveTime/L3ScalarTreeElapsedTime

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
        evt_dict["EVENT_TYPE"] = arr_dict["nTelArr"]
        if "mswArr" in arr_dict:
            evt_dict["MSW"] = arr_dict["mswArr"]
        if "mslArr" in arr_dict:
            evt_dict["MSL"] = arr_dict["mslArr"]
        # Filling Header info
        evt_dict["OBS_ID"] = runHeader.getRunNumber()
        evt_dict["DATE-OBS"] = Time(
            startTime.getString(), format="iso", scale="utc"
        ).to_value("fits")
        evt_dict["TIME-OBS"] = startTime.getString().split()[1]
        evt_dict["DATE-END"] = Time(
            endTime.getString(), format="iso", scale="utc"
        ).to_value("fits")
        evt_dict["TIME-END"] = endTime.getString().split()[1]
        evt_dict["DATE-AVG"] = time_avg.to_value("fits")
        evt_dict["TSTART"] = startTime_s
        evt_dict["TSTOP"] = endTime_s
        evt_dict["ONTIME"] = ontime_after_timecuts
        evt_dict["LIVETIME"] = (
            LiveTimeFraction * ontime_after_timecuts
        )  # True to suppress error warnings
        evt_dict["DEADC"] = LiveTimeFraction
        evt_dict["OBJECT"] = runHeader.getSourceId()
        evt_dict["RA_PNT"] = avRA
        evt_dict["DEC_PNT"] = avDec
        evt_dict["ALT_PNT"] = avAlt
        evt_dict["AZ_PNT"] = avAz
        evt_dict["RA_OBJ"] = np.rad2deg(runHeader.getSourceRA())
        evt_dict["DEC_OBJ"] = np.rad2deg(runHeader.getSourceDec())
        evt_dict["TELLIST"] = produceTelList(runHeader.fRunInfo.fConfigMask)
        evt_dict["N_TELS"] = runHeader.pfRunDetails.fTels

    if st6_configs is not None:
        split_configs = {opt.split()[0]: opt.split()[1] for opt in st6_configs}
        if "EA_ApplyEnergyCorrectionForExperimentalBias" in split_configs.keys():
            corr_EB = bool(split_configs["EA_ApplyEnergyCorrectionForExperimentalBias"])
    if corr_EB:
        offset = SkyCoord(avRA * units.deg, avDec * units.deg).separation(
            SkyCoord(arr_dict["raArr"] * units.deg, arr_dict["decArr"] * units.deg)
        )
        offset = offset.degree
        # This is a problem but I don't know if it's VEGAS or V2DL3
        # The azimuth and zenith need to be used in the previous event
        # So this means that the first event has no reference....
        eList = np.array([arr_dict["energyArr"][0]])
        eList = np.append(
            eList,
            energyBiasCorr(
                arr_dict["energyArr"],
                arr_dict["azArr"],
                arr_dict["altArr"],
                arr_dict["fNoise"],
                offset,
                effective_area_files[event_class_idx],
                irf_to_store,
                psf_king_params,
            )[1:],
        ).flatten()
        evt_dict["ENERGY"] = eList
    else:
        evt_dict["ENERGY"] = arr_dict["energyArr"]

    return (
        {
            "goodTimeStart": goodTimeStart,
            "goodTimeStop": goodTimeStop,
            "TSTART": startTime_s,
            "TSTOP": endTime_s,
        },
        {
            "azimuth": avAz,
            "zenith": (90.0 - avAlt),
            "noise": avNonNegativeNonSuppressedNoises,
        },
        returned_dicts,
    )


"""
Check event against the provided FoV upper and lower bounds

Parameters:
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
        frame="icrs",
        unit=(units.deg, units.deg),
    )

    pointing_position = SkyCoord(
        np.rad2deg(reco.fArrayTrackingRA_J2000_Rad),
        np.rad2deg(reco.fArrayTrackingDec_J2000_Rad),
        frame="icrs",
        unit=(units.deg, units.deg),
    )

    tel_sep = pointing_position.separation(event_skycoord).degree

    if tel_sep > fov_cut_upper or tel_sep < fov_cut_lower:
        return True, tel_sep

    return False, tel_sep


"""
Applies the Energy Bias Correction for Experimental Bias (might not be necessary
as gammapy already has the migration matrix which it uses for energy reconstruction
but including here just in case and making it easy to turn off)

Parameters:
    energy -- Stage 5 energy of the Events (in TeV)
    effective_area_files -- The Effective Area


"""
# from pprint import pprint


def energyBiasCorr(
    energy,
    azimuth,
    zenith,
    noise,
    offset,
    effective_area_file,
    irf_to_store,
    psf_king_params,
):
    logger.debug("Using Energy Bias Correction")
    # axis_dict = effective_area_file.axis_dict
    manager = effective_area_file.manager
    energyCorr = []

    for i in range(0, len(energy)):
        shift = (
            i - 1
        )  # this is the kludge I don't know why this needs to be done it does not make sense
        effectiveAreaParameters = ROOT.VAEASimpleParameterData()
        effectiveAreaParameters.fAzimuth = azimuth[shift]
        effectiveAreaParameters.fZenith = 90 - zenith[shift]
        effectiveAreaParameters.fNoise = noise[shift]
        # effectiveAreaParameters.fOffset = offset[i]
        effectiveAreaParameters = manager.getVectorParamsFromSimpleParameterData(
            effectiveAreaParameters
        )

        correction = manager.getCorrectionForExperimentalBias(
            effectiveAreaParameters, energy[i] * 1000
        )
        logger.debug(f"Correction = {correction}, Original E = {energy[i] * 1000}")
        energytemp = energy[i] / correction
        energyCorr.append(energytemp)

    return energyCorr

    # offset_index = axis_dict["AbsoluteOffset"]
    # __, ebias_dict, __ = get_irf_not_safe(
    #  manager,
    #  offset_index,
    #  azimuth,
    #  zenith,
    #  noise,
    #  irf_to_store["point-like"],
    #    psf_king=psf_king_params,
    # )

    # energyCorr = np.zeros(len(energy))
    # correction = 1.0
    # for i in range(len(energy)):
    #   e_near = np.argwhere(ebias_dict["ELow"] > energy[i])[0][0]
    #   offset_near = np.argmin(np.abs(offset_index - offset[i]))
    #  mig_near = np.argmax(ebias_dict["Data"][offset_near, :, e_near])
    #  correction = (
    #      (
    #          ebias_dict["MigrationHigh"][mig_near]
    #          - ebias_dict["MigrationLow"][mig_near]
    #       )      #       / 2
    #   ) + ebias_dict["MigrationLow"][mig_near]
    #    if correction < 1e-14:
    #        correction = 1.0  # This correction should never be negative or zero
    #    energyCorr[i] = (energy[i]) / correction
    # return energyCorr
