import numpy as np
import logging
from root_numpy import tree2array

# from pyV2DL3.vegas.util import produceTelList,decodeConfigMask
# from pyV2DL3.vegas.util import (getGTArray,
#                           getTimeCut,
#                           mergeTimeCut)

logger = logging.getLogger(__name__)

windowSizeForNoise = 7

def __fillEVENTS_not_safe__(edFileIO):
    evt_dict = {}

    # Load header ,array info and selected event tree ( vegas > v2.5.7)
    runHeader       = tree2array(edFileIO.Get("total_1/stereo/tRunSummary"))
    runNumber       = runHeader['runOn'][0]
    selectedEventsTree = tree2array(edFileIO.Get("run_{}/stereo/TreeWithEventsForCtools".format(runNumber)))
    # qStatsData = edFileIO.loadTheQStatsData()
    # pixelData = edFileIO.loadThePixelStatusData()
    # arrayInfo          = edFileIO.loadTheArrayInfo(0)
    # cuts = edFileIO.loadTheCutsInfo()

    # Start filling events
    avAlt = []
    avAz = []
    avRA = []
    avDec = []

    evNumArr = []
    timeArr = []
    raArr = []
    decArr = []
    azArr = []
    altArr = []
    energyArr = []
    nTelArr = []
    logger.debug("Start filling events ...")

    for ev in selectedEventsTree:
        evNumArr.append(selectedEventsTree['eventNumber'])
        # FIXME: TimeOfDay is NOT what we need!!
        # timeOfDay: double con suficiente precision. Segundos despu'es de MJDOn??
        timeArr.append(selectedEventsTree['timeOfDay'])
        raArr.append(selectedEventsTree['RA'])
        decArr.append(selectedEventsTree['DEC'])
        azArr.append(selectedEventsTree['Az'])
        altArr.append(selectedEventsTree['El'])
        energyArr.append(selectedEventsTree['Energy'])
        # nTelArr.append(ev.S.fImages)

        avAlt.append(ev.S.fArrayTrackingElevation_Deg)
        avAz.append(ev.S.fArrayTrackingAzimuth_Deg)
        avRA.append(ev.S.fArrayTrackingRA_J2000_Rad)
        avDec.append(ev.S.fArrayTrackingDec_J2000_Rad)

    avAlt = np.mean(avAlt)
    # Calculate average azimuth angle from average vector on a circle
    # https://en.wikipedia.org/wiki/Mean_of_circular_quantities
    avAz_deg = np.deg2rad(avAz)
    avAz = np.rad2deg(np.arctan2(np.sum(np.sin(avAz_deg)),np.sum(np.cos(avAz_deg))))
    avAz = avAz if avAz > 0 else avAz + 360

    avRA = np.rad2deg(np.mean(avRA))
    avDec = np.rad2deg(np.mean(avDec))
    # Filling Event List
    evt_dict['EVENT_ID'] = evNumArr
    evt_dict['TIME']     = timeArr
    evt_dict['RA']       = raArr
    evt_dict['DEC']      = decArr
    evt_dict['ALT']      = altArr
    evt_dict['AZ']       = azArr
    evt_dict['ENERGY']   = energyArr
    # evt_dict['EVENT_TYPE'] =nTelArr


    # Calculate Live Time
    startTime = runHeader.getStartTime()
    endTime = runHeader.getEndTime()
    
    startTime_s = float(startTime.getDayNS()) / 1e9
    endTime_s = float(endTime.getDayNS()) / 1e9
    startTime = runHeader.getStartTime()
    endTime = runHeader.getEndTime()
    # Get Time Cuts and build GTI start and stop time array
    for k in cuts:
        tmp =k.fCutsFileText
        tc = getTimeCut(k.fCutsFileText)
    
    goodTimeStart,goodTimeStop = getGTArray(startTime_s,endTime_s,mergeTimeCut(tc))
    real_live_time = np.sum(goodTimeStop - goodTimeStart)
    
    startTime_s = float(startTime.getDayNS()) / 1e9
    endTime_s = float(endTime.getDayNS()) / 1e9 

    # Filling Header info
    evt_dict['OBS_ID'] = runHeader.getRunNumber()
    evt_dict['DATE-OBS'] = startTime.getString().split()[0]
    evt_dict['TIME-OBS'] = startTime.getString().split()[1]
    evt_dict['DATE-END'] = endTime.getString().split()[0]
    evt_dict['TIME-END'] = endTime.getString().split()[1]
    evt_dict['TSTART']   = startTime_s
    evt_dict['TSTOP']    = endTime_s
    evt_dict['MJDREFI']  = int(startTime.getMJDInt())
    evt_dict['ONTIME']   = endTime_s - startTime_s
    evt_dict['LIVETIME'] = runHeader.getLiveTimeFrac()*real_live_time
    evt_dict['DEADC']    =  evt_dict['LIVETIME']/evt_dict['ONTIME']
    evt_dict['OBJECT']   = runHeader.getSourceId()
    evt_dict['RA_PNT']   = avRA
    evt_dict['DEC_PNT']   = avDec
    evt_dict['ALT_PNT']   = avAlt
    evt_dict['AZ_PNT']    = avAz 
    evt_dict['RA_OBJ']    = np.rad2deg(runHeader.getSourceRA()) 
    evt_dict['DEC_OBJ']    = np.rad2deg(runHeader.getSourceDec()) 
    evt_dict['TELLIST']    = produceTelList(runHeader.fRunInfo.fConfigMask) 
    evt_dict['N_TELS']    =runHeader.pfRunDetails.fTels 
    evt_dict['GEOLON']    = np.rad2deg(arrayInfo.longitudeRad())
    evt_dict['GEOLAT']    = np.rad2deg(arrayInfo.latitudeRad())
    evt_dict['ALTITUDE']  = arrayInfo.elevationM()

    avNoise = 0
    nTels = 0
    for telID in decodeConfigMask(runHeader.fRunInfo.fConfigMask):
        avNoise += qStatsData.getCameraAverageTraceVarTimeIndpt(telID-1, windowSizeForNoise, pixelData, arrayInfo)
        nTels += 1
    
    avNoise /= nTels
    return ({'goodTimeStart':goodTimeStart,'goodTimeStop':goodTimeStop,'TSTART':startTime_s,'TSTOP':endTime_s},
           {'azimuth':avAz,'zenith':(90. - avAlt),'noise':avNoise},
           evt_dict)
