import numpy as np
import logging

from pyV2DL3.vegas.util import produceTelList, decodeConfigMask
from pyV2DL3.vegas.util import getGTArray, getTimeCut, mergeTimeCut
from pyV2DL3.constant import VTS_REFERENCE_MJD, VTS_REFERENCE_LAT, VTS_REFERENCE_LON, VTS_REFERENCE_HEIGHT
from astropy.time import Time

logger = logging.getLogger(__name__)

windowSizeForNoise = 7


def __fillEVENTS_not_safe__(vegasFileIO):
    evt_dict = {}

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

    t_ref = Time(VTS_REFERENCE_MJD,format='mjd',scale='utc')
    seconds_from_reference_t0 =  (Time(mjd,format='mjd',scale='utc') - t_ref).sec
    startTime_s = startTime_s + seconds_from_reference_t0
    endTime_s   = endTime_s   + seconds_from_reference_t0
    
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
        # seconds since first light         
        time_relative_to_reference = (float(ev.S.fTime.getDayNS())/1e9 +
                                      seconds_from_reference_t0)
        evNumArr.append(ev.S.fArrayEventNum)
        timeArr.append(time_relative_to_reference)
        raArr.append(np.rad2deg(ev.S.fDirectionRA_J2000_Rad))
        decArr.append(np.rad2deg(ev.S.fDirectionDec_J2000_Rad))
        azArr.append(np.rad2deg(ev.S.fDirectionAzimuth_Rad))
        altArr.append(np.rad2deg(ev.S.fDirectionElevation_Rad))
        energyArr.append(ev.S.fEnergy_GeV / 1000.)
        nTelArr.append(ev.S.fImages)

        avAlt.append(ev.S.fArrayTrackingElevation_Deg)
        avAz.append(ev.S.fArrayTrackingAzimuth_Deg)
        avRA.append(ev.S.fArrayTrackingRA_J2000_Rad)
        avDec.append(ev.S.fArrayTrackingDec_J2000_Rad)

    avAlt = np.mean(avAlt)
    # Calculate average azimuth angle from average vector on a circle
    # https://en.wikipedia.org/wiki/Mean_of_circular_quantities
    avAz_rad = np.deg2rad(avAz)
    avAz = np.rad2deg(np.arctan2(np.sum(np.sin(avAz_rad)),np.sum(np.cos(avAz_rad))))
    avAz = avAz if avAz > 0 else avAz + 360

    avRA = np.rad2deg(np.mean(avRA))
    avDec = np.rad2deg(np.mean(avDec))
    # Filling Event List
    evt_dict['EVENT_ID'] = evNumArr
    evt_dict['TIME'] = timeArr
    evt_dict['RA'] = raArr
    evt_dict['DEC'] = decArr
    evt_dict['ALT'] = altArr
    evt_dict['AZ'] = azArr
    evt_dict['ENERGY'] = energyArr
    evt_dict['EVENT_TYPE'] = nTelArr

    # Get Time Cuts and build GTI start and stop time array
    # and calculate live time
    for k in cuts:
        tmp =k.fCutsFileText
        tc = getTimeCut(k.fCutsFileText)

    goodTimeStart, goodTimeStop = getGTArray(startTime_s,endTime_s,mergeTimeCut(tc))
    real_live_time = np.sum(np.array(goodTimeStop) - np.array(goodTimeStart))

    # Filling Header info
    evt_dict['OBS_ID'] = runHeader.getRunNumber()
    evt_dict['DATE-OBS'] = startTime.getString().split()[0]
    evt_dict['TIME-OBS'] = startTime.getString().split()[1]
    evt_dict['DATE-END'] = endTime.getString().split()[0]
    evt_dict['TIME-END'] = endTime.getString().split()[1]
    evt_dict['TSTART'] = startTime_s
    evt_dict['TSTOP'] = endTime_s
    evt_dict['ONTIME'] = endTime_s - startTime_s
    evt_dict['LIVETIME'] = runHeader.getLiveTimeFrac(True)*real_live_time # True to suppress error warnings
    evt_dict['DEADC'] = evt_dict['LIVETIME']/evt_dict['ONTIME']
    evt_dict['OBJECT'] = runHeader.getSourceId()
    evt_dict['RA_PNT'] = avRA
    evt_dict['DEC_PNT'] = avDec
    evt_dict['ALT_PNT'] = avAlt
    evt_dict['AZ_PNT'] = avAz
    evt_dict['RA_OBJ'] = np.rad2deg(runHeader.getSourceRA())
    evt_dict['DEC_OBJ'] = np.rad2deg(runHeader.getSourceDec())
    evt_dict['TELLIST'] = produceTelList(runHeader.fRunInfo.fConfigMask)
    evt_dict['N_TELS'] = runHeader.pfRunDetails.fTels

    avNoise = 0
    nTels = 0
    for telID in decodeConfigMask(runHeader.fRunInfo.fConfigMask):
        avNoise += qStatsData.getCameraAverageTraceVarTimeIndpt(telID-1, windowSizeForNoise, pixelData, arrayInfo)
        nTels += 1
    
    avNoise /= nTels
    return ({'goodTimeStart': goodTimeStart, 'goodTimeStop': goodTimeStop, 'TSTART': startTime_s, 'TSTOP': endTime_s},
            {'azimuth': avAz, 'zenith': (90. - avAlt), 'noise': avNoise},
            evt_dict)
