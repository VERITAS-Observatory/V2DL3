from astropy.io import fits
import numpy as np
import logging

from pyV2DL3.util import produceTelList,decodeConfigMask
from pyV2DL3.util import (getGTArray,
                          getTimeCut,
                          mergeTimeCut)

logger = logging.getLogger(__name__)

windowSizeForNoise = 7

def fillEVENTS(vegasFileIO):

    # Load header ,array info and selected event tree ( vegas > v2.5.7)
    runHeader          = vegasFileIO.loadTheRunHeader()
    selectedEventsTree = vegasFileIO.loadTheCutEventTree()    
    qStatsData = vegasFileIO.loadTheQStatsData()
    pixelData = vegasFileIO.loadThePixelStatusData()
    arrayInfo          = vegasFileIO.loadTheArrayInfo(0)
    cuts = vegasFileIO.loadTheCutsInfo()

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
    detXArr = []
    detYArr = []
    nTelArr = []
    logger.debug("Start filling events ...")
    
    for ev in selectedEventsTree:
        evNumArr.append(ev.S.fArrayEventNum)
        timeArr.append(float(ev.S.fTime.getDayNS())/1e9)
        raArr.append(np.rad2deg(ev.S.fDirectionRA_J2000_Rad))
        decArr.append(np.rad2deg(ev.S.fDirectionDec_J2000_Rad))
        azArr.append(np.rad2deg(ev.S.fDirectionAzimuth_Rad))
        altArr.append(np.rad2deg(ev.S.fDirectionElevation_Rad))
        energyArr.append(ev.S.fEnergy_GeV / 1000.)
        detXArr.append(ev.S.fDirectionXCamPlane_Deg)
        detYArr.append(ev.S.fDirectionYCamPlane_Deg)
        nTelArr.append(ev.S.fImages)
        
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
    
    logger.debug("Done!")    
    logger.debug("Create EVENT HDU")    

    # Create HDU
    hdu1 = fits.BinTableHDU.from_columns([
    fits.Column(name='EVENT_ID', format='1J', array=evNumArr), 
    fits.Column(name='TIME', format='1D', array=timeArr, unit="s"), 
    fits.Column(name='RA', format='1E', array=raArr, unit = "deg"), 
    fits.Column(name='DEC', format='1E', array=decArr, unit = "deg"), 
    fits.Column(name='ALT', format='1E', array=altArr, unit = "deg"), 
    fits.Column(name='AZ', format='1E', array=azArr, unit = "deg"), 
    fits.Column(name='ENERGY', format='1E', array=energyArr, unit = "TeV"), 
    # Remove to avoid confusion
    # fits.Column(name='DETX', format='1E', array=detXArr, unit = "deg"), 
    # fits.Column(name='DETY', format='1E', array=detYArr, unit = "deg"),
    fits.Column(name="EVENT_TYPE", format="1J", array=nTelArr)
    ])
    hdu1.name = "EVENTS"

    # Fill Header
    hdu1.header.set('OBS_ID  ', runHeader.getRunNumber(), 'Run Number')
    hdu1.header.set('TELESCOP', 'VERITAS', 'Data from VERITAS')
    
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

    hdu1.header.set('DATE-OBS',
                    startTime.getString().split()[0],
                    'start date (UTC) of obs yy-mm-dd')
    hdu1.header.set('TIME-OBS',
                    startTime.getString().split()[1],
                    'start time (UTC) of obs hh-mm-ss')
    hdu1.header.set('DATE-END',
                    endTime.getString().split()[0],
                    'end date (UTC) of obs yy-mm-dd')
    hdu1.header.set('TIME-END',
                    endTime.getString().split()[1],
                    'end time (UTC) of obs hh-mm-ss')
    
    hdu1.header.set('TSTART  ',
                    startTime_s,
                    'mission time of start of obs [s]')
    hdu1.header.set('TSTOP   ',
                    endTime_s,
                    'mission time of end of obs [s]')
    hdu1.header.set('MJDREFI ',
                    int(startTime.getMJDInt()), 'int part of reference MJD [days]')
    hdu1.header.set('MJDREFF ', 0., 'fractional part of reference MJD [days]')
    
    hdu1.header.set('TIMEUNIT', 's', 'time unit is seconds since MET start')
    hdu1.header.set('TIMESYS ', 'utc', 'time scale is UTC')
    hdu1.header.set('TIMEREF ', 'local', 'local time reference')
    
    hdu1.header.set('ONTIME  ', 
                    endTime_s - startTime_s,
                    'time on target (including deadtime)')
#    hdu1.header.set('LIVETIME', runHeader.pfRunDetails.fRunNominalLiveTimeSeconds,
#                    '(dead=ONTIME-LIVETIME) [s] ')
    # Correct live time for time cuts
    hdu1.header.set('LIVETIME', runHeader.getLiveTimeFrac()*real_live_time,
                    '(dead=ONTIME-LIVETIME) [s] ')

    hdu1.header.set('DEADC   ', runHeader.getLiveTimeFrac(),
                    'Average dead time correction (LIVETIME/ONTIME)')
    
    hdu1.header.set('OBJECT  ', runHeader.getSourceId(), 'observed object')
    
    hdu1.header.set('RA_PNT  ', avRA, 'observation position RA [deg]')
    hdu1.header.set('DEC_PNT ', avDec, 'observation position DEC [deg]')
    hdu1.header.set('ALT_PNT ', avAlt, 'average altitude of pointing [deg]')
    hdu1.header.set('AZ_PNT  ', avAz, 'average azimuth of pointing [deg]')
    
    hdu1.header.set('RA_OBJ  ',
                    np.rad2deg(runHeader.getSourceRA()),
                    'observation position RA [deg]')
    hdu1.header.set('DEC_OBJ ',
                    np.rad2deg(runHeader.getSourceDec()),
                    'observation position DEC [deg]')
    
    # get the list of telescopes that participate in the event
    hdu1.header.set('TELLIST',
                    produceTelList(runHeader.fRunInfo.fConfigMask),
                    'comma-separated list of tel IDs')
    hdu1.header.set('N_TELS', runHeader.pfRunDetails.fTels,
                    'number of telescopes in event list')
    
    # other info - weather? pointing mode
    
    hdu1.header.set('EUNIT   ', 'TeV', 'energy unit')
    hdu1.header.set('GEOLON  ', np.rad2deg(arrayInfo.longitudeRad()), 'longitude of array center [deg]')
    hdu1.header.set('GEOLAT  ', np.rad2deg(arrayInfo.latitudeRad()), 'latitude of array center [deg]')
    hdu1.header.set('ALTITUDE', arrayInfo.elevationM(), 'altitude of array center [m]')
    
    # What are these for? - May note be needed, leave out for now.
    # hdu1.header.set('DSTYP1', 'TIME    ', 'Data selection type')
    # hdu1.header.set('DSUNI1', 's       ', 'Data selection unit')
    # hdu1.header.set('DSVAL1', 'TABLE   ', 'Data selection value')
    # hdu1.header.set('DSREF1', ':GTI    ', 'Data selection reference')
    # hdu1.header.set('DSTYP2', 'POS(RA,DEC)', 'Data selection type')
    # hdu1.header.set('DSUNI2', 'deg     ', 'Data selection unit')
    # hdu1.header.set('DSVAL2', 'CIRCLE(83.63,22.01,5)', 'Data selection value')
    # hdu1.header.set('DSTYP3', 'ENERGY  ', 'Data selection type')
    # hdu1.header.set('DSUNI3', 'TeV     ', 'Data selection unit')
    # hdu1.header.set('DSVAL3', '0.05:100', 'Data selection value')
    # hdu1.header.set('NDSKEYS', '3       ', 'Number of data selections')

    # Calculate average noise
    avNoise = 0
    nTels = 0
    for telID in decodeConfigMask(runHeader.fRunInfo.fConfigMask):
        avNoise += qStatsData.getCameraAverageTraceVarTimeIndpt(telID-1, windowSizeForNoise, pixelData, arrayInfo)
        nTels += 1
    
    avNoise /= nTels
    return {'goodTimeStart':goodTimeStart,'goodTimeStop':goodTimeStop},{'azimuth':avAz,'zenith':(90. - avAlt),'noise':avNoise},hdu1
