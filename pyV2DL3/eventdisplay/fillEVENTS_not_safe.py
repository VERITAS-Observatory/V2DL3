import numpy as np
import logging
from pyV2DL3.eventdisplay.util import produceTelList
from root_numpy import tree2array
from astropy.time import Time

# EventDisplay imports
from ROOT import VEvndispRunParameter, VSkyCoordinatesUtilities, VAnaSumRunParameter


# from pyV2DL3.vegas.util import produceTelList,decodeConfigMask
# from pyV2DL3.vegas.util import (getGTArray,
#                           getTimeCut,
#                           mergeTimeCut)

logger = logging.getLogger(__name__)
windowSizeForNoise = 7


def __fillEVENTS_not_safe__(edFileIO):
    evt_dict = {}

    # FIXME: This should be taken from a common script (by VEGAS also)
    reference_mjd = 53402.0

    # Load required trees within the anasum file:
    runSummary = tree2array(edFileIO.Get("total_1/stereo/tRunSummary"))
    runNumber = runSummary['runOn'][0]
    telConfig = tree2array(edFileIO.Get("run_{}/stereo/telconfig".format(runNumber)))
    vAnaSumRunParameter = edFileIO.Get("run_{}/stereo/VAnaSumRunParameter".format(runNumber))
    runParametersV2 = edFileIO.Get("run_{}/stereo/runparameterV2".format(runNumber))
    selectedEventsTree = tree2array(edFileIO.Get("run_{}/stereo/TreeWithEventsForCtools".format(runNumber)))
    # qStatsData = edFileIO.loadTheQStatsData()
    # pixelData = edFileIO.loadThePixelStatusData()
    # arrayInfo          = edFileIO.loadTheArrayInfo(0)
    # cuts = edFileIO.loadTheCutsInfo()

    # Get start and stop time within the run.
    startDateTime = (runParametersV2.fDBRunStartTimeSQL).split(" ")
    stopDateTime = (runParametersV2.fDBRunStoppTimeSQL).split(" ")

    start_year, start_month, start_day = [int(k) for k in startDateTime[0].split("-")]
    stop_year, stop_month, stop_day = [int(k) for k in stopDateTime[0].split("-")]

    start_mjd = VSkyCoordinatesUtilities.getMJD(start_year, start_month, start_day)
    stop_mjd = VSkyCoordinatesUtilities.getMJD(stop_year, stop_month, stop_day)

    deadtime = vAnaSumRunParameter.fScalarDeadTimeFrac

    # Number of seconds between reference time and run MJD at 00:00:00:
    t_ref = Time(reference_mjd, format='mjd', scale='utc')
    seconds_from_reference = (Time(start_mjd, format='mjd', scale='utc') - t_ref).sec
    tstart_from_reference = (Time(runParametersV2.fDBRunStartTimeSQL, format='iso', scale='utc') - t_ref).sec
    tstop_from_reference = (Time(runParametersV2.fDBRunStoppTimeSQL, format='iso', scale='utc') - t_ref).sec

    # Start filling events
    avAlt = []
    avAz = []
    avRA = []
    avDec = []

    evNumArr = selectedEventsTree['eventNumber']
    # This should already have microsecond resolution if stored with double precision.
    timeArr = seconds_from_reference + selectedEventsTree['timeOfDay']
    raArr = selectedEventsTree['RA']
    decArr = selectedEventsTree['DEC']
    azArr = selectedEventsTree['Az']
    altArr = selectedEventsTree['El']
    energyArr = selectedEventsTree['Energy']
    # Not used for the moment by science tools.
    # nTelArr = selectedEventsTree['NImages']

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
    evt_dict['TIME'] = timeArr
    evt_dict['RA'] = raArr
    evt_dict['DEC'] = decArr
    evt_dict['ALT'] = altArr
    evt_dict['AZ'] = azArr
    evt_dict['ENERGY'] = energyArr
    # evt_dict['EVENT_TYPE'] =nTelArr

    # FIXME: Get Time Cuts and build GTI start and stop time array
    for k in cuts:
        tmp =k.fCutsFileText
        tc = getTimeCut(k.fCutsFileText)
    
    goodTimeStart,goodTimeStop = getGTArray(startTime_s,endTime_s,mergeTimeCut(tc))
    real_live_time = np.sum(goodTimeStop - goodTimeStart)
    
    startTime_s = float(startTime.getDayNS()) / 1e9
    endTime_s = float(endTime.getDayNS()) / 1e9 

    # Filling Header info
    evt_dict['OBS_ID'] = runNumber
    evt_dict['DATE-OBS'] = startDateTime[0]
    evt_dict['TIME-OBS'] = startDateTime[1]
    evt_dict['DATE-END'] = stopDateTime[0]
    evt_dict['TIME-END'] = stopDateTime[1]
    evt_dict['TSTART'] = tstart_from_reference
    evt_dict['TSTOP'] = tstop_from_reference
    evt_dict['MJDREFI'] = int(reference_mjd)
    evt_dict['ONTIME'] = tstop_from_reference - tstart_from_reference
    evt_dict['LIVETIME'] = (tstop_from_reference - tstart_from_reference) * (1 - deadtime)
    evt_dict['DEADC'] = 1 - deadtime
    evt_dict['OBJECT'] = runParametersV2.fTargetName
    evt_dict['RA_PNT'] = avRA
    evt_dict['DEC_PNT'] = avDec
    evt_dict['ALT_PNT'] = avAlt
    evt_dict['AZ_PNT'] = avAz
    evt_dict['RA_OBJ'] = runParametersV2.fTargetRA
    evt_dict['DEC_OBJ'] = runParametersV2.fTargetDec
    evt_dict['TELLIST'] = produceTelList(telConfig)
    evt_dict['N_TELS'] = len(telConfig['TelID'])
    evt_dict['GEOLON'] = np.rad2deg(arrayInfo.longitudeRad())
    evt_dict['GEOLAT'] = np.rad2deg(arrayInfo.latitudeRad())
    evt_dict['ALTITUDE'] = arrayInfo.elevationM()

    31.6747333333334, -110.9528

    avNoise = 0
    nTels = 0
    for telID in decodeConfigMask(runHeader.fRunInfo.fConfigMask):
        avNoise += qStatsData.getCameraAverageTraceVarTimeIndpt(telID-1, windowSizeForNoise, pixelData, arrayInfo)
        nTels += 1
    
    avNoise /= nTels
    return ({'goodTimeStart':goodTimeStart,'goodTimeStop':goodTimeStop,'TSTART':startTime_s,'TSTOP':endTime_s},
           {'azimuth':avAz,'zenith':(90. - avAlt),'noise':avNoise},
           evt_dict)
