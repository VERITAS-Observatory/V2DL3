import numpy as np
import logging
import uproot4
from pyV2DL3.eventdisplay.util import produce_tel_list
#from root_numpy import tree2array
from astropy.time import Time
from pyV2DL3.constant import VTS_REFERENCE_MJD, VTS_REFERENCE_LAT, VTS_REFERENCE_LON, VTS_REFERENCE_HEIGHT

logger = logging.getLogger(__name__)
windowSizeForNoise = 7


def __fillEVENTS_not_safe__(edFileIO):
    # EventDisplay imports
    #from ROOT import VEvndispRunParameter, VSkyCoordinatesUtilities, VAnaSumRunParameter

    evt_dict = {}

    #implementing changes to only use uproot

    # Load required trees within the anasum file:
    #filename = str(edFileIO.GetName())
    file = uproot4.open(edFileIO)

    #runSummary = tree2array(edFileIO.Get("total_1/stereo/tRunSummary"))
    runSummary = file['total_1/stereo/tRunSummary'].arrays(library='np')
    runNumber = runSummary['runOn'][0]

    #telConfig = tree2array(edFileIO.Get("run_{}/stereo/telconfig".format(runNumber)))
    telConfig = file['run_{}/stereo/telconfig'.format(runNumber)].arrays(library='np')
    #print (telConfig)
    #vAnaSumRunParameter = edFileIO.Get("run_{}/stereo/VAnaSumRunParameter".format(runNumber))
    #runParametersV2 = edFileIO.Get("run_{}/stereo/runparameterV2".format(runNumber))
    #selectedEventsTree = tree2array(edFileIO.Get("run_{}/stereo/TreeWithEventsForCtools".format(runNumber)))
    # qStatsData = edFileIO.loadTheQStatsData()
    # pixelData = edFileIO.loadThePixelStatusData()
    # arrayInfo          = edFileIO.loadTheArrayInfo(0)
    # cuts = edFileIO.loadTheCutsInfo()

    # Get start and stop time within the run.
    #startDateTime = (runParametersV2.fDBRunStartTimeSQL).split(" ")
    start_mjd = file['total_1/stereo/tRunSummary/MJDrunstart'].array(library='np')[0]
    #stopDateTime = (runParametersV2.fDBRunStoppTimeSQL).split(" ")
    stop_mjd = file['total_1/stereo/tRunSummary/MJDrunstop'].array(library='np')[0]
    print ("start,stop mjd", start_mjd,stop_mjd)
    #convert mjd to fits format
    t_start_fits = Time(start_mjd, format='mjd', scale='utc').to_value('fits')
    t_stop_fits = Time(stop_mjd, format='mjd', scale='utc').to_value('fits')
    print ("t_start_fits,stop fits",t_start_fits,t_stop_fits)
    #start_year, start_month, start_day = [int(k) for k in startDateTime[0].split("-")]
    #stop_year, stop_month, stop_day = [int(k) for k in stopDateTime[0].split("-")]

    #start_mjd = VSkyCoordinatesUtilities.getMJD(start_year, start_month, start_day)
    #stop_mjd = VSkyCoordinatesUtilities.getMJD(stop_year, stop_month, stop_day)

    #deadtime = vAnaSumRunParameter.fScalarDeadTimeFrac
    deadtime = file['total_1/stereo/tRunSummary/DeadTimeFracOn'].array(library='np')[0]

    # Number of seconds between reference time and run MJD at 00:00:00:
    t_ref = Time(VTS_REFERENCE_MJD, format='mjd', scale='utc')
    seconds_from_reference = (Time(start_mjd, format='mjd', scale='utc') - t_ref).sec

    #check further if acceptable using only MJDrunstart: 2012-10-13T10:52:08.106 MJDrunstop: 2012-10-13T11:12:09.895 before: runParametersV2.fDBRunStartTimeSQL: 2012-10-13 10:52:08 runParametersV2.fDBRunStoppTimeSQL: 2012-10-13 11:12:09
    #tstart_from_reference = (Time(runParametersV2.fDBRunStartTimeSQL, format='iso', scale='utc') - t_ref).sec
    #tstop_from_reference = (Time(runParametersV2.fDBRunStoppTimeSQL, format='iso', scale='utc') - t_ref).sec
    tstart_from_reference = (Time(start_mjd, format='mjd', scale='utc') - t_ref).sec
    tstop_from_reference = (Time(stop_mjd, format='mjd', scale='utc') - t_ref).sec
    #print ("tstart from reference:",Time(runParametersV2.fDBRunStartTimeSQL, format='iso', scale='utc'))
    #print ("MJDstart:",start_mjd)


    DL3EventTree = file['run_{}/stereo/DL3EventTree'.format(runNumber)].arrays(library='np')
    #evNumArr = selectedEventsTree['eventNumber']
    evNumArr = DL3EventTree['eventNumber']
    # This should already have microsecond resolution if stored with double precision.
    time_of_day = DL3EventTree['timeOfDay']
    timeArr = seconds_from_reference + time_of_day
    #timeArr = seconds_from_reference + selectedEventsTree['timeOfDay']

    #changing all to the new DL3EventTree
    raArr = DL3EventTree['RA']
    decArr = DL3EventTree['DEC']
    azArr = DL3EventTree['Az']
    altArr = DL3EventTree['El']
    #offset = DL3EventTree['Woff']
    energyArr = DL3EventTree['Energy']
    # Not used for the moment by science tools.
    nTelArr = DL3EventTree['NImages']

    avAlt = np.mean(altArr)
    # Calculate average azimuth angle from average vector on a circle
    # https://en.wikipedia.org/wiki/Mean_of_circular_quantities
    avAz_rad = np.deg2rad(azArr)
    avAz = np.rad2deg(np.arctan2(np.sum(np.sin(avAz_rad)),np.sum(np.cos(avAz_rad))))
    avAz = avAz if avAz > 0 else avAz + 360

    # RA and DEC already in degrees.
    avRA = np.mean(raArr)
    avDec = np.mean(decArr)

    # Filling Event List
    evt_dict['EVENT_ID'] = evNumArr
    evt_dict['TIME'] = timeArr
    evt_dict['RA'] = raArr
    evt_dict['DEC'] = decArr
    evt_dict['ALT'] = altArr
    evt_dict['AZ'] = azArr
    #evt_dict['OFFSET'] = offset #new
    evt_dict['ENERGY'] = energyArr
    evt_dict['EVENT_TYPE'] = nTelArr

    # FIXME: Get Time Cuts and build GTI start and stop time array
    # for k in cuts:
    #     tmp =k.fCutsFileText
    #     tc = getTimeCut(k.fCutsFileText)
    #
    # goodTimeStart,goodTimeStop = getGTArray(startTime_s,endTime_s,mergeTimeCut(tc))
    # real_live_time = np.sum(goodTimeStop - goodTimeStart)
    # startTime_s = float(startTime.getDayNS()) / 1e9
    # endTime_s = float(endTime.getDayNS()) / 1e9

    # Filling Header info
    evt_dict['OBS_ID'] = runNumber  # this does not allow for event type files
    #check if one datetime in fits format is enough
    #evt_dict['DATE-OBS'] = startDateTime[0]
    evt_dict['DATE-OBS']= t_start_fits
    #evt_dict['TIME-OBS'] = startDateTime[1]
    #evt_dict['DATE-END'] = stopDateTime[0]
    evt_dict['DATE-END'] = t_stop_fits
    #evt_dict['TIME-END'] = stopDateTime[1]
    evt_dict['TSTART'] = tstart_from_reference
    evt_dict['TSTOP'] = tstop_from_reference
    evt_dict['MJDREFI'] = int(VTS_REFERENCE_MJD)
    evt_dict['ONTIME'] = tstop_from_reference - tstart_from_reference
    evt_dict['LIVETIME'] = (tstop_from_reference - tstart_from_reference) * (1 - deadtime)
    evt_dict['DEADC'] = 1 - deadtime
    #evt_dict['OBJECT'] = runParametersV2.fTargetName
    evt_dict['OBJECT'] = runSummary['TargetName'][0]
    evt_dict['RA_PNT'] = avRA
    evt_dict['DEC_PNT'] = avDec
    evt_dict['ALT_PNT'] = avAlt
    evt_dict['AZ_PNT'] = avAz
    #evt_dict['RA_OBJ'] = runParametersV2.fTargetRA
    evt_dict['RA_OBJ'] = runSummary['TargetRAJ2000'][0]
    #evt_dict['DEC_OBJ'] = runParametersV2.fTargetDec
    evt_dict['DEC_OBJ'] = runSummary['TargetDecJ2000'][0]
    evt_dict['TELLIST'] = produce_tel_list(telConfig)
    evt_dict['N_TELS'] = len(telConfig['TelID'])
    evt_dict['GEOLON'] = VTS_REFERENCE_LON
    evt_dict['GEOLAT'] = VTS_REFERENCE_LAT
    evt_dict['ALTITUDE'] = VTS_REFERENCE_HEIGHT

    avNoise = runSummary['pedvarsOn'][0]

    # FIXME: For now we are not including any good time interval (GTI).
    # This should be improved in the future, reading the time masks.
    return ({'goodTimeStart': [tstart_from_reference], 'goodTimeStop': [tstop_from_reference],
             'TSTART': tstart_from_reference, 'TSTOP': tstop_from_reference},
            {'azimuth': avAz, 'zenith': (90. - avAlt), 'noise': avNoise},
            evt_dict)
