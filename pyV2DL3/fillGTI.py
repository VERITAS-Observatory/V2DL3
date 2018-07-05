from astropy.io import fits
import numpy as np
import logging
from pyV2DL3.util import (getGTArray,
                          getTimeCut,
                          mergeTimeCut)

logger = logging.getLogger(__name__)

def fillGTI(vegasFileIO):
    runHeader = vegasFileIO.loadTheRunHeader()    
    cuts = vegasFileIO.loadTheCutsInfo()
    startTime = runHeader.getStartTime()
    endTime = runHeader.getEndTime()

    startTime_s = float(startTime.getDayNS()) / 1e9
    endTime_s = float(endTime.getDayNS()) / 1e9
    # Get Time Cuts and build GTI start and stop time array
    for k in cuts:
        tmp =k.fCutsFileText
        tc = getTimeCut(k.fCutsFileText)
    goodTimeStart,goodTimeStop = getGTArray(startTime_s,endTime_s,mergeTimeCut(tc))


    hdu2 = fits.BinTableHDU.from_columns([
    fits.Column(name='START', format='1D', array=goodTimeStart), 
    fits.Column(name='STOP', format='1D', array=goodTimeStop)
])
    hdu2.name = "GTI"
    hdu2.header.set('TSTART',startTime_s,'start time same unit and system used in the rate table')
    hdu2.header.set('TSTOP',endTime_s,'stop time same unit and system used in the rate table')
    hdu2.header.set('TIMEZERO',startTime_s,'zero time same unit and system used in the rate table') # is this correct ? 
    hdu2.header.set('TTYPE1','START   ' ,' start of good time interval')
    hdu2.header.set('TTYPE2','STOP    ' ,' start of good time interval')
    hdu2.header.set('EXTNAME','GTI     ' ,' name: Good Time Intervals')
    return hdu2
