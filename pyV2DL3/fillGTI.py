from astropy.io import fits
import numpy as np
import logging

logger = logging.getLogger(__name__)

def fillGTI(vegasFileIO,goodTimeStart=None,goodTimeStop=None):
    runHeader = vegasFileIO.loadTheRunHeader()    
    cuts = vegasFileIO.loadTheCutsInfo()
    startTime = runHeader.getStartTime()
    endTime = runHeader.getEndTime()

    startTime_s = float(startTime.getDayNS()) / 1e9
    endTime_s = float(endTime.getDayNS()) / 1e9

    if(goodTimeStart is None):
        goodTimeStart = np.array([startTime_s])
    if(goodTimeStop is None):
        goodTimeStop = np.array([endTime_s])


    hdu2 = fits.BinTableHDU.from_columns([
    fits.Column(name='START', format='1D', array=goodTimeStart), 
    fits.Column(name='STOP', format='1D', array=goodTimeStop)
])
    hdu2.name = "GTI"
    hdu2.header.set('TSTART',startTime_s,'start time [s]')
    hdu2.header.set('TSTOP',endTime_s,'stop time same [s]')
    hdu2.header.set('TIMEZERO',startTime_s,'zero time [s]') # is this correct ? 
    hdu2.header.set('TTYPE1','START   ' ,' start of good time interval')
    hdu2.header.set('TTYPE2','STOP    ' ,' start of good time interval')
    hdu2.header.set('EXTNAME','GTI     ' ,' name: Good Time Intervals')
    return hdu2
