from astropy.io import fits
import numpy as np
import logging

logger = logging.getLogger(__name__)

def fillGTI(datasource,goodTimeStart=None,goodTimeStop=None):
    gti_dict = datasource.get_gti_data()
    goodTimeStart = gti_dict['goodTimeStart']     
    goodTimeStop = gti_dict['goodTimeStop']     
    startTime_s = gti_dict['TSTART']
    endTime_s   = gti_dict['TSTOP']
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
