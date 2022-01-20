from astropy.io import fits
import numpy as np
import logging
from pyV2DL3.addHDUClassKeyword import addHDUClassKeyword
from pyV2DL3.constant import (VERSION,VTS_REFERENCE_MJD,VTS_REFERENCE_LAT,
                              VTS_REFERENCE_LON,VTS_REFERENCE_HEIGHT,RADECSYS,EQUINOX)

logger = logging.getLogger(__name__)

def fillEVENTS(datasource, save_multiplicity=False):
    evt_dict = datasource.get_evt_data()
    logger.debug("Create EVENT HDU")    

    # Columns to be saved
    columns = [fits.Column(name='EVENT_ID', format='1K', array=evt_dict['EVENT_ID']),
               fits.Column(name='TIME', format='1D', array=evt_dict['TIME'], unit="s"),
               fits.Column(name='RA', format='1E', array=evt_dict['RA'], unit="deg"),
               fits.Column(name='DEC', format='1E', array=evt_dict['DEC'], unit="deg"),
               # Not necessary for gammapy and spec
               # fits.Column(name='ALT', format='1E', array=evt_dict['ALT'], unit = "deg"),
               # fits.Column(name='AZ', format='1E', array=evt_dict['AZ'], unit = "deg"),
               fits.Column(name='ENERGY', format='1E', array=evt_dict['ENERGY'], unit="TeV")
               ]
    try:
        columns.append(fits.Column('IS_GAMMA', format='1L', array=evt_dict['IS_GAMMA']))
        columns.append(fits.Column('BDT_SCORE', format='1E', array=evt_dict['BDT_SCORE']))
    except KeyError:
        pass

    # Add number of triggered telescope if necessary 
    if save_multiplicity:
        columns.append(fits.Column(name="EVENT_TYPE", format="1J", array=evt_dict['EVENT_TYPE']))

    # Create HDU
    hdu1 = fits.BinTableHDU.from_columns(columns)
    hdu1.name = "EVENTS"

    # Fill Standard HDUCLASS headers
    hdu1 = addHDUClassKeyword(hdu1, class1='EVENTS')

    # Fill Header
    hdu1.header.set('RADECSYS', RADECSYS, 'equatorial system type')
    hdu1.header.set('EQUINOX', EQUINOX, 'base equinox')
    hdu1.header.set('CREATOR', 'pyV2DL3 v{}::{}'.format(VERSION, datasource.get_source_name()))
    hdu1.header.set('ORIGIN', 'VERITAS Collaboration', 'Data from VERITAS')
    hdu1.header.set('TELESCOP', 'VERITAS')
    hdu1.header.set('INSTRUME', 'VERITAS')

    hdu1.header.set('OBS_ID  ', evt_dict['OBS_ID'], 'Run Number')

    hdu1.header.set('DATE-OBS',
                    evt_dict['DATE-OBS'],
                    'start date (UTC) of obs yy-mm-dd hh:mm:ss')
    #include time and date in DATE-OBS
    #hdu1.header.set('TIME-OBS',
    #                evt_dict['TIME-OBS'],
    #                'start time (UTC) of obs hh-mm-ss')
    hdu1.header.set('DATE-END',
                    evt_dict['DATE-END'],
                    'end date (UTC) of obs yy-mm-dd hh:mm:ss')
    #hdu1.header.set('TIME-END',
    #                evt_dict['TIME-END'],
    #                'end time (UTC) of obs hh-mm-ss')
    
    hdu1.header.set('TSTART  ',
                    evt_dict['TSTART'],
                    'mission time of start of obs [s]')
    hdu1.header.set('TSTOP   ',
                    evt_dict['TSTOP'],
                    'mission time of end of obs [s]')
    hdu1.header.set('MJDREFI ',
                    VTS_REFERENCE_MJD, 'int part of reference MJD [days]')
    hdu1.header.set('MJDREFF ', 0., 'fractional part of reference MJD [days]')
    
    hdu1.header.set('TIMEUNIT', 's', 'time unit is seconds since MET start')
    hdu1.header.set('TIMESYS ', 'utc', 'time scale is UTC')
    hdu1.header.set('TIMEREF ', 'local', 'local time reference')
    
    hdu1.header.set('ONTIME  ', 
                    evt_dict['ONTIME'],
                    'sum of good time intervals [s]')

    # Correct live time for time cuts
    hdu1.header.set('LIVETIME', evt_dict['LIVETIME'],
                    '(ontime * deadtime time correction) [s] ')

    hdu1.header.set('DEADC   ', evt_dict['DEADC'],
                    'Average dead time correction (LIVETIME/ONTIME)')
    
    hdu1.header.set('OBJECT  ', evt_dict['OBJECT'], 'observed object')
    hdu1.header.set('RA_OBJ  ',
                    evt_dict['RA_OBJ'],
                    'observation position RA [deg]')
    hdu1.header.set('DEC_OBJ ',
                    evt_dict['DEC_OBJ'],
                    'observation position DEC [deg]')

    hdu1.header.set('RA_PNT  ', evt_dict['RA_PNT'], 'observation position RA [deg]')
    hdu1.header.set('DEC_PNT ', evt_dict['DEC_PNT'], 'observation position DEC [deg]')
    hdu1.header.set('ALT_PNT ', evt_dict['ALT_PNT'], 'average altitude of pointing [deg]')
    hdu1.header.set('AZ_PNT  ', evt_dict['AZ_PNT'], 'average azimuth of pointing [deg]')

    # get the list of telescopes that participate in the event
    hdu1.header.set('TELLIST',
                    evt_dict['TELLIST'],
                    'comma-separated list of tel IDs')
    hdu1.header.set('N_TELS',evt_dict['N_TELS'],
                    'number of telescopes in event list')

    hdu1.header.set('EUNIT   ', 'TeV', 'energy unit')
    hdu1.header.set('GEOLON  ', VTS_REFERENCE_LON, 'longitude of array center [deg]')
    hdu1.header.set('GEOLAT  ', VTS_REFERENCE_LAT, 'latitude of array center [deg]')
    hdu1.header.set('ALTITUDE', VTS_REFERENCE_HEIGHT, 'altitude of array center [m]')

    try:
        hdu1.header.set('QUALITY', evt_dict['QUALITY'], 'Run quality flag based on VPM data used or not')
    except (KeyError):
        print("\033[1;31m The Keyword QUALITY not set in the EVENTS header")
        print(" Please make sure to use ED >= 486 \033[0;0m")

    # Calculate average noise
    return hdu1
