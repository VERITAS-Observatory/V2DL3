from astropy.io import fits
import numpy as np
import logging
from pyV2DL3.addHDUClassKeyword import addHDUClassKeyword 
logger = logging.getLogger(__name__)


def fill_response(datasource):
    response_dict = datasource.get_response_data()
    evt_dict = datasource.get_evt_data()
    irfs_to_store = response_dict['irfs_to_store']

    response_hdus = list()
    if irfs_to_store['point-like']:
        x = response_dict['EA']
        hdu3 = fits.BinTableHDU(data=x)
        hdu3.name = "EFFECTIVE AREA"
        # Fill Standard HDUCLASS keyword
        hdu3 = addHDUClassKeyword(hdu3, class1='RESPONSE', class2='EFF_AREA',
                                  class3='POINT-LIKE', class4='AEFF_2D')
        hdu3.header.set('OBS_ID', evt_dict['OBS_ID'],'Run Number')
        hdu3.header.set('TUNIT1 ', 'TeV', "")
        hdu3.header.set('TUNIT2 ', 'TeV', "")
        hdu3.header.set('TUNIT3 ', 'deg', "")
        hdu3.header.set('TUNIT4 ', 'deg', "")
        hdu3.header.set('TUNIT5 ', 'm2', "")

        hdu3.header.set('LO_THRES', response_dict['LO_THRES'],
                        'Low energy threshold of validity [TeV]')
        hdu3.header.set('HI_THRES', response_dict['HI_THRES'],
                        'High energy threshold of validity [TeV]')
        hdu3.header.set('RAD_MAX ', response_dict['RAD_MAX'], 'Direction cut applied [deg]')

        hdu3.header.set('CREF5', '(ENERG_LO:ENERG_HI,THETA_LO:THETA_HI)','')
        x = response_dict['MIGRATION']
        hdu4 = fits.BinTableHDU(data=x)
        hdu4.name = "ENERGY DISPERSION"
        # Fill Standard HDUCLASS keyword
        hdu4 = addHDUClassKeyword(hdu4, class1='RESPONSE',
                                        class2='EDISP',
                                        class3='POINT-LIKE',
                                        class4='EDISP_2D')

        hdu4.header.set('TUNIT1 ', 'TeV', "")
        hdu4.header.set('TUNIT2 ', 'TeV', "")
        hdu4.header.set('TUNIT5 ', 'deg', "")
        hdu4.header.set('TUNIT6 ', 'deg', "")
        # hdu3.header.set('TUNIT5 ', 'm^2', "")

        hdu4.header.set('LO_THRES', response_dict['LO_THRES'],
                        'Low energy threshold of validity [TeV]')
        hdu4.header.set('HI_THRES', response_dict['HI_THRES'],
                        'High energy threshold of validity [TeV]')
        hdu4.header.set('RAD_MAX ', response_dict['RAD_MAX'], 'Direction cut applied [deg]')
        #Axis order
        hdu4.header.set('CREF7', '(ETRUE_LO:ETRUE_HI,MIGRA_LO:MIGRA_HI,THETA_LO:THETA_HI)', '')
        response_hdus.append(hdu3)
        response_hdus.append(hdu4)
    if irfs_to_store['full-enclosure']:
        # TODO: Implement full-enclosure IRFs:
        print("Full-enclosure IRFs not yet implemented.")
    return response_hdus
