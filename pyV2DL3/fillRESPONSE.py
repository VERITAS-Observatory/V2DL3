from astropy.io import fits
import numpy as np
import logging
 
logger = logging.getLogger(__name__)

def fillRESPONSE(datasource):
    response_dict = datasource.get_response_data()
    x = response_dict['EA']    
    hdu3 = fits.BinTableHDU(data=x)
    hdu3.name = "EFFECTIVE AREA"
    
    hdu3.header.set('TUNIT1 ', 'TeV', "")
    hdu3.header.set('TUNIT2 ', 'TeV', "")
    hdu3.header.set('TUNIT3 ', 'deg', "")
    hdu3.header.set('TUNIT4 ', 'deg', "")
    hdu3.header.set('TUNIT5 ', 'm^2', "")
    
    hdu3.header.set('HDUCLASS', 'GADF',
                    'FITS file following the GADF data format.')
    hdu3.header.set('HDUCLAS1', 'RESPONSE', 'HDU class')
    hdu3.header.set('HDUCLAS2', 'EFF_AREA', 'HDU class')
    hdu3.header.set('HDUCLAS3', 'POINT-LIKE', 'HDU class')
    hdu3.header.set('HDUCLAS4', 'AEFF_2D', 'HDU class')
    hdu3.header.set('LO_THRES', response_dict['LO_THRES'],
                    'Low energy threshold of validity [TeV]')
    hdu3.header.set('HI_THRES', response_dict['HI_THRES'],
                    'High energy threshold of validity [TeV]')
    hdu3.header.set('RAD_MAX ', response_dict['RAD_MAX'], 'Direction cut applied [deg]')

    hdu3.header.set('CREF5','(ENERG_LO:ENERG_HI,THETA_LO:THETA_HI)','') 
    x = response_dict['MIGRATION']    
    hdu4 = fits.BinTableHDU(data=x)
    hdu4.name = "ENERGY DISPERSION"
    
    hdu4.header.set('TUNIT1 ', 'TeV', "")
    hdu4.header.set('TUNIT2 ', 'TeV', "")
    hdu4.header.set('TUNIT5 ', 'deg', "")
    hdu4.header.set('TUNIT6 ', 'deg', "")
    # hdu3.header.set('TUNIT5 ', 'm^2', "")
    
    hdu4.header.set('HDUCLASS', 'GADF',
                    'FITS file following the GADF data format.')
    hdu4.header.set('HDUCLAS1', 'RESPONSE', 'HDU class')
    hdu4.header.set('HDUCLAS2', 'EDISP', 'HDU class')
    hdu4.header.set('HDUCLAS3', 'POINT-LIKE', 'HDU class')
    hdu4.header.set('HDUCLAS4', 'EDISP_2D', 'HDU class')
    hdu4.header.set('LO_THRES', response_dict['LO_THRES'],
                    'Low energy threshold of validity [TeV]')
    hdu4.header.set('HI_THRES', response_dict['HI_THRES'],
                    'High energy threshold of validity [TeV]')
    hdu4.header.set('RAD_MAX ', response_dict['RAD_MAX'], 'Direction cut applied [deg]')
    #Axis order 
    hdu4.header.set('CREF7','(ETRUE_LO:ETRUE_HI,MIGRA_LO:MIGRA_HI,THETA_LO:THETA_HI)','') 
    return hdu3,hdu4
