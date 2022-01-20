from astropy.io import fits
import logging
from pyV2DL3.addHDUClassKeyword import addHDUClassKeyword
logger = logging.getLogger(__name__)


def fillRESPONSE(datasource):
    response_dict = datasource.get_response_data()
    evt_dict = datasource.get_evt_data()

    response_hdus = list()
    if datasource.__irf_to_store__['point-like']:

        # Effective area (Point-like)
        x = response_dict['EA']
        hdu_ea = fits.BinTableHDU(data=x)
        hdu_ea.name = "EFFECTIVE AREA"
        hdu_ea.header.set('TELESCOP ', 'VERITAS', "")
        hdu_ea.header.set('INSTRUME ', 'Epoch V6', "")

        # Fill Standard HDUCLASS keyword
        hdu_ea = addHDUClassKeyword(hdu_ea, class1='RESPONSE',
                                    class2='EFF_AREA',
                                    class3='POINT-LIKE', class4='AEFF_2D')
        hdu_ea.header.set('OBS_ID', evt_dict['OBS_ID'], 'Run Number')
        hdu_ea.header.set('TUNIT1 ', 'TeV', "")
        hdu_ea.header.set('TUNIT2 ', 'TeV', "")
        hdu_ea.header.set('TUNIT3 ', 'deg', "")
        hdu_ea.header.set('TUNIT4 ', 'deg', "")
        hdu_ea.header.set('TUNIT5 ', 'm2', "")

        hdu_ea.header.set('LO_THRES', response_dict['LO_THRES'],
                          'Low energy threshold of validity [TeV]')
        hdu_ea.header.set('HI_THRES', response_dict['HI_THRES'],
                          'High energy threshold of validity [TeV]')
        hdu_ea.header.set('RAD_MAX ', response_dict['RAD_MAX'],
                          'Direction cut applied [deg]')

        hdu_ea.header.set('CREF5', '(ENERG_LO:ENERG_HI,THETA_LO:THETA_HI)', '')
        # Energy dispersion (Point-like)
        x = response_dict['MIGRATION']
        hdu_edisp = fits.BinTableHDU(data=x)
        hdu_edisp.name = "ENERGY DISPERSION"
        # Fill Standard HDUCLASS keyword
        hdu_edisp = addHDUClassKeyword(hdu_edisp,
                                       class1='RESPONSE',
                                       class2='EDISP',
                                       class3='POINT-LIKE',
                                       class4='EDISP_2D')

        hdu_edisp.header.set('TUNIT1 ', 'TeV', "")
        hdu_edisp.header.set('TUNIT2 ', 'TeV', "")
        hdu_edisp.header.set('TUNIT3 ', '', "")
        hdu_edisp.header.set('TUNIT4 ', '', "")
        hdu_edisp.header.set('TUNIT5 ', 'deg', "")
        hdu_edisp.header.set('TUNIT6 ', 'deg', "")

        hdu_edisp.header.set('RAD_MAX ',
                             response_dict['RAD_MAX'],
                             'Direction cut applied [deg]')
        # Axis order
        hdu_edisp.header.set('CREF7',
                             '(ETRUE_LO:ETRUE_HI,MIGRA_LO:MIGRA_HI,THETA_LO:THETA_HI)',
                             '')
        response_hdus.append(hdu_ea)
        response_hdus.append(hdu_edisp)
    if datasource.__irf_to_store__['full-enclosure']:
        # Effective area (full-enclosure)
        x = response_dict['FULL_EA']
        hdu_fe_ea = fits.BinTableHDU(data=x)
        hdu_fe_ea.name = "EFFECTIVE AREA"
        hdu_fe_ea.header.set('TELESCOP ', 'VERITAS', "")
        hdu_fe_ea.header.set('INSTRUME ', 'Epoch V6', "")

        # Fill Standard HDUCLASS keyword
        hdu_fe_ea = addHDUClassKeyword(hdu_fe_ea,
                                       class1='RESPONSE',
                                       class2='EFF_AREA',
                                       class3='FULL-ENCLOSURE',
                                       class4='AEFF_2D')
        hdu_fe_ea.header.set('OBS_ID', evt_dict['OBS_ID'], 'Run Number')
        hdu_fe_ea.header.set('TUNIT1 ', 'TeV', "")
        hdu_fe_ea.header.set('TUNIT2 ', 'TeV', "")
        hdu_fe_ea.header.set('TUNIT3 ', 'deg', "")
        hdu_fe_ea.header.set('TUNIT4 ', 'deg', "")
        hdu_fe_ea.header.set('TUNIT5 ', 'm2', "")

        hdu_fe_ea.header.set('LO_THRES', response_dict['LO_THRES'],
                             'Low energy threshold of validity [TeV]')
        hdu_fe_ea.header.set('HI_THRES', response_dict['HI_THRES'],
                             'High energy threshold of validity [TeV]')
        hdu_fe_ea.header.set('CREF5',
                             '(ENERG_LO:ENERG_HI,THETA_LO:THETA_HI)',
                             '')
        # Energy dispersion (full-enclosure)
        x = response_dict['FULL_MIGRATION']
        hdu_fe_edisp = fits.BinTableHDU(data=x)
        hdu_fe_edisp.name = "ENERGY DISPERSION"
        # Fill Standard HDUCLASS keyword
        hdu_fe_edisp = addHDUClassKeyword(hdu_fe_edisp,
                                          class1='RESPONSE',
                                          class2='EDISP',
                                          class3='FULL-ENCLOSURE',
                                          class4='EDISP_2D')
        hdu_fe_edisp.header.set('TUNIT1 ', 'TeV', "")
        hdu_fe_edisp.header.set('TUNIT2 ', 'TeV', "")
        hdu_fe_edisp.header.set('TUNIT3 ', '', "")
        hdu_fe_edisp.header.set('TUNIT4 ', '', "")
        hdu_fe_edisp.header.set('TUNIT5 ', 'deg', "")
        hdu_fe_edisp.header.set('TUNIT6 ', 'deg', "")
        hdu_fe_edisp.header.set('CREF7',
                                '(ETRUE_LO:ETRUE_HI,MIGRA_LO:MIGRA_HI,THETA_LO:THETA_HI)',
                                '')
        # Direction dispersion (full-enclosure)
        x = response_dict['PSF']
        hdu_fe_psf = fits.BinTableHDU(data=x)
        hdu_fe_psf.name = "PSF"
        # Fill Standard HDUCLASS keyword
        hdu_fe_psf = addHDUClassKeyword(hdu_fe_psf,
                                        class1='RESPONSE',
                                        class2='PSF',
                                        class3='FULL-ENCLOSURE',
                                        class4='PSF_TABLE')

        hdu_fe_psf.header.set('TUNIT1 ', 'TeV', "")
        hdu_fe_psf.header.set('TUNIT2 ', 'TeV', "")
        hdu_fe_psf.header.set('TUNIT3 ', 'deg', "")
        hdu_fe_psf.header.set('TUNIT4 ', 'deg', "")
        hdu_fe_psf.header.set('TUNIT5 ', 'deg', "")
        hdu_fe_psf.header.set('TUNIT6 ', 'deg', "")
        hdu_fe_psf.header.set('TUNIT7 ', 'sr^-1', "")
        # Axis order
        hdu_fe_psf.header.set('CREF7',
                              '(ENERG_LO:ENERG_HI,THETA_LO:THETA_HI,RAD_LO:RAD_HI)',
                              '')
        response_hdus.append(hdu_fe_ea)
        response_hdus.append(hdu_fe_edisp)
        response_hdus.append(hdu_fe_psf)
    if(not (datasource.__irf_to_store__['point-like']
            or datasource.__irf_to_store__['full-enclosure'])):
        raise Exception('No IRF to store...')
    return response_hdus
