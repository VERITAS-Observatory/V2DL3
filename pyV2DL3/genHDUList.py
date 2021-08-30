from astropy.io import fits
import logging
from pyV2DL3.fillGTI import fillGTI
from pyV2DL3.fill_response import fill_response
from pyV2DL3.fillEVENTS import fillEVENTS
from pyV2DL3.vegas.VegasDataSource import VegasDataSource
from pyV2DL3.eventdisplay.EventDisplayDataSource import EventDisplayDataSource

logger = logging.getLogger(__name__)


def genPrimaryHDU():
    """
    Generate primary hdu
    """
    hdu0 = fits.PrimaryHDU()
    hdu0.header.set('TELESCOP',
                    'VERITAS',
                    'Telescope')
    hdu0.header.set('LICENSE ',
                    '',
                    'Copyright (c) 2018,The VERITAS Collaboration')
    hdu0.header['COMMENT'] = "FITS (Flexible Image Transport System) format is defined in 'Astronomy"
    hdu0.header['COMMENT'] = "and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H"
    return hdu0


def loadROOTFiles(data_file, effective_area_file, file_type='VEGAS'):
    if file_type == 'VEGAS':
        return VegasDataSource(data_file, effective_area_file)
    if file_type == 'ED':
        return EventDisplayDataSource(data_file, effective_area_file)
    else:
        raise Exception('File type not supported: {}'.format(file_type))


def genHDUlist(datasource, save_multiplicity=False):
    hdus = list()
    hdus.append(genPrimaryHDU())
    hdus.append(fillEVENTS(datasource, save_multiplicity=save_multiplicity))
    hdus.append(fillGTI(datasource))
    hdus.extend(fill_response(datasource))
    hdulist = fits.HDUList(hdus)
    return hdulist
