from astropy.io import fits
import logging
from pyV2DL3.load_vegas import VEGASStatus
# Load VEGAS
vegas = VEGASStatus()
vegas.loadVEGAS()
from ROOT import VARootIO
from pyV2DL3.fillGTI import fillGTI
from pyV2DL3.fillRESPONSE import fillRESPONSE
from pyV2DL3.fillEVENTS import fillEVENTS

logger = logging.getLogger(__name__)
def genPrimaryHDU():
    """
    Generate primary hdu
    """
    hdu0 = fits.PrimaryHDU()
    hdu0.header.set('TELESCOP', 
                    'VERITAS' ,
                    'Telescope')
    hdu0.header.set('LICENSE ', 
                            '', 
                    'Copyright (c) 2018,The VERITAS Collaboration')
    hdu0.header['COMMENT'] = "FITS (Flexible Image Transport System) format is defined in 'Astronomy"
    hdu0.header['COMMENT'] = "and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H"
    return hdu0

def loadROOTFiles(st5File,eaFile):
    return VARootIO(st5File, True),VARootIO(eaFile, True)

def genHDUlist(vegasFileIO,effectiveAreaIO):
    hdu0 = genPrimaryHDU() 
    config,hdu1 = fillEVENTS(vegasFileIO)
    hdu2       =  fillGTI(vegasFileIO)
    hdu3,hdu4  = fillRESPONSE(effectiveAreaIO,offset=0.5,**config) 
    hdulist = fits.HDUList([hdu0, hdu1, hdu2, hdu3, hdu4])
    return hdulist
