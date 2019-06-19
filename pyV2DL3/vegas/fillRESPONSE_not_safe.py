from ctypes import c_float
import numpy as np
import logging
from root_numpy import hist2array
from pyV2DL3.vegas.util import decodeConfigMask,getThetaSquareCut
from pyV2DL3.vegas.irfloader import  IRFLoader
import ROOT
 
logger = logging.getLogger(__name__)


def __fillRESPONSE_not_safe__(effectiveAreaIO, azimuth, zenith, noise, offset,irf_to_store={}):

    response_dict = {}
    effectiveAreaIO.loadTheRootFile()
    irfloader = IRFLoader(effectiveAreaIO, pointlike=irf_to_store['point-like'])
    ea_final_data, ebias_final_data, abias_final_data = irfloader.getIRF(azimuth, zenith, noise)
    minEnergy,maxEnergy = irfloader.getSafeEnergy(azimuth, zenith, noise)
    response_dict['LO_THRES'] = minEnergy
    response_dict['HI_THRES'] = maxEnergy

    # Point-like
    if irf_to_store['point-like']:
        response_dict['EA'] = ea_final_data
        response_dict['MIGRATION'] = ebias_final_data

        # Load the theta squared cut
        logger.debug('Getting Theta2 cut from EA file')
        cuts = effectiveAreaIO.loadTheCutsInfo()
        for k in cuts:
            theta2cut = getThetaSquareCut(k.fCutsFileText)
        logger.debug('Theta2 cut is {:.2f}'.format(theta2cut))
        response_dict['RAD_MAX'] = np.sqrt(theta2cut)

    # Full-enclosure
    elif irf_to_store['full-enclosure']:
        response_dict['FULL_EA'] = ea_final_data
        response_dict['FULL_MIGRATION'] = ebias_final_data
        response_dict['PSF'] = abias_final_data
    else:
        raise ValueError("IRF requested should be point-like or full-enclosure")

    return response_dict



