import logging

import numpy as np

from pyV2DL3.vegas.irfloader import IRFLoader
from pyV2DL3.vegas.util import getThetaSquareCut

logger = logging.getLogger(__name__)


def __fillRESPONSE_not_safe__(
    effectiveAreaIO, azimuth, zenith, noise, irf_to_store=None
):

    if irf_to_store is None:
        irf_to_store = {}

    response_dict = {}
    effectiveAreaIO.loadTheRootFile()
    irfloader = IRFLoader(effectiveAreaIO, pointlike=irf_to_store["point-like"])
    ea_final_data, ebias_final_data, abias_final_data = irfloader.getIRF(
        azimuth, zenith, noise
    )
    minEnergy, maxEnergy = irfloader.getSafeEnergy(azimuth, zenith, noise)
    response_dict["LO_THRES"] = minEnergy
    response_dict["HI_THRES"] = maxEnergy

    # Point-like
    if irf_to_store["point-like"]:
        response_dict["EA"] = ea_final_data
        response_dict["MIGRATION"] = ebias_final_data

        # Load the theta squared cut
        logger.debug("Getting Theta2 cut from EA file")
        cuts = effectiveAreaIO.loadTheCutsInfo()
        for k in cuts:
            theta2cut = getThetaSquareCut(k.fCutsFileText)
        logger.debug(f"Theta2 cut is {theta2cut:.2f}")
        response_dict["RAD_MAX"] = np.sqrt(theta2cut)

    # Full-enclosure
    elif irf_to_store["full-enclosure"]:
        response_dict["FULL_EA"] = ea_final_data
        response_dict["FULL_MIGRATION"] = ebias_final_data
        response_dict["PSF"] = abias_final_data
    else:
        raise ValueError("IRF requested should be point-like or full-enclosure")

    return response_dict
