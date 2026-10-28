import logging

import numpy as np

from pyV2DL3.vegas.irfloader import getIRF

logger = logging.getLogger(__name__)


def __fillRESPONSE_not_safe__(
    effective_area_file,
    azimuth,
    zenith,
    noise,
    irf_to_store,
    psf_king_params=None,
    st6_configs=None,
):
    response_dict = {}
    ea_final_data, ebias_final_data, abias_final_data = getIRF(
        azimuth,
        zenith,
        noise,
        effective_area_file,
        irf_to_store["point-like"],
        psf_king_params=psf_king_params,
    )

    manager = effective_area_file.manager

    lo_thres = 0.0
    hi_thres = np.inf

    for data in (ea_final_data, ebias_final_data, abias_final_data):
        lo_thres = max(lo_thres, np.min(data["ENERG_LO"]))
        hi_thres = min(hi_thres, np.max(data["ENERG_HI"]))

    if lo_thres >= hi_thres:
        raise ValueError(
            f"No overlapping energy range across IRFs: LO_THRES={lo_thres}, HI_THRES={hi_thres}"
        )

    response_dict["LO_THRES"] = lo_thres
    response_dict["HI_THRES"] = hi_thres

    logger.debug(f"RAD_MAX = {np.sqrt(effective_area_file.theta_square_upper)}")
    # Point-like
    if irf_to_store["point-like"]:
        response_dict["EA"] = ea_final_data
        response_dict["MIGRATION"] = ebias_final_data
        response_dict["RAD_MAX"] = np.sqrt(effective_area_file.theta_square_upper)

    # Full-enclosure
    elif irf_to_store["full-enclosure"]:
        response_dict["FULL_EA"] = ea_final_data
        response_dict["FULL_MIGRATION"] = ebias_final_data
        response_dict["PSF"] = abias_final_data
    else:
        raise ValueError("IRF requested should be point-like or full-enclosure")

    return response_dict
