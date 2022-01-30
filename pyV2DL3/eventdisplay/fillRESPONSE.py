import logging

import numpy as np
import uproot

from pyV2DL3.eventdisplay.IrfInterpolator import IrfInterpolator
from pyV2DL3.eventdisplay.util import bin_centers_to_edges

logger = logging.getLogger(__name__)


def find_energy_range(log_energy_tev):
    """Find min and max of energy axis

    """
    energy_low = np.power(
        10, log_energy_tev - (log_energy_tev[1] - log_energy_tev[0]) / 2.0
    )
    energy_high = np.power(
        10, log_energy_tev + (log_energy_tev[1] - log_energy_tev[0]) / 2.0
    )
    return energy_low, energy_high


def print_logging_info(irf_to_store, camera_offsets, pedvar, zenith):
    """Print information of parameter space to access

    """

    str_info = "Extracting "
    if irf_to_store["point-like"]:
        str_info += "Point-like IRFs "
    else:
        str_info += "Full-enclosure IRFs "
    str_info += "for zenith: {1:.1f} deg, pedvar: {0:.1f}".format(pedvar, zenith)
    logging.info(str_info)
    str_woff = "\tcamera offset: "
    for w in camera_offsets.tolist():
        str_woff += " {0:.2f},".format(w)
    logging.info(str_woff)


def check_parameter_range(par, par_irf, par_name):
    """Check that coordinates are in range of provided IRF

    """

    logging.info("\t{0} range of a given IRF: {1:.1f} - {2:.1f}"
                 .format(par_name, np.min(par_irf), np.max(par_irf)))
    if np.all(par_irf < par) or np.all(par_irf > par):
        raise ValueError("Coordinate not inside IRF {0} range"
                         .format(par_name))


def find_camera_offsets(camera_offsets):
    """Find camera offsets, depending on  availability in the effective area file.

    """

    if len(camera_offsets) == 1:
        # Many times, just IRFs for 0.5 deg are available.
        # Assume that offset for the whole camera.
        logger.debug(
            "IMPORTANT: Only one camera offset bin " +
            + "({} deg) simulated within the effective area file selected."
            .format(camera_offsets[0])
        )
        logger.debug(
            "IMPORTANT: Setting the IRFs of that given camera \
                     offset value to the whole camera"
        )
        return [0.0, 10.0], [0.0, 10.0]

    # Note in the camera offset _low and _high may refer
    # to the simulated "points", and
    # not to actual bins.
    return camera_offsets, camera_offsets


def fill_effective_area(
                     irf_name,
                     irf_interpolator,
                     camera_offsets, pedvar, zenith, offset,
                     theta_low, theta_high):
    """Effective areas

    """

    irf_interpolator.set_irf(irf_name)

    ea_final = []

    # Loop over offsets and store
    for offset in camera_offsets:
        eff_area, axis = irf_interpolator.interpolate([pedvar, zenith, offset])
        ea_final.append(np.array(eff_area))

    # Always same axis values in loop, therefore calculate afterwards
    energy_low, energy_high = find_energy_range(axis[0])

    x = np.array(
        [(energy_low, energy_high, theta_low, theta_high, ea_final)],
        dtype=[
            ("ENERG_LO", ">f4", np.shape(energy_low)),
            ("ENERG_HI", ">f4", np.shape(energy_high)),
            ("THETA_LO", ">f4", np.shape(theta_low)),
            ("THETA_HI", ">f4", np.shape(theta_high)),
            ("EFFAREA", ">f4", np.shape(ea_final)),
        ],
    )

    return x, min(energy_low), max(energy_high)


def fill_energy_migration(
                     irf_name,
                     irf_interpolator,
                     camera_offsets, pedvar, zenith, offset,
                     theta_low, theta_high):
    """Energy migration matrix

    """

    irf_interpolator.set_irf(irf_name)
    ac_final = []

    for offset in camera_offsets:
        bias, axis = irf_interpolator.interpolate([pedvar, zenith, offset])

        _, eLow, eHigh = bin_centers_to_edges(axis[0])
        _, bLow, bHigh = bin_centers_to_edges(axis[1], False)

        ac = []

        for aa in bias.transpose():
            ab = aa / np.sum(aa * (bHigh - bLow)) if np.sum(aa) > 0 else aa
            try:
                ac = np.vstack((ac, ab))
            except ValueError:
                ac = ab

        ac = ac.transpose()
        ac_final.append(ac)

    return np.array([
        (eLow, eHigh, bLow, bHigh, theta_low, theta_high, ac_final)],
        dtype=[
            ("ENERG_LO", ">f4", (len(eLow),)),
            ("ENERG_HI", ">f4", (len(eHigh),)),
            ("MIGRA_LO", ">f4", (len(bLow),)),
            ("MIGRA_HI", ">f4", (len(bHigh),)),
            ("THETA_LO", ">f4", (len(theta_low),)),
            ("THETA_HI", ">f4", (len(theta_high),)),
            ("MATRIX", ">f4", (np.shape(ac_final))),
        ],
    )


def fill_direction_migration(
                     irf_interpolator,
                     camera_offsets, pedvar, zenith, offset,
                     theta_low, theta_high):
    """Direction dispersion (for full-enclosure IRFs)

    """

    irf_interpolator.set_irf("hAngularLogDiffEmc_2D")

    rpsf_final = []

    for offset in camera_offsets:

        direction_diff, axis = irf_interpolator.interpolate([pedvar,
                                                             zenith,
                                                             offset])

        _, eLow, eHigh = bin_centers_to_edges(axis[0])
        rad_edges, rLow, rHigh = bin_centers_to_edges(axis[1])

        # Normalize rpsf by solid angle
        rad_width_deg = np.diff(np.power(10, rad_edges))
        e_sum = np.sum(direction_diff * 2 * rad_width_deg[:, np.newaxis]
                       * np.pi * np.power(10, axis[1])[:, np.newaxis], axis=0)
        normsum = np.divide(((180 / np.pi) ** 2), e_sum, where=e_sum != 0)

        rpsf = direction_diff * normsum
        rpsf_final.append(np.nan_to_num(rpsf).cumsum(axis=0))

    # PSF (3-dim with axes: psf[rad_index, offset_index, energy_index]
    rpsf_final = np.swapaxes(rpsf_final, 0, 1)

    return np.array(
        [(eLow, eHigh, theta_low, theta_high, rLow, rHigh, rpsf_final)],
        dtype=[
            ("ENERG_LO", ">f4", (np.shape(eLow))),
            ("ENERG_HI", ">f4", (np.shape(eHigh))),
            ("THETA_LO", ">f4", (np.shape(theta_low))),
            ("THETA_HI", ">f4", (np.shape(theta_high))),
            ("RAD_LO", ">f4", (np.shape(rLow))),
            ("RAD_HI", ">f4", (np.shape(rHigh))),
            ("RPSF", ">f4", (np.shape(rpsf_final))),
        ],
    )


def __fillRESPONSE__(
    edFileIO, effectiveArea,
    azimuth, zenith, pedvar, offset,
    irf_to_store=None
):
    if irf_to_store is None:
        irf_to_store = {}
    
    response_dict = {}

    # IRF interpolator
    irf_interpolator = IrfInterpolator(effectiveArea, azimuth)

    # Extract camera offsets available from the effective areas file.
    fast_eff_area = uproot.open(effectiveArea)["fEffArea"]
    camera_offsets = np.unique(
        np.round(fast_eff_area["Woff"].array(library="np"), decimals=2)
    )
    zeniths_irf = np.unique(
        np.round(fast_eff_area["ze"].array(library="np"), decimals=0)
    )
    pedvar_irf = np.unique(
        np.round(fast_eff_area["pedvar"].array(library="np"), decimals=2)
    )

    print_logging_info(irf_to_store, camera_offsets, pedvar, zenith)

    check_parameter_range(zenith, zeniths_irf, 'zenith')
    check_parameter_range(pedvar, pedvar_irf, 'pedvar')
    theta_low, theta_high = find_camera_offsets(camera_offsets)

    if irf_to_store["point-like"]:

        # Effective area (full-enclosure)
        response_dict["EA"], response_dict["LO_THRES"], response_dict["HI_THRES"] = fill_effective_area(
                                             "eff",
                                             irf_interpolator,
                                             camera_offsets,
                                             pedvar,
                                             zenith,
                                             offset,
                                             theta_low, theta_high
                                      )

        # Get RAD_MAX; cuts don't depend on energy/wobble
        file = uproot.open(edFileIO)
        runSummary = file["total_1/stereo/tRunSummary"].arrays(library="np")
        theta2cut = runSummary["Theta2Max"][0]
        response_dict["RAD_MAX"] = np.sqrt(theta2cut)

        # Energy dispersion (point-like)
        response_dict["MIGRATION"] = fill_energy_migration(
                                             "hEsysMCRelative2D",
                                             irf_interpolator,
                                             camera_offsets,
                                             pedvar,
                                             zenith,
                                             offset,
                                             theta_low, theta_high
                                      )

    elif irf_to_store["full-enclosure"]:

        # require multiple offsets for full enclosure
        if len(camera_offsets) <= 1:
            logger.warning(
                "IRF used for interpolation should be "
                "defined for several offsets for"
                "Full-Enclosure conversion"
            )

        # Effective area (full-enclosure)
        response_dict["FULL_EA"], response_dict["LO_THRES"], response_dict["HI_THRES"] = fill_effective_area(
                                             "effNoTh2",
                                             irf_interpolator,
                                             camera_offsets,
                                             pedvar,
                                             zenith,
                                             offset,
                                             theta_low, theta_high
                                      )

        # Energy dispersion (full-enclosure)
        response_dict["FULL_MIGRATION"] = fill_energy_migration(
                                             "hEsysMCRelative2DNoDirectionCut",
                                             irf_interpolator,
                                             camera_offsets,
                                             pedvar,
                                             zenith,
                                             offset,
                                             theta_low, theta_high
                                      )

        # Direction dispersion (for full-enclosure IRFs)
        response_dict["PSF"] = fill_direction_migration(
                                        irf_interpolator,
                                        camera_offsets,
                                        pedvar,
                                        zenith,
                                        offset,
                                        theta_low, theta_high
                                     )

    return response_dict
