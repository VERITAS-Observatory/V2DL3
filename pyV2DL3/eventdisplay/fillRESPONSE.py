import click
import logging

import numpy as np
import uproot

from pyV2DL3.eventdisplay.IrfInterpolator import IrfInterpolator
from pyV2DL3.eventdisplay.util import bin_centers_to_edges

logger = logging.getLogger(__name__)


def find_energy_range(log_energy_tev):
    """Find min and max of energy axis"""
    energy_low = np.power(
        10, log_energy_tev - (log_energy_tev[1] - log_energy_tev[0]) / 2.0
    )
    energy_high = np.power(
        10, log_energy_tev + (log_energy_tev[1] - log_energy_tev[0]) / 2.0
    )
    return energy_low, energy_high


def print_logging_info(irf_to_store, camera_offsets, pedvar, zenith):
    """Print information of parameter space to access"""

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


def check_parameter_range(par, irf_stored_par, par_name, **kwargs):
    """Check that coordinates are in range of provided IRF and whether extrapolation is to be done
       0. checks if command line parameter force_extrapolation is given. If given,
          the extrapolation will happen when parameter is outside IRF range. If parameter is
          within IRF range, it works as normal. Default is False.
       1. Further checks for fuzzy boundary (parameter close to boundary value).
          If fuzzy boundary is within a given tolerance then IRF is interpolated for
          at boundary value. Default is 0.0 tolerance.
    """

    logging.info(
        "\t{0} range of a given IRF: {1:.2f} - {2:.2f}".format(
            par_name, np.min(irf_stored_par), np.max(irf_stored_par)
        )
    )

    if kwargs.get("use_click", True):
        clk = click.get_current_context()
        tolerance = clk.params["fuzzy_boundary"]
        extrapolation = clk.params["force_extrapolation"]
    else:
        tolerance = kwargs.get("fuzzy_boundary", 0.0)
        extrapolation = kwargs.get("force_extrapolation", False)

    if np.all(irf_stored_par < par) or np.all(irf_stored_par > par):
        if extrapolation:
            logging.warning("IRF extrapolation allowed for coordinate not inside IRF {0} range".format(par_name))
        elif tolerance > 0.0:
            if check_fuzzy_boundary(par, np.max(irf_stored_par), tolerance):
                par = np.max(irf_stored_par)
            elif check_fuzzy_boundary(par, np.min(irf_stored_par), tolerance):
                par = np.min(irf_stored_par)
            else:
                raise ValueError("Tolerance not calculated for coordinate {0}".format(par_name))
        else:
            raise ValueError(
                "Coordinate not inside IRF {0} range! Try using --fuzzy_boundary".format(par_name)
            )
    return par


def check_fuzzy_boundary(par, boundary, tolerance):
    """" Checks if the parameter value is within the given tolerance.
    tolerance parameter is defined as ratio of absolute difference
    between boundary and par to the boundary.

    Parameters
    ----------
    par: parameter of given run, it can be pedvar, zenith or camera offset
    boundary: lower or upper boundary value of stored IRF
    tolerance: allowed value of --fuzzy_boundary command line argument

    Returns
    -------
    Boolean. Default is False. True if tolerance is within given allowed value.
    If boundary zero then also returns default False.

    """
    if boundary == 0:
        return False

    if boundary > 0:
        fuzzy_diff = np.abs(boundary - par) / boundary
        if fuzzy_diff < tolerance:
            logging.warning(
                "Coordinate tolerance is {0:0.3f} and is within {1:0.3f}".format(
                    fuzzy_diff, tolerance
                )
            )
            return True
        else:
            raise ValueError(
                "Coordinate tolerance is {0:0.3f} and is outside {1:0.3f}".format(
                    fuzzy_diff, tolerance
                )
            )

    return False


def find_camera_offsets(camera_offsets):
    """Find camera offsets, depending on  availability in the effective area file."""

    if len(camera_offsets) == 1:
        # Many times, just IRFs for 0.5 deg are available.
        # Assume that offset for the whole camera.
        logger.debug(
            "IMPORTANT: Only one camera offset bin "
            + "({} deg) simulated within the effective area file selected.".format(
                camera_offsets[0]
            )
        )
        logger.debug(
            "IMPORTANT: Setting the IRFs of that given camera offset value to the whole camera"
        )
        return [0.0, 10.0], [0.0, 10.0]

    # Note in the camera offset _low and _high may refer
    # to the simulated "points", and
    # not to actual bins.
    return camera_offsets, camera_offsets


def fill_effective_area(
        irf_name,
        irf_interpolator,
        camera_offsets,
        pedvar,
        zenith,
        theta_low,
        theta_high,
        **kwargs
):
    """Effective areas"""

    irf_interpolator.set_irf(irf_name, **kwargs)

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
        camera_offsets,
        pedvar,
        zenith,
        theta_low,
        theta_high,
        **kwargs
):
    """Energy migration matrix"""

    irf_interpolator.set_irf(irf_name, **kwargs)
    ac_final = []

    for offset in camera_offsets:
        bias, axis = irf_interpolator.interpolate([pedvar, zenith, offset])

        _, e_low, e_high = bin_centers_to_edges(axis[0])
        _, b_low, b_high = bin_centers_to_edges(axis[1], False)

        ac = []

        for aa in bias.transpose():
            ab = aa / np.sum(aa * (b_high - b_low)) if np.sum(aa) > 0 else aa
            try:
                ac = np.vstack((ac, ab))
            except ValueError:
                ac = ab

        ac = ac.transpose()
        ac_final.append(ac)

    return np.array(
        [(e_low, e_high, b_low, b_high, theta_low, theta_high, ac_final)],
        dtype=[
            ("ENERG_LO", ">f4", (len(e_low),)),
            ("ENERG_HI", ">f4", (len(e_high),)),
            ("MIGRA_LO", ">f4", (len(b_low),)),
            ("MIGRA_HI", ">f4", (len(b_high),)),
            ("THETA_LO", ">f4", (len(theta_low),)),
            ("THETA_HI", ">f4", (len(theta_high),)),
            ("MATRIX", ">f4", (np.shape(ac_final))),
        ],
    )


def fill_direction_migration(
        irf_interpolator, camera_offsets, pedvar, zenith, theta_low, theta_high, **kwargs
):
    """Direction dispersion (for full-enclosure IRFs)"""

    irf_interpolator.set_irf("hAngularLogDiffEmc_2D", **kwargs)

    rpsf_final = []
    rpsf_test = []
    test_psf = False  # use PSF distribution from IRFs by default

    for offset in camera_offsets:

        # direction diff (rad, energy),
        direction_diff, axis = irf_interpolator.interpolate([pedvar, zenith, offset])

        # energy axis from ~ 0.1 - 100 TeV
        energy_axis_index_lb = np.searchsorted(np.power(10, axis[0]), 0.1)
        energy_axis_index_ub = np.searchsorted(np.power(10, axis[0]), 100) - len(axis[0])

        axis[0] = axis[0][energy_axis_index_lb:energy_axis_index_ub]
        _, e_low, e_high = bin_centers_to_edges(axis[0])

        # generate psf data from halfnorm pdf
        if test_psf:
            from scipy.stats import halfnorm

            # interpolation test
            rad_edges, r_low, r_high = bin_centers_to_edges(np.linspace(0, 10, 4000), logaxis=False)

            # use linspace instead of rad_edges
            rad_width_deg = np.diff(rad_edges)

            x = np.linspace(0, 10, 4000)
            sigma = 0.5
            scale = sigma * np.sqrt(1 - 2 / np.pi)
            y = halfnorm.pdf(x, loc=0, scale=scale)
            cumsum = (2 * np.pi * rad_width_deg * y * (r_low + r_high) / 2).cumsum()
            normed = y / cumsum.max() * ((180 / np.pi) ** 2)
            normed = np.nan_to_num(normed)

            # PSF should be normed (deg**2 / sr), test should give 3200:
            # values = 2 * np.pi * rad_width_deg * normed * (r_low + r_high) / 2
            # print("PSF normed? ( ≈ 3282 (deg**2 / sr))", values.cumsum().max())

            y = np.array(normed)
            test = np.repeat(y[np.newaxis, ...], len(axis[0]), axis=0)
            rpsf_test.append(test)

        else:
            direction_diff = direction_diff[:, energy_axis_index_lb:energy_axis_index_ub]

            # Using rad**2 bins to normalize, dN/dlog(rad) ~ rad*dN/d(rad)
            rad_edges, r_low, r_high = bin_centers_to_edges(axis[1], logaxis=True)

            rad_width_deg = np.diff(np.power(10, rad_edges))
            # this step makes sure all arrays have the same dimensions, rad_width_deg and the central rad values are
            # repeated by the length of the energy axis.
            norm = np.sum(direction_diff * np.repeat(rad_width_deg[..., np.newaxis], len(axis[0]), axis=1)
                          / np.repeat(((r_low + r_high) / 2)[..., np.newaxis], len(axis[0]), axis=1), axis=0)
            norm = norm * 2 * np.pi
            direction_diff = direction_diff / (
                np.repeat(((r_low + r_high) / 2)[..., np.newaxis], len(axis[0]), axis=1) ** 2)
            normed = direction_diff / norm * ((180 / np.pi) ** 2)
            rpsf_final.append(np.nan_to_num(normed))

    # PSF (3-dim with axes: psf[rad_index, offset_index, energy_index]
    if test_psf:
        rpsf_test = np.swapaxes(rpsf_test, 0, 1)
        rpsf_final = np.swapaxes(rpsf_test, 0, 2)
    else:
        rpsf_final = np.swapaxes(rpsf_final, 0, 1)

    return np.array(
        [(e_low, e_high, theta_low, theta_high, r_low, r_high, rpsf_final)],
        dtype=[
            ("ENERG_LO", ">f4", (np.shape(e_low))),
            ("ENERG_HI", ">f4", (np.shape(e_high))),
            ("THETA_LO", ">f4", (np.shape(theta_low))),
            ("THETA_HI", ">f4", (np.shape(theta_high))),
            ("RAD_LO", ">f4", (np.shape(r_low))),
            ("RAD_HI", ">f4", (np.shape(r_high))),
            ("RPSF", ">f4", (np.shape(rpsf_final))),
        ],
    )


def __fill_response__(
        ed_file_io, effective_area, azimuth, zenith, pedvar, irf_to_store=None, **kwargs
):
    if irf_to_store is None:
        irf_to_store = {}

    response_dict = {}

    # IRF interpolator
    irf_interpolator = IrfInterpolator(effective_area, azimuth)

    # Extract camera offsets available from the effective areas file.
    fast_eff_area = uproot.open(effective_area)["fEffAreaH2F"]
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

    zenith = check_parameter_range(zenith, zeniths_irf, "zenith", **kwargs)
    pedvar = check_parameter_range(pedvar, pedvar_irf, "pedvar", **kwargs)
    theta_low, theta_high = find_camera_offsets(camera_offsets)

    if irf_to_store["point-like"]:

        # Effective area (full-enclosure)
        (
            response_dict["EA"],
            response_dict["LO_THRES"],
            response_dict["HI_THRES"],
        ) = fill_effective_area(
            "eff",
            irf_interpolator,
            camera_offsets,
            pedvar,
            zenith,
            theta_low,
            theta_high,
            **kwargs
        )

        # Get RAD_MAX; cuts don't depend on energy/wobble
        file = uproot.open(ed_file_io)
        run_summary = file["total_1/stereo/tRunSummary"].arrays(library="np")
        theta2cut = run_summary["Theta2Max"][0]
        response_dict["RAD_MAX"] = np.sqrt(theta2cut)

        # Energy dispersion (point-like)
        response_dict["MIGRATION"] = fill_energy_migration(
            "hEsysMCRelative2D",
            irf_interpolator,
            camera_offsets,
            pedvar,
            zenith,
            theta_low,
            theta_high,
            **kwargs
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
        (
            response_dict["FULL_EA"],
            response_dict["LO_THRES"],
            response_dict["HI_THRES"],
        ) = fill_effective_area(
            "effNoTh2",
            irf_interpolator,
            camera_offsets,
            pedvar,
            zenith,
            theta_low,
            theta_high,
            **kwargs
        )

        # Energy dispersion (full-enclosure)
        response_dict["FULL_MIGRATION"] = fill_energy_migration(
            "hEsysMCRelative2DNoDirectionCut",
            irf_interpolator,
            camera_offsets,
            pedvar,
            zenith,
            theta_low,
            theta_high,
            **kwargs
        )

        # Direction dispersion (for full-enclosure IRFs)
        response_dict["PSF"] = fill_direction_migration(
            irf_interpolator,
            camera_offsets,
            pedvar,
            zenith,
            theta_low,
            theta_high,
            **kwargs
        )

    return response_dict
