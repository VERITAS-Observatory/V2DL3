import logging
import numpy as np
from pyV2DL3.eventdisplay.IrfInterpolator import IrfInterpolator
from pyV2DL3.eventdisplay.util import bin_centers_to_edges
import uproot

logger = logging.getLogger(__name__)


def __fillRESPONSE__(
    edFileIO, effectiveArea, azimuth, zenith, pedvar, offset, irf_to_store={}
):
    response_dict = {}

    # EventDisplay IRF interpolator object
    irf_interpolator = IrfInterpolator(effectiveArea, azimuth)

    # Extract the camera offsets simulated within the effective areas file.
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
    # check that coordinates are in range of provided IRF

    print("\tzenith range of a given IRF:", np.min(zeniths_irf), "-", np.max(zeniths_irf))
    print("\tpedvar range of a given IRF:", np.min(pedvar_irf), "-", np.max(pedvar_irf))

    if np.all(zeniths_irf < zenith) or np.all(zeniths_irf > zenith):
        raise ValueError("Coordinate not inside IRF zenith range")
    if np.all(pedvar_irf < pedvar) or np.all(pedvar_irf > pedvar):
        raise ValueError("Coordinate not inside IRF pedvar range")

    # Check the camera offset bins available in the effective area file.
    theta_low = []
    theta_high = []
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
            "IMPORTANT: Setting the IRFs of that given camera \
                     offset value to the whole camera"
        )
        theta_low = [0.0, 10.0]
        theta_high = [0.0, 10.0]
    if len(camera_offsets) > 1:
        # Note in the camera offset _low and _high may refer
        # to the simulated "points", and
        # not to actual bins.
        theta_low = camera_offsets
        theta_high = camera_offsets

    if irf_to_store["point-like"]:
        print("Point-like IRF: ")
        print("\tcamera offset: ", camera_offsets, " deg")
        print("\tpedvar: %.1f" % pedvar, end=" ")
        print(", zenith: %.1f deg" % zenith)
        #
        # Interpolate effective area  (point-like)
        #
        irf_interpolator.set_irf("eff")

        ea_final = []

        # Loop over offsets
        for offset in camera_offsets:
            eff_area, axis = irf_interpolator.interpolate([pedvar, zenith, offset])
            ea_final.append(np.array(eff_area))

        # Always same axis in loop
        log_energy_tev = axis[0]
        energy_low = np.power(
            10, log_energy_tev - (log_energy_tev[1] - log_energy_tev[0]) / 2.0
        )
        energy_high = np.power(
            10, log_energy_tev + (log_energy_tev[1] - log_energy_tev[0]) / 2.0
        )

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
        response_dict["EA"] = x
        response_dict["LO_THRES"] = min(energy_low)
        response_dict["HI_THRES"] = max(energy_high)

        # Get RAD_MAX; cuts don't depend on energy/wobble
        file = uproot.open(edFileIO)
        runSummary = file["total_1/stereo/tRunSummary"].arrays(library="np")
        theta2cut = runSummary["Theta2Max"][0]
        response_dict["RAD_MAX"] = np.sqrt(theta2cut)
        #
        # Energy dispersion (point-like)
        #
        irf_interpolator.set_irf("hEsysMCRelative2D")
        ac_final = []
        for offset in camera_offsets:
            bias, axis = irf_interpolator.interpolate([pedvar, zenith, offset])

            energy_edges = bin_centers_to_edges(axis[0])
            bias_edges = bin_centers_to_edges(axis[1])

            eLow = np.power(10, [energy_edges[:-1]])[0]
            eHigh = np.power(10, [energy_edges[1:]])[0]

            bLow = np.array([bias_edges[:-1]])[0]
            bHigh = np.array([bias_edges[1:]])[0]

            ac = []
            for aa in bias.transpose():
                if np.sum(aa) > 0:
                    ab = aa / np.sum(aa * (bHigh - bLow))
                else:
                    ab = aa
                try:
                    ac = np.vstack((ac, ab))
                except Exception:
                    ac = ab
            ac = ac.transpose()
            ac_final.append(ac)

        x = np.array(
            [(eLow, eHigh, bLow, bHigh, theta_low, theta_high, ac_final)],
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
        response_dict["MIGRATION"] = x
        response_dict["RAD_MAX"] = np.sqrt(theta2cut)
        print("IRF interpolation done")

    if irf_to_store["full-enclosure"]:
        print("Full-enclosure: ")
        print("\tcamera offset: ", camera_offsets)
        print("\tpedvar: %.1f" % pedvar, end=" ")
        print(", zenith: %.1f deg" % zenith)

        # check if IRF contains multiple offsets:
        if len(camera_offsets) <= 1:
            logger.warning(
                "IRF used for interpolation should be "
                "defined for several offsets for"
                "Full-Enclosure conversion"
            )
        #
        # Interpolate effective area (full-enclosure)
        #
        irf_interpolator.set_irf("effNoTh2")

        ea_final = []

        # Loop over offsets and store
        for offset in camera_offsets:
            eff_area, axis = irf_interpolator.interpolate([pedvar, zenith, offset])

            y = np.array(eff_area)
            ea = y  # [y, y]
            ea_final.append(ea)

        # Always same axis values in loop, therefore calculate afterwards
        log_energy_tev = axis[0]
        energy_low = np.power(
            10, log_energy_tev - (log_energy_tev[1] - log_energy_tev[0]) / 2.0
        )
        energy_high = np.power(
            10, log_energy_tev + (log_energy_tev[1] - log_energy_tev[0]) / 2.0
        )

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
        response_dict["LO_THRES"] = min(energy_low)
        response_dict["HI_THRES"] = max(energy_high)
        response_dict["FULL_EA"] = x
        #
        # Energy dispersion (full-enclosure)
        #
        irf_interpolator.set_irf("hEsysMCRelative2DNoDirectionCut")
        ac_final = []

        for offset in camera_offsets:
            bias, axis = irf_interpolator.interpolate([pedvar, zenith, offset])

            energy_edges = bin_centers_to_edges(axis[0])
            bias_edges = bin_centers_to_edges(axis[1])

            eLow = np.power(10, [energy_edges[:-1]])[0]
            eHigh = np.power(10, [energy_edges[1:]])[0]

            bLow = np.array([bias_edges[:-1]])[0]
            bHigh = np.array([bias_edges[1:]])[0]

            ac = []

            for aa in bias.transpose():
                if np.sum(aa) > 0:
                    ab = aa / np.sum(aa * (bHigh - bLow))
                else:
                    ab = aa
                try:
                    ac = np.vstack((ac, ab))
                except Exception:
                    ac = ab

            ac = ac.transpose()
            ac_final.append(ac)

        x = np.array(
            [(eLow, eHigh, bLow, bHigh, theta_low, theta_high, ac_final)],
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

        response_dict["FULL_MIGRATION"] = x
        #
        # Direction dispersion (for full-enclosure IRFs)
        #
        irf_interpolator.set_irf("hAngularLogDiffEmc_2D")

        rpsf_final = []

        for offset in camera_offsets:

            direction_diff, axis = irf_interpolator.interpolate([pedvar, zenith, offset])
                                                                
            energy_edges = bin_centers_to_edges(axis[0])
            rad_edges = bin_centers_to_edges(axis[1])

            eLow = np.power(10, [energy_edges[:-1]])[0]
            eHigh = np.power(10, [energy_edges[1:]])[0]

            rLow = np.power(10, [rad_edges[:-1]])[0]
            rHigh = np.power(10, [rad_edges[1:]])[0]

            # Normalize rpsf by solid angle
            rad_width_deg = np.diff(np.power(10, rad_edges))
            e_sum = np.sum(direction_diff * 2 * rad_width_deg[:, np.newaxis]
                           * np.pi * np.power(10, axis[1])[:, np.newaxis], axis=0)
            normsum = np.divide(((180 / np.pi) ** 2), e_sum, where=e_sum != 0)

            rpsf = direction_diff * normsum
            rpsf_final.append(np.nan_to_num(rpsf).cumsum(axis=0))

        # PSF (3-dim with axes: psf[rad_index, offset_index, energy_index]
        rpsf_final = np.swapaxes(rpsf_final, 0, 1)

        x = np.array(
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
        response_dict["PSF"] = x
    return response_dict
