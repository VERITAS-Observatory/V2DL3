import logging

from astropy.io import fits

from pyV2DL3.addHDUClassKeyword import addHDUClassKeyword


def fill_bintablehdu(
    irf_name, hdu_name, class_2, class_3, class_4, response_dict, evt_dict, epoch_str
):
    """fill bin table HDU for IRF"""

    hdu = fits.BinTableHDU(data=response_dict[irf_name])
    hdu.name = hdu_name
    hdu = addHDUClassKeyword(
        hdu, class1="RESPONSE", class2=class_2, class3=class_3, class4=class_4
    )

    hdu.header.set("TELESCOP ", "VERITAS", "")
    hdu.header.set("INSTRUME ", epoch_str, "")
    hdu.header.set("OBS_ID", evt_dict["OBS_ID"], "Run Number")
    hdu.header.set("TUNIT1 ", "TeV", "")
    hdu.header.set("TUNIT2 ", "TeV", "")
    # effective areas
    if hdu_name == "EFFECTIVE AREA":
        hdu.header.set("TUNIT3 ", "deg", "")
        hdu.header.set("TUNIT4 ", "deg", "")
        hdu.header.set("TUNIT5 ", "m2", "")

        hdu.header.set(
            "LO_THRES",
            response_dict["LO_THRES"],
            "Low energy threshold of validity [TeV]",
        )
        hdu.header.set(
            "HI_THRES",
            response_dict["HI_THRES"],
            "High energy threshold of validity [TeV]",
        )
        hdu.header.set("CREF5", "(ENERG_LO:ENERG_HI,THETA_LO:THETA_HI)", "")
    # energy migration matrices
    elif hdu_name == "ENERGY DISPERSION":
        hdu.header.set("TUNIT3 ", "", "")
        hdu.header.set("TUNIT4 ", "", "")
        hdu.header.set("TUNIT5 ", "deg", "")
        hdu.header.set("TUNIT6 ", "deg", "")
        hdu.header.set(
            "CREF7", "(ETRUE_LO:ETRUE_HI,MIGRA_LO:MIGRA_HI,THETA_LO:THETA_HI)", ""
        )
    # PSF
    elif hdu_name == "PSF":
        # PSF king function
        if class_4 == "PSF_KING":
            hdu.header.set('TUNIT1 ', 'TeV', "")  # ENERGY_LO
            hdu.header.set('TUNIT2 ', 'TeV', "")  # ENERGY_HI
            hdu.header.set('TUNIT3 ', 'deg', "")  # THETA_LO
            hdu.header.set('TUNIT4 ', 'deg', "")  # THETA_HI
            hdu.header.set('TUNIT5 ', "", "")     # GAMMA
            hdu.header.set('TUNIT6 ', 'deg', "")  # SIGMA
            hdu.header.set('CREF6', '(ENERG_LO:ENERG_HI,THETA_LO:THETA_HI,GAMMA:SIGMA)', '')

        else:
            # PSF table
            hdu.header.set("TUNIT3 ", "deg", "")
            hdu.header.set("TUNIT4 ", "deg", "")
            hdu.header.set("TUNIT5 ", "deg", "")
            hdu.header.set("TUNIT6 ", "deg", "")
            hdu.header.set("TUNIT7 ", "sr^-1", "")
            hdu.header.set(
                "CREF7", "(ENERG_LO:ENERG_HI,THETA_LO:THETA_HI,RAD_LO:RAD_HI)", ""
            )

    # point-like IRFs
    if class_3 == "POINT-LIKE":
        hdu.header.set(
            "RAD_MAX ", response_dict["RAD_MAX"], "Direction cut applied [deg]"
        )

    return hdu


def fillRESPONSE(datasource, instrument_epoch=None, event_class_index=None):
    response_dict = datasource.get_response_data()
    evt_dict = datasource.get_evt_data()

    if event_class_index is not None:
        response_dict = response_dict[event_class_index]
        evt_dict = evt_dict[event_class_index]

    epoch_str = "VERITAS"
    if instrument_epoch:
        epoch_str = "Epoch " + instrument_epoch
    logging.debug("FITS header instrument epoch set to: {0}".format(epoch_str))

    response_hdus = []
    if datasource.__irf_to_store__["point-like"]:

        response_hdus.append(
            fill_bintablehdu(
                "EA",
                "EFFECTIVE AREA",
                "EFF_AREA",
                "POINT-LIKE",
                "AEFF_2D",
                response_dict,
                evt_dict,
                epoch_str,
            )
        )
        response_hdus.append(
            fill_bintablehdu(
                "MIGRATION",
                "ENERGY DISPERSION",
                "EDISP",
                "POINT-LIKE",
                "EDISP_2D",
                response_dict,
                evt_dict,
                epoch_str,
            )
        )

    elif datasource.__irf_to_store__["full-enclosure"]:

        response_hdus.append(
            fill_bintablehdu(
                "FULL_EA",
                "EFFECTIVE AREA",
                "EFF_AREA",
                "FULL-ENCLOSURE",
                "AEFF_2D",
                response_dict,
                evt_dict,
                epoch_str,
            )
        )
        response_hdus.append(
            fill_bintablehdu(
                "FULL_MIGRATION",
                "ENERGY DISPERSION",
                "EDISP",
                "FULL-ENCLOSURE",
                "EDISP_2D",
                response_dict,
                evt_dict,
                epoch_str,
            )
        )
        # PSF King format if king function was used
        if "psf-king" in datasource.__irf_to_store__:
            if datasource.__irf_to_store__["psf-king"]:
                response_hdus.append(
                    fill_bintablehdu(
                        "PSF",
                        "PSF",
                        "PSF",
                        "FULL-ENCLOSURE",
                        "PSF_KING",
                        response_dict,
                        evt_dict,
                        epoch_str,
                    )
                )
        # Else PSF table
        else:
            response_hdus.append(
                fill_bintablehdu(
                    "PSF",
                    "PSF",
                    "PSF",
                    "FULL-ENCLOSURE",
                    "PSF_TABLE",
                    response_dict,
                    evt_dict,
                    epoch_str,
                )
            )

    else:
        raise Exception("No IRF to store...")

    return response_hdus
