import logging

from astropy.io import fits

from pyV2DL3.fillEVENTS import fillEVENTS
from pyV2DL3.fillGTI import fillGTI
from pyV2DL3.fillRESPONSE import fillRESPONSE

logger = logging.getLogger(__name__)


def genPrimaryHDU():
    """Generate primary hdu"""

    hdu0 = fits.PrimaryHDU()
    hdu0.header.set("TELESCOP", "VERITAS", "Telescope")
    hdu0.header.set("LICENSE ", "", "")
    hdu0.header[
        "COMMENT"
    ] = "FITS (Flexible Image Transport System) \
                              format is defined in 'Astronomy"
    hdu0.header[
        "COMMENT"
    ] = "and Astrophysics', volume 376, page 359; \
                              bibcode: 2001A&A...376..359H"
    return hdu0


def loadROOTFiles(data_file, effective_area_file, file_type="VEGAS"):
    if file_type == "VEGAS":
        from pyV2DL3.vegas.VegasDataSource import VegasDataSource

        return VegasDataSource(data_file, effective_area_file)
    if file_type != "ED":
        raise Exception("File type not supported: {}".format(file_type))
    from pyV2DL3.eventdisplay.EventDisplayDataSource import EventDisplayDataSource

    return EventDisplayDataSource(data_file, effective_area_file)


def genHDUlist(datasource, save_multiplicity=False, instrument_epoch=None, trans_finder=False):
    hdus = [
        genPrimaryHDU(),
        fillEVENTS(
            datasource,
            save_multiplicity=save_multiplicity,
            instrument_epoch=instrument_epoch,
            trans_finder=trans_finder
        ),
        fillGTI(datasource),
    ]
    hdus.extend(fillRESPONSE(datasource, instrument_epoch))
    return fits.HDUList(hdus)
