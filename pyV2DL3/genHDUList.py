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


def loadROOTFiles(data_file, effective_area_file, file_type="VEGAS",
                  event_classes=None,
                  user_cut_file=None,
                  reco_type=1,
                  store_msw_msl=False,
                  ):

    if effective_area_file is None and event_classes is None:
        raise Exception("Running V2DL3 without effective area file(s) is currently unsupported."
                        "\nRemove this exception if you are implementing this behavior.")
    
    if file_type == "VEGAS":
        from pyV2DL3.vegas.VegasDataSource import VegasDataSource
        return VegasDataSource(data_file, effective_area_file, 
                               event_classes=event_classes, 
                               user_cut_file=user_cut_file, 
                               reco_type=reco_type, 
                               store_msw_msl=store_msw_msl,
                               )

    if file_type != "ED":
        raise Exception("File type not supported: {}".format(file_type))

    if effective_area_file is not None:
        from pyV2DL3.eventdisplay.EventDisplayDataSource import EventDisplayDataSource
        return EventDisplayDataSource(data_file, effective_area_file)
    else:
         raise Exception("EventDisplay requires an effective area file."
                        "\nRunning EventDisplay with event classes is unsupported.")


def genHDUlist(datasource, save_multiplicity=False, instrument_epoch=None, event_class_idx=0):
    hdus = [
        genPrimaryHDU(),
        fillEVENTS(
            datasource,
            save_multiplicity=save_multiplicity,
            instrument_epoch=instrument_epoch,
            event_class_idx=event_class_idx,
        ),
        fillGTI(datasource),
    ]
    hdus.extend(fillRESPONSE(datasource, instrument_epoch, event_class_idx=event_class_idx))

    return fits.HDUList(hdus)
