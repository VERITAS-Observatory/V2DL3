"""
Provide tools to access DB FITS files.

"""

import logging

import astropy.io.registry
import numpy as np
from astropy.table import Table

logger = logging.getLogger(__name__)


def read_db_fits_file(db_fits_file):
    """
    Read DB FITS file and return a dictionary with the DQM table

    """

    db_dict = {}
    if db_fits_file is None:
        return db_dict

    logger.info("Reading DB FITS file: %s", db_fits_file)
    try:
        db_table = Table.read(db_fits_file, hdu="DQM")
        db_dict = {key: value[0] for key, value in db_table.items()}
    except FileNotFoundError:
        logger.error("DB FITS file not found: %s", db_fits_file)
        raise
    except astropy.io.registry.base.IORegistryError:
        logger.error("DB FITS file is not a FITS file: %s", db_fits_file)
        raise
    except KeyError:
        logger.error("DB FITS file does not contain DQM table: %s", db_fits_file)
        raise

    ensure_nan_instead_masked_arrays(db_dict)

    return db_dict


def ensure_nan_instead_masked_arrays(db_dict):
    """
    Replace masked arrays with NaNs
    (FITS headers do not allow for NaN values)

    """

    for key, value in db_dict.items():
        if isinstance(value, np.ma.core.MaskedConstant):
            db_dict[key] = None
