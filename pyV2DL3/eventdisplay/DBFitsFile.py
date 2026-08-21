"""
Provide tools to access DB FITS files.

"""

import logging

import astropy.io.registry
import numpy as np
from astropy.table import Table

from pyV2DL3.fillEVENTS import non_standard_hdu_keys_and_comments

logger = logging.getLogger(__name__)


# DQM columns that are intentionally exposed as optional EVENTS headers.  All
# other DQM columns are ignored so that an arbitrary database column cannot
# become part of the conversion state by accident.
SUPPORTED_DB_COLUMNS = frozenset(non_standard_hdu_keys_and_comments())

# The DQM row identifier has had several spellings in database exports.  The
# first matching name below is preferred, while matching itself is
# case/underscore insensitive.
RUN_ID_COLUMNS = (
    "runNumber",
    "run_number",
    "run_id",
    "run",
    "obs_id",
    "OBS_ID",
)


def _normalise_name(name):
    """Normalize FITS column names for matching known fields."""

    return "".join(
        character.lower() for character in str(name) if character.isalnum()
    )


def _normalise_run_id(value):
    """Return a comparable scalar run identifier."""

    if np.ma.is_masked(value) or value is None:
        return None

    try:
        numeric_value = float(value)
    except (TypeError, ValueError):
        return str(value).strip()

    if numeric_value.is_integer():
        return int(numeric_value)
    return numeric_value


def _scalar_value(value):
    """Convert FITS scalar values to ordinary Python values."""

    if np.ma.is_masked(value):
        return None
    if isinstance(value, np.generic):
        return value.item()
    return value


def _find_run_id_column(column_names):
    """Find the preferred run identifier column in a DQM table."""

    by_normalised_name = {
        _normalise_name(name): name for name in column_names
    }
    for candidate in RUN_ID_COLUMNS:
        column = by_normalised_name.get(_normalise_name(candidate))
        if column is not None:
            return column
    return None


def read_db_fits_file(db_fits_file, run_number=None, protected_keys=()):
    """
    Read the DQM metadata for one run.

    A database file is only safe to use when its DQM row is associated with
    the run being converted.  Therefore a run number is required whenever a
    database file is supplied, and exactly one matching row must exist.
    ``protected_keys`` contains metadata already derived from the anasum file;
    a DQM column matching one of those keys is rejected instead of replacing
    it.

    """

    db_dict = {}
    if db_fits_file is None:
        return db_dict
    if run_number is None:
        raise ValueError("run_number is required when reading a DB FITS file")

    logger.info("Reading DB FITS file: %s", db_fits_file)
    try:
        db_table = Table.read(db_fits_file, hdu="DQM")
    except FileNotFoundError:
        logger.error("DB FITS file not found: %s", db_fits_file)
        raise
    except astropy.io.registry.base.IORegistryError:
        logger.error("DB FITS file is not a FITS file: %s", db_fits_file)
        raise
    except KeyError:
        logger.error(
            "DB FITS file does not contain DQM table: %s", db_fits_file
        )
        raise

    run_id_column = _find_run_id_column(db_table.colnames)
    if run_id_column is None:
        raise ValueError(
            "DQM table in {} has no supported run identifier column; "
            "expected one of {}".format(db_fits_file, RUN_ID_COLUMNS)
        )

    expected_run = _normalise_run_id(run_number)
    matching_rows = [
        index
        for index, value in enumerate(db_table[run_id_column])
        if _normalise_run_id(value) == expected_run
    ]
    if len(matching_rows) != 1:
        raise ValueError(
            "DQM table in {} has {} rows for run {} "
            "(expected exactly one)".format(
                db_fits_file, len(matching_rows), run_number
            )
        )

    row = db_table[matching_rows[0]]
    supported_columns = {
        _normalise_name(column): column for column in SUPPORTED_DB_COLUMNS
    }
    protected_columns = {_normalise_name(key) for key in protected_keys}
    run_id_name = _normalise_name(run_id_column)

    for column in db_table.colnames:
        normalised_column = _normalise_name(column)
        if normalised_column == run_id_name:
            continue
        if normalised_column in protected_columns:
            raise ValueError(
                "DQM column {!r} would overwrite core "
                "Eventdisplay metadata".format(column)
            )
        if normalised_column in protected_columns:
            raise ValueError(
                "DQM column {!r} would overwrite core "
                "Eventdisplay metadata".format(column)
            )
        canonical_name = supported_columns.get(normalised_column)
        if canonical_name is not None:
            db_dict[canonical_name] = _scalar_value(row[column])
        else:
            logger.debug("Ignoring unsupported DQM column %s", column)

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
