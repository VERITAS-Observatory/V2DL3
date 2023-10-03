import logging
import os

import astropy.units as u
from astropy.io import fits
from astropy.io.fits import table_to_hdu
from astropy.table import Table, vstack

from pyV2DL3.addHDUClassKeyword import addHDUClassKeyword

logger = logging.getLogger(__name__)

hdu_class_type = {
    ("EVENTS", None): ("events", "events"),
    ("GTI", None): ("gti", "gti"),
    ("RESPONSE", "EFF_AREA"): ("aeff", "aeff_2d"),
    ("RESPONSE", "EDISP"): ("edisp", "edisp_2d"),
    ("RESPONSE", "PSF"): ("psf", None),
    ("RESPONSE", "BKG"): ("bkg", None),
    ("RESPONSE", "RAD_MAX"): ("rad_max", "rad_max_2d"),
}


class NoFitsFileError(Exception):
    pass


def get_hdu_type_and_class(header: fits.header):
    hdu_key = tuple([header.get("HDUCLAS1", None), header.get("HDUCLAS2", None)])
    hdu_type, hdu_class = hdu_class_type[hdu_key]

    if hdu_class is None:
        hdu_class = header["HDUCLAS4"].lower()

    return hdu_type, hdu_class


def gen_hdu_index(filelist, index_file_dir="./"):
    """create HDU index"""

    hdu_tables = []
    # loop through the files
    for _file in filelist:
        # Get relative path from the index file output dir to
        # fits files.
        _rel_path = os.path.relpath(_file, start=index_file_dir)
        _filename = os.path.basename(_rel_path)
        _path = os.path.dirname(_rel_path)

        try:
            with fits.open(_file) as dl3_hdu:
                obsid = dl3_hdu[1].header["OBS_ID"]
                n_hdus = len(dl3_hdu) - 1
                obs_id = [obsid] * n_hdus

                hdu_type_name = []
                hdu_type = []
                hdu_name = []
                file_dir = [_path] * n_hdus
                file_name = [_filename] * n_hdus

                for hdu in dl3_hdu[1:]:
                    type_, class_ = get_hdu_type_and_class(hdu.header)
                    hdu_type_name.append(type_)
                    hdu_type.append(class_)
                    hdu_name.append(hdu.name)

                t = Table(
                    [obs_id, hdu_type_name, hdu_type, file_dir, file_name, hdu_name],
                    names=(
                        "OBS_ID",
                        "HDU_TYPE",
                        "HDU_CLASS",
                        "FILE_DIR",
                        "FILE_NAME",
                        "HDU_NAME",
                    ),
                    dtype=(">i8", "S6", "S10", "S40", "S54", "S20"),
                )

                hdu_tables.append(t)

        except FileNotFoundError:
            logger.warning("{} does not exist. Skipped!".format(_file))

    if not hdu_tables:
        raise NoFitsFileError("No fits file found in the list.")

    hdu_table = vstack(hdu_tables)
    hdu_table = table_to_hdu(hdu_table)
    hdu_table.name = "HDU_INDEX"
    hdu_table = addHDUClassKeyword(hdu_table, "INDEX", class2="HDU")

    return hdu_table


def get_unit_string_from_comment(comment_string):
    """
    Return unit string from FITS comment

    Examples: 'average pointing azimuth [deg]'

    """
    _unit_string = None
    bstart = comment_string.find("[")
    bstopp = comment_string.find("]")
    if bstart != -1 and bstopp != -1:
        _unit_string = comment_string[bstart + 1: bstopp]
        # this is a hack to fix a bug in the VERITAS DL3 files
        if _unit_string == "muA":
            return "uA"
    return _unit_string


def gen_obs_index(filelist, index_file_dir="./", dqm_header=False):
    names = [
        "OBS_ID",
        "RA_PNT",
        "DEC_PNT",
        "ZEN_PNT",
        "ALT_PNT",
        "AZ_PNT",
        "ONTIME",
        "LIVETIME",
        "DEADC",
        "TSTART",
        "TSTOP",
        "N_TELS",
        "TELLIST",
        "OBJECT",
        "RA_OBJ",
        "DEC_OBJ",
        "DATE-OBS",
        "DATE-AVG",
        "DATE-END",
        "NSBLEVEL",
    ]
    data_type = [
        ">i8",
        ">f4",
        ">f4",
        ">f4",
        ">f4",
        ">f4",
        ">f4",
        ">f4",
        ">f4",
        ">f4",
        ">f4",
        ">i8",
        "S20",
        "S20",
        ">f4",
        ">f4",
        "S20",
        "S20",
        "S20",
        ">f4",
    ]
    if dqm_header:
        names, data_type = _add_auxiliary_headers(names, data_type)
    _table_units = {}
    _table_data = {n: [] for n in names}
    missing_keys = set()

    for _file in filelist:
        logger.debug("Reading obs header from %s ...", _file)
        try:
            with fits.open(_file) as dl3_hdu:
                # get values and units from fits header entries
                for key, value in _table_data.items():
                    if not _fill_table_data(
                        key, value, data_type[names.index(key)], dl3_hdu[1].header, dqm_header
                    ):
                        logger.debug(
                            "Keyword %s not found in %s when building observation index file",
                            key, _file
                        )
                        missing_keys.add(key)
                        continue
                    _add_table_units(key, _table_units, dl3_hdu[1].header)
        except FileNotFoundError:
            logger.warning("%s does not exist. Skipped!", _file)
            continue

    for key in missing_keys:
        key_idx = names.index(key)
        del names[key_idx]
        del data_type[key_idx]
        del _table_data[key]

    for key, value in _table_data.items():
        value, _table_units[key] = _check_unit_consistency(
            key, value, _table_units[key]
        )

    obs_table = Table(_table_data, names=names, dtype=data_type)

    for key, value in _table_units.items():
        if value:
            obs_table[key].unit = value

    if len(obs_table) == 0:
        raise NoFitsFileError("No fits file found in the list.")
    obs_table = vstack(obs_table)

    obs_table.meta["MJDREFI"] = dl3_hdu[1].header["MJDREFI"]
    obs_table.meta["MJDREFF"] = dl3_hdu[1].header["MJDREFF"]
    obs_table.meta["TIMEUNIT"] = dl3_hdu[1].header["TIMEUNIT"]
    obs_table.meta["TIMESYS"] = dl3_hdu[1].header["TIMESYS"]
    obs_table.meta["TIMEREF"] = dl3_hdu[1].header["TIMEREF"]
    obs_table.meta["ALTITUDE"] = dl3_hdu[1].header["ALTITUDE"]
    obs_table.meta["GEOLAT"] = dl3_hdu[1].header["GEOLAT"]
    obs_table.meta["GEOLON"] = dl3_hdu[1].header["GEOLON"]

    obs_table = table_to_hdu(obs_table)
    obs_table.name = "OBS_INDEX"
    obs_table = addHDUClassKeyword(obs_table, "INDEX", class2="OBS")

    return obs_table


def create_obs_hdu_index_file(
        filelist,
        index_file_dir="./",
        hdu_index_file="hdu-index.fits.gz",
        obs_index_file="obs-index.fits.gz",
        dqm_header=False,
):
    """Create Observation Index File and HDU index file

    For each directory tree, two files should be present:

    **obs-index.fits.gz**
    defined
    in http://gamma-astro-data-formats.readthedocs.io/en/latest/data_storage/obs_index/index.html

    **hdu-index.fits.gz**
    defined in
    http://gamma-astro-data-formats.readthedocs.io/en/latest/data_storage/hdu_index/index.html

    This function will create the necessary data format, starting from the
    path that contains the DL3 converted fits file.

    Parameters
    ----------
    filelist : list
        list of VERITAS DL3 files.

    index_file_dir : path
        directory to save the index files.

    hdu_index_file : filename
        HDU index file name to be written

    obs_index_file : filename
        Observation index file name to be written

    """

    logger.info("Filling HDU index")
    hdu_table = gen_hdu_index(filelist, index_file_dir)
    logger.debug("Writing {} ...".format(hdu_index_file))
    hdu_table.writeto(f"{index_file_dir}/{hdu_index_file}", overwrite=True)

    logger.info("Filling Obs index")
    obs_table = gen_obs_index(filelist, index_file_dir, dqm_header)
    logger.debug("Writing {} ...".format(obs_index_file))
    obs_table.writeto(f"{index_file_dir}/{obs_index_file}", overwrite=True)


def _add_auxiliary_headers(names, data_type):
    """
    Add auxiliary header entries (mostly DQM related)

    """

    aux_header_list = [
        "PED_VAR ", "QUALITY ", "RUNTYPE ", "OBSMODE ", "RUNSTAT ", "WEATHER ",
        "CONFIG  ", "TRIGCFG ", "DATACAT ", "DQMSTAT ", "DQMREAS ", "DQMMASK ",
        "VPMCFG  ", "LIGHTLEV", "L3RATE  ", "L3RATESD", "CURRMEAN", "CURRSTD ",
        "CURRMED ", "WINDSPE ", "WINDMAX ", "WINDMIN ", "WINDDIR ", "AIRTEMP ",
        "RELHUMID", "FIRMEAN0", "FIRMEAN1", "FIRMEAN3", "FIRSTD0 ",
        "FIRSTD1 ", "FIRSTD3 ", "FIRCORM0", "FIRCORM1", "FIRCORM3",
    ]

    aux_data_type = [
        ">f4", ">i8", "S20", "S20", "S20", "S20",
        ">i8", "S20", "S20", "S20", "S20", ">i8",
        ">i8", "S20", ">f4", ">f4", ">f4", ">f4",
        ">f4", ">f4", ">f4", ">f4", ">f4", ">f4",
        ">f4", ">f4", ">f4", ">f4", ">f4",
        ">f4", ">f4", ">f4", ">f4", ">f4",
    ]

    names.extend(aux_header_list)
    data_type.extend(aux_data_type)

    return names, data_type


def _default_null_value(data_type):
    """
    Return null (none value) for a given data_type

    """

    if data_type == '>f4':
        return -9999.0
    if data_type == '>i8':
        return -9999
    return ""


def _check_unit_consistency(key, value, units):
    """
    Use wherever possible SI units

    """
    if key is not None and key.find("WIND") != -1:
        mph = u.imperial.mile / u.hour
        return [v * mph.to(u.km / u.h) if v is not None else None for v in value], "km/h"
    return value, units


def _add_table_units(key, _table_units, header):
    """
    Fill unit string for table

    """
    _tmp_key = "ALT_PNT" if key == "ZEN_PNT" else key
    try:
        _table_units[key] = get_unit_string_from_comment(
            header.comments[_tmp_key]
        )
    except KeyError:
        if key not in _table_units:
            _table_units[key] = None


def _fill_table_data(key, value, data_type, header, dqm_header):
    """
    Fill table data

    """

    if key == "ZEN_PNT":
        value.append(90.0 - float(header["ALT_PNT"]))
        return True

    try:
        if header[key] == 'NULL' or header[key] is None:
            value.append(_default_null_value(data_type))
        else:
            value.append(header[key])
    except KeyError:
        if dqm_header:
            value.append(_default_null_value(data_type))
        else:
            return False

    return True
