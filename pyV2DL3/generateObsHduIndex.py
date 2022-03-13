import logging
import os

from astropy.io import fits
from astropy.io.fits import table_to_hdu
from astropy.table import Table
from astropy.table import vstack

from pyV2DL3.addHDUClassKeyword import addHDUClassKeyword

hdu_class_type = {
    ("EVENTS", None): ("events", "events"),
    ("GTI", None): ("gti", "gti"),
    ("RESPONSE", "EFF_AREA"): ("aeff", "aeff_2d"),
    ("RESPONSE", "EDISP"): ("edisp", "edisp_2d"),
    ("RESPONSE", "PSF"): ("psf", None),
    ("RESPONSE", "BKG"): ("bkg", None),
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

        if not os.path.exists(_file):
            logging.warning("{} does not exist. Skipped!".format(_file))
            continue

        # open the fits file
        dl3_hdu = fits.open(_file)

        # informations to be stored
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
    if not hdu_tables:
        raise NoFitsFileError("No fits file found in the list.")

    hdu_table = vstack(hdu_tables)
    hdu_table = table_to_hdu(hdu_table)
    hdu_table.name = "HDU_INDEX"
    hdu_table = addHDUClassKeyword(hdu_table, "INDEX", class2="HDU")

    return hdu_table


def get_unit_string_from_comment(comment_string):
    """return unit string from FITS comment"""
    bstart = comment_string.find("[")
    bstopp = comment_string.find("]")
    if bstart != -1 and bstopp != -1:
        return comment_string[bstart + 1 : bstopp]
    return None


def gen_obs_index(filelist, index_file_dir="./"):
    names = (
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
        "DATE-END",
    )
    dtype = (
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
    )
    _tableunits = {}
    _tabledata = {n: [] for n in names}
    # loop through the files
    for _file in filelist:
        # fits files.
        if not os.path.exists(_file):
            logging.warning("{} does not exist. Skipped!".format(_file))
            continue
        dl3_hdu = fits.open(_file)
        # get values and units from fits header entries
        for key, value in _tabledata.items():
            if key == "ZEN_PNT":
                value.append(90.0 - float(dl3_hdu[1].header["ALT_PNT"]))
                _tableunits[key] = get_unit_string_from_comment(
                    dl3_hdu[1].header.comments["ALT_PNT"]
                )
            else:
                value.append(dl3_hdu[1].header[key])
                _tableunits[key] = get_unit_string_from_comment(
                    dl3_hdu[1].header.comments[key]
                )

    obs_table = Table(_tabledata, names=names, dtype=dtype)

    for key, value in _tableunits.items():
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

    hdu_table = gen_hdu_index(filelist, index_file_dir)
    logging.debug("Writing {} ...".format(hdu_index_file))
    hdu_table.writeto(f"{index_file_dir}/{hdu_index_file}", overwrite=True)

    obs_table = gen_obs_index(filelist, index_file_dir)
    logging.debug("Writing {} ...".format(obs_index_file))
    obs_table.writeto(f"{index_file_dir}/{obs_index_file}", overwrite=True)
