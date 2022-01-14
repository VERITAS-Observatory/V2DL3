import logging
import os
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.io.fits import table_to_hdu
from pyV2DL3.addHDUClassKeyword import addHDUClassKeyword
logger = logging.getLogger(__name__)


class NoFitsFileError(Exception):
    pass


def get_hdu_type_and_class(header):
    class1 = header['HDUCLAS1']
    if(class1 == 'EVENTS'):
        return 'events', 'events'
    elif(class1 == 'GTI'):
        return 'gti', 'gti'
    elif(class1 == 'RESPONSE'):
        class2 = header['HDUCLAS2']
        if(class2 == 'EFF_AREA'):
            return 'aeff', 'aeff_2d'
        elif(class2 == 'EDISP'):
            return 'edisp', 'edisp_2d'
        elif(class2 == 'PSF'):
            class4 = header['HDUCLAS4']
            return 'psf', class4.lower()
        elif(class2 == 'BKG'):
            class4 = header['HDUCLAS4']
            return 'bkg', class4.lower()


def gen_hdu_index(filelist, index_file_dir='./'):
    """create HDU index"""

    hdu_tables = []
    # loop through the files
    for _file in filelist:
        # Get relative path from the index file output dir to
        # fits files.
        _rel_path = os.path.relpath(_file, start=index_file_dir)
        _filename = os.path.basename(_rel_path)
        _path = os.path.dirname(_rel_path)

        if(not os.path.exists(_file)):
            logger.warning('{} does not exist. Skipped!'.format(_file))
            continue
        # open the fits file
        dl3_hdu = fits.open(_file)
        # informations to be stored
        obs_id = []
        hdu_type_name = []
        hdu_type = []
        hdu_name = []
        file_dir = []
        file_name = []
        obsid = dl3_hdu[1].header['OBS_ID']
        for hdu in dl3_hdu[1:]:
            obs_id.append(obsid)
            type_, class_ = get_hdu_type_and_class(hdu.header)
            hdu_type_name.append(type_)
            hdu_type.append(class_)
            hdu_name.append(hdu.name)
            file_dir.append(_path)
            file_name.append(_filename)

        t = Table([obs_id, hdu_type_name, hdu_type, file_dir, file_name, hdu_name],
                  names=('OBS_ID', 'HDU_TYPE', 'HDU_CLASS', 'FILE_DIR', 'FILE_NAME', 'HDU_NAME'),
                  dtype=('>i8', 'S6', 'S10', 'S40', 'S54', 'S20')
                  )

        hdu_tables.append(t)
    if len(hdu_tables) == 0:
        raise NoFitsFileError('No fits file found in the list.')

    hdu_table = vstack(hdu_tables)
    hdu_table = table_to_hdu(hdu_table)
    hdu_table.name = 'HDU_INDEX'
    hdu_table = addHDUClassKeyword(hdu_table, 'INDEX', class2='HDU')

    return hdu_table


def get_unit_string_from_comment(comment_string):
    """return unit string from FITS comment"""
    bstart = comment_string.find('[')
    bstopp = comment_string.find(']')
    if bstart != -1 and bstopp != -1:
        return comment_string[bstart+1:bstopp]
    return None


def gen_obs_index(filelist, index_file_dir='./'):
    names = (
        'OBS_ID', 'RA_PNT', 'DEC_PNT', 'ZEN_PNT', 'ALT_PNT', 'AZ_PNT', 'ONTIME', 'LIVETIME', 'DEADC', 'TSTART',
        'TSTOP', 'N_TELS', 'TELLIST', 'OBJECT', 'RA_OBJ', 'DEC_OBJ', 'DATE-OBS', 'DATE-END')
    dtype = (
        '>i8', '>f4', '>f4', '>f4', '>f4', '>f4', '>f4', '>f4', '>f4', '>f4', '>f4', '>i8', 'S20', 'S20',
        '>f4', '>f4', 'S20', 'S20')
    _tabledata = {}
    _tableunits = {}
    for n in names:
        _tabledata[n] = []

    # loop through the files
    for _file in filelist:
        # Get relative path from the index file output dir to
        # fits files.
        _rel_path = os.path.relpath(_file, start=index_file_dir)
        _filename = os.path.basename(_rel_path)
        _path = os.path.dirname(_rel_path)
        if(not os.path.exists(_file)):
            logger.warning('{} does not exist. Skipped!'.format(_file))
            continue
        dl3_hdu = fits.open(_file)
        # get values and units from fits header entries
        for key, value in _tabledata.items():
            if key is 'ZEN_PNT':
                value.append(90.-float(dl3_hdu[1].header['ALT_PNT']))
                _tableunits[key] = get_unit_string_from_comment(dl3_hdu[1].header.comments['ALT_PNT'])
            else:
                value.append(dl3_hdu[1].header[key])
                _tableunits[key] = get_unit_string_from_comment(dl3_hdu[1].header.comments[key])

    obs_table = Table(
       _tabledata,
       names=names,
       dtype=dtype)

    for key, value in _tableunits.items():
        if value:
            obs_table[key].unit = value

    if(len(obs_table) == 0):
        raise NoFitsFileError('No fits file found in the list.')
    obs_table = vstack(obs_table)

    obs_table.meta['MJDREFI'] = dl3_hdu[1].header['MJDREFI']
    obs_table.meta['MJDREFF'] = dl3_hdu[1].header['MJDREFF']
    obs_table.meta['TIMEUNIT'] = dl3_hdu[1].header['TIMEUNIT']
    obs_table.meta['TIMESYS'] = dl3_hdu[1].header['TIMESYS']
    obs_table.meta['TIMEREF'] = dl3_hdu[1].header['TIMEREF']
    obs_table.meta['ALTITUDE'] = dl3_hdu[1].header['ALTITUDE']
    obs_table.meta['GEOLAT'] = dl3_hdu[1].header['GEOLAT']
    obs_table.meta['GEOLON'] = dl3_hdu[1].header['GEOLON']

    obs_table = table_to_hdu(obs_table)
    obs_table.name = 'OBS_INDEX'
    obs_table = addHDUClassKeyword(obs_table, 'INDEX',
                                   class2='OBS')

    return obs_table


def create_obs_hdu_index_file(filelist, index_file_dir='./',
                              hdu_index_file='hdu-index.fits.gz',
                              obs_index_file='obs-index.fits.gz'):
    """Create Observation Index File and HDU index file

    For each directory tree, two files should be present:
    **obs-index.fits.gz**
    (defined in http://gamma-astro-data-formats.readthedocs.io/en/latest/data_storage/obs_index/index.html)
    **hdu-index.fits.gz**
    (defined in http://gamma-astro-data-formats.readthedocs.io/en/latest/data_storage/hdu_index/index.html)

    This function will create the necessary data format, starting from the path that contains the DL3
    converted fits file.

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
    logger.debug('Writing {} ...'.format(hdu_index_file))
    hdu_table.writeto('{}/{}'.format(index_file_dir, hdu_index_file), overwrite=True)

    obs_table = gen_obs_index(filelist, index_file_dir)
    logger.debug('Writing {} ...'.format(obs_index_file))
    obs_table.writeto('{}/{}'.format(index_file_dir, obs_index_file), overwrite=True)
