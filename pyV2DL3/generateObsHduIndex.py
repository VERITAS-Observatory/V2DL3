import logging
import os
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.io.fits import table_to_hdu
from pyV2DL3.addHDUClassKeyword import addHDUClassKeyword
logger = logging.getLogger(__name__)


class NoFitsFileError(Exception):
    pass


def gen_hdu_index(filelist, index_file_dir='./'):
    # create the hdu-index.fits.gz
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
        obs_id = 4 * [dl3_hdu[1].header['OBS_ID']]
        hdu_type_name = ['gti', 'events', 'aeff', 'edisp']
        hdu_type = ['gti', 'events', 'aeff_2d', 'edisp_2d']
        file_dir = 4 * [_path]
        file_name = 4 * [_filename]
        hdu_name = ['GTI', 'EVENTS', 'EFFECTIVE AREA', 'ENERGY DISPERSION']

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
    hdu_table = addHDUClassKeyword(hdu_table,'INDEX', class2='HDU')

    return hdu_table


def gen_obs_index(filelist,index_file_dir='./'):
     # empty lists with the quantities we want
    obs_id = []
    ra_pnt = []
    dec_pnt = []
    zen_pnt = []
    alt_pnt = []
    az_pnt = []
    ontime = []
    livetime = []
    deadc = []
    tstart = []
    tstop = []
    N_TELS = []
    TELLIST = []

    # loop through the files
    for _file in filelist:
        # Get relative path from the index file output dir to
        # fits files.
        _rel_path = os.path.relpath(_file,start=index_file_dir)
        _filename = os.path.basename(_rel_path)
        _path = os.path.dirname(_rel_path)
        if(not os.path.exists(_file)):
            logger.warning('{} does not exist. Skipped!'.format(_file))
            continue
        dl3_hdu = fits.open(_file)
        # let's fill all of them
        obs_id.append(dl3_hdu[1].header['OBS_ID'])
        ra_pnt.append(dl3_hdu[1].header['RA_PNT'])
        dec_pnt.append(dl3_hdu[1].header['DEC_PNT'])
        zen_pnt.append(90 - float(dl3_hdu[1].header['ALT_PNT']))
        alt_pnt.append(dl3_hdu[1].header['ALT_PNT'])
        az_pnt.append(dl3_hdu[1].header['AZ_PNT'])
        ontime.append(dl3_hdu[1].header['ONTIME'])
        livetime.append(dl3_hdu[1].header['LIVETIME'])
        deadc.append(dl3_hdu[1].header['DEADC'])
        tstart.append(dl3_hdu[1].header['TSTART'])
        tstop.append(dl3_hdu[1].header['TSTOP'])
        N_TELS.append(4)
        TELLIST.append(dl3_hdu[1].header['TELLIST'])

    obs_table = Table(
        [obs_id, ra_pnt, dec_pnt, zen_pnt, alt_pnt, az_pnt, ontime, livetime, deadc, tstart, tstop, N_TELS, TELLIST],
        names=(
            'OBS_ID', 'RA_PNT', 'DEC_PNT', 'ZEN_PNT', 'ALT_PNT', 'AZ_PNT', 'ONTIME', 'LIVETIME', 'DEADC', 'TSTART',
            'TSTOP',
            'N_TELS', 'TELLIST'),
        dtype=('>i8', '>f4', '>f4', '>f4', '>f4', '>f4', '>f4', '>f4', '>f4', '>f4', '>f4', '>i8', 'S20')
    )
    #obs_table = Table(
    #    [obs_id, ra_pnt, dec_pnt, tstart, tstop],
    #    names=('OBS_ID', 'RA_PNT', 'DEC_PNT', 'TSTART','TSTOP'),
    #    dtype=('>i8', '>f4', '>f4', '>f4', '>f4')
    #)

    # Set units
    obs_table['RA_PNT'].unit  = 'deg'
    obs_table['DEC_PNT'].unit = 'deg'

    obs_table['ZEN_PNT'].unit = 'deg'
    obs_table['ALT_PNT'].unit = 'deg'
    obs_table['AZ_PNT'].unit = 'deg'
    obs_table['ONTIME'].unit = 's'
    obs_table['LIVETIME'].unit = 's'

    obs_table['TSTART'].unit  = 's'
    obs_table['TSTOP'].unit   = 's'
    if(len(obs_table) ==0):
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
    obs_table = addHDUClassKeyword(obs_table,'INDEX',
                                        class2='OBS')

    return obs_table
  

def create_obs_hdu_index_file(filelist, index_file_dir='./',
                              hdu_index_file= 'hdu-index.fits.gz',
                              obs_index_file= 'obs-index.fits.gz'):
    """Before we start to work with gammapy we will stick to the format established for the data.
    What the pointLikeDL3 is generating is in the following format

    >>> from astropy.io import fits
    >>> hdu = fits.open('MAGIC_data/Output_ST0307_pointLikeDL3.fits')
    >>> hdu.info()
    Filename: MAGIC_data/Output_ST0307_pointLikeDL3.fits
    No.    Name         Type      Cards   Dimensions   Format
            0  PRIMARY     PrimaryHDU       8   ()
            1  EVENTS      BinTableHDU     64   17001R x 7C   [1E, 1E, 1E, 1E, 1E, 1E, 1E]
            2  GTI         BinTableHDU     15   1R x 2C   [1D, 1D]
            3  EFFECTIVE AREA  BinTableHDU     36   1R x 5C   [150E, 150E, 2E, 2E, 300E]
            4  ENERGY DISPERSION  BinTableHDU     39   1R x 7C   [130E, 130E, 80E, 80E, 2E, 2E, 20800E]
            5  BACKGROUND  BinTableHDU     36   1R x 5C   [30E, 30E, 2E, 2E, 60E]

    As you can see in MAGIC we store event list and IRFs in the same FITS file

    A two-level index file scheme is used in the IACT community to allow an arbitrary folder structures
    For each directory tree, two files should be present:
    **obs-index.fits.gz**
    (defined in http://gamma-astro-data-formats.readthedocs.io/en/latest/data_storage/obs_index/index.html)
    **hdu-index.fits.gz**
    (defined in http://gamma-astro-data-formats.readthedocs.io/en/latest/data_storage/hdu_index/index.html)

    obs-index contains the single run informations e.g.
    (OBS_ID, RA_PNT, DEC_PNT, ZEN_PNT, ALT_PNT)
    while hdu-index contains the informations about the locations of the other HDU (header data units)
    necessary to the analysis e.g. A_eff, E_disp and so on...
        http://gamma-astro-data-formats.readthedocs.io/en/latest/

    This function will create the necessary data format, starting from the path that contains the DL3
    converted fits file.

    Parameters
    ----------
    filelist : list
        list of VERITAS dl3 files.

    index_file_dir : path
        directory to save the index files.     

    Example
    -------
    >>> create_obs_hdu_index(['data/magic/run05029747_05029748/run05029747/20131004_05029747_CrabNebula.fits',
                              'data/magic/run05029747_05029748/run05029748/20131004_05029748_CrabNebula.fits'])
    """

    hdu_table = gen_hdu_index(filelist, index_file_dir)
    logger.debug('Writing {} ...'.format(hdu_index_file))
    hdu_table.writeto('{}/{}'.format(index_file_dir, hdu_index_file), overwrite=True)

    obs_table = gen_obs_index(filelist,index_file_dir)
    logger.debug('Writing {} ...'.format(obs_index_file))
    obs_table.writeto('{}/{}'.format(index_file_dir, obs_index_file), overwrite=True)


