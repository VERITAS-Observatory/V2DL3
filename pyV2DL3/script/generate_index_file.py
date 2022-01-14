import click
import logging
import os
from pyV2DL3.generateObsHduIndex import create_obs_hdu_index_file


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--folder_location', '-f', nargs=1, type=click.Path(exists=True),
              help='Path containing DL3 files.')
@click.option('--index_file_dir', '-i', nargs=1, type=click.Path(exists=True),
              help='Output path to write the index file.')
@click.option('--hdu_index_file', '-hdu', nargs=1, 
              default='hdu-index.fits.gz',
              help='HDU index file (default:hdu-index.fits.gz)')
@click.option('--obs_index_file', '-obs', nargs=1, 
              default='obs-index.fits.gz',
              help='Observation index file (default:obs-index.fits.gz)')
@click.option('--recreate', '-r', is_flag=True)
@click.option('--debug', '-d', is_flag=True)
def cli(folder_location, index_file_dir, 
        hdu_index_file, obs_index_file,
        recreate, debug):
    """Command line tool for generating index file from a set of DL3 files

    \b
    Index files are created from all DL3 given in folder.

    """
    if folder_location is None:
        click.secho("No DL3 fits files path specified. " +
                    "Assume the current " +
                    "folder is the location of DL3 files.",
                    fg='yellow')
        folder_location = os.getcwd()
    if index_file_dir is None:
        click.secho("No output path specified. " +
                    "Assume the current folder is the output path.",
                    fg='yellow')
        index_file_dir = os.getcwd()

    if debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    if not recreate:
        if os.path.exists('{}/obs-index.fits.gz'.format(index_file_dir)):
            logging.info('Existing {}/obs-index.fits.gz'.format(index_file_dir))
            logging.info('Please remove before continuing or use recreate option -r')
            return
        if os.path.exists('{}/hdu-index.fits.gz'.format(index_file_dir)):
            logging.info('Existing {}/hdu-index.fits.gz'.format(index_file_dir))
            logging.info('Please remove before continuing or use recreate option -r')
            return

    logging.debug(
        "Start by searching all DL3 files in:\n{}".format(folder_location))

    __fits_files = [_file[:-1] for _file in list(os.popen(f'ls {folder_location}/*.fits*'))]
    fits_files = []
    for f in __fits_files:
        if f.find(obs_index_file) == -1 and \
            f.find(hdu_index_file) == -1:
            fits_files.append(f)

    if len(fits_files) == 0:
        logging.error('No fits files found')
        return

    logging.info('Found the following fits files:')
    for f in fits_files:
        logging.info(' -> {0}'.format(f))

    logging.info('Generating index files {}/obs-index.fits.gz and {}/hdu-index.fits.gz'.format(index_file_dir,
                                                                                               index_file_dir))
    create_obs_hdu_index_file(fits_files,
                              index_file_dir,
                              hdu_index_file,
                              obs_index_file)


if __name__ == '__main__':
    cli()
