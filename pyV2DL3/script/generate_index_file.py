import click
import logging
import os

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--folder_location', '-f', nargs=1, type=click.Path(exists=True),
              help='Path containing DL3 files.')
@click.option('--index_file_dir', '-f', nargs=1, type=click.Path(exists=True),
              help='Output path to write the index file.')
@click.option('--debug', '-d', is_flag=True)
def cli(folder_location, index_file_dir, debug):
    """Command line tool for generating index file from a set of DL3 files

    \b
    Just input the folder containing the DL3 files.

    """
    if folder_location is None:
        click.secho("No fits files path was specified. We assume the current " +
                    "folder is the location of all DL3 files.",
                    fg='yellow')
        folder_location = os.getcwd()
    if index_file_dir is None:
        click.secho("No output path was specified. We assume the current folder is the output path.",
                    fg='yellow')
        index_file_dir = os.getcwd()

    if debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    logging.debug("Start by searching all DL3 files withon the folder:\n{}".format(folder_location))

    from pyV2DL3.generateObsHduIndex import create_obs_hdu_index_file

    fits_files = [_file[:-1] for _file in list(os.popen(f'ls {folder_location}/*.fits'))]

    logging.info('Found the following fits files:', fits_files)

    logging.info('Generating index files {}/obs-index.fits.gz and {}/hdu-index.fits.gz'.format(index_file_dir,
                                                                                               index_file_dir))
    create_obs_hdu_index_file(fits_files, index_file_dir)


if __name__ == '__main__':
    cli()
