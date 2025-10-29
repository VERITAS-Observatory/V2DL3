"""
Generate index files required for the analysis with science tools
like gammapy or ctools.

The index contain each contain the following information:

Observation index table describing the observation (mainly the pointing
direction and observation start and stop time.

The HDU index table contains the information file locations and directories
of each DL3 file.

"""

import glob
import logging
import os

import click

from pyV2DL3.generateObsHduIndex import create_obs_hdu_index_file
from pyV2DL3.version import __version__

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(f'pyV2DL3 version {__version__}')
    ctx.exit()


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--version', '-v', is_flag=True, callback=print_version,
              expose_value=False, is_eager=True,
              help='Show version and exit.')
@click.option(
    "--folder_location",
    "-f",
    nargs=1,
    type=click.Path(exists=True),
    help="Path containing DL3 files.",
)
@click.option(
    "--index_file_dir",
    "-i",
    nargs=1,
    type=click.Path(exists=True),
    help="Output path to write the index file.",
)
@click.option(
    "--hdu_index_file",
    "-hdu",
    nargs=1,
    default="hdu-index.fits.gz",
    help="HDU index file (default:hdu-index.fits.gz)",
)
@click.option(
    "--obs_index_file",
    "-obs",
    nargs=1,
    default="obs-index.fits.gz",
    help="Observation index file (default:obs-index.fits.gz)",
)
@click.option(
    "--recreate",
    "-r",
    is_flag=True,
    help="Recreate index files even if they already exist.",
)
@click.option("--debug", "-d", is_flag=True)
@click.option(
    "--dqm_header",
    "-dqm",
    is_flag=True,
    help="Add DQM header to the HDU index file.",
)
def cli(
    folder_location, index_file_dir, hdu_index_file, obs_index_file, recreate, debug, dqm_header
):
    """Command line tool for generating index file from a set of DL3 files

    \b
    Index files are created from all DL3 fits files found in the given folder.

    """
    if folder_location is None:
        click.secho(
            "No DL3 fits files path specified. "
            "Assume the current "
            "folder is the location of DL3 files.",
            fg="yellow",
        )
        folder_location = os.getcwd()
    if index_file_dir is None:
        click.secho(
            "No output path specified. "
            "Assume the current folder is the output path.",
            fg="yellow",
        )
        index_file_dir = os.getcwd()

    if debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    logging.info("Generating HUD and OBS index files")

    if not recreate:
        if os.path.exists(f"{index_file_dir}/obs-index.fits.gz"):
            logging.info("Existing %s/obs-index.fits.gz", index_file_dir)
            logging.info("Remove before continuing or use recreate option -r")
            return
        if os.path.exists(f"{index_file_dir}/hdu-index.fits.gz"):
            logging.info("Existing %s/hdu-index.fits.gz", index_file_dir)
            logging.info("Remove before continuing or use recreate option -r")
            return

    logging.debug("Start by searching all DL3 files in:\n%s", folder_location)

    __fits_files = glob.glob(f"{folder_location}/[0-9]*.fits*")
    if len(__fits_files) == 0:
        logging.info("No FITS files found, trying Eventdisplay-style DL3 archive folder.")
        __fits_files = glob.glob(f"{folder_location}/[0-9]*/[0-9]*.fits*")
    fits_files = [
        f
        for f in __fits_files
        if f.find(obs_index_file) == -1 and f.find(hdu_index_file) == -1
    ]
    if not fits_files:
        logging.error("No fits files found")
        return
    fits_files.sort()

    logging.debug("Found the following fits files:")
    for _file in fits_files:
        logging.debug(" -> %s", _file)

    logging.info("Found %d fits files", len(fits_files))
    create_obs_hdu_index_file(
        fits_files, index_file_dir, hdu_index_file, obs_index_file, dqm_header
    )


if __name__ == "__main__":
    cli()
