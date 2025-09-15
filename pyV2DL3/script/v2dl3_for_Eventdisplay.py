import logging
import os

import click

from pyV2DL3.genHDUList import genHDUlist, loadROOTFiles
from pyV2DL3.version import __version__

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
IRF_AXIS = ["zenith", "pedvar"]


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
    "--file_pair",
    "-f",
    nargs=2,
    type=click.Path(exists=True),
    help="anasum file (<file 1>) and the corresponding effective area (<file 2>).",
)
@click.option(
    "--full-enclosure",
    is_flag=True,
    help="Store full-enclosure IRFs (no direction cut applied)",
)
@click.option(
    "--point-like", is_flag=True, help="Store point-like IRFs (direction cut applied)"
)
@click.argument("output", metavar="<output file>")
@click.option("--debug", "-d", is_flag=True, help="Set log level to debug")
@click.option(
    "--logfile", "-l", nargs=1, default=None, help="Store log information to file"
)
@click.option(
    "--instrument_epoch",
    "-i",
    nargs=1,
    default=None,
    help="Instrument epoch to be stored in EVENTS header",
)
@click.option(
    "--save_multiplicity",
    "-m",
    is_flag=True,
    help="Save telescope multiplicity into event list",
)
@click.option(
    "--filename_to_obsid",
    "-I",
    is_flag=True,
    help="Override OBS_ID with output filename",
)
@click.option(
    "--evt_filter",
    type=click.Path(exists=True),
    help="Load condition to filter events form json or yaml file.",
)
@click.option(
    "--force_extrapolation",
    is_flag=True,
    help="IRF is extrapolated when parameter is found to be outside IRF range",
)
@click.option(
    "--fuzzy_boundary",
    multiple=True,
    nargs=2,
    type=(click.Choice(IRF_AXIS), click.FLOAT),
    default=None,
    help="Parameter outside IRF range but within a given tolerance is interpolated\
at boundary value. tolerance = ratio of absolute difference between boundary and parameter\
value to boundary. Given for each IRF axes (zenith, pedvar) as key, value pair.",
)
@click.option(
    "--db_fits_file",
    nargs=1,
    type=click.Path(exists=True),
    help="FITS file containing the database tables (including DQM table).",
    default=None,
)
@click.option(
    "--interpolator_name",
    type=click.Choice(["KNeighborsRegressor", "RegularGridInterpolator"]),
    help="Name of the interpolator to be used for IRF interpolation",
    default="KNeighborsRegressor",
)
def cli(
    file_pair,
    full_enclosure,
    point_like,
    output,
    debug,
    logfile,
    instrument_epoch,
    save_multiplicity,
    filename_to_obsid,
    evt_filter,
    force_extrapolation,
    fuzzy_boundary,
    db_fits_file,
    interpolator_name
):
    """Convert Eventdisplay anasum files and corresponding IRFs to DL3"""
    if len(file_pair) == 0:
        click.echo(cli.get_help(click.Context(cli)))
        raise click.Abort()

    if debug:
        logging.basicConfig(
            format="%(levelname)s:v2dl3: %(message)s",
            level=logging.DEBUG,
            filename=logfile,
        )
    else:
        logging.basicConfig(
            format="%(levelname)s:v2dl3: %(message)s",
            level=logging.INFO,
            filename=logfile,
        )
    logging.debug("logging level %s", logging.getLevelName(logging.getLogger().level))

    # default: point like IRFs
    if not full_enclosure and not point_like:
        point_like = True
        full_enclosure = False
    irfs_to_store = {"full-enclosure": full_enclosure, "point-like": point_like}
    logging.info("Converting to DL3 (full-enclosure:  %s, point-like: %s)",
                 full_enclosure, point_like)

    anasum_str, ea_str = file_pair
    logging.info("Eventdisplay: anasum file %s", anasum_str)
    logging.info("Effective: area file %s", ea_str)
    logging.info("IRF axes force extrapolation: %s", force_extrapolation)
    if fuzzy_boundary is not None:
        for key, value in fuzzy_boundary:
            logging.info("Fuzzy boundary setting for %s axis: %s", key, value)
    logging.info("IRF interpolator name: %s", interpolator_name)
    logging.info("Database FITS file: %s", db_fits_file)

    datasource = loadROOTFiles(anasum_str, ea_str, "Eventdisplay")
    datasource.set_irfs_to_store(irfs_to_store)
    datasource.fill_data(
        evt_filter=evt_filter,
        db_fits_file=db_fits_file,
    )
    hdulist = genHDUlist(
        datasource,
        save_multiplicity=save_multiplicity,
        instrument_epoch=instrument_epoch,
    )
    fname_base = os.path.splitext(os.path.basename(output))[0]
    if filename_to_obsid:
        logging.info(
            "Overwriting OBS_ID=%d with OBS_ID=%d",
            hdulist[1].header['OBS_ID'], fname_base
        )
        hdulist[1].header["OBS_ID"] = fname_base
    hdulist.writeto(output, overwrite=True)
    logging.info("FITS output written to %s", output)


if __name__ == "__main__":
    cli()
