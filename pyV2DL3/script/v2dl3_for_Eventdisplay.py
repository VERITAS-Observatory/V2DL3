import click
import logging
import os
from pyV2DL3.genHDUList import genHDUlist
from pyV2DL3.genHDUList import loadROOTFiles

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--file_pair', '-f', nargs=2, type=click.Path(exists=True),
              help='A anasum file (<file 1>) and \
              the corresponding effective area (<file 2>).')
@click.option('--gen_index_file', '-g', is_flag=True,
              help='Generate hdu and observation index list files. \
              Only have effect in file list mode.')
@click.option('--save_multiplicity', '-m', is_flag=True,
              help='Save telescope multiplicity into event list')
@click.option('--filename_to_obsid', '-I', is_flag=True,
              help='Override OBS_ID with output filename')
@click.option('--full-enclosure', is_flag=True,
              help='Store full-enclosure IRFs (no direction cut applied)')
@click.option('--point-like', is_flag=True,
              help='Store point-like IRFs (direction cut applied)')
@click.option('--debug', '-d', is_flag=True)
@click.option('--verbose', '-v', is_flag=True, help='Print root output')
@click.option('--evt_filter', type=click.Path(exists=True),
              help='Load condition to filter events form json or yaml file.')
@click.argument('output', metavar='<output>')
def cli(file_pair, gen_index_file, save_multiplicity,
        filename_to_obsid, full_enclosure, point_like,
        debug, verbose, output, evt_filter):
    """Tool for converting Eventdisplay anasum files and corresponding IRFs to DL3

    """
    if len(file_pair) == 0:
        click.echo(cli.get_help(click.Context(cli)))
        raise click.Abort()

    if debug:
        logging.basicConfig(level=logging.DEBUG)
        print("Logging level DEBUG")
    else:
        logging.basicConfig(level=logging.INFO)
        print("Logging level INFO")

    # By default we will only store point-like IRFs.
    if not full_enclosure and not point_like:
        point_like = True
        full_enclosure = False
    irfs_to_store = {'full-enclosure': full_enclosure,
                     'point-like': point_like}

    anasum_str, ea_str = file_pair
    datasource = loadROOTFiles(anasum_str, ea_str, 'ED')
    datasource.set_irfs_to_store(irfs_to_store)
    datasource.fill_data(evt_filter=evt_filter)
    hdulist = genHDUlist(datasource, save_multiplicity=save_multiplicity)
    fname_base = os.path.splitext(os.path.basename(output))[0]
    if filename_to_obsid:
        logging.info('Overwriting OBS_ID={0} with OBS_ID={1}'.format(
                     hdulist[1].header['OBS_ID'], fname_base))
        hdulist[1].header['OBS_ID'] = fname_base
    hdulist.writeto(output, overwrite=True)


if __name__ == '__main__':
    cli()
