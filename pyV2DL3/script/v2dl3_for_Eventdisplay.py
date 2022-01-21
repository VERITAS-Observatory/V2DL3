import click
import logging
import os
from pyV2DL3.genHDUList import genHDUlist
from pyV2DL3.genHDUList import loadROOTFiles

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--file_pair', '-f', nargs=2, type=click.Path(exists=True),
              help='anasum file (<file 1>) and \
              the corresponding effective area (<file 2>).')
@click.option('--full-enclosure', is_flag=True,
              help='Store full-enclosure IRFs (no direction cut applied)')
@click.option('--point-like', is_flag=True,
              help='Store point-like IRFs (direction cut applied)')
@click.argument('output', metavar='<output>')
@click.option('--debug', '-d', is_flag=True,
              help='Set log level to debug')
@click.option('--logfile', '-l', nargs=1,
              default=None, help='Store log information to file')
@click.option('--instrument_epoch', '-i', nargs=1, default=None,
              help='Instrument epoch to be stored in EVENTS header')
@click.option('--save_multiplicity', '-m', is_flag=True,
              help='Save telescope multiplicity into event list')
@click.option('--filename_to_obsid', '-I', is_flag=True,
              help='Override OBS_ID with output filename')
@click.option('--evt_filter', type=click.Path(exists=True),
              help='Load condition to filter events form json or yaml file.')
def cli(file_pair, full_enclosure, point_like,
        output, debug, logfile, instrument_epoch,
        save_multiplicity, filename_to_obsid, evt_filter):
    """Tool for converting Eventdisplay anasum files and corresponding IRFs to DL3

    """
    if len(file_pair) == 0:
        click.echo(cli.get_help(click.Context(cli)))
        raise click.Abort()

    if debug:
        logging.basicConfig(format='%(levelname)s:v2dl3: %(message)s',
                            level=logging.DEBUG,
                            filename=logfile)
    else:
        logging.basicConfig(format='%(levelname)s:v2dl3: %(message)s',
                            level=logging.INFO,
                            filename=logfile)
    logging.debug('logging level {0}'
                  .format(logging.getLevelName(logging.getLogger().level)))

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
    hdulist = genHDUlist(datasource,
                         save_multiplicity=save_multiplicity,
                         instrument_epoch=instrument_epoch)
    fname_base = os.path.splitext(os.path.basename(output))[0]
    if filename_to_obsid:
        logging.info('Overwriting OBS_ID={0} with OBS_ID={1}'.format(
                     hdulist[1].header['OBS_ID'], fname_base))
        hdulist[1].header['OBS_ID'] = fname_base
    hdulist.writeto(output, overwrite=True)


if __name__ == '__main__':
    cli()
