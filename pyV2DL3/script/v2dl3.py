import click
import logging
import os


def runlist_to_file_pair(rl_dict):
    eas = rl_dict['EA']
    st5s = rl_dict['RUNLIST']
    file_pair = []
    for k in st5s.keys():
        ea = eas[k][0]
        for f in st5s[k]:
            file_pair.append((f, ea))
    return file_pair


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--file_pair', '-f', nargs=2, type=click.Path(exists=True),
              help='A stage5 or anasum file (<file 1>) and \
              the corresponding effective area (<file 2>).')
@click.option('--runlist', '-l', nargs=1,
              type=click.Path(exists=True),
              help='Stage6 runlist')
@click.option('--gen_index_file', '-g', is_flag=True,
              help='Generate hdu and observation index list files. \
              Only have effect in file list mode.')
@click.option('--save_multiplicity', '-m', is_flag=True,
              help='Save telescope multiplicity into event list')
@click.option('--ed', '-e', is_flag=True, help='Eventdisplay mode')
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
def cli(file_pair, runlist, gen_index_file, save_multiplicity,
        ed, filename_to_obsid, full_enclosure, point_like,
        debug, verbose, output, evt_filter):
    """Tool for converting VEGAS stage5 or Eventdisplay anasum files to DL3

    \b
    There are two modes:
        1) Single file mode
            When --file_pair is invoked, the path to the stage5/anasum file
            and the corresponding effective area should be provided.
            The <output> argument is then the resulting fits file name.
        2) File list mode
            When using the option --runlist, the path to a stage6 runlist
            should be used.  The <output> is then the directory to which
            the fits files will be saved to.

    Note: One one mode can be used at a time.
    """
    if len(file_pair) == 0 and runlist is None:
        click.echo(cli.get_help(click.Context(cli)))
        raise click.Abort()
    if len(file_pair) > 0 and runlist is not None:
        click.echo(cli.get_help(click.Context(cli)))
        click.secho("Only one file source can be used.", fg='yellow')
        raise click.Abort()

    if debug:
        logging.basicConfig(format='%(levelname)s:v2dl3: %(message)s',
                            level=logging.DEBUG)
        print("Logging level DEBUG")
    else:
        logging.basicConfig(format='%(levelname)s:v2dl3: %(message)s',
                            level=logging.INFO)
        print("Logging level INFO")

    logging.debug("Start importing ROOT")
    from pyV2DL3.genHDUList import genHDUlist
    from pyV2DL3.genHDUList import loadROOTFiles
    from pyV2DL3.root_lib_util import cpp_print_context

    # Store in a dict the IRFs to be stored within a file.
    # By default we will only store point-like IRFs.
    if not full_enclosure and not point_like:
        point_like = True
        full_enclosure = False
    irfs_to_store = {'full-enclosure': full_enclosure,
                     'point-like': point_like}

    if len(file_pair) > 0:
        st5_str, ea_str = file_pair
        if ed or st5_str.find('.anasum.root') >= 0:
            datasource = loadROOTFiles(st5_str, ea_str, 'ED')
        else:
            datasource = loadROOTFiles(st5_str, ea_str, 'VEGAS')
        datasource.set_irfs_to_store(irfs_to_store)
        with cpp_print_context(verbose=verbose):
            datasource.fill_data(evt_filter=evt_filter)
        hdulist = genHDUlist(datasource, save_multiplicity=save_multiplicity)
        fname_base = os.path.splitext(os.path.basename(output))[0]
        if filename_to_obsid:
            logging.info('Overwriting OBS_ID={0} with OBS_ID={1}'.format(
                         hdulist[1].header['OBS_ID'], fname_base))
            hdulist[1].header['OBS_ID'] = fname_base
        hdulist.writeto(output, overwrite=True)
    else:
        from pyV2DL3.generateObsHduIndex import create_obs_hdu_index_file
        from pyV2DL3.parseSt6RunList import parseRunlistStrs
        from pyV2DL3.parseSt6RunList import RunlistParsingError
        from pyV2DL3.parseSt6RunList import RunlistValidationError
        from pyV2DL3.parseSt6RunList import validateRunlist

        with open(runlist) as f:
            lines = f.readlines()
        try:
            rl_dict = parseRunlistStrs(lines)
        except RunlistParsingError as e:
            click.secho(str(e), fg='red')
            raise click.Abort()
        try:
            validateRunlist(rl_dict)
        except RunlistValidationError as e:
            click.secho(str(e), fg='red')
            raise click.Abort()
        if not os.path.exists(output):
            os.makedirs(output)
        elif os.path.isfile(output):
            click.secho("{} already exists as a file. \
                        <output> needs to be a directory for runlist mode."
                        .format(output),
                        fg='yellow')
            raise click.Abort()

        file_pairs = runlist_to_file_pair(rl_dict)
        flist = []
        for st5_str, ea_str in file_pairs:
            logging.info('Processing file: {}'.format(st5_str))
            logging.debug('Stage5 file:{}, EA file:{}'.format(st5_str, ea_str))
            fname_base = os.path.splitext(os.path.basename(st5_str))[0]
            if ed or st5_str.find('.anasum.root') >= 0:
                datasource = loadROOTFiles(st5_str, ea_str, 'ED')
            else:
                datasource = loadROOTFiles(st5_str, ea_str, 'VEGAS')

            datasource.set_irfs_to_store(irfs_to_store)
            with cpp_print_context(verbose=verbose):
                datasource.fill_data()
            hdulist = genHDUlist(datasource,
                                 save_multiplicity=save_multiplicity)
            if filename_to_obsid:
                logging.info('Overwriting OBS_ID={0} with OBS_ID={1}'
                             .format(hdulist[1].header['OBS_ID'], fname_base))
                hdulist[1].header['OBS_ID'] = fname_base
            hdulist.writeto('{}/{}.fits'.format(output, fname_base),
                            overwrite=True)
            flist.append('{}/{}.fits'.format(output, fname_base))
            # Generate hdu obs index file
        if gen_index_file:
            logging.info('Generating index files {}/obs-index.fits.gz \
                          and {}/hdu-index.fits.gz'.format(output, output))
            create_obs_hdu_index_file(flist, output)


if __name__ == '__main__':
    cli()
