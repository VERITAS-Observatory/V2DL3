import click
import logging
import os

def runlist2FP(rl_dict):
    eas  = rl_dict['EA']
    st5s = rl_dict['RUNLIST']
    file_pair = []
    for k in st5s.keys():
        ea = eas[k][0]
        for f in st5s[k]:
            file_pair.append((f,ea))
    return file_pair

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--file_pair','-f',nargs=2,type=click.Path(exists=True),
             help='A stage5 file (<file 1>) and the corresponding effective area (<file 2>).')
@click.option('--runlist','-l',nargs=1,type=click.Path(exists=True),help='Stage6 runlist')
@click.option('--gen_index_file','-g',is_flag=True,
              help='Generate hdu and observation index list files. Only have effect in file list mode.')
@click.option('--save_multiplicity','-m',is_flag=True,
              help='Save telescope multiplicity into event list')
@click.option('--ed','-e',is_flag=True,help='ED mode')
@click.option('--debug','-d',is_flag=True)
@click.option('--verbose','-v',is_flag=True,help='Print root output')
@click.argument('output',metavar='<output>')

def cli(file_pair,runlist,gen_index_file,
        save_multiplicity,ed,debug,verbose,output):
    """Command line tool for converting stage5 file to DL3

    \b
    There are two modes:
        1) Single file mode
            When --file_pair is invoked, the path to the stage5 file and the
            corresponding effective area should be provided. The <output> argument
            is then the resulting fits file name.
        2) File list mode
            When using the option --runlist, the path to a stage6 runlist should be used.
            The <output> is then the directory to which the fits files will be saved to.

    Note: One one mode can be used at a time.
    """
    if((len(file_pair) == 0) and (runlist is None)):
        click.echo(cli.get_help(click.Context(cli)) )
        raise click.Abort()
    if((len(file_pair) > 0) and (runlist is not None)):
        click.echo(cli.get_help(click.Context(cli)) )
        click.secho("Only one file source can be used.",fg='yellow')
        raise click.Abort()

    if(debug):
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    logging.debug("Start importing ROOT")
    from pyV2DL3.genHDUList import loadROOTFiles,genHDUlist
    from pyV2DL3.root_lib_util import CppPrintContext
    from pyV2DL3.parseSt6RunList import (parseRunlistStrs,validateRunlist,
                                         RunlistValidationError,RunlistParsingError)
    from pyV2DL3.generateObsHduIndex import create_obs_hdu_index_file

    if(len(file_pair) > 0):
        st5_str,ea_str = file_pair
        if(ed):
            datasource = loadROOTFiles(st5_str,ea_str,'ED')
        else:
            datasource = loadROOTFiles(st5_str,ea_str,'VEGAS')

        with CppPrintContext(verbose=verbose):
           datasource.fill_data()
        hdulist = genHDUlist(datasource,save_multiplicity=save_multiplicity)
        hdulist.writeto(output, overwrite=True)
    else:
        with open(runlist) as f:
            lines = f.readlines()
        try:
            rl_dict = parseRunlistStrs(lines)
        except RunlistParsingError as e:
            click.secho(str(e),fg='red')
            raise click.Abort()
        try:
            validateRunlist(rl_dict)
        except RunlistValidationError as e:
            click.secho(str(e),fg='red')
            raise click.Abort()
        if(not os.path.exists(output)):
            os.makedirs(output)
        elif(os.path.isfile(output)):
            click.secho("{} already exists as a file. <output> needs to be a directory for runlist mode.".format(output),fg='yellow')
            raise click.Abort()

        file_pairs = runlist2FP(rl_dict)
        flist = []
        for st5_str,ea_str in file_pairs:
           logging.info('Processing file: {}'.format(st5_str))
           logging.debug('Stage5 file:{}, EA file:{}'.format(st5_str,ea_str))
           fname_base = os.path.splitext(os.path.basename(st5_str))[0]
           if(ed):
               datasource = loadROOTFiles(st5_str,ea_str,'ED')
           else:
               datasource = loadROOTFiles(st5_str,ea_str,'VEGAS')

           with CppPrintContext(verbose=verbose):
              datasource.fill_data()
           hdulist = genHDUlist(datasource,save_multiplicity=save_multiplicity)
           hdulist.writeto('{}/{}.fits'.format(output,fname_base), overwrite=True)
           flist.append('{}/{}.fits'.format(output,fname_base))
           # Generate hdu obs index file
        if(gen_index_file):
           logging.info('Generating index files {}/obs-index.fits.gz and {}/hdu-index.fits.gz'.format(output,output))
           create_obs_hdu_index_file(flist,output)

if __name__ == '__main__':
    cli()
