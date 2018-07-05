import click
import logging
import os


@click.command()
@click.option('--file_pair','-f',nargs=2,type=click.Tuple([click.Path(exists=True),click.Path(exists=True)]),
             help='A stage5 file (<file 1>) and the corresponding effective area (<file 2>).')
@click.option('--runlist','-l',nargs=1,type=click.Path(exists=True),help='Stage6 runlist')
@click.option('--debug','-d',is_flag=True)
@click.option('--verbose','-v',is_flag=True,help='Print root output')
@click.argument('output',metavar='<output>')
def cli(file_pair,runlist,debug,verbose,output):
    """Command line tool for converting stage5 file to DL3

    \b 
    There are two modes:
        1) Single file mode
            When --fifle_par is invoked, the path to the stage5 file and the 
            corresponding effective area should be provided. The <output> argument
            is then the resulting fits file name.
        2) File list mode 
            When using the option --runlist, the path to a stage6 runlist should be used.
            The <output> is then the directory to which the fits while will be saved to.

    Note: One one mode can be used at a time.
    """
    if((len(file_pair) == 0) and (runlist is None)):
        click.echo(cli.get_help(click.Context(cli)) )
        raise click.Abort()    
    if((len(file_pair) > 0) and (runlist is not None)):
        click.echo(cli.get_help(click.Context(cli)) )
        click.echo("Only one file source can be used.")
        raise click.Abort()    

    if(debug):
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    logging.debug("Start importing VEGAS/ROOT") 
    from pyV2DL3.genHDUList import loadROOTFiles,genHDUlist
    from pyV2DL3.load_vegas import CppPrintContext    
    

    if(len(file_pair) > 0):
        st5,ea = loadROOTFiles(file_pair[0],file_pair[1])
        with CppPrintContext(verbose=verbose): 
            hdulist = genHDUlist(st5,ea)
        hdulist.writeto(output, overwrite=True)        
    else:
        click.echo("Not implemented yet!!!!")


if __name__ == '__main__':
    cli()
