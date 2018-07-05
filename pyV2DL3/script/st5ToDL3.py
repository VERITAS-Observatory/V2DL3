import click
import logging
import os



@click.command()
@click.option('--file_pair','-f',nargs=2,type=click.Tuple([click.Path(exists=True),click.Path(exists=True)]))
@click.option('--runlist','-l',nargs=1,type=click.Path(exists=True))
@click.option('--debug','-d',is_flag=True)
@click.option('--verbose','-v',is_flag=True)
@click.argument('output',type=str)
def cli(file_pair,runlist,debug,verbose,output):

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
