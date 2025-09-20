#!/usr/bin/python
#
# compare two fits files and write
# differences into a file
import click
from astropy.io import fits

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--file_pair",
    "-f",
    nargs=2,
    type=click.Path(exists=True),
    help="pair of fits files to be compared",
)
@click.option("--diff_file", "-d", nargs=1, help="summary of differences")
def cli(file_pair, diff_file):
    """Difference report of two FITS files"""

    if not file_pair or len(file_pair) == 0:
        click.echo(cli.get_help(click.Context(cli)))
        raise click.Abort()

    file1, file2 = file_pair

    fd = fits.FITSDiff(file1, file2, rtol=1.e-4, ignore_keywords=["CREATOR"])
    fd.report(diff_file, overwrite=True)


if __name__ == "__main__":
    cli()
