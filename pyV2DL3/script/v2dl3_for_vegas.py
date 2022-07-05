import logging
import os

import click

from pyV2DL3 import EventClass


# In event_class_mode, The user can provide multiple EAs per tag to define the event classes for that group.
def runlist_to_file_pair(rl_dict, event_class_mode=False):
    eas = rl_dict["EA"]
    st5s = rl_dict["RUNLIST"]
    file_pair = []
    # If using multiple EAs for event classes
    if event_class_mode:
        for k in st5s.keys():
            event_classes = []
            for ea in eas[k]:
                event_classes.append(EventClass(ea))
            if len(event_classes) == 0:
                raise Exception("No EA filenames defined for runlist tag: " + k)
            for f in st5s[k]:
                file_pair.append((f, event_classes))
    else:
        for k in st5s.keys():
            ea = eas[k][0]
            for f in st5s[k]:
                file_pair.append((f, ea))
    return file_pair


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--file_pair",
    "-f",
    nargs=2,
    type=click.Path(exists=True),
    help="A stage5 file (<file 1>) and "
    "the corresponding effective area (<file 2>).",
)
@click.option(
    "--runlist", "-l", nargs=1, type=click.Path(exists=True), help="Stage6 runlist"
)
@click.option(
    '--event_class_mode',
    '-ec',
    is_flag=True,
    help="Use EA(s) of the same runlist ID to define event class(es) for that runlist ID.",
)
@click.option(
    "--gen_index_file",
    "-g",
    is_flag=True,
    help="Generate hdu and observation index list files. "
    "Only have effect in file list mode.",
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
    "--full-enclosure",
    is_flag=True,
    help="Store full-enclosure IRFs (no direction cut applied)",
)
@click.option(
    "--point-like", is_flag=True, help="Store point-like IRFs (direction cut applied)"
)
@click.option("--debug", "-d", is_flag=True)
@click.option("--verbose", "-v", is_flag=True, help="Print root output")
@click.argument("output", metavar="<output>")
def cli(
    file_pair,
    runlist,
    event_class_mode,
    gen_index_file,
    save_multiplicity,
    filename_to_obsid,
    full_enclosure,
    point_like,
    debug,
    verbose,
    output,
):
    """Tool for converting VEGAS stage5 files to DL3

    \b
    There are two modes:
        1) Single file mode
            When --file_pair is invoked, the path to the stage5 file
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
    if len(file_pair) > 0:
        if runlist is not None:
            click.echo(cli.get_help(click.Context(cli)))
            click.secho("Only one file source can be used.", fg="yellow")
            raise click.Abort()
        if event_class_mode:
            click.echo(cli.get_help(click.Context(cli)))
            click.secho("Event class mode requires runlist", fg="yellow")
            raise click.Abort()

    if debug:
        logging.basicConfig(
            format="%(levelname)s:v2dl3: %(message)s", level=logging.DEBUG
        )
        print("Logging level DEBUG")
    else:
        logging.basicConfig(
            format="%(levelname)s:v2dl3: %(message)s", level=logging.INFO
        )
        print("Logging level INFO")

    logging.debug("Start importing ROOT")
    from pyV2DL3.genHDUList import genHDUlist
    from pyV2DL3.genHDUList import loadROOTFiles
    from pyV2DL3.vegas.root_lib_util import cpp_print_context

    # Store in a dict the IRFs to be stored within a file.
    # By default we will only store point-like IRFs.
    if not full_enclosure and not point_like:
        point_like = True
        full_enclosure = False
    irfs_to_store = {"full-enclosure": full_enclosure, "point-like": point_like}

    if len(file_pair) > 0:
        st5_str, ea_str = file_pair
        datasource = loadROOTFiles(st5_str, ea_str, "VEGAS")
        datasource.set_irfs_to_store(irfs_to_store)
        with cpp_print_context(verbose=verbose):
            datasource.fill_data()
        hdulist = genHDUlist(datasource, save_multiplicity=save_multiplicity)
        fname_base = os.path.splitext(os.path.basename(output))[0]
        if filename_to_obsid:
            logging.info(
                f"Overwriting OBS_ID={hdulist[1].header['OBS_ID']} with OBS_ID={fname_base}"
            )
            hdulist[1].header["OBS_ID"] = fname_base
        hdulist.writeto(output, overwrite=True)
    else:
        from pyV2DL3.generateObsHduIndex import create_obs_hdu_index_file
        from pyV2DL3.vegas.parseSt6RunList import parseRunlistStrs
        from pyV2DL3.vegas.parseSt6RunList import RunlistParsingError
        from pyV2DL3.vegas.parseSt6RunList import RunlistValidationError
        from pyV2DL3.vegas.parseSt6RunList import validateRunlist

        with open(runlist) as f:
            lines = f.readlines()
        try:
            rl_dict = parseRunlistStrs(lines)
        except RunlistParsingError as e:
            click.secho(str(e), fg="red")
            raise click.Abort()
        try:
            validateRunlist(rl_dict, event_class_mode=event_class_mode)
        except RunlistValidationError as e:
            click.secho(str(e), fg="red")
            raise click.Abort()
        if not os.path.exists(output):
            os.makedirs(output)
        elif os.path.isfile(output):
            click.secho(
                f"{output} already exists as a file. "
                "<output> needs to be a directory for runlist mode.",
                fg="yellow",
            )
            raise click.Abort()

        file_pairs = runlist_to_file_pair(rl_dict, event_class_mode=event_class_mode)
        flist = []
        for st5_str, ea_str in file_pairs:
            logging.info(f"Processing file: {st5_str}")
            event_classes = None
            # Reassign vars if using event classes
            if event_class_mode:
                event_classes = ea_str
                ea_str = None
                logging.debug(f"Stage5 file:{st5_str}, Event classes:{event_classes}")
            else:
                logging.debug(f"Stage5 file:{st5_str}, EA file:{ea_str}")
            fname_base = os.path.splitext(os.path.basename(st5_str))[0]
            datasource = loadROOTFiles(st5_str, ea_str, "VEGAS",
                                       event_classes=event_classes,
                                       )

            datasource.set_irfs_to_store(irfs_to_store)
            with cpp_print_context(verbose=verbose):
                datasource.fill_data()
            hdulist = genHDUlist(datasource, save_multiplicity=save_multiplicity)
            if filename_to_obsid:
                logging.info(
                    f"Overwriting OBS_ID={hdulist[1].header['OBS_ID']} with OBS_ID={fname_base}"
                )
                hdulist[1].header["OBS_ID"] = fname_base
            hdulist.writeto(f"{output}/{fname_base}.fits", overwrite=True)
            flist.append(f"{output}/{fname_base}.fits")
            # Generate hdu obs index file
        if gen_index_file:
            logging.info(
                f"Generating index files {output}/obs-index.fits.gz "
                f"and {output}/hdu-index.fits.gz"
            )
            create_obs_hdu_index_file(flist, output)


if __name__ == "__main__":
    cli()
