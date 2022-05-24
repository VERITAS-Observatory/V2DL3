import logging
import os

import click

from pyV2DL3.EventClass import EventClass

# In event_class_mode, The user can provide multiple EAs per tag to define the event classes for that group.
def runlist_to_file_pair(rl_dict, event_class_mode=False, irfs_to_store=None, override_cuts_validation=False):
    eas = rl_dict["EA"]
    st5s = rl_dict["RUNLIST"]
    file_pair = []
    # Using multiple EAs for event classes
    if event_class_mode:
        for k in st5s.keys():
            event_classes = []
            for ea in eas[k]:
                    event_classes.append(EventClass(ea, irfs_to_store, override_cuts_validation=override_cuts_validation))
            if len(event_classes) == 0:
                raise Exception("No EA filenames defined for runlist tag: " + k)
            for f in st5s[k]:
                file_pair.append((f, event_classes))
    # Default single EA behavior
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
    help="A stage5 or anasum file (<file 1>) and "
    "the corresponding effective area (<file 2>).",
)
@click.option(
    "--runlist",
    "-l",
    nargs=1,
    type=click.Path(exists=True),
    help="Stage6 runlist"
)
@click.option(
    '--event_class_mode',
    '-ec',
    is_flag=True,
    help="Treat EAs from the runlist as different event classes of the same run.",
)
@click.option(
    '--_4dmlm',
    '-4d',
    nargs=2,
    help=("Presets -l, -e, -k, -r 2, --full-enclosure. "
          "ARGS: <runlist> <king params file>")
)
@click.option(
    '--king_function',
    '-k',
    nargs=1,
    type=click.Path(exists=True),
    help='list of PSF variables'
)
@click.option(
    '--event_cut_file',
    '-c',
    nargs=1,
    type=click.Path(exists=True),
    help='File defining spatial exclusion regions to cut'
)
@click.option(
    "--reconstruction_type",
    "-r",
    nargs=1,
    default=0,
    type=click.INT,
    help="1: Standard (.S) - default, 2: ITM (.M3D)"
)
@click.option(
    "--store_msw_msl",
    is_flag=True,
    help="Append MSW and MSL columns to event tables."
    "This is done automatically if using MSW-based event classes."
)
@click.option(
    "--override_cuts_validation",
    "-o",
    is_flag=True,
    help="Override EventClass cuts validations. Use responsibly."
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
@click.option("--ed", "-e", is_flag=True, help="Eventdisplay mode")
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
@click.option(
    "--evt_filter",
    type=click.Path(exists=True),
    help="Load condition to filter events form json or yaml file.",
)
@click.argument("output", metavar="<output>")
def cli(
    file_pair,
    runlist,
    event_class_mode,
    _4dmlm,
    king_function,
    event_cut_file,
    reconstruction_type,
    store_msw_msl,
    override_cuts_validation,
    gen_index_file,
    save_multiplicity,
    ed,
    filename_to_obsid,
    full_enclosure,
    point_like,
    debug,
    verbose,
    output,
    evt_filter,
):
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

    # In older `click`, options with narg >= 2 returned len 0 tuples when not chosen.
    # In newer `click`, they return None to conform with other `option` behaviors.
    if file_pair is not None:
        if len(file_pair) == 0:
            file_pair = None
    if _4dmlm is not None:
        if len(_4dmlm) == 0:
            _4dmlm = None
        else:
            # Set 4dmlm presets
            runlist, king_function = _4dmlm
            event_class_mode = True
            if point_like:
                click.echo(cli.get_help(click.Context(cli)))
                click.secho(
                    "4DMLM uses full-enclosure IRFs. --point-like not supported.", fg="yellow")
                raise click.Abort()
            full_enclosure = True
            if reconstruction_type == 0:
                reconstruction_type = 2

    if file_pair is None and runlist is None:
        click.echo(cli.get_help(click.Context(cli)))
        raise click.Abort()
    if file_pair is not None and runlist is not None:
        click.echo(cli.get_help(click.Context(cli)))
        click.secho("Only one file source can be used.", fg="yellow")
        raise click.Abort()
    if file_pair is not None and event_class_mode:
        click.echo(cli.get_help(click.Context(cli)))
        click.secho("Event_class_mode requires a runlist", fg="yellow")
        raise click.Abort()
    if point_like and full_enclosure:
        click.echo(cli.get_help(click.Context(cli)))
        click.secho("Choose point-like OR full-enclosure (both were selected)", fg="yellow")
        raise click.Abort()

    # Store in a dict the IRFs to be stored within a file.
    # By default we will only store point-like IRFs.
    if not full_enclosure and not point_like:
        point_like = True
        full_enclosure = False
    irfs_to_store = {"full-enclosure": full_enclosure,
                     "point-like": point_like}

    if king_function is not None:
        if point_like:
            click.echo(cli.get_help(click.Context(cli)))
            click.secho(
                "Point-like king function is currently unsupported", fg="yellow")
            raise click.Abort()
        irfs_to_store["psf-king"] = True
        irfs_to_store["psf-king-filename"] = king_function
        # Need EC mode even if not using event classes in order to load MSW ranges for king params
        # EC mode with a single EA per group will behave the same as not using EC mode
        event_class_mode = True
    else:
        irfs_to_store["psf-king"] = False

    # Use .S if not specified
    if reconstruction_type == 0:
        reconstruction_type = 1

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

    if file_pair is not None:
        st5_str, ea_str = file_pair
        if ed or st5_str.find(".anasum.root") >= 0:
            datasource = loadROOTFiles(st5_str, ea_str, "ED", reco_type=reconstruction_type)
        else:
            datasource = loadROOTFiles(st5_str, ea_str, "VEGAS", reco_type=reconstruction_type)
        datasource.set_irfs_to_store(irfs_to_store)
        with cpp_print_context(verbose=verbose):
            datasource.fill_data(evt_filter=evt_filter)
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

        file_pairs = runlist_to_file_pair(rl_dict, 
                                          event_class_mode=event_class_mode,
                                          irfs_to_store=irfs_to_store, 
                                          override_cuts_validation=override_cuts_validation)
        flist = []
        for st5_str, ea_str in file_pairs:
            logging.info(f"Processing file: {st5_str}")
            event_classes = None
            # Rearrange var names for control flow if using event classes
            if event_class_mode:
                event_classes = ea_str
                ea_str = None
                logging.debug(f"Stage5 file:{st5_str}, Event classes:{event_classes}")
            else:
                logging.debug(f"Stage5 file:{st5_str}, EA file:{ea_str}")
            fname_base = os.path.splitext(os.path.basename(st5_str))[0]
            file_type = "ED" if (ed or st5_str.find(".anasum.root") >= 0) else "VEGAS"
            datasource = loadROOTFiles( st5_str, ea_str, file_type, 
                                        event_classes=event_classes,
                                        user_cut_file=event_cut_file, 
                                        reco_type=reconstruction_type, 
                                        store_msw_msl=store_msw_msl,
                                        )

            datasource.set_irfs_to_store(irfs_to_store)
            with cpp_print_context(verbose=verbose):
                datasource.fill_data()
            output_path = f"{output}/{fname_base}"
            # This will be 1 when not using event classes
            eclass_count = len(datasource.get_evt_data())
            if eclass_count < 1:
                raise Exception("No event data found")
            for i in range(0, eclass_count):
                # If multiple event classes, create a subdirectory and append eclass # to filename
                if eclass_count > 1:
                    if not os.path.exists(output_path):
                        os.makedirs(output_path)
                    stage_idx = fname_base.find(".")
                    # Splice our '_ec#' identifier into the filename just before the first '.'
                    if(stage_idx > -1):
                        eclass_fname = fname_base[:stage_idx] + "_ec" + str(i) + fname_base[stage_idx:]
                    # If no '.' found, append to end
                    else:
                        eclass_fname = fname_base + "_ec" + str(i)
                    final_output_path = output_path + "/" + eclass_fname
                else:
                    final_output_path = output_path
                hdulist = genHDUlist(datasource, save_multiplicity=save_multiplicity, event_class_idx=i)
                if filename_to_obsid:
                    logging.info(
                        f"Overwriting OBS_ID={hdulist[1].header['OBS_ID']} with OBS_ID={fname_base}"
                    )
                    hdulist[1].header["OBS_ID"] = fname_base
                final_output_path += ".fits"
                hdulist.writeto(final_output_path, overwrite=True)
                flist.append(final_output_path)

        # Generate hdu obs index file
        if gen_index_file:
            logging.info(
                f"Generating index files {output}/obs-index.fits.gz "
                f"and {output}/hdu-index.fits.gz"
            )
            create_obs_hdu_index_file(flist, output, psf_king=irfs_to_store["psf-king"])


if __name__ == "__main__":
    cli()
