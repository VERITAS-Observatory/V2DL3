import logging
import os

import click

from pyV2DL3.generateObsHduIndex import create_obs_hdu_index_file
from pyV2DL3.vegas.EffectiveAreaFile import EffectiveAreaFile


def runlist_to_file_pair(rl_dict):
    eas = rl_dict["EA"]
    st5s = rl_dict["RUNLIST"]
    file_pair = []
    for k in st5s.keys():
        ea_files = []
        for ea in eas[k]:
            ea_files.append(EffectiveAreaFile(ea))
        if len(ea_files) == 0:
            raise Exception("No EA filenames defined for runlist tag: " + k)
        for f in st5s[k]:
            file_pair.append((f, ea_files))

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
    help="Use EA(s) of the same runlist ID to define event class(es) for that runlist ID. "
         + "Event classes sort events to separate fits files based on EA cuts parameters. "
         + "Uses only MSW intervals for now. See EffectiveAreaFile.py to extend."
)
@click.option(
    "--reconstruction_type",
    "-r",
    nargs=1,
    default=1,
    type=click.INT,
    help="1: Standard (.S) - default, 2: ITM (.M3D)"
)
@click.option(
    "--no_fov_cut",
    "-nf",
    is_flag=True,
    help="Disable automatic event cuts according to the EA file's FoVCut parameters"
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
    "--save_msw_msl",
    is_flag=True,
    help="Append MSW and MSL columns to event tables."
    "This is done automatically if using MSW-based event classes."
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
    reconstruction_type,
    no_fov_cut,
    gen_index_file,
    save_multiplicity,
    save_msw_msl,
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

    # File pair mode
    if len(file_pair) > 0:
        st5_str, ea_files = file_pair
        datasource = loadROOTFiles(st5_str, None, "VEGAS",
                                   bypass_fov_cut=no_fov_cut,
                                   ea_files=ea_files,
                                   event_class_mode=event_class_mode,
                                   reco_type=reconstruction_type,
                                   save_msw_msl=save_msw_msl,
                                   )
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
    # Runlist mode
    else:
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

        file_pairs = runlist_to_file_pair(rl_dict)
        flist = []
        for st5_str, ea_files in file_pairs:
            logging.info(f"Processing file: {st5_str}")
            logging.debug(f"Stage5 file:{st5_str}, Event classes:{ea_files}")
            fname_base = os.path.splitext(os.path.basename(st5_str))[0]
            datasource = loadROOTFiles(st5_str, None, "VEGAS",
                                       bypass_fov_cut=no_fov_cut,
                                       ea_files=ea_files,
                                       event_class_mode=event_class_mode,
                                       reco_type=reconstruction_type,
                                       save_msw_msl=save_msw_msl,
                                       )

            datasource.set_irfs_to_store(irfs_to_store)
            with cpp_print_context(verbose=verbose):
                datasource.fill_data()

            # Prepare output paths
            output_path = os.path.join(output, fname_base)
            # This is length 1 when not using event class mode
            num_event_groups = len(datasource.get_evt_data())
            if num_event_groups < 1:
                raise Exception("No event data found")
            for i in range(0, num_event_groups):
                # Make event class subdirectories if there is more than one event group in the VegasDataSource
                if num_event_groups > 1:
                    output_path = make_eclass_path(output, fname_base, i)

                # Write out the fits files
                hdulist = genHDUlist(datasource, save_multiplicity=save_multiplicity, event_class_idx=i)
                if filename_to_obsid:
                    logging.info(
                        f"Overwriting OBS_ID={hdulist[1].header['OBS_ID']} with OBS_ID={fname_base}"
                    )
                    hdulist[1].header["OBS_ID"] = fname_base
                output_path += ".fits"
                hdulist.writeto(output_path, overwrite=True)
                flist.append(output_path)

        if gen_index_file:
            gen_index_files(flist, output, eclass_count=num_event_groups)


"""
Generates the index files for a list of files.

Arguments:
    flist         --  List of .fits filepaths
    output        --  Destination to write the generated index files.
    eclass_count  --  Number of event classes for this batch
"""


def gen_index_files(flist, output, eclass_count=1):
    # Generate master index files
    logging.info(
        f"Generating index files {output}/obs-index.fits.gz "
        f"and {output}/hdu-index.fits.gz"
    )
    create_obs_hdu_index_file(flist, output)

    # Generate index files per event class
    if eclass_count > 1:
        for i in range(0, eclass_count):
            eclass_flist = []
            eclass_output_dir = os.path.join(output, "ec" + str(i))
            for fname in os.listdir(eclass_output_dir):
                eclass_flist.append(os.path.join(eclass_output_dir, fname))
            logging.info(
                f"Generating index files {eclass_output_dir}/obs-index.fits.gz "
                f"and {eclass_output_dir}/hdu-index.fits.gz"
            )
            create_obs_hdu_index_file(eclass_flist, eclass_output_dir)


"""
Sorts output files to subdirectories and appends "_ec#" to their filename
according to event class.

Arguments:
    output       --  Output directory
    fname_base   --  Name of the fits file (without '.fits')
    eclass_idx   --  The event class # that this file belongs to

Returns:
    New output path (excluding '.fits') as a string
"""


def make_eclass_path(output, fname_base, eclass_idx):
    output_path = os.path.join(output, "ec" + str(eclass_idx))
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    stage_idx = fname_base.find(".")
    # Splice an '_ec#' identifier into the filename just before the first '.'
    if(stage_idx > -1):
        eclass_fname = fname_base[:stage_idx] + "_ec" + str(eclass_idx) + fname_base[stage_idx:]
    # If no '.' found, append to end
    else:
        eclass_fname = fname_base + "_ec" + str(eclass_idx)

    return os.path.join(output_path, eclass_fname)


if __name__ == "__main__":
    cli()
