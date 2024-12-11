import logging
import os

import click

from pyV2DL3.generateObsHduIndex import create_obs_hdu_index_file

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--file_pair",
    "-f",
    nargs=2,
    type=click.Path(exists=True),
    help="A stage5 file (<file 1>) and " "the corresponding effective area (<file 2>).",
)
@click.option(
    "--runlist", "-l", nargs=1, type=click.Path(exists=True), help="Stage6 runlist"
)
@click.option(
    "--event_class_mode",
    "-ec",
    is_flag=True,
    help="Use EA file(s) of the same runlist ID to define event class(es) for that runlist ID. "
    + "Event classes sort events to separate fits files based on EA file cuts parameters. "
    + "See README_vegas.md for more information.",
)
@click.option(
    "--psf_king",
    "-k",
    nargs=1,
    type=click.Path(exists=True),
    help="Provide a file of King PSF parameter values. See README_vegas.md for more information",
)
@click.option(
    "--reconstruction_type",
    "-r",
    nargs=1,
    default=1,
    type=click.INT,
    help="1: Standard (.S) - default, 2: ITM (.M3D)",
)
@click.option(
    "--no_fov_cut",
    "-nf",
    is_flag=True,
    help="Disable automatic event cuts according to the EA file's FoVCut parameters",
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
    "This is done automatically if using MSW-based event classes.",
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
@click.option(
    "--energy-bias-corr", is_flag=True, help="Apply the energy bias correction"
)
@click.option("--debug", "-d", is_flag=True)
@click.option("--verbose", "-v", is_flag=True, help="Print root output")
@click.argument("output", metavar="<output>")
def cli(
    file_pair,
    runlist,
    event_class_mode,
    psf_king,
    reconstruction_type,
    no_fov_cut,
    gen_index_file,
    save_multiplicity,
    save_msw_msl,
    filename_to_obsid,
    full_enclosure,
    point_like,
    energy_bias_corr,
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
            The <output> argument is then the directory to which
            the fits files will be saved to.
        2) File list mode
            When using the option --runlist, the path to a stage6 runlist
            should be used.  The <output> is then the directory to which
            the fits files will be saved to.

    Note: One one mode can be used at a time.
    """

    # Before click 8+, options with narg >= 2 returned len 0 tuples when not chosen.
    # Both should be supported as many existing setups for VEGAS are unable to upgrade to click 8+
    if file_pair is not None:
        if len(file_pair) == 0:
            file_pair = None

    if file_pair is None and runlist is None:
        click.echo(cli.get_help(click.Context(cli)))
        raise click.Abort()
    if file_pair is not None:
        if runlist is not None:
            click.echo(cli.get_help(click.Context(cli)))
            click.secho("Only one file source can be used.", fg="yellow")
            raise click.Abort()
        if event_class_mode:
            click.echo(cli.get_help(click.Context(cli)))
            click.secho("Event class mode requires runlist", fg="yellow")
            raise click.Abort()

    if psf_king is not None and not full_enclosure:
        click.echo(cli.get_help(click.Context(cli)))
        click.secho(
            "PSF king function should be used for full-enclosure analysis", fg="yellow"
        )
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
    from pyV2DL3.genHDUList import genHDUlist, loadROOTFiles
    from pyV2DL3.vegas.root_lib_util import (  # Lundy- I don't like that this is where these are imported
        cpp_print_context,
    )

    # Store in a dict the IRFs to be stored within a file.
    # By default we will only store point-like IRFs.
    if not full_enclosure and not point_like:
        point_like = True
        full_enclosure = False
    irfs_to_store = {"full-enclosure": full_enclosure, "point-like": point_like}

    # These will be passed to VegasDataSource
    datasource_kwargs = {
        "bypass_fov_cut": no_fov_cut,
        "event_class_mode": event_class_mode,
        "reco_type": reconstruction_type,
        "save_msw_msl": save_msw_msl,
        "corr_EB_params" : energy_bias_corr,
    }

    if psf_king is not None:
        from pyV2DL3.vegas.util import load_psf_king_parameters

        psf_king_params = load_psf_king_parameters(psf_king)
        datasource_kwargs["psf_king_params"] = psf_king_params
        irfs_to_store["psf-king"] = True
    """
    The main loop for all of the files to call and store the contents
    into the fits files. This function is in a weird spot because
    of the imports of the previous functions.
    """

    def processFilePair(
        st5_str,
        ea_files,
        flist,
        failed_list,
        output,
        save_multiplicity,
        filename_to_obsid,
        datasource_kwargs,
        irfs_to_store,
        st6_config,
        verbose,
    ):
        logging.info(f"Processing file: {st5_str}")
        logging.debug(f"Stage5 file:{st5_str}, Event classes:{ea_files}")
        fname_base = os.path.splitext(os.path.basename(st5_str))[0]
        datasource_kwargs.update({"st6_configs" : st6_config})
        datasource = loadROOTFiles(st5_str, ea_files, "VEGAS", **datasource_kwargs)
        datasource.set_irfs_to_store(irfs_to_store)
        with cpp_print_context(verbose=verbose):
            try:
                datasource.fill_data()
            except Exception as e:
                logging.info("Exception encountered in " + st5_str + ":")
                logging.info(e)
                # We don't want one run's problem to stop the entire batch
                logging.info("Skipping " + st5_str)
                failed_list[st5_str] = e
                return flist, failed_list, 1

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
            hdulist = genHDUlist(
                datasource, save_multiplicity=save_multiplicity, event_class_idx=i
            )
            if filename_to_obsid:
                logging.info(
                    f"Overwriting OBS_ID={hdulist[1].header['OBS_ID']} with OBS_ID={fname_base}"
                )
                hdulist[1].header["OBS_ID"] = fname_base
            output_path += ".fits"
            hdulist.writeto(output_path, overwrite=True)
            flist.append(output_path)

        return flist, failed_list, num_event_groups

    flist = []
    failed_list = {}
    # File pair mode
    if file_pair is not None:
        from pyV2DL3.vegas.EffectiveAreaFile import EffectiveAreaFile

        if not os.path.exists(output):
            os.makedirs(output)
        st5_str, ea_file = file_pair
        ea_file = EffectiveAreaFile(ea_file)
        flist, failed_list, num_event_groups = processFilePair(
            st5_str,
            ea_file,
            flist,
            failed_list,
            output,
            save_multiplicity,
            filename_to_obsid,
            datasource_kwargs,
            irfs_to_store,
            None,
            verbose,
        )
    # Runlist mode
    else:
        file_pairs = runlist_to_file_pairs(runlist, event_class_mode, output)
        for st5_str, ea_file, st6_config in file_pairs:

            flist, failed_list, num_event_groups = processFilePair(
                st5_str,
                ea_file,
                flist,
                failed_list,
                output,
                save_multiplicity,
                filename_to_obsid,
                datasource_kwargs,
                irfs_to_store,
                st6_config,
                verbose,
            )

    # Clean up and report files that didn't work

    if gen_index_file and len(flist) > 0:
        gen_index_files(flist, output, eclass_count=num_event_groups)

    logging.info("Processing complete.")
    if len(failed_list) > 0:
        logging.info("V2DL3 was unable to process the following files:")
        for key in failed_list:
            logging.info(key + ": " + str(failed_list[key]))


"""
Generates the index files for a list of files.

Parameters:
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

Parameters:
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
    if stage_idx > -1:
        eclass_fname = (
            fname_base[:stage_idx] + "_ec" + str(eclass_idx) + fname_base[stage_idx:]
        )
    # If no '.' found, append to end
    else:
        eclass_fname = fname_base + "_ec" + str(eclass_idx)

    return os.path.join(output_path, eclass_fname)


"""
Creates `file_pair` tuples of runs and effective area files.

Parameters:
    rl_dict  --  runlist dict made by parseRunlistStrs()

Returns:
    List of tuples (stage5 filename, EffectiveAreaFile)
"""


def runlist_to_file_pairs(runlist, event_class_mode, output):
    from pyV2DL3.vegas.parseSt6RunList import (
        RunlistParsingError,
        RunlistValidationError,
        parseRunlistStrs,
        validateRunlist,
    )

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
            f"{output} already exists as a file. " "<output> needs to be a directory.",
            fg="yellow",
        )
        raise click.Abort()

    # This object imports ROOT, so it should imported after click's CLI is allowed to run
    from pyV2DL3.vegas.EffectiveAreaFile import EffectiveAreaFile

    eas = rl_dict["EA"]
    st5s = rl_dict["RUNLIST"]
    configs = None
    if "CONFIG" in rl_dict.keys():
        configs = rl_dict["CONFIG"]

    file_pairs = []
    for runlist_id in st5s.keys():
        if len(eas[runlist_id]) < 1:
            raise Exception("No EA filenames defined for runlist tag: " + runlist_id)

        ea_files = [EffectiveAreaFile(ea) for ea in eas[runlist_id]]
        if "CONFIG" in rl_dict.keys():
            file_pairs.extend([(st5_file, ea_files, configs[runlist_id]) for st5_file in st5s[runlist_id]])
        else:
            file_pairs.extend(
                [(st5_file, ea_files, configs) for st5_file in st5s[runlist_id]]
            )

    return file_pairs


if __name__ == "__main__":
    cli()
