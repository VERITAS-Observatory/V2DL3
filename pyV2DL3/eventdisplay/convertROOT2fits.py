import os
import glob
from pathlib import Path


def convertROOT2fits(files, eff, **kwargs):
    from pyV2DL3.genHDUList import genHDUlist, loadROOTFiles
    from pyV2DL3 import generateObsHduIndex

    if type(files) == str:
        files = [files]
    else:
        files = glob.glob(files + "/*anasum.root")

    full_enclosure = kwargs.pop("full_enclosure", True)
    point_like = kwargs.pop("point_like", True)
    instrument_epoch = kwargs.pop("instrument_epoch", None)
    save_multiplicity = kwargs.pop("save_multiplicity", False)
    filename_to_obsid = kwargs.pop("filename_to_obsid", True)
    evt_filter = kwargs.pop("evt_filter", None)

    if evt_filter is not None:
        evt_filter = Path(evt_filter)

    force_extrapolation = kwargs.get("force_extrapolation", False)
    fuzzy_boundary = kwargs.get("fuzzy_boundary", 0.0)

    if not (full_enclosure) and not (point_like):
        point_like = True
        full_enclosure = False

    irfs_to_store = {"full-enclosure": full_enclosure, "point-like": point_like}

    for file in files:
        datasource = loadROOTFiles(Path(file), Path(eff), "ED")
        datasource.set_irfs_to_store(irfs_to_store)

        datasource.fill_data(
            evt_filter=evt_filter,
            use_click=False,
            force_extrapolation=force_extrapolation,
            fuzzy_boundary=fuzzy_boundary,
            **kwargs)

        hdulist = genHDUlist(
            datasource,
            save_multiplicity=save_multiplicity,
            instrument_epoch=instrument_epoch,
        )

        fname_base = os.path.basename(file)
        obs_id = int(fname_base.split(".")[0])

        if filename_to_obsid:
            hdulist[1].header["OBS_ID"] = obs_id

        output = kwargs.get("output", file.replace(".root", ".fits"))

        if ".fits" not in output:
            output += ".fits"

        hdulist.writeto(output, overwrite=True)

    datadir = str(Path(file).absolute().parent)
    filelist = glob.glob(f"{datadir}/*anasum.fit*")

    generateObsHduIndex.create_obs_hdu_index_file(filelist, index_file_dir=datadir)
