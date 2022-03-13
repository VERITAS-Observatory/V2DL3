from gammapy.data import DataStore
from gammapy.data import HDUIndexTable
from gammapy.data import ObservationTable
import os
from pyV2DL3.generateObsHduIndex import gen_hdu_index
from pyV2DL3.generateObsHduIndex import gen_obs_index


def getDSfromList(flist):
    hdu_index = gen_hdu_index(flist, os.path.realpath("./"))
    obs_index = gen_obs_index(flist, os.path.realpath("./"))
    hdu_index.meta["BASE_DIR"] = os.path.realpath("./")
    return DataStore(
        hdu_table=HDUIndexTable(hdu_index),
        obs_table=ObservationTable(obs_index),
    )
