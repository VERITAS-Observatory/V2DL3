from gammapy.data import HDUIndexTable ,ObservationTable,DataStore
from pyV2DL3.generateObsHduIndex import gen_hdu_index,gen_obs_index
import os


def getDSfromList(flist):
    hdu_index  = gen_hdu_index(flist,os.path.realpath('./'))
    obs_index  = gen_obs_index(flist,os.path.realpath('./'))
    hdu_index.meta['BASE_DIR']  = os.path.realpath('./')
    data_store = DataStore(hdu_table=HDUIndexTable(hdu_index),obs_table=ObservationTable(obs_index))
    return data_store    
