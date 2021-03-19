import numpy as np
import uproot4
import pandas as pd
from pyV2DL3.eventdisplay.IrfInterpolator import IrfInterpolator

def __importer__(etv):
    file = uproot4.open(etv)
    runSummary = file['total_1/stereo/tRunSummary'].arrays(library='np')
    runSummary = pd.DataFrame.from_dict(runSummary)
    runNumber = runSummary['runOn'][0]
    DL3EventTree = file['run_{}/stereo/DL3EventTree'.format(runNumber)].arrays(library='np')
    dataframe = pd.DataFrame.from_dict(DL3EventTree)
    return dataframe, runSummary

def is_sorted(array):
    value = all(array[i] <= array[i + 1] for i in range(len(array) - 1))
    if value == False:
        print("Time Array not in ascending order")
        return (value);
    return (value)

def get_halftime_bin(time_of_day):
    # check if sorted
    is_sorted(time_of_day)
    x = (time_of_day[-1] - time_of_day[0]) / 2
    # get elapsed time since run started
    t_elapsed = time_of_day - time_of_day[0]
    # End if resulting split run is shorter than 60s.
    if t_elapsed.max() < 300:
        return;
    return (next(i for i, y in enumerate(t_elapsed) if y >= x))

def counted(f):
    def wrapped(*args, **kwargs):
        wrapped.calls += 1
        return f(*args, **kwargs)
    wrapped.calls = 0
    return wrapped

# decorator: split_half.calls gives times split_half was called
@counted
def split_half(df):
    # takes single pd.Dataframe or a list of pd.DataFrames and splits according to halftime
    if (type(df) == pd.core.frame.DataFrame):
        idx = get_halftime_bin(df['timeOfDay'].to_numpy())
        splits = [df.iloc[:idx, :], df.iloc[idx:, :]]
        return splits
    elif (type(df) == list and type(df[0]) == pd.core.frame.DataFrame):
        split_list = []
        for frames in df:
            idx = get_halftime_bin(frames['timeOfDay'].to_numpy())
            splits = [frames.iloc[:idx, :], frames.iloc[idx:, :]]
            split_list.extend(splits)
        return split_list


def split_counter():
    # for each time split_half is called, the number of split files is doubled
    count = 2 ** split_half.calls
    if count > 8:
        return ('Reached splitting limit')
    else:
        return count

def check_dict_length():
    [len(x) for x in DL3EventTree.values()]


def systematics_checks(effectiveArea, df, runSummary, **kwargs):

    ##should accept the effective area file and an arbitrary number of split df to check parallel
    irf_file = effectiveArea
    fast_eff_area = uproot4.open(irf_file)['fEffArea']
    
    # run systematics checks for each split df and append boolean value to this list:
    tolerance_split = []

    if (type(df) == pd.core.frame.DataFrame):
        df = [df]
    for df_split in df:

        azimuth, az_max, az_min = df_split['Az'].agg(['mean', 'max', 'min'])
        elevation, el_max, el_min = df_split['El'].agg(['mean', 'max', 'min'])
        zenith, zenith_max, zenith_min = 90 - elevation, 90 - el_max, 90 - el_min
        noise = runSummary['pedvarsOn'][0]

        irf_interpolator = IrfInterpolator(irf_file, azimuth)
        irf_interpolator.set_irf('effNoTh2')

        offset = 0.5  # Get this also from AnasumFile
        eff_area, axis0 = irf_interpolator.interpolate([noise, zenith, offset])
        
        eff_area_00, axis00 = irf_interpolator.interpolate([noise, zenith_min, offset])
        eff_area_01, axis01 = irf_interpolator.interpolate([noise, zenith_max, offset])    

        change_start = (eff_area_00 - eff_area) / eff_area
        change_end = (eff_area_01 - eff_area) / eff_area
        change_ave_s0 = []
        change_ave_e0 = []
        change_ave_s1 = []
        change_ave_e1 = []
        
        #define standard parameters here for check
        Emin, Emax, Tolerance = get_ethreshold_at_zenith(zenith,'moderate'), 30, 0.05 #cut type needs to be read from RunSummary

        # Checking the Effective area changes below and above 1 TeV
        for i in range(axis00[0].size):
            if Emin <= 10.0 ** axis00[0][i] < 1.0:
                change_ave_s0.append(change_start[i])
                change_ave_e0.append(change_end[i])

            if 1.0 <= 10.0 ** axis00[0][i] <= Emax:
                change_ave_s1.append(change_start[i])
                change_ave_e1.append(change_end[i])

        change = np.array(
            [np.mean(change_ave_s0), np.mean(change_ave_e0), np.mean(change_ave_s1), np.mean(change_ave_e1)])
        print('Average change below and above 1 TeV:', change)

        ind = np.abs(change) < Tolerance
        is_all_true = np.all((ind == True))

        tolerance_split.append(is_all_true)
        print(tolerance_split)

    is_all_true = np.all((tolerance_split == True))
    return is_all_true


def __splitter__(effectiveArea, etv):
    dataframe, runSummary = __importer__(etv)
    # repeat systematic checks and split
    while not systematics_checks(effectiveArea, dataframe, runSummary):
        dataframe = split_half(dataframe)
    counts = split_half.calls

    return dataframe, counts


def get_ethreshold_at_zenith(zenith,cut):
    #These are the values from the Figure 2 of ICRC 2015
    #These numbers are kept just for reference
    Zenith = np.array([11.5, 19.0, 28.5, 39.0, 48.7,  58.6])
    Energy_threshold_soft     = np.array([133.63, 140.0, 184.54, 254.54, 458.18, 668.18])
    Energy_threshold_moderate = np.array([222.72, 216.36,260.90, 375.45, 649.09, 1050.0]) 
    Energy_threshold_hard     = np.array([311.81, 318.18,400.90, 598.18, 1030.90, 1501.81])
    
    #Following are the fit parameter for energy-threshold vs zenith with polynomial of order 6
    soft_para = np.array([-3.85861543e-07,  4.14010537e-05, -2.42103742e-04, -1.12552414e-01,
                         5.03169246e+00, -7.91225611e+01,  5.46075976e+02])
    moderate_para = np.array([-2.56357442e-07,  2.79090148e-05, -1.61483213e-04, -7.46857615e-02,
                              3.48114049e+00, -5.83724064e+01,  5.45013373e+02])
    hard_para = np.array([-4.72081854e-07,  4.88485410e-05, -1.99392547e-04, -1.28812793e-01,
                          5.73595442e+00, -9.14165628e+01,  7.95182812e+02])
    
    if cut=='soft':
        poly = np.poly1d(soft_para)
    elif cut=='moderate':    
        poly = np.poly1d(moderate_para)
    elif cut=='hard':
        poly = np.poly1d(hard_para)
    else:
        print ('cut value invalid')
    
    ethreshold = poly(zenith)*1.0e-3
    
    return ethreshold
