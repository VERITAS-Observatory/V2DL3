import numpy as np
import uproot4
import pandas as pd
from pyV2DL3.eventdisplay.IrfInterpolator import IrfInterpolator

def __importer__(etv):
    file = uproot4.open(etv)
    runSummary = file['total_1/stereo/tRunSummary'].arrays(library='np')
    runNumber = runSummary['runOn'][0]
    DL3EventTree = file['run_{}/stereo/DL3EventTree'.format(runNumber)].arrays(library='np')
    dataframe = pd.DataFrame.from_dict(DL3EventTree)
    return dataframe

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

def systematics_checks(effectiveArea, df, **kwargs):
    # define standard parameters here for check
    Emin, Emax, Tolerance = 0.1, 30, 0.1

    ##should accept the effective area file and an arbitrary number of split df to check parallel
    azimuth = np.mean(df['Az'])
    az_min = np.min(df['Az'])
    az_max = np.max(df['Az'])
    print('azimuth', azimuth, az_min, az_max)

    elevation = np.mean(df['El'])
    ele_min = np.min(df['El'])
    ele_max = np.max(df['El'])
    zenith = 90 - elevation
    zenith_min = 90 - ele_max
    zenith_max = 90 - ele_min
    print('zenith', zenith, zenith_min, zenith_max)

    # tRunSummary = file['total_1/stereo/tRunSummary'].arrays(library='np')
    # noise = np.mean(tRunSummary['pedvarsOn'])
    noise = 7.252038728534893  # Get this also from AnssumFile
    print('noise', noise)

    irf_file = effectiveArea
    fast_eff_area = uproot4.open(irf_file)['fEffArea']
    irf_interpolator = IrfInterpolator(irf_file, azimuth)
    irf_interpolator.set_irf('effNoTh2')

    offset = 0.5  # Get this also from AnasumFile
    eff_area, axis0 = irf_interpolator.interpolate([noise, zenith, offset])

    irf_interpolator_azmin = IrfInterpolator(irf_file, az_min)
    irf_interpolator_azmin.set_irf('effNoTh2')
    eff_area_00, axis00 = irf_interpolator_azmin.interpolate([noise, zenith_min, offset])
    eff_area_01, axis01 = irf_interpolator_azmin.interpolate([noise, zenith_max, offset])

    irf_interpolator_azmax = IrfInterpolator(irf_file, az_max)
    irf_interpolator_azmax.set_irf('effNoTh2')
    eff_area_10, axis10 = irf_interpolator_azmax.interpolate([noise, zenith_min, offset])
    eff_area_11, axis11 = irf_interpolator_azmax.interpolate([noise, zenith_max, offset])

    change_start = (eff_area_00 - eff_area) / eff_area
    change_end = (eff_area_11 - eff_area) / eff_area
    change_ave_s0 = []
    change_ave_e0 = []
    change_ave_s1 = []
    change_ave_e1 = []

    # Checking the Effective area changes below and above 1 TeV
    for i in range(axis00[0].size):
        if Emin <= 10.0 ** axis00[0][i] < 1.0:
            # print (axis00[0][i], change_start[i])
            change_ave_s0.append(change_start[i])
            change_ave_e0.append(change_end[i])

        if 1.0 <= 10.0 ** axis00[0][i] <= Emax:
            # print (axis00[0][i], change_start[i])
            change_ave_s1.append(change_start[i])
            change_ave_e1.append(change_end[i])

    change = np.array([np.mean(change_ave_s0), np.mean(change_ave_e0), np.mean(change_ave_s1), np.mean(change_ave_e1)])
    print('Average change below and above 1 TeV:', change)

    ind = change[np.abs(change) > Tolerance]
    if ind.size == 0:
        return 0
    else:
        return 1

def __splitter__(effectiveArea, etv):
    dataframe = __importer__(etv)
    #repeat systematic checks and split
    while systematics_checks(effectiveArea, dataframe):
        dataframe = split_half(dataframe)
        counts = split_half.calls

    return dataframe, counts
