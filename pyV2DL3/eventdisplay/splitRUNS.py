import numpy as np
import uproot4
import pandas as pd

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

#decorator: split_half.calls gives times split_half was called
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
    #for each time split_half is called, the number of split files is doubled
    count = 2**split_half.calls
    if count > 8:
            return ('Reached splitting limit')
    else:
        return count

def check_dict_length():
    [len(x) for x in DL3EventTree.values()]

def systematics_checks(effectiveArea,df,**kwargs):
    #define standard parameters here for check
    Emax, Emin, Tolerance = 0.1, 30, 0.1

    ##should accept the effective area file and an arbitrary number of split df to check parallel

    return #boolean

def __splitter__(dataframe):
    #repeat systematic_checks and split
    while systematics_checks(dataframe):
        dataframe = split_half(dataframe)

    return dataframe
