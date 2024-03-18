import uproot as up


def get_epoch_effective_area(anasum_file, run):
    '''
    Reads and returns the epoch from the mscw log and the effective area file name
    from the anasum log stored within a anasum.root file.

    Parameters
    ----------
    file: str
        path and filename of ANASUM file
    run: int
        VERITAS run number
    Returns
    -------
    Tuple
        epoch and effective area
    '''
    file = up.open(anasum_file)
    data_list = file['anasumLog;1'].members['fLines']._data
    sub = "reading effective areas from"
    string = str([s for s in data_list if sub in s][0])
    string = string.replace(r"reading effective areas from ", '')
    effective_area = string.strip(' \n\t')
    effective_area = effective_area[effective_area.find("effArea"):]
    data_list = file[f'run_{run};1/stereo/mscwTableLog;1'].members['fLines']._data
    sub = "Evaluating instrument epoch"
    string = str([s for s in data_list if sub in s][0])
    string = string.split(sep=",")
    epoch = str(string[1].replace(r" is: ", '').replace(r")", ''))
    return epoch, effective_area
