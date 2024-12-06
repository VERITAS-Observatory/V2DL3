import sys

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


def main():
    '''
    Main function to execute the script.

    Usage: python query_anasum_runparameters.py <anasum_file> <run_number>
    '''
    if len(sys.argv) != 3:
        print("Usage: python query_anasum_runparameters.py <anasum_file> <run_number>")
        sys.exit(1)

    anasum_file = sys.argv[1]
    run_number = int(sys.argv[2])

    epoch, effective_area = get_epoch_effective_area(anasum_file, run_number)
    print("Epoch:", epoch)
    print("Effective Area:", effective_area)


if __name__ == "__main__":
    main()
