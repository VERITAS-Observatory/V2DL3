"""
Query anasum file for run parameters (epoch and effective area file name).

Reads log information stored in anasum.root file and extracts run parameters.

"""
import sys

import uproot as up


def get_epoch_effective_area(anasum_file, run):
    """Return epoch and effective area from an anasum.root file."""
    f = up.open(anasum_file)

    # Effective area
    lines = f['anasumLog;1'].members['fLines']._data
    eff_line = next((s for s in lines if "reading effective areas from" in s), None)
    if eff_line is None:
        raise ValueError("Could not find 'reading effective areas from' in anasumLog;1")
    eff = eff_line.split("reading effective areas from", 1)[1].strip()
    effective_area = eff[eff.find("effArea"):]

    # Epoch
    lines = f[f'run_{run};1/stereo/mscwTableLog;1'].members['fLines']._data
    epoch_line = next((s for s in lines if "Evaluating instrument epoch" in s), None)
    if epoch_line is None:
        raise ValueError("Could not find 'Evaluating instrument epoch' in mscwTableLog;1")
    epoch = epoch_line.split("is:", 1)[1].split(")")[0].strip()

    return epoch, effective_area


def main():
    """Query anasum file for run parameters (epoch and effective area file name)."""

    if len(sys.argv) != 3:
        print("Usage: v2dl3-eventdisplay-query-runparameters <anasum_file> <run_number>")
        sys.exit(1)

    anasum_file = sys.argv[1]
    run_number = int(sys.argv[2])

    epoch, effective_area = get_epoch_effective_area(anasum_file, run_number)
    print("Epoch:", epoch)
    print("Effective Area:", effective_area)


if __name__ == "__main__":
    main()
