"""
Query anasum file for run parameters (epoch and effective area file name).

Reads log information stored in anasum.root file and extracts run parameters.

"""
import sys

import uproot as up

from pyV2DL3.eventdisplay.util import get_root_log_lines


def get_epoch_effective_area(anasum_file, run):
    """Return epoch and effective area from an anasum.root file."""
    with up.open(anasum_file) as file:
        try:
            effective_lines = get_root_log_lines(file["anasumLog;1"])
            epoch_lines = get_root_log_lines(
                file[f"run_{run};1/stereo/mscwTableLog;1"]
            )
        except (KeyError, up.exceptions.KeyInFileError, ValueError) as error:
            raise ValueError(
                f"Could not decode Eventdisplay logs in {anasum_file} for run {run}: {error}"
            ) from error

    effective_matches = [
        line for line in effective_lines if "reading effective areas from" in line
    ]
    if len(effective_matches) != 1:
        raise ValueError(
            f"Expected one effective-area line in {anasum_file}, found "
            f"{len(effective_matches)}"
        )
    eff = effective_matches[0].split("reading effective areas from", 1)[1].strip()
    marker = eff.find("effArea")
    if marker < 0:
        raise ValueError(f"Could not parse effective-area path in {anasum_file}")
    effective_area = eff[marker:]

    epoch_matches = [line for line in epoch_lines if "Evaluating instrument epoch" in line]
    if len(epoch_matches) != 1:
        raise ValueError(
            f"Expected one instrument-epoch line in {anasum_file} for run {run}, "
            f"found {len(epoch_matches)}"
        )
    try:
        epoch = epoch_matches[0].split("is:", 1)[1].split(")")[0].strip()
    except IndexError as error:
        raise ValueError(
            f"Could not parse instrument epoch in {anasum_file} for run {run}"
        ) from error

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
