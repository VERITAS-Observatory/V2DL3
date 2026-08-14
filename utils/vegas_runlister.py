import argparse
import os
from pathlib import Path

parser = argparse.ArgumentParser(
    description="""
    Create a runlist for v2dl3-vegas.

    - The runlist will group all files under the same runlist ID. Make multiple runlists for multiple stage5/EA groupings.
    - Input arguments may be called multiple times to add multiple directories or files.
    - All files will be added as absolute paths unless you provide the `--relative` flag.
    - Example: python vegas_runlister -rd STAGE5_DIR1 -rd STAGE5_DIR2 -ea EA1 -ea EA2
    """,
    formatter_class=argparse.RawTextHelpFormatter,
)

parser.add_argument(
    "output",
    type=str,
    help="Resulting runlist filepath. Enter just the filename to output to working directory.",
)
parser.add_argument(
    "-rd",
    "--run_dir",
    type=str,
    action="append",
    help="Directory of stage5 runs. Adds the path of every .root file in the directory to the runlist.",
)
parser.add_argument(
    "-r",
    "--run_file",
    type=str,
    action="append",
    help="Add the provided stage5 run file to the runlist",
)
parser.add_argument(
    "-ed",
    "--ea_dir",
    action="append",
    type=str,
    help="Directory of effective area files. Adds the path of every .root file in the directory to the runlist.\n"
    "NOTE: multiple EAs in a single runlist ID is for event class mode.",
)
parser.add_argument(
    "-e",
    "--ea_file",
    action="append",
    type=str,
    help="Adds the provided EA file to the runlist.",
)
parser.add_argument(
    "--relative",
    action="store_true",
    help="Add the files' relative paths (default is absolute paths)",
)
parser.add_argument(
    "--no_prompt", action="store_true", help="Answer yes to all confirmation prompts"
)

args = parser.parse_args()

# At least 1 stage5 input is required
if args.run_file is None and args.run_dir is None:
    raise Exception("No stage5 input arguments provided")


def root_files_from_paths(paths, relative=False):
    filepaths = []

    for filename in paths:
        # Confirm file exists
        if not os.path.exists(filename):
            raise Exception(str(filename) + " not found.")

        # Confirm .root extension
        if not filename.lower().endswith(".root"):
            raise Exception(str(filename) + " is not a .root file.")

        # Append filepath
        if relative:
            filepaths.append(os.path.relpath(filename))
        else:
            filepaths.append(Path(filename).absolute())

    return filepaths


def root_files_from_dirs(dirs, relative=False):
    filepaths = []

    for dir in dirs:
        if relative:
            dir = os.path.relpath(dir)
        else:
            dir = Path(dir).absolute()

        root_file_found = False

        # Append .root files in dir
        for filename in os.listdir(dir):
            if filename.lower().endswith(".root"):
                # Append filepath
                filepaths.append(os.path.join(dir, filename))

                root_file_found = True

        # Exit if no .root files in dir
        if not root_file_found:
            raise Exception(str(dir) + " does not contain any .root files.")

    return filepaths


if (
    (args.ea_dir is not None and len(args.ea_dir) > 1)
    or (args.run_dir is not None and len(args.run_dir) > 1)
) and args.run_file is not None:
    raise Exception("Cannot use run files and multiple run groups, use directories.")

stage5_paths = []
if args.run_file is not None:
    stage5_paths += root_files_from_paths(args.run_file, relative=args.relative)

if args.run_dir is not None:
    if len(args.run_dir) == 1:
        stage5_paths += root_files_from_dirs(args.run_dir, relative=args.relative)
    else:
        for dirs in args.run_dir:
            stage5_paths.append(root_files_from_dirs([dirs], relative=args.relative))


# Warn user and prompt confirmation if no EAs are provided
if args.ea_dir is None and args.ea_file is None:
    if not args.no_prompt:
        response = input(
            "No effective area files provided as arguments. Enter 'y' to make runlist without EA files.\n"
        )
        if response != "y" and response != "'y'":
            exit(0)
        print("Continuing without effective area files...")

# Collect effective area filepaths
ea_paths = []
if args.ea_file is not None:
    ea_paths += root_files_from_paths(args.ea_file, relative=args.relative)

if args.ea_dir is not None:
    ea_paths += root_files_from_dirs(args.ea_dir, relative=args.relative)

# Check runlist overwrite
if os.path.exists(args.output):
    if not args.no_prompt:
        response = input(
            str(args.output) + " already exists. Enter 'y' to overwrite.\n"
        )
        if response != "y" and response != "'y'":
            exit(0)

# Check Multiple Run Groups
if len(args.run_dir) != len(ea_paths):
    raise Exception("Mismatch between length of EA and Stage5 directory list. ")

# Write runlist
with open(args.output, "w") as runlist:
    if len(args.run_dir) > 1:
        for group, rundirs in enumerate(args.run_dir):
            # Write stage5 paths
            runlist.write(f"[RUNLIST ID: {group}]\n")
            for path in stage5_paths[group]:
                runlist.write(str(path) + "\n")
            runlist.write(f"[/RUNLIST ID: {group}]\n")

            # Write effective area filepaths
            runlist.write(f"[EA ID: {group}]\n")
            runlist.write(str(ea_paths[group]) + "\n")
            runlist.write(f"[/EA ID: {group}]\n")
    else:
        # Write stage5 paths
        runlist.write("[RUNLIST ID: 0]\n")
        for path in stage5_paths:
            runlist.write(str(path) + "\n")
        runlist.write("[/RUNLIST ID: 0]\n")
        # Write effective area filepaths
        runlist.write("[EA ID: 0]\n")
        for path in ea_paths:
            runlist.write(str(path) + "\n")
        runlist.write("[/EA ID: 0]\n")


print("Runlist written to " + str(args.output))
