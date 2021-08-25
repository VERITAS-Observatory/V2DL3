"""
Submit v2dl3 EventDisplay jobs to batch farm. Uses input from script that contains induvidual commands per line.
The required script is created using the script ANALYSIS.anasum_parallel_from_runlist_v2dl3.sh.
"""

import os
import subprocess
import click

import pyV2DL3
from datetime import date

# get the directory of the run_loop.py file
dirname = os.path.dirname(pyV2DL3.__file__)
# get today's data for the log dir
today = date.today().strftime("%y%m%d")

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--v2dl3_script',
              type=click.Path(exists=True),
              help='Script that contains individual commands per line that can also be used'
                   ' in stand alone mode with the v2dl3 converter'
              )
@click.option('--conda_env',
              default="v2dl3",
              help='Name of the conda environment. (Default: "v2dl3")'
              )
@click.option('--conda_exe',
              default="$CONDA_EXE",
              help='Path to the conda executable. Change if required.'
              )
@click.option('--rootsys',
              default="$ROOTSYS",
              help='Path of the ROOTSYS. Change if required.'
              )
@click.option('--add_option',
              default="",
              help='Option to add when running v2dl3.'
              )
def cli(v2dl3_script, conda_env, conda_exe, rootsys, add_option):
    with open(v2dl3_script, "r") as script:
        commands = []
        for line in script.read().splitlines():
            if line.startswith("v2dl3"):
                commands.append(line)
    if len(commands) == 0:
        raise ValueError("No commands found.")

    conda_exe = os.path.expandvars(conda_exe)
    rootsys = os.path.expandvars(rootsys)
    logdir = os.path.expandvars(f"$VERITAS_USER_LOG_DIR/v2dl3/{today}")
    qsub_cmd = f"qsub -t 1"

    for ii, command in enumerate(commands):
        runnumber = command.split("/")[-1].split(".")[0]

        if add_option:
            command_split = command.split(None)
            command_split.insert(1, add_option)
            command = " ".join(command_split)

        logs = f"-j y -o {logdir}/{runnumber}.log"
        sub_cmd = f"{qsub_cmd} {logs} {dirname}/script/helper/qsub_convert.sh " \
                  f"'{command}' '{conda_exe}' '{conda_env}' '{rootsys}'"

        subprocess.call(sub_cmd, shell=True)

if __name__ == '__main__':
    cli()