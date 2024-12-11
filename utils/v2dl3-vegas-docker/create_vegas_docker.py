#!/usr/bin/env python

from __future__ import print_function

import argparse
import os
import subprocess
import sys
from argparse import RawTextHelpFormatter
from shutil import copyfile


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# parse arguments
parser = argparse.ArgumentParser(
    description="""Create and push a VEGAS docker image.

Files will be copied and the docker file will be built in your current directory.
The docker image will be exported and saved in .tar.gz form.

You will need to either clone the VDB and VBF repository into the current directory,
or you will be prompted to download these repositories.

This script assumes that you have ssh keys setup for the VERITAS github.""",
    formatter_class=RawTextHelpFormatter,
)

parser.add_argument(
    "commit",
    type=str,
    help="""[String] The commit of VEGAS to be built in the image. Can be specified by tag (preferred), branch, or hash.""",
)
parser.add_argument(
    "--recipe",
    type=str,
    default="Dockerfile_vegas_v2dl3_x64",
    help="""Dockerfile to use for building. Default is `Dockerfile_vegas_v2dl3_x64`""",
)
parser.add_argument(
    "--numproc",
    type=int,
    default=4,
    help="[Int] Number of processors to use in cmake build",
)
parser.add_argument(
    "--suppress_build_image",
    action="store_true",
    help="Suppress building of the docker image.",
)
parser.add_argument(
    "--export_tgz", action="store_true", help="Export image to .tar.gz."
)
parser.add_argument(
    "--tag",
    type=str,
    help="""[String] Specify the tag to give the docker image.  If not specified, the commit name will be used""",
)
parser.add_argument(
    "--builddir",
    default=os.getcwd(),
    help="""[String] The directory in which to copy files and build image  (default is current directory)""",
)
parser.add_argument(
    "--vegasdir",
    default="vegas",
    help="""[String] The directory in which VEGAS is cloned (default is current directory)""",
)

args = parser.parse_args()

print("Building with {0:d} processors in cmake".format(args.numproc))

# for accessing VERITAS source code
veruser = "veritas"
vw = None

# retrieve VDB/VBF tars
for lib in ["VBF", "VDB"]:
    # git clone VBF/VDB
    if not os.path.isdir(lib):
        cmd = ["git", "clone", "-b", "c++17", f"git@github.com:VERITAS-Observatory/{lib}.git", lib]
        pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        out = pipe.communicate()

    else:
        # try: input = raw_input
        # except NameError: pass
        choice = input(f"{lib} exists. pull updates? y or n\n")
        if choice.lower() == "y":
            cmd = ["git", "-C", lib, "pull", "--rebase"]
            pipe = subprocess.Popen(cmd, stderr=subprocess.PIPE)
            err = pipe.communicate()[0]
            print(err)

# download ROOT
rootTGZ = "root_v6.26.14.source.tar.gz"
rootPath = "https://root.cern.ch/download/"

if not os.path.isfile(rootTGZ):
    # process takes string or list of args
    cmd = ["wget", ("{0:s}/{1:s}".format(rootPath, rootTGZ))]

    # create subprocess, direct output to vars
    pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = pipe.communicate()
    pipe.wait()
    print("out: ")
    print(out)
    print("err: ")
    print(err)
else:
    print("archive file {0:s} already exists".format(rootTGZ))

# git clone VEGAS
if not os.path.isdir(args.vegasdir):
    cmd = ["git", "clone", "git@github.com:VERITAS-Observatory/VEGAS.git", "vegas"]
    pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    out = pipe.communicate()

else:
    # try: input = raw_input
    # except NameError: pass
    choice = input("vegas exists. pull updates? y or n\n")
    if choice.lower() == "y":
        cmd = ["git", "-C", "vegas", "pull", "--rebase"]
        pipe = subprocess.Popen(cmd, stderr=subprocess.PIPE)
        err = pipe.communicate()[0]
        print(err)


cmd = ["git", "-C", "vegas", "checkout", args.commit]
pipe = subprocess.Popen(cmd, stderr=subprocess.PIPE)
err = pipe.communicate()[0]
if err:
    print("there were errors!:")
    print(err)

# produce Dockerfile for build
Dockerfile = "Dockerfile"

dirname = os.path.dirname(__file__)
base_dockerfile = os.path.join(dirname, args.recipe)

with open(Dockerfile, "w") as new_file:
    with open(base_dockerfile) as old_file:
        for line in old_file:
            if line.startswith("RUN cmake --build . -- -j8"):
                new_file.write("RUN cmake --build . -- -j{0:d}\n".format(args.numproc))
            else:
                new_file.write(line)

# create dockerignore
dockerignore = os.path.join(dirname, ".dockerignore")
if os.path.isfile(dockerignore):
    copyfile(dockerignore, ".dockerignore")


# build docker image
commit_name = args.commit.replace("/", "__")
if not args.tag:
    tag = commit_name
else:
    tag = args.tag

fulltag = "vegas:{0:s}".format(tag)
print("Tag name is {0:s}".format(tag))
print("Full tag is {0:s}".format(fulltag))


# check for Docker file
if not os.path.isfile("Dockerfile"):
    print("Dockerfile is missing! must be in a directory with a Dockerfile")
    raise IOError
    sys.exit(1)

if not args.suppress_build_image:
    # first we need to make sure that there is not a tag of the same name
    cmd = ["docker", "images", "-q", ("{0:s}".format(fulltag))]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    imageID = proc.communicate()[0]

    if len(imageID) > 0:
        choice = input(
            ("An image name {0:s} exists. Overwrite? [y/n] : ".format(fulltag))
        )
        if choice.lower() == "y":
            cmd = ["docker", "rmi", ("{0:s}".format(fulltag))]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            ret = proc.communicate()
            if proc.returncode != 0:
                print("failure to delete docker image!!")
                sys.exit(1)
        else:
            print("Since image of name {0:s} already exists".format(fulltag))
            print("and you dont want to overwrite it we are exiting.")
            sys.exit(1)

    # build docker image
    cmd = [
        "docker",
        "build",
        "-t",
        ("{0:s}".format(fulltag)),
        args.builddir,
    ]
    proc = subprocess.Popen(cmd)
    proc.wait()
    out, err = proc.communicate()
    if proc.returncode != 0:
        print("failure to build docker image!!")
        sys.exit(1)

if args.export_tgz:
    cmd = [
        "docker",
        "save",
        "-o",
        ("VEGAS-{0:s}.tar".format(tag)),
        ("{0:s}".format(fulltag)),
    ]
    proc = subprocess.Popen(cmd, stderr=subprocess.PIPE)
    err = proc.communicate()
    proc.wait()
    if proc.returncode != 0:
        print(err)
        eprint("Save to tar failed!")

    cmd = ["gzip", "-f", ("VEGAS-{0:s}.tar".format(tag))]
    proc = subprocess.Popen(cmd, stderr=subprocess.PIPE)
    err = proc.communicate()
    proc.wait()
    if proc.returncode != 0:
        print(err)
        eprint("GZip of tar failed!")
    cmd = ["chmod", "g+x", ("VEGAS-{0:s}.tar.gz".format(tag))]
    proc = subprocess.Popen(cmd, stderr=subprocess.PIPE)
    err = proc.communicate()
