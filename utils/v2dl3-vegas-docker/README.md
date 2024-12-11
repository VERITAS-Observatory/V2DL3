The contents of this directory are modified and streamlined from the VEGAS docker shifter, optimized for v2dl3-vegas CI testing. It is also suitable as a means to run v2dl3-vegas yourself without the heavy prerequisite software setup.

It is recommended to install gammapy somewhere else which allows for the latest versions of gammapy and python.

# Recipes

#### Dockerfile_vegas_v2dl3_x64
- Ubuntu 18.04
- ROOT 6.13.08
- VBF 0.3.4
- VDB 4.3.2
- VEGAS version chosen by user (script arg)
- latest V2DL3 main branch installed to base environment

#### Dockerfile_vegas_v2dl3_arm64
- Same as above but for ARM architecture

#### Dockerfile_vegas_ci
- Image used by v2dl3-vegas.yml CI workflow
- Does not install V2DL3 (is installed by the workflow)

# Instructions for Building

This assumes that you have a git ssh key set up to download VEGAS from the VERITAS repository.

Make and use a build sub-directory
```Bash
 mkdir BuildDir
 cd BuildDir
```

You may copy VBF and VDB into the above directory. If not, it will download them using wget after prompting you for the password.

Run the Docker build script where <branchName> is the name of the branch you want to build (e.g. `v2_5_9_rc1`, which the original v2dl3-vegas CI image uses)
```Bash
 python /path/to/create_vegas_docker.py <branchName>
```
Options that can be passed to create_vegas_docker.py
```Bash

branchName [String] The commit of VEGAS to be built in the image. Can be specified by tag (preferred), branch, or hash.

  -h, --help            show this help message and exit

  --recipe              Dockerfile to use for building. Default is `Dockerfile_vegas_v2dl3_x64`

  --numproc NUMPROC     [Int] Number of processors to use in cmake build

  --suppress_build_Image
                        Suppress building of the docker image.

  --export_tgz          Export image to .tar.gz

  --tag TAG             [String] Specify the tag to give the docker image.  If not specified, the commit name will be used

  --builddir BUILDDIR   [String] The directory in which to copy files and build image  (default is current directory)
```
