# V2DL3 - VERITAS (VEGAS and Eventdisplay) to DL3 Converter

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://github.com/VERITAS-Observatory/V2DL3/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/138582622.svg)](https://zenodo.org/badge/latestdoi/138582622)
[![PyPI version](https://badge.fury.io/py/v2dl3.svg)](https://badge.fury.io/py/v2dl3)

V2DL3 is a tool to convert [VERITAS](https://veritas.sao.arizona.edu/) data products to DL3 FITS format, allowing to use e.g. the [gammapy science tools](https://gammapy.org/) for  the high-level analysis.

DL3 files include event lists, instrument response functions (IRFs) and observation index tables.
The V2DL3 converter can be used to convert point-like and full-enclosure IRFs.
The FITS output follows the data formats for gamma-ray astronomy as defined in open [gamma-astro-data-formats](https://github.com/open-gamma-ray-astro/gamma-astro-data-formats) (GADF) repository.

The V2DL3 project tries to share as many tools as possible between VEGAS and [Eventdisplay](https://github.com/VERITAS-Observatory/EventDisplay_v4), especially those used for writing the FITS files.

Two steps are required to convert VERITAS data products to DL3 FITS format and use them with gammapy.
Each of these steps are covered by one of the following tools:

- converter of event lists and instrument response functions to DL3 (`v2dl3-vegas` for VEGAS, `v2dl3-eventdisplay` for Eventdisplay)
- `v2dl3-generate-index-file` tool to generate observation index tables

## User installation

Users of V2DL3 can install the package with the following steps:

```bash
python3 -m venv v2dl3
source v2dl3/bin/activate
pip install v2dl3
```

All scripts (`v2dl3-vegas`, `v2dl3-eventdisplay`, `v2dl3-generate-index-file`) are now available as callable command line tools.

## V2DL3 for VEGAS

- VEGAS version >= 2.5.7
- Requirements are listed in the ```environment-vegas.yml``` file.
- Alternatively, a script which builds a Docker image with the latest V2DL3 and the prerequisite software for v2dl3-vegas is available. See *utils/v2dl3-vegas-docker/README.md*

### Developer Installation (VEGAS)

Use the [conda package manager](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) to install the dependencies:

```bash
conda env create -f environment-vegas.yml
```

The environment ```v2dl3-vegas``` will be created and can be activated with:

```bash
conda activate v2dl3-vegas
```

Install now pyV2DL3:

```bash
pip install .
```

### Docker recipe

To use a Docker image with v2dl3-vegas pre-installed, see *utils/v2dl3-vegas-docker/README.md*

### The command line tool v2dl3 with VEGAS

Run `v2dl3-vegas --help` to see all options.

See *README_vegas.md* for more information on newer v2dl3-vegas features such as ITM reconstruction, full-enclosure, and event classes.

Make sure you have ROOT with pyROOT enabled and VEGAS(>=v2.5.7) installed to proceed.
Now, lets create the DL3 fits files from the stage 5 files in the ```./VEGAS/``` folder.

#### One file at a time

To convert a single stage 5 file to DL3 fits you need to provide the path to the stage 5 file as well as the corresponding effective area file using the flag ```-f```. The last argument is the name of the output directory .Beware that the file names for the outputs are inferred from the root file name (xxx.root -> xxx.fits)

```bash
v2dl3-vegas -f ./VEGAS/54809.med.ED.050.St5_Stereo.root ./VEGAS/EA_na21stan_medPoint_050_ED_GRISU.root ./test
```

#### Generate from a VEGAS stage6 run list

You can also provide a stage6 run list to the command line tool. In this case the last argument is the folder where all the output DL3 files will be saved. Beware that the file names for the outputs are inferred from the root file name (xxx.root -> xxx.fits)

```bash
v2dl3-vegas -l ./runlist.txt  ./test
```

Run lists may be generated via a utility script.

```bash
python utils/vegas_runlister.py --help
```

## V2DL3 for EventDisplay

The pip installation as discussed above is recommended for all users.

### Developer Installation (Eventdisplay)

Install dependencies and activate the environment using the [conda package manager](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html):

```bash
conda env create -f environment-eventdisplay.yml
conda activate v2dl3Eventdisplay
pip install -e .
```

### Converting Eventdisplay data products to DL3

Run `v2dl3-eventdisplay --help` to see all options.

Convert an anasum output file to DL3.
The following input is required:

- anasum file for a given run
- effective area file for the corresponding cut applied during the preparation of the anasum file (DL3 version)

Example for point-like analysis:

```bash
v2dl3-eventdisplay \
    -f 54809.anasum.root [Effective Area File] \
     ./output_dir/54809.anasum.fits
```

Example for full-enclosure analysis:

```bash
v2dl3-eventdisplay \
     --full-enclosure \
    -f 64080.anasum.root [Effective Area File] \
     ./output_dir/64080.anasum.fits
```

Runs with observational parameters (i.e., zenith, night sky background) outside but close to corresponding IRF axes range can be converted with the one of the following two command line parameters:

- `--fuzzy_boundary tolerance`: This option interpolates the IRF at the boundary value if the run parameter value is within the given tolerance. The tolerance is define as the ratio of absolute difference between boundary and run parameter value to boundary. This option is preferable over `--force_extrapolation`.
- `--force_extrapolation`: This option extrapolates linearly the IRF at the run parameter value. Use this option with a caution since the extrapolation is applied even for run parameter values very far from the corresponding IRF axes range.

Recommended options is: `--fuzzy_boundary zenith 0.05 --fuzzy_boundary pedvar 0.5`.
This takes into account that extrapolation of the IRF zenith axis is applied to very large zenith angles only, where shower properties changes significantly with small changes in zenith angle.

Further options are:

- `--interpolator_name [KNeighborsRegressor|RegularGridInterpolator]`: select the interpolator used to interpolate the IRFs values. The default is `KNeighborsRegressor`, which allows to robustly interpolate from a irregular grid of IRFS.
- `--save_multiplicity`: write telescope multiplicity as `EVENTTYPE` keyword.
- `--instrument_epoch [epoch]`: write instrument epoch as `INSTRUME` keyword.
- `--db_fits_file [db_fits_file]`: copy run information (e.g., weather, L3 rates) from DB fits file into header information. Requires access to DB fits files.
- `--evt_filter [filter_file]`: apply an event filter defined in a yaml or json file (examples below).

### Eventdisplay CI

Test files for the Eventdisplay CI can be generated using the manual workflow of the [v2dl3-eventdisplay.yml](.github/workflows/v2dl3-eventdisplay.yml) GitHub action.

#### Event filter examples

Selects all events with `IsGamma == 1` (relevant only for anasum files produced with the `WRITE_ALL_EVENTS` set to 1).

```yaml
IsGamma: 1
```

Select all events with energies between 1 and 10 TeV:

```yaml
ENERGY: [1, 10]
```

## Data storage and generating index files

Generate observation index and HDU tables for DL3 data storage are required to use with *gammapy* in for reading and analysis of the generated DL3 data.
This steps is independent of VEGAS or Eventdisplay.
The two index files are generated with the tool `v2dl3-generate-index-file` (check all options with `v2dl3-generate-index-file --help`).

The tables are described on the [GADF website](https://gamma-astro-data-formats.readthedocs.io/en/v0.2/data_storage/index.html):

- [Observation index table](https://gamma-astro-data-formats.readthedocs.io/en/v0.2/data_storage/obs_index/index.html)
- [HDU index table](https://gamma-astro-data-formats.readthedocs.io/en/v0.2/data_storage/hdu_index/index.html)

## Contributing to V2DL3 and Developing Code

Your contribution is welcome!

A few remarks when contributing code:

- goal is to keep as much common code for converting from VEGAS or Eventdisplay data products
- put package specific code into the [pyV2DL3/vegas](pyV2DL3/vegas) and [pyV2DL3/eventdisplay](pyV2DL3/eventdisplay) directories. As different environments are used for both packages, do not put any imports to vegas/eventdisplay in modules in pyV2DL3

To ensure readability, we try follow the Python [PEP8](https://www.python.org/dev/peps/pep-0008/) style guide.
Best to install the [precommit](https://pre-commit.com/) tool to check the code style before committing.

Functions and classes should contain a docstring with a short description.

Unit tests are encouraged and are available for few cases at this point. Unit tests are in the tests directory and can be called using [pytest](http://docs.pytest.org/).

Use the [python logging system](https://docs.python.org/3/howto/logging.html) instead of the ‘print()’ function to output text. This allows to pipe all output into a log file and for different logging levels (INFO, DEBUG, …).
