# V2DL3 - VERITAS (VEGAS and Eventdisplay) to DL3 Converter.

V2DL3 is a tool to convert VERITAS data products to DL3 FITS format, allowing to use e.g. the [Gammapy science tools](https://gammapy.org/) for analysis. 

The converter can be used to convert point-like and full-enclosure IRFs. 
The FITS output follows the data formats for gamma-ray astronomy as defined in open [gamma-astro-data-formats](https://github.com/open-gamma-ray-astro/gamma-astro-data-formats) (GADF) repository.

The projects tries to share as many tools as possible between VEGAS and Eventdisplay, especially those used for writing the FITS files.

The two main tools required to convert VERITAS data products to DL3 FITS format and use them with gammapy are:
- converter to DL3 (`v2dl3-vegas` for VEGAS, `v2dl3_for_Eventdisplay.py` for Eventdisplay)
- tool to generate observation index tables

For contributors: please note the section for developers below.

---
# V2DL3 for VEGAS

* vegas version >= 2.5.7
* requirements are listed in the ```environment-vegas.yml``` file.

## Installation

Use the [conda package manager](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) to install the dependenies:
```
conda env create -f environment-vegas.yml
```
The environment ```v2dl3-vegas``` will be created and can be activated with:

```
conda activate v2dl3-vegas
```

Install now pyV2DL3:
```
pip install .
```

## The commandline tool v2dl3 with VEGAS

Run `v2dl3-vegas --help` to see all options.

Make sure you have ROOT with pyROOT enabled and VEGAS(>=v2.5.7) installed to proceed.
Now, lets create the DL3 fits files from the stage 5 files in the ```./VEGAS/``` folder. 

### One file at a time

To convert a single stage 5 file to DL3 fits you need to provide the path to the stage 5 file as well as the corresponding effective area file using the flag ```-f```. The last argument is the name of the ouput DL3 file.

```
v2dl3-vegas -f ./VEGAS/54809.med.ED.050.St5_Stereo.root ./VEGAS/EA_na21stan_medPoint_050_ED_GRISU.root ./test.fits
```

### Generate from a VEGAS stage6 runlist

You can also provide a stage6 runlist to the command line tool. In this case the last argument is the folder where all the output DL3 files will be saved. Beware that the file names for the outputs are inferred from the root file name (xxx.root -> xxx.fits)

```
v2dl3-vegas -l ./runlist.txt  ./test
```

---

# V2DL3 for EventDisplay

- use Eventdisplay version >= 487
- recipes for Docker containers are available from the [V2DL3-Docker-recipes](https://github.com/Eventdisplay/V2DL3-Docker-recipes) repository.

## Installation

Use the [conda package manager](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) to install the dependenies:
```
conda env create -f environment-eventdisplay.yml
```

Activate the environment ```v2dl3ED``` and set `PYTHONPATH`:

```
conda activate v2dl3ED
export PYTHONPATH=$PYTHONPATH:"${PWD}"
```

Note that no pip is required for using the v2dl3 tool with Eventdisplay.

## Usage of commandline tool v2dl3

Run `python pyV2DL3/script/v2dl3_for_Eventdisplay.py --help` to see all options.

Convert an anasum output file to DL3.
The following input is required:
- anasum file for a given run
- effective area file for the corresponding cut applied during the preparation of the anasum file (DL3 version)

Example for point-like analysis:
```
python pyV2DL3/script/v2dl3_for_Eventdisplay.py -f 54809.anasum.root [Effective Area File] ./outputdir/54809.anasum.fits
```
Example for full-enclosure analysis:
```
python pyV2DL3/script/v2dl3_for_Eventdisplay.py --full-enclosure -f 64080.anasum.root [Effective Area File] ./outputdir/64080.anasum.fits
```

The run having their observational parameters (zenith, night sky background) outside but close to corresponding IRF axes range can be run with the one of the following two commandline parameters: 

- `--force_extrapolation`: This option extrapolates linearly the IRF at the run parameter value. Use this option with a caution since the exptrapolation happens even for run parameter values very far from the corresponding IRF axes range.

- `--fuzzy_boundary tolerance`: This option interpolates the IRF at the boundary value if the run parameter value is within the given tolerance. The tolerance is define as the ratio of absolute difference between boundary and run parameter value to boundary. This option is preferable over `--force_extrapolation`. 

---
# Data storage and generating index files

Two index files are required for DL3-type analysis and can be generated with the tool `generate_index_file.py`.

The tables are descriped on the [GADF website](https://gamma-astro-data-formats.readthedocs.io/en/v0.2/data_storage/index.html):
- [Observation index table](https://gamma-astro-data-formats.readthedocs.io/en/v0.2/data_storage/obs_index/index.html)
- [HDU index table](https://gamma-astro-data-formats.readthedocs.io/en/v0.2/data_storage/hdu_index/index.html)

To use `generate_index_file.py`, run:
- `generate_index_file --help` when using VEGAS
- `python pyV2DL3/script/generate_index_file.py --help` when using Eventdisplay 

---

# Contributing and Developing Code

A few remarks when contributing code:
- goal is to keep as much common code for converting from VEGAS or Eventdisplay data products
- put package specific code into the [pyV2DL3/vegas](pyV2DL3/vegas) and [pyV2DL3/eventdisplay](pyV2DL3/eventdisplay) directories. As different environments are used for both packages, do not put any imports to vegas/eventdisplay in modules in pyV2DL3

To ensure readability, we try follow the Python [PEP8](https://www.python.org/dev/peps/pep-0008/) style guide. 

Functions and classes should contain a docstring with a short description.

Unit tests are encouraged and are available for few cases at this point. Unit tests are in the tests directory and can be called using [pytest](http://docs.pytest.org/). 

Use the [python logging system](https://docs.python.org/3/howto/logging.html) instead of the ‘print()’ function to output text. This allows to pipe all output into a log file and for different logging levels (INFO, DEBUG, …).

---

##### Multi file processing

To convert many runs at once with different Effective Area files there is a anasum script [ANALYSIS.anasum_parallel_from_runlist.sh](https://github.com/VERITAS-Observatory/Eventdisplay_AnalysisScripts_VTS/blob/main/scripts/ANALYSIS.anasum_parallel_from_runlist.sh), that can be used to create a ``` v2dl3_for_runlist_from_EDxxxx-anasum.sh ``` script. This script then contains one line for each processed file in the formatting as shown above in the point-like case. Here, xxxx is the Eventdisply version (for eg. v487).

Then in your bash run 
```
./v2dl3_for_runlist_from_EDxxxx-anasum.sh
```
to create the fits files one after another. 

Alternatively, you can submit one job for each entry of this script using
```bash
v2dl3_qsub --v2dl3_script <script>
```
where `<script>` is the script that was written out by `ANALYSIS.anasum_parallel_from_runlist_v2dl3.sh`.
`v2dl3_qsub` has the following options:
 - `--conda_env` name of the conda environment. Default is `v2dl3` 
 - `--conda_exe` path to the conda executable. Only needed if `$CONDA_EXE` is not set.
 - `--rootsys` path to rootsys. Only needed if `$ROOTSYS` is not set
 - `--add_option` allows to add further options to v2dl3. (e.g. `--add_option '--evt_filter /path/to/file.yaml'`)

---
**TEXT BELOW REQUIRES REVIEW**

#### Filter events
Using --evt_filter option, you can filter which events are written to the fits file. The argument takes the path of a 
yaml file that stores conditions. E.g. to select only events between 0.5 and 1.0 TeV:
```yaml
Energy: [0.5, 1.0]
```

---
### Git pushing
If two people have used the same notebook at the same time it gets a bit nasty with a merge due to differences in outputs and cell run counts.  To overcome this I have followed the instructions in http://timstaley.co.uk/posts/making-git-and-jupyter-notebooks-play-nice/

Specifically, this requires that you have jq (https://stedolan.github.io/jq/) which should be easy enough to install, I get it through brew on my mac (I'm not sure what happens if you don't have jq installed - maybe you will find out!).

The following files have been edited to allow for this (in this repository directory, if you want you can set some of this up globally).

.git/config
```
[core]
attributesfile = .gitattributes

[filter "nbstrip_full"]
clean = "jq --indent 1 \
        '(.cells[] | select(has(\"outputs\")) | .outputs) = []  \
        | (.cells[] | select(has(\"execution_count\")) | .execution_count) = null  \
        | .metadata = {\"language_info\": {\"name\": \"python\", \"pygments_lexer\": \"ipython3\"}} \
        | .cells[].metadata = {} \
        '"
smudge = cat
required = true
```

.gitattributes
```
*.ipynb filter=nbstrip_full
```
