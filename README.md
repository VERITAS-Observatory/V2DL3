# V2DL3
VERITAS (VEGAS and Eventdisplay) to DL3 Converter

Contact:
	Ralph Bird (ralph.bird.1@gmail.com)
	Tarek Hassan (tarek.hassan@desy.de)
	Tony Lin (xyxlin@gmail.com)
	Tobias Kleiner (tobias.kleiner@desy.de)
        
## Aim

This repository is for the code that will be used to convert VERITAS data into DL3 format. It can be used to convert point-like IRFs (as included in the joint-Crab paper), as well as full-enclosure IRFs.

The project follows the most recent DL3 format definition from the [open gamma-ray astro data formats repository](https://github.com/open-gamma-ray-astro/gamma-astro-data-formats).

---
# pyV2DL3 

The python package for converting stage5/anasum files to the DL3 FITS format. Other than useful functions that can be called, a commandline tool `v2dl3` comes with the package.

### Requirements
The requirements are listed in the ```environment.yml``` file.

#### VEGAS

* vegas version >= 2.5.7

#### EventDisplay

* The converter does not depend on EventDisplay. However, make sure that the EventDisplay ansum stage runs with version >= v485 to include the DL3EventTree in the ansum file.

### Install pyV2DL3

To install the needed python dependencies, use of conda is recommended. The necessary python environment can be created from ```environment.yml```.

Just run:
```
conda env create -f environment.yml
```
and the environment ```v2dl3``` will be created. After activating the environment (`conda activate v2dl3`), install pyV2DL3 as follows:

```
pip install .
```
### Usage of commandline tool v2dl3

```
Usage: v2dl3 [OPTIONS] <output>

  Command line tool for converting stage5/anasum file to DL3

  There are two modes:
      1) Single file mode
          When --file_pair is invoked, the path to the stage5/anasum file and the
          corresponding effective area should be provided. The <output> argument
          is then the resulting fits file name.
      2) File list mode
          When using the option --runlist, the path to a stage6/anasum runlist should be used.
          The <output> is then the directory to which the fits files will be saved to.

  Note: One one mode can be used at a time.

Options:
  -f, --file_pair PATH...  A stage5/anasum file (<file 1>) and the corresponding
                           effective area (<file 2>).
  -l, --runlist PATH       Stage6/anasum runlist
  -g, --gen_index_file     Generate hdu and observation index list files. Only
                           have effect in file list mode.
  -m, --save_multiplicity  Save telescope multiplicity into event list
  -e, --ed                 Eventdisplay mode
  -d, --debug
  -v, --verbose            Print root output
  --help                   Show this message and exit.
```

---
### Examples

### VEGAS

Make sure you have ROOT with pyROOT enabled and VEGAS(>=v2.5.7) installed to proceed.
Now, lets create the DL3 fits files from the stage 5 files in the ```./VEGAS/``` folder. 

##### One file at a time

To convert a single stage 5 file to DL3 fits you need to provide the path to the stage 5 file as well as the corresponding effective area file using the flag ```-f```. The last argument is the name of the ouput DL3 file.


```
v2dl3 -f ./VEGAS/54809.med.ED.050.St5_Stereo.root ./VEGAS/EA_na21stan_medPoint_050_ED_GRISU.root ./test.fits
```

##### Generate from a VEGAS stage6 runlist

You can also provide a stage6 runlist to the command line tool. In this case the last argument is the folder where all the output DL3 files will be saved. Beware that the file names for the outputs are inferred from the root file name (xxx.root -> xxx.fits)

```
v2dl3 -l ./runlist.txt  ./test
```

### EventDisplay

ROOT is installed directly into the environment from conda-forge following the above steps using the ```environment.yml```.
Now, lets create the DL3 fits files from the anasum files in the ```./eventDisplay/``` folder. 

##### One file at a time

To convert an anasum file to DL3 you need to provide the path to the anasum file as well as the corresponding effective area file using the flag ```-f```. The only difference with VEGAS is that you need to add the '--ed' flag. The last argument is the name of the ouput DL3 file.


```
v2dl3 --ed -f ./eventDisplay/54809.anasum.root [Effective Area File] ./eventDisplay/54809.anasum.fits
```

##### Full-enclosure

For full-enclosure IRFs you need to pass the additional flag --full-enclosure and be sure to provide the proper effective area files:
```
v2dl3 -ed --full-enclosure -f 64080.anasum.root $VERITAS_EVNDISP_AUX_DIR/EffectiveAreas/effArea-v485-auxv01-CARE_June2020-Cut-NTel2-PointSource-Hard-TMVA-BDT-GEO-V6_2012_2013a-ATM62-T1234.root ./FITS/64080.anasum.fits
```

##### Multi file processing

To convert many runs at once with different Effective Area files there is a modified anasum script here ( ``` VERITAS-Observatory/Eventdisplay_AnalysisScripts_VTS/scripts/ANALYSIS.anasum_parallel_from_runlist_v2dl3.sh ``` ), that can be used to create a ``` v2dl3_for_runlist_from_ED485-anasum.sh ``` script. This script then contains one line for each processed file in the formatting as shown above in the full-enclosure case. 
Then in your bash run 
```
./v2dl3_for_runlist_from_ED485-anasum.sh
```
to create the fits files one after another. Alternatively you can use this script to submit each line to your batch farm to run in parallel.

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
