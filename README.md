# V2DL3
VERITAS (VEGAS and ED) to DL3 Converter

Contact:
	Ralph Bird (ralph.bird.1@gmail.com)
	Tarek Hassan (tarek.hassan@desy.de)
	Tony Lin (xyxlin@gmail.com)
        
## Aim

This repository is for the code that will be used to convert VERITAS data into DL3 format.  Initially it is only being developed for point like irfs

---
# pyV2DL3 

The python package for converting stage5 files to DL3 fits file. Other than useful functions that can be called, a commandline tool `st5ToDL3` comes with the package.

### Requirements

1. vegas version >= 2.5.7
2. pyROOT
3. click
4. astropy
5. numpy
6. gammapy

### Install pyV2DL3

To insatll all the needed python dependencies, use of conda is recomended. The necessary python environment can be created from ```environment.yml```.

Just run:
```
conda env create -f environment.yml
```
and an environment named V2DL3 will be created. After activating the environment just install pyV2DL3 as follows.

```
pip install .
```
### Usage of commandline tool st5ToDL3

```
Usage: st5ToDL3 [OPTIONS] <output>

  Command line tool for converting stage5 file to DL3

  There are two modes:
      1) Single file mode
          When --file_pair is invoked, the path to the stage5 file and the
          corresponding effective area should be provided. The <output> argument
          is then the resulting fits file name.
      2) File list mode
          When using the option --runlist, the path to a stage6 runlist should be used.
          The <output> is then the directory to which the fits files will be saved to.

  Note: One one mode can be used at a time.

Options:
  -f, --file_pair PATH...  A stage5 file (<file 1>) and the corresponding
                           effective area (<file 2>).
  -l, --runlist PATH       Stage6 runlist
  -g, --gen_index_file     Generate hdu and observation index list files. Only
                           have effect in file list mode.
  -m, --save_multiplicity  Save telescope multiplicity into event list
  -e, --ed                 ED mode
  -d, --debug
  -v, --verbose            Print root output
  --help                   Show this message and exit.
```


### Simple Example

Only VEGAS is supported for now, make sure you have ROOT with pyROOT enabled and VEGAS(>=v2.5.7) installed to proceed.
Now, lets create the DL3 fits files from the stage 5 files in the ```./VEGAS/``` folder. 

#### One file at a time

To convert a single stage 5 file to DL3 fits you need to provide the path to the stage 5 file as well as the corresponding effective area file using the flag ```-f```. The last argument is the name of the ouput DL3 file.


```
st5ToDL3 -f ./VEGAS/54809.med.ED.050.St5_Stereo.root ./VEGAS/EA_na21stan_medPoint_050_ED_GRISU.root ./test.fits
```

#### Generate from a VEGAS stage6 runlist

You can also provide a stage6 runlist to the command line tool. In this case the last argument is the folder where all the output DL3 files will be saved. Beware that the file names for the outputs are inferred from the root file name (xxx.root -> xxx.fits)

```
st5ToDL3 -l ./runlist.txt  ./test
```

---
## Git pushing
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
