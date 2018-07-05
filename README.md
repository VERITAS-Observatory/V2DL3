# V2DL3
Simple VERITAS (VEGAS and ED) to DL3 Converter

Contact: Ralph Bird (ralph.bird.1@gmail.com) for details.

## Aim

This repository is for the code that will be used to convert VERITAS data into DL3 format.  Initially it is only being developed for point like irfs


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

---
# pyV2DL3 -- package the code

The python package for converting stage5 files to DL3 fits file. This is basically just a cleared up and packaged version of  the code written in the Jupyter notebook. Other than useful functions that can be called, a commandline tool `st5ToDL3` comes with the package.

### Requirements


1. vegas version >= 2.5.7
2. pyROOT
3. click
4. astropy
5. numpy

### Install pyV2DL3

```
pip install .
```

### Usage of commandline tool st5ToDL3

```
Usage: st5ToDL3 [OPTIONS] <output>

  Command line tool for converting stage5 file to DL3

  There are two modes:
      1) Single file mode
          When --fifle_par is invoked, the path to the stage5 file and the
          corresponding effective area should be provided. The <output> argument
          is then the resulting fits file name.
      2) File list mode
          When using the option --runlist, the path to a stage6 runlist should be used.
          The <output> is then the directory to which the fits while will be saved to.

  Note: One one mode can be used at a time.

Options:
  -f, --file_pair <PATH PATH>...  A stage5 file (<file 1>) and the
                                  corresponding effective area (<file 2>).
  -l, --runlist PATH              Stage6 runlist
  -d, --debug
  -v, --verbose                   Print root output
  --help                          Show this message and exit.
```

*Note: The runlist mode hasn't been fully implemented yet.  
