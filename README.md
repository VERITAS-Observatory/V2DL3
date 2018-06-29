# V2DL3
Simple VERITAS (VEGAS and ED) to DL3 Converter

Contact: Ralph Bird (ralph.bird.1@gmail.com) for details.

## Aim

This repository is for the code that will be used to 


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
