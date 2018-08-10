# eventDisplay

## Installation

The process we followed to have a working root + ED + python environment was:

* Install a conda environment with the required versions/packages

`conda create --name=v2dl3 python=2.7 click astropy numpy jupyter ipython matplotlib`

* Install ROOT v5-34-32, making sure it uses the python executable/includes/libraries from the conda environment. 
Also make sure that the packages required by ED are set (TMVA, MINUIT2, etc...).

* Install root_numpy with the activated environment. We used pip for that:

`pip install root_numpy`


## Contacts 

* Tarek Hassan (tarek.hassan@desy.de)
* Mireia Nievas (mireia.nievas-rosillo@desy.de)


