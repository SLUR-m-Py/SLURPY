# SLURPY

## Setting up Python Environment

Slurpy was developed using anaconda (python v 3.10.13). 
We recommend using conda to manage the python environment needed by slurpy.
Below are commands needed to set up the "bioenv" for running slurpy. 

```
## Make a python environment named bioenv 
conda create -n bioenv 

## Activate the environment
conda activate bioenv 

## Use conda to install needed mods
conda install numpy pandas matplotlib seaborn pysam samtools bwa macs2
```

If the above installation command hangs, we recommend removing macs2 from the list of libraries and trying again. Then installing macs2 via pip.

```
pip install macs2
```

## Installation
Downloading the repository as a .zip archive is easiest. For developers a clone command with git:

"""
git clone https://github.com/SLUR-m-Py/SLURPY.git
"""

Once slurpy is downloaded (and expanded), change the current directory to the local SLURPY directory and modify the python files as executables. 

"""
cd ./SLURPY
chmod +x *.py 
"""

### Checking environment 
Once within the slrupy directory, run the environment checking script, "modcheck.py".

"""
## Activate the bio env 
conda activate bioenv 
python modcheck.py
"""

If the environment was created successfully the script will run to completion and print the following:

"""
INFO: Modules loaded correctly.
INFO: We are ready to run slurpy!
"""

## Running slurpy


## Dependencies

SLURM 
samtools v 1.15.1 (or greater)

### Common core commands:
cat 
rm
echo 