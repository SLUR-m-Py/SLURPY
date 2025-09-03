#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## --------------------------------- v 1.0.0 ------------------------------------ ## 
##   MYECHO: A script for checking mods needed for SLUR(M)-py and echoing statments ... statments
## ------------------------------------------------------------------------------ ##
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## Import sys 
import sys

## Ftn for commenting out command lines for debuging
def debuglines(intxt:list) -> list: 
    ## Initilizse new lines and counter
    newlines = []
    ## Iterate thru the input txt lines 
    for l in intxt:
        ## If the first chracter is already a comment like #SBATCH
        if (l[0] == '#'): 
            newlines.append(l)
        ## If it is an echo statment, leave it as is 
        elif (l.split(' ')[0]=='echo'):
            newlines.append('sleep 10\n'+l)
        elif (l.split(' ')[0]=='myecho.py'):
            newlines.append('sleep 10\n'+l)
        else: ## Othewise, comment out the lines 
            newlines.append('##'+l)
    ## Return the commented out lines 
    return newlines 

## Ftn to write to file
def writetofile(inpath:str,intxt:list,debug:bool,mode='w') -> str:
    """Opens a file to write lines to file."""
    ## Modify the input text lines if in debug mode
    intxt = debuglines(intxt) if debug else intxt
    ## Open the input path and write out to file 
    with open(inpath,mode) as ofile:
        ofile.writelines(intxt)
    ## Return the path
    return inpath

## List the needed mods
"""
numpy
pandas
matplotlib
seaborn 
"""
## List of standard libraries with python 3
"""
glob
os
subprocess
argparse
time
sys
json
multiprocessing
datetime
biopython
dask
"""
## Needed stand alones 
"""
samtools
bwa
macs3
"""
## Conda install command 
"""
conda install numpy pandas matplotlib seaborn ipython pysam samtools macs3 bwa 
"""
## If the conda install dosn't work for macs2, try:
"""
pip install macs3
"""
## ----------------------------------------------- MAIN EXECUTABLE --------------------------------------------------- ## 
## If the library is called as an executable
if __name__ == "__main__":
    ## If we are doing a check
    if (len(sys.argv) > 1):
        ## Set the output message and output file 
        outmessage, outfile = ' '.join(sys.argv[1:-1]), sys.argv[-1]
        ## Write to file
        writetofile(outfile,outmessage,debug=False,mode='a')
    else: ## Otherwise check the mods
        ## Load in date and time
        from datetime import datetime
        ## Load in glob
        from glob import glob 
        ## Load in SeqIO
        from Bio import SeqIO
        ## Bring in basename, get file size, and path-exists 
        from os.path import isfile, basename, exists, getsize 
        ## Bring in make dirs
        from os import makedirs, remove, getcwd
        ## bring in matplot lib 
        from matplotlib import pyplot as plt
        ## Load pandas, numpy and seaborn 
        import pandas as pd, numpy as np, seaborn as sns 
        ## Load in system, json, time, argparese and subprocess 
        import sys, json, time, argparse, subprocess
        ## Loadin checksam
        from biotools import checksam
        ## Bring in dask 
        import dask.dataframe as dd 
        ## Check samtools
        assert checksam(), 'ERROR: The detected version of samtools is not greater than or equal to v 1.15.1!\nPlease update samtools and try again.\n'
        ## Print to screen that this ran
        print("INFO: Modules loaded correctly.\nINFO: We are ready to run slurpy!")
## End of file 