## Function for the SLUR(M)-py pipeline
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      DEFAULT SLURPY PIPELINE FUNCTIONS 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      MODULE LOADING 
## Load in glob
from glob import glob 
## Load in ftns from os. path
from os.path import exists, getsize, basename, splitext
## Bring in sub process
import subprocess, pandas as pd
## Load in remove
from os import remove

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      SLURPY FUNCTIONS 
## Ftn for makign a list of zipped things
def listzip(a,b):
    """Forms a list out of items zipped in a and b."""
    ## Return the list of zipped a and b
    return list(zip(a,b))

## Ftn for making a dicitonary for zipped lists
def dictzip(a,b) -> dict:
    """Makes a dictionary from zipped items in lists a and b."""
    ## Return the zipped dicts
    return dict(zip(a,b))

## Write a ftn for printing given a condition
def ifprint(message,bool):
    """Prints message given the boolean state."""
    ## Print the message
    print(message) if bool else None
    ## Return the boolean
    return bool 
    
## Ftn for returning a sorted glob
def sortglob(wildcard:str) -> list:
    """Retuns the sorted glob of input wild card."""
    ## Returns the sorted glob of input
    return sorted(glob(wildcard))

## Turns an input into a list 
def makelist(input):
    """Transforms an input into a list."""
    ## Returns list 
    return input if (type(input) is list) else [input]

## Write ftn for checking file is a real deal!
def fileexists(filepath:str) -> bool:
    """Checks the existance and size of an input file (foodstuffs)."""
    ## Via os, check if the path to "food" exists and meets our size threshold
    return exists(filepath) and (getsize(filepath)>0)

## Ftn for reseting file
def reset(infiles:list) -> list:
    """Removes the files given with in the intput infiles."""
    ## Return the none type
    return [(remove(f) if exists(f) else None) for f in infiles]

## Ftn for calling sub-process
def submitsbatch(pathtoscript:str,returnid=True) -> int | list[str]:
    """Checks output from sub-process submission of input path-to-script assuming it is an sbatch command. Returns the job id."""
    ## call sub process
    process_lines = subprocess.check_output(pathtoscript, shell=True).decode("utf-8").split('\n')
    ## returns the job id or the lines in process 
    return int(process_lines[0].split(' ')[-1]) if returnid else process_lines

## Ftn for getting basename with no extension
def basenoext(filename:str) -> str:
    """Returns the basename of an input file with no file extension."""
    ## Returns the basename of in .sh file with no .sh
    return splitext(basename(filename))[0]

## Ftn for formating an sbatch text
def sbatch(nameojob:str, cpus:int, cwd:str, report:str, partition=None, nodes=1, tasks=1, runtime='200:00:00',nice=10**7, nodelist=None, memory=None) -> list:
    """Generates the sbatch settings for a script with a give input jobname, cpu count, working directory, and others."""
    ## Gather the extension
    jobext = nameojob.split('.')[-1] if nameojob else 'sh'
    ## Set the runner
    runner = 'python' if (jobext == 'py') else 'bash'
    ## return the formated sbatch txt
    settings = [f'#!/usr/bin/env {runner}\n',                                                 ##      The shebang
                 '#SBATCH --job-name=%s\n'%basenoext(nameojob) if nameojob else '\n',         ##      Name of job 
                 '#SBATCH -o %s\n'%report,                                                    ##      Set the debug dir and report anme 
                 '#SBATCH --nodes=%s\n'%str(nodes),                                           ##      Number of nodes
                 '#SBATCH --ntasks-per-node=%s\n'%str(tasks),                                 ##      Tasks per node
                 '#SBATCH --cpus-per-task=%s\n'%str(cpus),                                    ##      Number of cpus
                 '#SBATCH --time=%s\n'%runtime,                                               ##      Set nice parameter
                 '#SBATCH --nice=%s\n'%str(nice)]                                             ##      The allowed run time of the job 
    ## Add the current working dir                                                            ##  
    settings = settings + ['#SBATCH --chdir=%s\n'%cwd] if cwd else settings                   ##      Set the current workign dir
    ## Add the partition                                                                      ##      
    settings = settings + ['#SBATCH --partition=%s\n'%partition] if partition else settings   ##      The partitions
    ## Add nodes if they were passed
    settings = settings + ['#SBATCH --nodelist=%s\n'%','.join(nodelist)] if nodelist else settings ## List of nodes to run on
    ## Add the max memory
    settings = settings + ['#SBATCH --mem=%s'%memory] if memory else settings
    ## return settings 
    return settings

## Write dict on extensions
extdict = dictzip(['tsv','fai','bed','txt','csv'],['\t','\t','\t',' ',','])

## Ftn set setp
def setsep(ext:str) -> str:
    """Given input extension returns the deliminator."""
    ## Return the extension 
    return extdict[ext]

## Ftn for reading in table 
def readtable(inpath:str,header=None,names=None) -> pd.DataFrame:
    """Reads in a table give its extension."""
    ## Return the table 
    return pd.read_csv(inpath,sep=setsep(inpath.split('.')[-1]),header=header,names=names)

## Ftn for readin an ann file from bwa index
def readann(inpath:str) -> list:
    """Parses an input reference .ann file from bwa index to generate chromosome dataframe of contig names and lengths."""
    ## open and read in the lines
    with open(inpath,'r') as inhandle:
        ## Take up to the third column of each line
        newlines = [l.split(' ')[:3] for l in inhandle]
        inhandle.close()
    ## Gather chromosomes
    contigs = [l[1] for l in newlines[1::2]]
    lengths = [int(l[1]) for l in newlines[2::2]]
    ## Format and return a dataframe
    return pd.DataFrame(list(zip(contigs,lengths)))

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
## EOF 