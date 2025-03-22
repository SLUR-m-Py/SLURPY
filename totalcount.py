#!/usr/bin/env python
#SBATCH --job-name=total.count          ## Name of job
#SBATCH --output=%x.%j.out              ## Name stdout
#SBATCH --error=%x.%j.err               ## Name stderr
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=1               ## Number of tasks to be launched
#SBATCH --partition=mpi                 ## Set the partition
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
#######################################
##         FASTP TOTAL COUNT        ##
#######################################
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## ----------------------------------------- LOAD IN MODULES ------------------------------------------ ## 
## Load time and sys mods 
import json, sys 
## Load in current wd 
from os import getcwd
## append path path
sys.path.append(getcwd()) 

## Load in sortg glob ftn 
from defaults import sortglob

## Get the total data counts 
def gettotal(data) -> int: 
    return int(data['summary']['before_filtering']['total_reads'])/2

## Set common name list 
def common_name(str1, str2):
    common = ""
    index2 = 0
    for char1 in str1:
        while index2 < len(str2):
            if char1 == str2[index2]:
                common += char1
                index2 += 1
                break
            index2 += 1
    return common

## Ftn for getting json files
def getjson(inpath:str) -> list:
    return sortglob(f'{inpath}/0.fastp.*.json')

## Set ftn for counting sums
def sumcounts(inpaths:list) -> int:
    ## Set counts 
    counts = []
    ## Iterate over fastp paths 
    for filename_in in inpaths:
        ## Load in json report
        with open(filename_in,'r') as file:
            data = json.load(file)
            file.close()
        ## Calc total pair reads 
        counts.append(gettotal(data))
    ## Sum the counts 
    return sum(counts)

## If the script is envoked 
if __name__ == "__main__":
    ## ------------------------------------------ Input Variables ------------------------------------------ ##
    ## Gather the file name, the start time stamp from input and calculate the end time stamp 
    new_name       = sys.argv[1]
    directory_path = sys.argv[2]

    ## Correct directory path
    directory_path = directory_path[:-1] if (directory_path[-1] == '/') else directory_path

    ## Set sjon paths 
    fastp_paths = getjson(directory_path)

    ## Calc totlas
    total_count = sumcounts(fastp_paths)

    ## Reset new name if none was passed
    new_name = directory_path.split('/diag')[0].split('/')[-1] if ((new_name.lower() =='none') or (new_name is None)) else new_name

    ## Format output file name
    out_file_name = directory_path + f'/{new_name}.counts.csv'
    #print(out_file_name)

    ## Open and write to file, Format new lines to print to file 
    with open(out_file_name,'w') as outf:
        outf.writelines([f'{new_name}, {len(fastp_paths)}, {total_count}\n'])
        outf.close()
## End of file 