#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
filterdesc = "Filters an input SAM file for unmapped, placed, and mitochondiral mapping reads."
## ----------------------------- v 0.0.0 ------------------------------ ## 
##   Pre (Hi-C) Filter   
## -------------------------------------------------------------------- ##
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## ----------------------------------- MODULE LOADING ------------------------------------ ##
## Bring in column names of sam file 
from pysamtools import samnames, issam, writeset

## Bring in ftns and options from slurpy
from defaults import mito, M_help

## ------------------------ FUNCTION and VARIABLE FILTERING ------------------------------ ## 
## Set help strings
S_help = "Path to input .sam file from bwa."

## ---------------------------- MAIN SCRIPT & ARGUMENT PARSING -------------------------- ## 
## If the script is envoked
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse

    ## Make the parse
    parser = argparse.ArgumentParser(description = filterdesc)
    ## Add the required arguments
    parser.add_argument("-s", dest="S", type=str,  required=True,  help=S_help, metavar='./path/to/input.sam')
    ## Add optional variables 
    parser.add_argument("-M", dest="M", type=str,  required=False, help=M_help, metavar='n', default=mito)
    ## Set the paresed values as inputs
    inputs = parser.parse_args()

    ## Parse inputs
    inpath, mito = inputs.S, inputs.M

    ## Load in dask dataframe
    import dask.dataframe as dd 

    ## Ftn for loading in sam file with dask 
    def dataloader(inpath):
        """Given an input path, loads a sam file as a dask object"""
        ## Check that this is a sam file 
        assert issam(inpath), "ERROR: The input file is not a .sam file!"
        ## Return the dask obj 
        return dd.read_csv(inpath, sep='\t', comment='@', usecols=range(len(samnames)), names=samnames)

    ## Set the path head
    path_head = inpath.split('.sam')[0]

    ## Load in sam
    thesam = dataloader(inpath)

    ## Define mito mapping file
    mito_path = path_head + f'.{mito}.txt'

    ## Gather the set of mito
    mito_reads = set(thesam[(thesam.Rname==mito) | (thesam.Rnext==mito)].Qname.compute())

    ## Write the mito mapping reads to file
    writeset(mito_path,mito_reads)

    ## Set path to the unmapped reads
    unmp_path = path_head + f'.unmapped.txt'

    ## Gather the set of reads with at least one read unmapped
    unmp_reads = set(thesam[(thesam.Rname=='*') | (thesam.Rnext=='*') | (thesam.Cigar=='*')].Qname.compute())

    ## Write the un-mapped reads set to file
    writeset(unmp_path,unmp_reads)

    ## Filter those names from above and remove from data
    savesam = thesam[~((thesam.Qname.isin(unmp_reads)) | (thesam.Qname.isin(mito_reads)))]

    ## Set name of output path
    save_path = path_head + f'.genomic.txt'

    ## Save the filtered, smaller dataframe
    savesam.to_csv(save_path,header=False,index=False,mode='w',single_file=True,compute=True,sep='\t')
## End of file 