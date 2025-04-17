#!/usr/bin/env python
#SBATCH --job-name=mergetoshort         ## Name of job
#SBATCH --output=%x.%j.out              ## Name stdout
#SBATCH --error=%x.%j.err               ## Name stderr
#SBATCH --partition=mpi                 ## Size of partitions
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=12              ## Number of tasks to be launched
##################################################
##     Transform Merged_nodups to Short
##################################################
"""
Written by:

    Cullen Roth (croth@lanl.gov)

Genomics and Bioanalytics (B-GEN), Los Alamos National Laboratory, Los Alamos, NM 87545
"""
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      IMPORTING MODS, FUNCTIONS, and VARIABLES
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## Load in pandas and numpy and others 
import pandas as pd
## Bring in exists from os 
from os.path import exists
## from os load in remove
from os import remove

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      SETING DEFAULTS AND HELP MESSAGES 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## Set the description
descrit = "Parse merged_nodups.txt.gz to short_nodups.txt.gz"

## Set the defaults
chunks  = 500000
hicsep  = ' '
mapq    = 30

## Set help messages 
i_help = "Path to an input Hi-C text file."
M_help = "Minimum mapping quality. Default: %s"%mapq
O_help = "The output path to save out result of random downsampling. The default behavior is to save in the same directory as the input file."
Z_help = "Number of rows (default: %s) loaded into pandas at a time. WARNING: while increasing could speed up pipeline it could also cause memeory issues."%chunks
F_help = "Boolean flag to force overwrite of output file." 

## JUICER PARAMETERS
"""
https://github.com/aidenlab/juicer/wiki/Pre#long-format
Juicer short format:

        'str1','chr1','pos1','frag1','str2','chr2','pos2','frag2'

str:  the strand of the contact, 0 for positive orientation anything else negtive
chr:  the chormosome of the contact
pos:  the position of the contact
frag: the fragment, if using dummy var must be different for the pair 
"""
## Set Juicer columns and data types 
juicer_cols  = [ 'Str1','Chr1','Pos1','Frag1','Str2','Chr2','Pos2','Frag2','Mapq1','Cigar1','Seq1','Mapq2','Cigar2','Seq2','Qname1','Qname2','Space']
juicer_types = [   int,  str,   int,    int,   int,    str,  int,   int,    int,     str,    str,    int,    str,    str,    str,      str,    str ]
short_cols   = [ 'Str1','Chr1','Pos1','Frag1','Str2','Chr2','Pos2','Frag2','Mapq1','Mapq2']

## zip the dict
juicer_dict = dict(zip(juicer_cols,juicer_types))

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      MAIN SCRIPT & ARGUMENT PARSING 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## If the script is called 
if __name__ == "__main__":
    ## Bring in argparse, os and set parser
    import argparse

    ## Make the parse
    parser = argparse.ArgumentParser(description = descrit)

    ## Add the required arguments
    parser.add_argument("-i", "--input-path",  dest="i", type=str,   required=True,  help = i_help, metavar = "./path/to/input.txt"    )

    ## Add optional arguments
    parser.add_argument("-M", "--mapq",        dest="M", type=int,   required=False, help = M_help, metavar = 'q',  default=mapq       )
    parser.add_argument("-O", "--output-path", dest="O", type=str,   required=False, help = O_help, metavar = './path/to/saveout.txt', default= None)
    parser.add_argument("-Z", "--chunk-size",  dest="Z", type=int,   required=False, help = Z_help, metavar = 'n',  default = chunks   )
    parser.add_argument("-F", "--force",       dest="F", action = 'store_true',      help = F_help)

    ## Set the parsed values as inputs
    inputs = parser.parse_args()

    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ##      PASS ARGUMENTS & SET VARIABLES 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Format input arguments
    input_path  = inputs.i       ## Input path
    output_path = inputs.O       ## Output path 
    c_size      = inputs.Z       ## Set the chunk size 
    force       = inputs.F       ## Force overwrite 
    mapq        = inputs.M       ## Set mapq

    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ##      INITIALIZATIONS
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Check if the input path exists 
    assert exists(input_path),      "ERROR: A non-existent input path -- %s -- was passed!"%input_path 
    ## Reset output path to default if none was provided
    output_path = output_path if output_path else input_path.replace('merged','short').split('.gz')[0]
    #print(output_path)
    ## Remove output path if it exsits and forcing was passed 
    remove(output_path) if (exists(output_path) and force) else None 
    ## Check that we are not about to overwright 
    assert not exists(output_path), "ERROR: output file -- %s -- already exists!\nINFO: Remove file and try running this script again."%output_path

    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ##      CHUNKING
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Initiate row counts 
    total_rows = 0
    valid_rows = 0
    ## Open the chunks with a with statment
    with pd.read_csv(input_path, sep=hicsep, names=juicer_cols, dtype=juicer_dict, chunksize=c_size, usecols=short_cols, header=None) as datachunks:
        ## Iterate over the data chunks 
        for i,chunk in enumerate(datachunks):
            ## Add to total rows
            total_rows += chunk.shape[0]
            ## Take the valid
            valid = chunk[(chunk.Frag1!=chunk.Frag2) & (chunk.Mapq1>=mapq) & (chunk.Mapq2>=mapq)][short_cols[:-2]]
            ## Cound valid
            valid_rows += valid.shape[0]
            ## Save out valid
            valid.to_csv(output_path,sep=hicsep,index=False,header=False,mode='a' if i else 'w')
        ## Close chunks 
        datachunks.close()
    ## Print counts
    print('Total rows: %s'%total_rows)
    print('Valid rows: %s'%valid_rows)
## End of file     