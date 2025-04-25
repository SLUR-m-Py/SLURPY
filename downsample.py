#!/usr/bin/env python
#SBATCH --job-name=downsample           ## Name of job
#SBATCH --output=%x.%j.out              ## Name stdout
#SBATCH --error=%x.%j.err               ## Name stderr
#SBATCH --partition=mpi                 ## Size of partitions
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=12              ## Number of tasks to be launched
##################################################
##          Hi-C Random Down-sampling
##################################################
"""
Written by:

    Sasha Bacot (sbacot@lanl.gov)

Edits by: 

    Cullen Roth (croth@lanl.gov)

Genomics and Bioanalytics (B-GEN), Los Alamos National Laboratory, Los Alamos, NM 87545
"""
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      IMPORTING MODS, FUNCTIONS, and VARIABLES
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## Load in pandas and numpy and others 
import pandas as pd, numpy as np, gzip 
## Bring in exists from os 
from os.path import exists
## from os load in remove
from os import remove

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      SETING DEFAULTS AND HELP MESSAGES 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## Set the description
downsample_desc = "Randomly downsamples rows from an input .txt file representing Hi-C contacts."

## Set the defaults
ranseed = 720
chunks  = 500000
ST      = 'store_true'

## Set help messages 
i_help = "Path to an input Hi-C text file."
n_help = "Number (or percent) of rows to randomly downsample. Float values less than 1 are interpreted as a fraction of rows (Hi-C contacts) to downsample."
O_help = "The output path to save out result of random downsampling. The default behavior is to save in the same directory as the input file."
R_help = "Set integer used as random seed in numpy (default: %s)"%ranseed
Z_help = "Number of rows (default: %s) loaded into pandas at a time. WARNING: while increasing could speed up pipeline it could also cause memeory issues."%chunks
F_help = "Boolean flag to force overwrite of output file." 
I_help = "Samples only the intra-chromosome contacts."
S_hlep = "Boolean flag indicating input is in short format."

## Set some default input variables used in these functions 
## Such as the input column names, file seperator, header, and replacement value 
hiccols     = ['Strand1', 'Chr1', 'Pos1', 'Dum1', 'Strand2', 'Chr2', 'Pos2', 'Dum2']
hicsep      = ' '   ## space deliminator 
hicheader   = None  ## No header 
replacement = False ## Do not change this, sampling with replacement is not supported in this script 
debugc      = False ## flag debug 
verbose     = False ## flag for verbosity setting 

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      FUNCTION DEFINING
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## Ftn for gzipped file 
def isgzip(inpath) -> bool: 
    ## Return the check 
    return inpath.split('.')[-1] == 'gz'

## Ftn for counting 
def countlines(inpath) -> list:
    ## Iniate counters 
    intra, inter_ix = 0, []
    ## open the ifle
    with (gzip.open(inpath,'rt') if isgzip(inpath) else open(inpath,'rt')) as infile:
        ## Iterate
        for i,l in enumerate(infile):
            line = l.split(hicsep)
            ## If is an inter chrom contact add to inter
            if line[1] != line[5]:
                inter_ix.append(i)
                if (len(inter_ix) < debugc) and verbose:
                    print(line)
                    print(line[1],line[5])
            else: ## Otherwise add to intra 
                intra += 1
    ## Return coutns
    return intra, inter_ix, i+1

## Ftn for reseting sample size
def getsamplesize(nsample:float,nlines:int) -> int:
    """Calculates the number of rows to subsample from given the total number of lines.""" 
    ## Reset the number of lines to randomly downsample if a fraction/floating point number was given as input
    return int(np.ceil(nlines * nsample) if (nsample < 1) else nsample)

## Ftn for precalculating random index
def getrandom(nlines:int,samplesize:int,intraonly=np.array([])) -> np.ndarray[int]:
    """Calculates an array of random integers given the total number of lines (nlines) and the sample size."""
    ## If intra only was passed remove from full index and return random numpy array from nlines, given the sample size 
    return np.random.choice(np.setdiff1d(np.arange(nlines),intraonly), samplesize, replace = replacement)

## Print counts 
def printcounts(infile,total,intras,inters) -> None:
    print(f'File: {infile}, Total: {total}, Intra-contacts: {intras}, Inter-contacts: {inters}')
    pass 

## Define the main function to write out sub sampled input file 
def writesubsample(infile:str, outfile:str, nsample:float, sampleintra:bool, csize:int, saveheader=False, saveindex=False) -> bool: 
    """Loads in Hi-C contacts from an input .txt file (via chunking) and randomly downsamples a set of rows, saving to an output file."""
    ## Count the total number of lines in file
    intra_count, inter_ix, linecount = countlines(infile)
    ## Print the counts
    printcounts(infile,linecount,intra_count,len(inter_ix))
    ## Calculate the smaple size and create  the random indices
    indices = getrandom(linecount,getsamplesize(nsample,linecount),intraonly=inter_ix) if sampleintra else getrandom(linecount,getsamplesize(nsample,linecount)) 

    ## Open the chunks with a with statment
    with pd.read_csv(infile, sep = hicsep, header = hicheader, names = hiccols, chunksize=csize) as datachunks:
        ## Iterate over the data chunks 
        for i,chunk in enumerate(datachunks):

            ## Gather the index of the chunk that are within the randomly selected list
            boolindex = chunk.index.isin(indices)
            # Write out the sub smaple of the chunk if there is a true value in the boolean index, ie. sum > 0
            chunk[boolindex].to_csv(outfile, sep=hicsep, index=saveindex, header=saveheader, mode='a' if i else 'w') if sum(boolindex) else None 

        ## Close the data chunks 
        datachunks.close()
    ## Return the boolean that the outfile exists 
    return exists(outfile)

## Write ftn for splitting on .txt 
def splitontxt(infile:str) -> str:
     """Splits and input string (representing path to a .txt file) on the .txt extension."""
     ## split on the .txt and return the first item in list 
     return infile.split('.txt')[0] 

## Write function to reset output path
def resetoutput(infile:str,outfile:str,nsmaple:float) -> str:
     """Resets the output file name given the input if no output file name is defined."""
     ## Set nsmaple to an integfer if greater than one 
     nsample = nsmaple if (nsmaple < 1) else int(nsmaple)
     ## Return the new output name
     return outfile if outfile else splitontxt(infile) + '.random_sample_%s.txt'%nsample

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      MAIN SCRIPT & ARGUMENT PARSING 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## If the script is called 
if __name__ == "__main__":
    ## Bring in argparse, os and set parser
    import argparse

    ## Make the parse
    parser = argparse.ArgumentParser(description = downsample_desc)

    ## Add the required arguments
    parser.add_argument("-i", "--input-path",  dest="i", type=str,   required=True,  help = i_help, metavar = "./path/to/input.txt"    )
    parser.add_argument("-n", "--downsample",  dest="n", type=float, required=True,  help = n_help, metavar = 'n'                      )

    ## Add optional arguments
    parser.add_argument("-O", "--output-path", dest="O", type=str,   required=False, help = O_help, metavar = './path/to/saveout.txt', default= None)
    parser.add_argument("-R", "--random-seed", dest="R", type=int,   required=False, help = R_help, metavar = 'n',  default = ranseed  ) 
    parser.add_argument("-Z", "--chunk-size",  dest="Z", type=int,   required=False, help = Z_help, metavar = 'n',  default = chunks   )

    ## Add boolean var
    parser.add_argument("-F", "--force",       dest="F", help = F_help, action = ST)
    parser.add_argument("-I", "--intra-only",  dest="I", help = I_help, action = ST)

    ## Set the parsed values as inputs
    inputs = parser.parse_args()

    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ##      PASS ARGUMENTS & SET VARIABLES 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Format input arguments
    input_path  = inputs.i       ## Input path
    n_sample    = inputs.n       ## Amount of downsampling
    output_path = inputs.O       ## Output path 
    r_seed      = inputs.R       ## Random seed
    c_size      = inputs.Z       ## Set the chunk size 
    force       = inputs.F       ## Force overwrite 
    intra       = inputs.I       ## Force intra-chrom sampling only 

    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ##      INITIALIZATIONS
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Set randome seed 
    np.random.seed(r_seed)
    
    ## Check if the input path exists 
    assert exists(input_path),      "ERROR: A non-existent input path -- %s -- was passed!"%input_path 
    ## Reset output path to default if none was provided
    output_path = resetoutput(input_path,output_path,n_sample)
    ## Remove output path if it exsits and forcing was passed 
    remove(output_path) if (exists(output_path) and force) else None 
    ## Check that we are not about to overwright 
    assert not exists(output_path), "ERROR: output file -- %s -- already exists!\nINFO: Remove file and try running this script again."%output_path

    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ##      PERFORM DOWNSAMPLING
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## If we are in short format, write out sub smaple of input file asserting it exists 
    assert writesubsample(input_path, output_path, n_sample, intra, c_size)
## End of file     