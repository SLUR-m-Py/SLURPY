#!/usr/bin/env python
#SBATCH --job-name=downsample           ## Name of job
#SBATCH --output=%x.%j.out              ## Name stdout
#SBATCH --error=%x.%j.err               ## Name stderr
#SBATCH --partition=mpi                 ## Size of partitions
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=12              ## Number of tasks to be launched
#####################################
##    Hi-C Random Down-sampling
#####################################
"""
Written by:

    Sasha Bacot (sbacot@lanl.gov) and Cullen Roth (croth@lanl.gov)

Genomics and Bioanalytics (B-GEN), Los Alamos National Laboratory, Los Alamos, NM 87545
"""
## ---------------------------------------------------------------------------- ##
##      IMPORTING MODS, FUNCTIONS, and VARIABLES
## ---------------------------------------------------------------------------- ##
## Load in pandas and numpy and others 
import pandas as pd, numpy as np, gzip 
## Bring in exists from os 
from os.path import exists
## from os load in remove
from os import remove

## ---------------------------------------------------------------------------- ##
##      SETING DEFAULTS AND HELP MESSAGES 
## ---------------------------------------------------------------------------- ##
## Set the description
downsample_desc = "Randomly downsamples rows from an input .txt file representing Hi-C contacts."

## Set the defaults
ranseed = 720
chunks  = 500000
ST      = 'store_true'

## Set some default input variables used in these functions, such as the input column names and file seperator
hiccols     = ['Strand1', 'Chr1', 'Pos1', 'Dum1', 'Strand2', 'Chr2', 'Pos2', 'Dum2']
hicsep      = ' '   ## space deliminator 

## ---------------------------------------------------------------------------- ##
##      FUNCTION DEFINING
## ---------------------------------------------------------------------------- ##
## Ftn for gzipped file 
def isgzip(inpath) -> bool: 
    ## Return the check 
    return inpath.split('.')[-1] == 'gz'

## Ftn for counting 
def countshort(inpath:str) -> tuple[int]:
    ## Iniate counters 
    inter_count = 0
    ## open the ifle
    with (gzip.open(inpath,'rt') if isgzip(inpath) else open(inpath,'rt')) as infile:
        ## Iterate
        for i,l in enumerate(infile):
            line = l.split(hicsep)
            ## If is an inter chrom contact add to inter
            if line[1] != line[5]:
                inter_count += 1
    ## Return counts
    return (i + 1) - inter_count, inter_count

## Ftn for counting 
def countbedpe(inpath:str) -> tuple[int]:
    ## Iniate counters 
    inter_count = 0
    ## open the ifle
    with (gzip.open(inpath,'rt') if isgzip(inpath) else open(inpath,'rt')) as infile:
        ## Iterate thru lines of file 
        for i,l in enumerate(infile):
            if not i: ## Assumes first line is header, gather index of "Inter"
                inter_ix = l.split(hicsep).index('Inter')
            ## When we have a first row representing the hic contact
            if i: ## If is an inter chrom contact add to inter
                inter_count += int(l.split(hicsep)[inter_ix])
    ## Return counts
    return i - inter_count, inter_count

## Get a random index
def randomix(total:int,samplesize:int):
    return np.random.choice(range(total),samplesize,replace=False)

## Ftn for getting random dhcun
def getchunk(df:pd.DataFrame,randomix:np.array,runningcount:int) -> pd.DataFrame:
    ## Reset the intra index
    df.reset_index(drop=True,inplace=True)
    ## Add the intra running count
    df.index = df.index + runningcount
    ## Get the random set of intra by new index
    return df[(df.index.isin(randomix))]
    
## Define the main function to write out sub sampled input bedpe file 
def subbedpe(infile:str, outfile:str, inter_sample:int, intra_sample:int, inter_total:int, intra_total:int, chunksize=chunks) -> bool: 
    """Loads in Hi-C contacts from an input bedpe file from SLUR(M)-py (via chunking) and randomly downsamples a set of rows, saving to an output file."""
    ## Gather the intra and inter ix to sample
    intra_ix = randomix(intra_total,intra_sample)
    inter_ix = randomix(inter_total,inter_sample)
    ## Iniate running count of intra and inter contacts
    runing_intra = 0 
    runing_inter = 0
    
    ## Open the chunks with a with statment
    with pd.read_csv(infile, sep = hicsep, chunksize=chunksize) as datachunks:
        ## Iterate over the data chunks 
        for i,chunk in enumerate(datachunks):

            ## Split the intra contacts from the chunk and gather the random chunk 
            intra_chunk = chunk[(chunk.Inter==0)]
            intra_randm = getchunk(intra_chunk,intra_ix,runing_intra)
            
            ## Split the inter contacts from the chunk and gather the randomd chunk 
            inter_chunk = chunk[(chunk.Inter==1)]
            inter_randm = getchunk(inter_chunk,inter_ix,runing_inter)

            ## Set the new chunk and write out the sub smaple of the chunks  
            pd.concat([intra_randm,inter_randm]).to_csv(outfile, sep=hicsep, index=False, header=not i, mode='a' if i else 'w') 

            ## Add to the running countof both inter and intra contacts
            runing_intra += intra_chunk.shape[0]
            runing_inter += inter_chunk.shape[0]

        ## Close the data chunks 
        datachunks.close()
    ## Return the boolean that the outfile exists 
    return exists(outfile)

## Define the main function to write out sub sampled input short file 
def subshort(infile:str, outfile:str, inter_sample:int, intra_sample:int, inter_total:int, intra_total:int, chunksize=chunks) -> bool: 
    """Loads in Hi-C contacts from an input short file from Juicer (via chunking) and randomly downsamples a set of rows, saving to an output file."""
    ## Gather the intra and inter ix to sample
    intra_ix = randomix(intra_total,intra_sample)
    inter_ix = randomix(inter_total,inter_sample)
    ## Iniate running count of intra and inter contacts
    runing_intra = 0 
    runing_inter = 0

    ## Open the chunks with a with statment
    with pd.read_csv(infile, sep = hicsep, names=hiccols, header=None, chunksize=chunksize) as datachunks:
        ## Iterate over the data chunks 
        for i,chunk in enumerate(datachunks):

            ## Split the intra contacts from the chunk and gather the random chunk 
            intra_chunk = chunk[(chunk.Chr1==chunk.Chr2)]
            intra_randm = getchunk(intra_chunk,intra_ix,runing_intra)
            
            ## Split the inter contacts from the chunk and gather the randomd chunk 
            inter_chunk = chunk[(chunk.Chr1!=chunk.Chr2)]
            inter_randm = getchunk(inter_chunk,inter_ix,runing_inter)

            ## Set the new chunk and write out the sub smaple of the chunks  
            pd.concat([intra_randm,inter_randm]).to_csv(outfile, sep=hicsep, index=False, header=False, mode='a' if i else 'w') 

            ## Add to the running countof both inter and intra contacts
            runing_intra += intra_chunk.shape[0]
            runing_inter += inter_chunk.shape[0]

        ## Close the data chunks 
        datachunks.close()
    ## Return the boolean that the outfile exists 
    return exists(outfile)

## Write ftn for splitting on .txt 
def splitontxt(infile:str) -> str:
     """Splits and input string (representing path to a .txt file) on the .txt extension."""
     ## split on the .txt and return the first item in list 
     return infile.split('.bedpe')[0] 

## Write function to reset output path
def resetoutput(infile:str,outfile:str,nsmaple:int) -> str:
     """Resets the output file name given the input if no output file name is defined."""
     ## Set nsmaple to an integfer if greater than one 
     nsample = nsmaple if (nsmaple < 1) else int(nsmaple)
     ## Return the new output name
     return outfile if outfile else splitontxt(infile) + '.random_sample_%s.bedpe'%nsample

## Set help messages 
b_help = "Path to an input Hi-C text file from SLUR(M)-py in bedpe format."
i_help = "Number of intra-chromosomal contacts (within chromosomes) to randomly sample without replacement."
I_help = "The total number of known intra-chromosomal contacts. Not needed but passing this variable will speed up script execution"
x_help = "Number of inter-chromosomal contacts (within chromosomes) to randomly sample without replacement."
X_help = "The total number of known inter-chromosomal contacts. Not needed but passing this variable will speed up script execution"
p_help = "Percent (0 to 1) of intra-chromosomal contacts (within chromosomes) to randomly sample without replacement. Float values less than 1 are interpreted as a fraction of intra-chromosomal Hi-C contacts to downsample."
f_help = "Percent (0 to 1) of inter-chromosomal contacts (between chromosomes) to randomly sample without replacement. Float values less than 1 are interpreted as a fraction of inter-chromosomal Hi-C contacts to downsample."
R_help = "Set integer used as random seed in numpy for random sampling (default: %s)"%ranseed
Z_help = "Number of rows (default: %s) loaded into pandas at a time. WARNING: while increasing could speed up pipeline it could also cause memeory issues."%chunks
O_help = "The output path to save out result of random downsampling. The default behavior is to save in the same directory as the input file."
F_help = "Boolean flag to force overwrite of output file." 
S_help = "Boolean flag indicating input is in juicers short format."

## ---------------------------------------------------------------------------- ##
##      MAIN SCRIPT & ARGUMENT PARSING 
## ---------------------------------------------------------------------------- ##
## If the script is called 
if __name__ == "__main__":
    ## Bring in argparse, os and set parser
    import argparse

    ## Make the parse
    parser = argparse.ArgumentParser(description = downsample_desc)

    ## Add the required arguments
    parser.add_argument("-b", "--bedpe-path",    dest="b", type=str,    required=True,   help = b_help, metavar = "./path/to/input.txt")
    
    ## Add optional arguments
    parser.add_argument("-i", "--intra-sample",  dest="i", type=int,    required=False,  help = i_help, default = 0)
    parser.add_argument("-I", "--total-intra",   dest="I", type=int,    required=False,  help = I_help, default = 0)
    parser.add_argument("-x", "--inter-sample",  dest="x", type=int,    required=False,  help = x_help, default = 0)
    parser.add_argument("-X", "--total-inter",   dest="X", type=int,    required=False,  help = X_help, default = 0)
    parser.add_argument("-p", "--intra-portion", dest="p", type=float,  required=False,  help = p_help, default = 0)
    parser.add_argument("-f", "--inter-portion", dest="f", type=float,  required=False,  help = f_help, default = 0)
    parser.add_argument("-R", "--random-seed",   dest="R", type=int,    required=False,  help = R_help, default = ranseed) 
    parser.add_argument("-Z", "--chunk-size",    dest="Z", type=int,    required=False,  help = Z_help, default = chunks)
    parser.add_argument("-O", "--output-path",   dest="O", type=str,    required=False,  help = O_help, metavar = './path/to/saveout.txt', default= None)

    ## Add boolean vars
    parser.add_argument("-F", "--force", dest="F", help = F_help, action = ST)
    parser.add_argument("-S", "--short", dest="S", help = S_help, action = ST)

    ## Set the parsed values as inputs
    inputs = parser.parse_args()

    ## ---------------------------------------------------------------------------- ##
    ##      PASS ARGUMENTS & SET VARIABLES 
    ## ---------------------------------------------------------------------------- ##
    ## Format input arguments
    bedpe_inpath = inputs.b      ## Input path
    intra_sample = inputs.i      ## N intra contacts to sample 
    intra_count  = inputs.I      ## Total N of intra counts
    inter_sample = inputs.x      ## N inter contacts to sample 
    inter_count  = inputs.X      ## Total N of inter counts 
    intra_frac   = inputs.p      ## Fraction of intra contacts to sample
    inter_frac   = inputs.f      ## Fraction of inter contacts to sample     
    ran_seed     = inputs.R      ## Random seed
    chunks       = inputs.Z      ## Set chunk size 
    output_path  = inputs.O      ## Output path 
    force        = inputs.F      ## Force overwrite?
    isshort      = inputs.S      ## Short format?

    ## ---------------------------------------------------------------------------- ##
    ##      INITIALIZATIONS
    ## ---------------------------------------------------------------------------- ##
    ## Set randome seed 
    np.random.seed(ran_seed)
    ## Check if the input path exists 
    assert exists(bedpe_inpath), "ERROR: A non-existent input path -- %s -- was passed!"%bedpe_inpath

    ## Gather the total counts, if non i.e. zero sum was passed
    if not (intra_count+inter_count):
        if isshort:
            intra_count,inter_count = countshort(bedpe_inpath)
        else:
            intra_count,inter_count = countbedpe(bedpe_inpath)

    ## print the totals
    print(f'INFO: Total inter-contacts: {inter_count}')
    print(f'INFO: Total intra-contacts: {intra_count}')

    ## Set the number of inter- and intra- contacts to sample was passed
    inter_sample = inter_sample if inter_sample else int(inter_frac*inter_count)
    intra_sample = intra_sample if intra_sample else int(intra_frac*intra_count)

    ## print the totals
    print(f'INFO: Total inter-contacts to sample: {inter_sample} ( {100*round(inter_sample/inter_count,4)} %)')
    print(f'INFO: Total intra-contacts to sample: {intra_sample} ( {100*round(intra_sample/intra_count,4)} %)')

    ## Check our work
    assert inter_sample + intra_sample, "ERROR: Unable to calculate the number of contacts to subsample! Please check inputs."

    ## Set the n sample
    n_sample = inter_sample + intra_sample
    ## Reset output path to default if none was provided
    output_path = resetoutput(bedpe_inpath,output_path,n_sample)
    ## Remove output path if it exsits and forcing was passed 
    remove(output_path) if (exists(output_path) and force) else None 
    ## Check that we are not about to overwright 
    assert not exists(output_path), "ERROR: output file -- %s -- already exists!\nINFO: Remove file and try running this script again."%output_path

    ## ---------------------------------------------------------------------------- ##
    ##      PERFORM DOWNSAMPLING
    ## ---------------------------------------------------------------------------- ##
    ## If we are in short format, write out sub smaple of input file asserting it exists 
    if isshort:
        assert subshort(bedpe_inpath, output_path, inter_sample, intra_sample, inter_count, intra_count, chunksize=chunks)
    else:
        assert subbedpe(bedpe_inpath, output_path, inter_sample, intra_sample, inter_count, intra_count, chunksize=chunks)
    ## Print to file
    print("Finished subsampleing input file: %s"%bedpe_inpath)
## End of file    