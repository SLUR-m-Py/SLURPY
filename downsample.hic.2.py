#!/usr/bin/env python
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

## Load in z help and chunk sizes from defaults 
from defaults import Z_help, chunks, isgzip

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      SETING DEFAULTS AND HELP MESSAGES 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## Set the description
downsample_desc = "Randomly downsamples rows from an input .txt file representing Hi-C contacts."

## Set the defaults
ranseed = 720

## Set help messages 
I_help = "Path to an input Hi-C text file."
O_help = "The output path to save out result of random downsampling. The default behavior is to save in the same directory as the input file."
N_help = "Number (or percent) of rows to randomly downsample. Float values less than 1 are interpreted as a fraction of rows (Hi-C contacts) to downsample."
S_help = "Set integer used as random seed in numpy (default: %s)"%ranseed

## Set some default input variables used in these functions 
## Such as the input column names, file seperator, header, and replacement value 
hiccols     = ['Strand1', 'Chr1', 'Pos1', 'Dum1', 'Strand2', 'Chr2', 'Pos2', 'Dum2']
hicsep      = ' '
hicheader   = None 
replacement = False ## Do not change this, sampling with replacement is not supported in this script 

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      FUNCTION DEFINING
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## Ftn for counting lines of file 
def countlines(inputfile:str) -> int:
    """Counts the number of lines from an input file. Input file can be gzipped."""
    ## If the file is gizpped, open with gzip open, else just open 
    with (gzip.open(inputfile, mode = 'rb') if isgzip(inputfile) else open(inputfile, mode = 'rb')) as infile:
        for i, line in enumerate(infile):
                pass
        ## Close the file 
        infile.close()
    ## Return the line count + one, since it is zero indexed
    return i + 1 

## Ftn for reseting sample size
def getsamplesize(nsample:float,nlines:int):
    """Calculates the number of rows to subsample from given the total number of lines.""" 
    ## Reset the number of lines to randomly downsample if a fraction/floating point number was given as input
    return int(np.ceil(nlines * nsample) if (nsample < 1) else nsample)

## Ftn for precalculating random index
def getrandom(nlines:int,samplesize:int):
    """Calculates an array of random integers given the total number of lines (nlines) and the sample size."""
    ## REturn random numpy array from nlines, given the sample size 
    return np.random.choice(nlines, samplesize, replace = replacement)

## Define the main function to write out sub sampled input file 
def writesubsample(inputfile:str, outputfile:str, nsample:float, buffersize=10000, saveheader=False, saveindex=False) -> bool: 
    """Loads in Hi-C contacts from an input .txt file (via chunking) and randomly downsamples a set of rows, saving to an output file."""
    ## Count the total number of lines in file
    linecount = countlines(inputfile)

    ## Get the sample
    nsample = getsamplesize(nsample,linecount)

    ## Calc those to remove
    toremove = linecount - nsample
    print('Found %s lines, removing %s'%(linecount,toremove))

    ## Calculate the smaple size and create the random indices to remove from file 
    dropindices = getrandom(linecount,toremove)

    ## Initate add lines 
    addlines = []
    ## Write out a gzipped file 
    with (gzip.open(inputfile, mode = 'rb') if isgzip(inputfile) else open(inputfile, mode = 'rb')) as infile:
        with open(outputfile, 'wb') as outfile:
            ## Write all lines with an index not in the dropindex 
            for i, line in enumerate(infile):
                if (i in dropindices):
                     pass 
                else: 
                    outfile.write(line)
                     #addlines.append(line)

                #if len(addlines) > buffersize:
                #    outfile.writelines(addlines)

                    ## Reset add lines 
                #    addlines = []
                #else:
                #    pass 
            
            ## Write out the last chunk
            #outfile.writelines(addlines) if len(addlines) else None 

        ## Closet the files 
        infile.close()
        outfile.close()

    ## Return the boolean that the outfile exists 
    return exists(outputfile)

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
    return outfile if outfile else (splitontxt(infile) + '.random_sample_%s.txt'%nsample)
    #tmpfile = outfile if outfile else (splitontxt(infile) + '.random_sample_%s.txt'%nsample)
    #return tmpfile if isgzip(tmpfile) else tmpfile + '.gz'
    
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
    parser.add_argument("-I", "--input-path",  dest="I", type=str,   required=True,  help = I_help, metavar = "./path/to/input.txt"    )
    parser.add_argument("-N", "--downsample",  dest="N", type=float, required=True,  help = N_help, metavar = 'n'                      )

    ## Add optional arguments
    parser.add_argument("-O", "--output-path", dest="O", type=str,   required=False, help = O_help, metavar = './path/to/saveout.txt', default= None)
    parser.add_argument("-S", "--random-seed", dest="S", type=int,   required=False, help = S_help, metavar = 'n',  default = ranseed  ) 
    parser.add_argument("-Z", "--chunk-size",  dest="Z", type=int,   required=False, help = Z_help, metavar = 'n',  default = chunks   )

    ## Set the parsed values as inputs
    inputs = parser.parse_args()

    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ##      PASS ARGUMENTS & SET VARIABLES 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Format input arguments
    input_path  = inputs.I       ## Input path
    output_path = inputs.O       ## Output path 
    n_sample    = inputs.N       ## Amount of downsampling
    r_seed      = inputs.S       ## Random seed
    c_size      = inputs.Z       ## Set the chunk size 

    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ##      INITIALIZATIONS
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Set randome seed 
    np.random.seed(r_seed)
    
    ## Check if the input path exists 
    assert exists(input_path),      "ERROR: A non-existent input path -- %s -- was passed!"%input_path 

    ## Reset output path to default if none was provided
    output_path = resetoutput(input_path,output_path,n_sample)
 
    ## Check that we are not about to overwright 
    assert not exists(output_path), "ERROR: output file -- %s -- already exists!\nINFO: Remove file and try running this script again."%output_path

    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ##      PERFORM DOWNSAMPLING
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Write out sub smaple of input file asserting it exists 
    assert writesubsample(input_path, output_path, n_sample, c_size)
## End of file     