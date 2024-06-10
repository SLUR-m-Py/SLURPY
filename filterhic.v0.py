#!/usr/bin/env python
#################################################
##           Hi-C Filter (v 0.0.0)             ##
#################################################
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## Set version
version = '0.0.0'

## Set the descriptions
filter_desc = "Hi-C FILTER (v %s): Filters an input juicer file (space deliminated) from an Hi-C experiment on ligation errors, self-circles, and dangling ends."%version

## ----------------------------------- MODULE LOADING ------------------------------------ ##
## load in ftns from pysam tools
from pysamtools import writeset, byseq, loadref

## Bring in ftns from slurpy 
from defaults import dictzip, circle_dist, error_dist, hicsep, C_help, E_help, L_help, r_help, juicer_cols, juicer_types, reset

## Ftn for saving dataframe 
def savetohic(indf,outpath,m='w',h=True):
    """Saves out the input dataframe (df) to output path."""
    ## Save out the passed df to the outpath as a sinfle file with header 
    indf.to_csv(outpath,sep=hicsep,mode=m,single_file=True,compute=True,header=h,index=False)
    pass 

## Ftn for defninging restriciton sites 
def returnsite(enzyme):
    """Returns the restriction site sequences and dangaling sequences based on input Hi-C library enzyme."""
    ## If we are working with the arima kit (which we are)
    if (enzyme.lower() == 'arima'):
        restsites, dangsites = ['GATC', 'GAATC', 'GAGTC', 'GATTC', 'GACTC'], ['GATC', 'AATC', 'AGTC','ATTC', 'ACTC']
    ## The MobI , DpnII, or Sau3AI enzymes
    elif enzyme.lower() in ['mobi', 'dpnii','sau3ai']:
        restsites, dangsites = ['GATC'], ['GATC']
    ## The HindIII enzyme
    elif (enzyme.lower() == 'hindiii'):
        restsites,dangsites = ['AAGCTT'],['AGCT']
    else: ## otherwise return none
        restsites, dangsites = None, None 
    ## Return the sites
    return restsites, dangsites

## Ftn for checking dangling end 
def dangcheck(seq, isforward, endsequence):
    """Checks a pysam read for a dangling sequence at the start or end of read (given its strand orientation)."""
    ## Check if the forward read stars with the restriction sequence
    if isforward and seq.upper().startswith(endsequence):
        dc = 1
    ## Check if the reverse, end  read ends with the reverse of the restriction sequence
    elif not isforward and seq.upper().endswith(endsequence[::-1]):
        dc = 2
    else: ## Otherwise return zero count 
        dc = 0
    ## Return the check 
    return dc 

## Write ftn for checking dangling ends 
def checkends(s1,s2,dends):
    """Sums scores of checks for dangling ends (dends) for first (s1) and second (s2) strands."""
    ## Return the summed scores of chekcs for each strand 
    return sum([dangcheck(s1,1,d) for d in dends]) + sum([dangcheck(s2,0,d) for d in dends])

##  Ftn to iterate over the restriction sites and check if they are within a given sequence
def restcheck(chrom,left,right,reference,restsites):
    """Checks the number or matching restrcition sites in a given fragment of the genome."""
    ## Set the maximum rest count
    max_rs = 0
    for r in restsites:
        ## Gather the sequence of the fragment formed by the reads
        chromseq = str(byseq(reference,chrom)[left+len(r):right-len(r)+1].seq)
        ## Count the nubmer of rest sties 
        count_rs = len(chromseq.split(r)) - 1
        ## Split the chromseq by the rest site and cound
        max_rs = count_rs if count_rs > max_rs else max_rs 
    ## Return the max count of the restriciton s
    return max_rs 

## Set input help message
T_help = "Path to input juicer long format (+ T-length and orientation) file."

## Set the filt columns
filt_cols = juicer_cols + ['Chrn1','Chrn2','Tlen','Orientation']
filt_type = juicer_types + [  int,   int,    int,   str]

## Zip the column names and types 
filt_dict = dictzip(filt_cols,filt_type)

## ----------------------------------------------- MAIN EXECUTABLE -------------------------------------------------- ##
## if the script is envoked
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse

    ## Make the parse
    parser = argparse.ArgumentParser(description = filter_desc)

    ## Add the required arguments
    parser.add_argument("-t", dest="T", type=str,  required=True,  help=T_help, metavar='./path/to/input.txt'       )
    parser.add_argument("-r", dest="R", type=str,  required=True,  help=r_help, metavar='./path/to/ref.fai'         )
    parser.add_argument("-l", dest="L", type=str,  required=True,  help=L_help, metavar='Armia'                     )
    
    ## Add optional variables 
    parser.add_argument("-C", dest="C", type=int,  required=False, help=C_help, metavar='n',    default=circle_dist )
    parser.add_argument("-E", dest="E", type=int,  required=False, help=E_help, metavar='n',    default=error_dist  )

    ## Set the paresed values as inputs
    inputs = parser.parse_args()

    ## Set required and optional input variables 
    inpath      = inputs.T  ## Input path 
    refpath     = inputs.R  ## reference path 
    library     = inputs.L  ## Library used to construct hi-c 
    circle_dist = inputs.C  ## Linear distance to look for circles
    error_dist  = inputs.E  ## Error distance

    ## Bring in needed mods 
    import dask.dataframe as dd

    ## Ftn for loading data
    def hicloader(toload):
        """Loads in a hic file (space delminated) in juicer long format with dask."""
        ## Loads in the dataframe with dask 
        return dd.read_csv(toload,sep=hicsep,usecols=filt_cols,dtype=filt_dict)

    ## Load in the refernece 
    reference = loadref(refpath)

    ## Load in restriction sites and dangling ends 
    rest_site, dang_ends = returnsite(library)

    ## Gather the file head 
    filehead = inpath.split('.txt')[0]

    ## Set the file ends and the file names from input
    new_ends = ['dangling','selfcircle','errors','tohic']
    newfiles = [filehead+f'.{n}.txt' for n in new_ends]

    ## Format the dictionary of file names 
    filedict = dictzip(new_ends,newfiles)

    ## Remove previous attempts 
    k = reset(newfiles)

    ## Load in the hic data with dask 
    df = hicloader(inpath)

    ## Initilize lists of errors, self circles, and reads with dangling ends 
    errors, self_circles, read_w_dangends = [], [], []

    ## IF we are checking for errors at a given distance and have rest sites 
    if error_dist and rest_site:
        ## Gather the outward, intra-chormosom contacts 
        pointing_intra = df[(df.Chr1==df.Chr2) & 
                            (df.Tlen<error_dist) & 
                            (df.Orientation.isin(['Left','Right']))][['Qname1','Chr1','Pos1','Pos2']].compute()
        
        ## Iterate over the restriction sites and count the occurance of each rest site 
        er = pointing_intra.apply(lambda x: restcheck(x['Chr1'], x['Pos1'], x['Pos2'], reference, rest_site), axis=1)
        ## Append the qnames to list
        errors = set(pointing_intra[(er==0)].Qname1)
        
        ## Write the set out to file
        k = writeset(filedict['errors'],errors,mode='w') if len(errors) else None

    ## IF we are checking for self circles and have rest sites 
    if circle_dist and rest_site:
        ## Gather the outward contacts by the first chromosome
        outward_intra = df[(df.Chr1==df.Chr2) & 
                            (df.Tlen<circle_dist) & 
                            (df.Orientation=='Outward')][['Qname1','Chr1','Pos1','Pos2']].compute()
        
        ## Iterate over the restriction sites and count the occurance of each rest site 
        sc = outward_intra.apply(lambda x: restcheck(x['Chr1'], x['Pos1'], x['Pos2'], reference, rest_site), axis=1)
        ## Append the qnames to list
        self_circles = set(outward_intra[(sc==0)].Qname1)
        
        ## Write the set out to file
        k = writeset(filedict['selfcircle'],self_circles,mode='w') if len(self_circles) else None

    ## If the dang ends is defined 
    if dang_ends:
        ## Gather the inward contacts by the first chromosome 
        inward_orientd = df[(df.Orientation=='Inward')][['Qname1','Seq1','Seq2']].compute()

        ## Do we have dang ends?
        have_dang_ends = inward_orientd.apply(lambda x: checkends(x['Seq1'], x['Seq2'], dang_ends),axis=1)
        ## Make the reads into a set 
        read_w_dangends = set(inward_orientd[(have_dang_ends>0)].Qname1)

        ## Write the set out to file
        k = writeset(filedict['dangling'],read_w_dangends,mode='w') if len(read_w_dangends) else None

    ## Remove the recoreds with Qnames in any of our above sets and save out the filtered records
    savetohic(df[~((df.Qname1.isin(read_w_dangends)) | (df.Qname1.isin(self_circles)) | (df.Qname1.isin(errors)))],filedict['tohic'])
## End of file 