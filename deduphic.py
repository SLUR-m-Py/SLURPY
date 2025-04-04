#!/usr/bin/env python
########################################################################################
##      Concat, sort, and deduplicate bedpe file representing Hi-C contacts           ##
########################################################################################
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## Set the descriptions
desc = "Concats, sorts, and removes duplicates from input bedpe files representing Hi-C contacts from a paired-end Hi-C experiment."

## ----------------------------------- MODULE LOADING ------------------------------------ ##
## Bring in pandas
import dask.dataframe as dd, numpy as np 
## Bring in params
from parameters import ST, hicsep
## Bring in ftn from defaults
from directories import bedtmpdir

## ---------------------------------- VARIABLE SETTING ------------------------------------ ##
## Set drop and sorting by columns 
drop_by = ['Chrn1','Chrn2','Pos1','Pos2','Seqrev1','Seqrev2','Mcount1','Mcount2']
sort_by = drop_by[:-4]
"""
Why?
Chrn1 and Chrn2:       Duplicates must map to the same pair of chromosomes / contigs.
Pos1 and Pos2:         Duplicates must have the same left most starting coordinates.
Mcount1 and Mcount2:   The number of "M" in the cigar string, duplicates have the same number.
Seqrev1 and Seqrev2:   Integers signifying the strand, 0 for positive, anything else for negative. 
"""

## Set help messages
sort_help = 'Flag to left right sort input rows by chromosome and position.'
dup_help  = 'Flag to drop duplicate rows. Duplicates are found by read positions and sequencing length of fragment.'
D_help    = 'Output file name and path to save duplciate read pairs to (in bedpe format).'
O_help    = 'Output file name and path to save results to (in bedpe format).'
B_help    = 'The ending pattern of bedpe files wanted as input into script.'

## Set run local
runlocal = False

## Ftn for returning the match coutn in a cigar str
def matchcount(cig:str) -> int:
    return cig.count('M')

## -------------------------------------- MAIN EXECUTABLE -------------------------------------------------- ##
## if the script is envoked
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse
    ## Make the parse
    parser = argparse.ArgumentParser(description = desc)
    ## Add the required arguments
    parser.add_argument("-b",      dest="B",     type=str,  required=True,    help=B_help ) 
    parser.add_argument("-o",      dest="O",     type=str,  required=True,    help=O_help )
    parser.add_argument("-d",      dest="D",     type=str,  required=False,   help=D_help )

    ## Add boolean vars 
    parser.add_argument("--sort",  dest="sort",   help = sort_help,  action = ST)
    parser.add_argument("--dedup", dest="dup",    help = dup_help,   action = ST)
   
    ## Set the paresed values as inputs
    inputs = parser.parse_args()

    ## Set inputs 
    filebackend = inputs.B      ## Set the file backend 
    output_path = inputs.O      ## Output path   
    dedupe_path = inputs.D      ## Deduplication out file
    sorting     = inputs.sort   ## Are we sorting
    deduplicate = inputs.dup    ## Are we deduplicating

    ## Check input
    assert filebackend.split('.')[-1] == 'bedpe', "ERROR: The given extension ending -- %s -- is not a bedpe file!"%filebackend

    ## Set wild card ofr input paths to dask dataframes 
    input_paths = f'{bedtmpdir}/*.{filebackend}' if not runlocal else f'*.{filebackend}'

    ## Load in bedpe file
    bedpe = dd.read_csv(input_paths,sep=hicsep)

    ## Preset duplicate counts
    interdup_counts = 0
    intradup_counts = 0

    ## Set up if statements, if we are BOTH deduplicateing and soritng our inputs 
    if sorting and deduplicate:
        ## Gather the counts and number of uniq pos1 and 2
        (pos1,pc1), (pos2,pc2) = np.unique(bedpe.Pos1,return_counts=True), np.unique(bedpe.Pos2,return_counts=True)
        ## Gather duplicated pos
        dpos1,dpos2 = pos1[(pc1>1)], pos2[(pc2>1)]
        
        ## Compute the possible duplicates
        pos_dups = bedpe[(bedpe.Pos1.isin(dpos1) & bedpe.Pos2.isin(dpos2))].compute()
        ## Calculate the mcounds, apply match count 
        pos_dups['Mcount1'] = pos_dups.Cigar1.apply(matchcount)
        pos_dups['Mcount2'] = pos_dups.Cigar2.apply(matchcount)

        ## Drop duplicates
        pos_uniq   = pos_dups.drop_duplicates(drop_by)
        ## Gather duplicate hits by read names 
        duplicates = pos_dups[~(pos_dups.Qname1.isin(pos_uniq.Qname1))] 
        ## Count the duplicates
        duplicate_count = duplicates.Pos1.count()
        ## If we have duplicates
        if duplicate_count:
            ## Calculate the inter duplicates
            interdup_count = duplicates.Inter.sum()
            ## Calc the intra dup count
            intradup_count = duplicate_count - interdup_count
            ## Update counts
            interdup_counts += interdup_count
            intradup_counts += intradup_count
            ## Save the duplicates to csv 
            duplicates.to_csv(dedupe_path,header=True,index=False,sep=hicsep) 

        ## Drop duplicates,sort values, and save out a csv file 
        bedpe[~(bedpe.Qname1.isin(duplicates.Qname1))].sort_values(sort_by).to_csv(output_path,index=False,header=True,single_file=True,sep=hicsep)

    ## Just sorting and not deduplciating 
    elif sorting and (not deduplicate):
        bedpe.sort_values(sort_by).to_csv(output_path,index=False,header=True,single_file=True,sep=hicsep)
    ## Otherwise just concat inputs 
    else: ## Do nothing
        print("WARNING: A combination of arguments left us doing nothing for this call of deduphic.py.")
    
    ## Format new names and print counts
    new_names = [  'Interdups',     'Intradups']
    new_count = [interdup_counts, intradup_counts]
    ## Iterate thru and print the counts to log 
    [print('INFO: %s\t%s'%(a,b)) for a,b in zip(new_names,new_count)]
    ## Print to log
    print("Finished concatonation, sorting, and deduplicating bedpe files ending with %s"%filebackend)
## End of file 