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
import dask.dataframe as dd
## Bring in params
from .defaults.parameters import hicsep
## Bring in ftn from defaults
from .defaults.defaults import bamtmpdir

## ---------------------------------- VARIABLE SETTING ------------------------------------ ##
## Set drop and sorting by columns 
drop_by = ['Chrn1','Chrn2','Pos1','Pos2','Len1','Len2']

"""
Why?
Chrn1 and Chrn2: Duplicates must map to the same pair of chromosomes / contigs
Pos1 and Pos2:   Duplicates must have the same left most starting coordinates
Len1 and Len2:   The lengths of the mapped fragments must be identical to be duplicates
"""

## Set help messages
sort_help = 'Flag to left right sort input rows by chromosome and position.'
dup_help  = 'Flag to drop duplicate rows. Duplicates are found by read positions and sequencing length of fragment.'
D_help    = 'Output file name and path to save duplciate read pairs to (in bedpe format).'
O_help    = 'Output file name and path to save results to (in bedpe format).'
B_help    = 'The ending pattern of bedpe files wanted as input into script.'

## Set run local
runlocal = False
## -------------------------------------- MAIN EXECUTABLE -------------------------------------------------- ##
## if the script is envoked
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse
    ## Make the parse
    parser = argparse.ArgumentParser(description = desc)
    ## Add the required arguments
    parser.add_argument("-b",      dest="B",  type=str,  required=True,    help=B_help ) 
    parser.add_argument("-o",      dest="O",  type=str,  required=True,    help=O_help )
    parser.add_argument("-d",      dest="D",  type=str,  required=False,   help=D_help )
    parser.add_argument("--sort",  dest="sort",   help = sort_help,  action = 'store_true')
    parser.add_argument("--dedup", dest="dup",    help = dup_help,   action = 'store_true')
   
    ## Set the paresed values as inputs
    inputs = parser.parse_args()

    filebackend = inputs.B
    output_path = inputs.O
    dedupe_path = inputs.D 
    sorting     = inputs.sort
    deduplicate = inputs.dup

    ## Check input
    assert filebackend.split('.')[-1] == 'bedpe', "ERROR: The given extension ending -- %s -- is not a bedpe file!"%filebackend

    ## Set wild card ofr input paths to dask dataframes 
    input_paths = f'{bamtmpdir}/*.{filebackend}' if not runlocal else f'*.{filebackend}'

    ## Load in bedpe file, calc starting number of rows 
    bedpe = dd.read_csv(input_paths,sep=hicsep)
    start_count = bedpe.Qname1.count().compute()

    ## Print to screen
    print("INFO: %s rows at start."%start_count)
    
    ## Set up if statements, if we are BOTH deduplicateing and soritng our inputs 
    if sorting and deduplicate:
        ## Drop duplicates 
        df = bedpe.drop_duplicates(drop_by).sort_values(drop_by[:-2])
        ## Calculate end count of rows and print 
        end_count = df.Qname1.count().compute()
        ## Calc n removed 
        nremoved = int(start_count - end_count)

        ## Print to file 
        print("INFO: %s rows at end."%end_count)
        print("INFO: %s rows removed."%(nremoved))

        ## Gather the duplicate reads and save to csv 
        bedpe[~(bedpe.Qname1.isin(df.Qname1))].to_csv(dedupe_path,index=False,sep=hicsep,single_file=True) if nremoved else None 

    ## Just sorting and not deduplciateing 
    elif sorting and (not deduplicate):
        df = bedpe.sort_values(drop_by[:-3])
    ## Otherwise just concat inputs 
    else:
        df = bedpe

    ## Save out valid non dedupe counts 
    df.to_csv(output_path,index=False,single_file=True,sep=hicsep)

    ## Print to log
    print("Finished concatonation, sorting, and deduplicating bedpe files ending with %s"%filebackend)
## End of file 