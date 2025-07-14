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
## Bring in params
from parameters import ST, hicsep, save_help
## Bring in ftn from defaults
from directories import bedtmpdir
## Load in from defaults
from defaults import sortglob
## Set debuging 
debuging = False 

## ---------------------------------- VARIABLE SETTING ------------------------------------ ##
## Set drop and sorting by columns 
drop_by = ['Chrn1','Chrn2','Pos1','Pos2','Seqrev1','Seqrev2','Mcount1','Mcount2']
sort_by = drop_by[:4]
"""
Why?
Chrn1 and Chrn2:       Duplicates must map to the same pair of chromosomes / contigs.
Pos1 and Pos2:         Duplicates must have the same left most starting coordinates.
Mcount1 and Mcount2:   The number of "M" in the cigar string, duplicates have the same number.
Seqrev1 and Seqrev2:   Integers signifying the strand, 0 for positive, anything else for negative. 
"""

## Set help messages
sort_help = 'Flag to left right sort input rows by chromosome and position.'
dup_help  = 'Flag to drop duplicate rows. Duplicates are found by read positions and sequencing signal of fragment.'
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
    parser.add_argument("-b",           dest="B",     type=str,  required=True,    help=B_help ) 
    parser.add_argument("-o",           dest="O",     type=str,  required=True,    help=O_help )
    parser.add_argument("-d",           dest="D",     type=str,  required=False,   help=D_help )

    ## Add boolean vars 
    parser.add_argument("--sort",       dest="sort",   help = sort_help,  action = ST)
    parser.add_argument("--dedup",      dest="dup",    help = dup_help,   action = ST)
    parser.add_argument("--save-dups",  dest="save",   help = save_help,  action = ST)
   
    ## Set the paresed values as inputs
    inputs = parser.parse_args()

    ## Bring in pandas
    import dask.dataframe as dd

    ## Set inputs 
    filebackend = inputs.B      ## Set the file backend 
    output_path = inputs.O      ## Output path   
    dedupe_path = inputs.D      ## Deduplication out file
    sorting     = inputs.sort   ## Are we sorting
    deduplicate = inputs.dup    ## Are we deduplicating
    keep_dups   = inputs.save   ## Flag to save out the duplicates

    ## Set wild card ofr input paths to dask dataframes 
    input_wc =  f'*{filebackend}' if runlocal else f'{bedtmpdir}/*{filebackend}'
    input_paths = sortglob(input_wc)
    ## Check our work 
    assert len(input_paths), "ERROR: No input files found."
    print(input_paths) if debuging else None 

    ## Load in bedpe file
    bedpe = dd.concat([dd.read_parquet(input_path) for input_path in input_paths])

    ## Preset duplicate counts
    interdup_counts = 0
    intradup_counts = 0

    ## Set up if statements, if we are BOTH deduplicateing and soritng our inputs 
    if sorting and deduplicate:
        ## Drop the non-unique rows via dask dataframes, and sort and save out csv
        deduped = bedpe.drop_duplicates(drop_by)
        deduped.sort_values(sort_by).to_csv(output_path,index=False,header=True,single_file=True,sep=hicsep)

        if keep_dups:
            ## Gather the duplicates and count them 
            duplicates = bedpe[(~bedpe.Qname1.isin(deduped.Qname1.compute().tolist()))]
            duplicate_count = duplicates.Pos1.count().compute()

            ## If we have duplicates
            if duplicate_count:
                ## Calculate the inter duplicates
                interdup_counts = duplicates.Inter.sum().compute()
                ## Calc the intra dup count
                intradup_counts = duplicate_count - interdup_counts
                ## Save the duplicates to csv 
                duplicates.to_csv(dedupe_path,index=False,header=True,single_file=True,sep=hicsep)
        else:
            ## Calc totals 
            Total_counts = bedpe.Chrn1.count().compute()
            Inter_counts = bedpe.Inter.sum().compute()
            Intra_counts = Total_counts - Inter_counts

            ## Calc new totals 
            Total_dedup = deduped.Chrn1.count().compute()
            Inter_dedup = deduped.Inter.sum().compute()
            Intra_dedup = Total_dedup - Inter_dedup

            ## Calc the removed
            interdup_counts = Inter_counts - Inter_dedup 
            intradup_counts = Intra_counts - Intra_dedup

    ## Just sorting and not deduplciating 
    elif sorting and (not deduplicate):
        bedpe.sort_values(sort_by).to_csv(output_path,index=False,header=True,single_file=True,sep=hicsep)
    ## Otherwise just concat inputs 
    else: ## Do nothing
        print("WARNING: A combination of arguments left us doing nothing for this call of deduphic.py.")
    
    ## Format new names and print counts
    new_names = ['InterDuplicates', 'IntraDuplicates']
    new_count = [interdup_counts, intradup_counts]
    ## Iterate thru and print the counts to log 
    [print('INFO: %s\t%s'%(a,b)) for a,b in zip(new_names,new_count)]
    ## Print to log
    print("Finished concatonation, sorting, and deduplicating bedpe files ending with %s"%filebackend)
## End of file 