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
from parameters import ST, hicsep, save_help, Z_help, chunksize
## Bring in ftn from defaults
from directories import bedtmpdir
## Load in from defaults
from defaults import sortglob, remove
## Set debuging 
debuging = False 
## Load in subprocess 
import subprocess, pandas as pd 
## Bring in conactonation ftn 
from pandacat import concatonation

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

## Bring in dask
#import dask.dataframe as dd, pandas as pd  

import pandas as pd 

## Define ftn for laoding in parquets with dask
#def loadparquets(inpaths:list,usecols = []):
#    ## Check inputs
#    assert len(inpaths) and (type(inpaths) == list), "ERROR: A list with paths was not provided for deduplication!"
#    ## initilize list of dask dataframes
#    to_cat = []
#    ## iterate over input paths
#    for p in inpaths:
#        try: ## Attempt ot cat the file 
#            to_cat.append(dd.read_parquet(p,columns=usecols) if len(usecols) else dd.read_parquet(p))
#        ## Soem seem to be empty
#        except ValueError:
#            print("INFO: One of the given parquet files did not load and might be empty: %s"%p)
#    ## Return the dask concatonated list
#    return dd.concat(to_cat)

## Def for concatonatin files

## -------------------------------------- MAIN EXECUTABLE -------------------------------------------------- ##
## if the script is envoked
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse
    ## Make the parse
    parser = argparse.ArgumentParser(description = desc)
    ## Add the required arguments
    parser.add_argument("-b", dest="B", type=str, required=True,  help=B_help ) 
    parser.add_argument("-o", dest="O", type=str, required=True,  help=O_help )
    parser.add_argument("-z", dest="Z", type=int, required=False, help=Z_help,   default=chunksize )
    parser.add_argument("-d", dest="D", type=str, required=False, help=D_help,   default=None      )

    ## Add boolean vars 
    parser.add_argument("--save-dups",  dest="save",   help = save_help,  action = ST)
   
    ## Set the paresed values as inputs
    inputs = parser.parse_args()

    ## Set inputs 
    filebackend = inputs.B          ## Set the file backend 
    output_path = inputs.O          ## Output path   
    chunksize   = inputs.Z          ## Chunksize for loading in rows 
    dedupe_path = inputs.D          ## Deduplication out file
    keep_dups   = inputs.save       ## Flag to save out the duplicates

    ## Was a path passed 
    deduplicate = True if dedupe_path else False  ## Are we deduplicating

    ## Set wild card ofr input paths to dask dataframes 
    input_wc =  f'*{filebackend}' if runlocal else f'{bedtmpdir}/*{filebackend}'
    ## Bring in paths by wild card
    input_paths = sortglob(input_wc)
    ## Check our work 
    print(input_paths) if debuging else None 

    ## Preset duplicate counts
    interdup_counts = 0
    intradup_counts = 0
    Total_counts    = 0 

    ## Load the first of the file  #bedpe = loadparquets(input_paths)
    bedpe = pd.read_csv(input_paths[0],sep=hicsep)
    ## Gather the column space
    column_names = bedpe.columns.tolist()
    ## set the index of the positions 
    chr1_ix = column_names.index('Chrn1') + 1
    chr2_ix = column_names.index('Chrn2') + 1
    pos1_ix = column_names.index('Pos1')  + 1
    pos2_ix = column_names.index('Pos2')  + 1

    ## Set the sorted ouput path 
    sorted_path = (output_path + '.sort.tmp') if deduplicate else output_path
    tmp_path    = (output_path + '.tmp') 

    ## Write the column space to file 
    with open(sorted_path,'w') as outfile:
        outfile.writelines(' '.join(column_names)+'\n')
        outfile.close()
    ## Sort the input parquet files #bedpe.sort_values(sort_by).to_csv(sorted_path,index=False,header=True,single_file=True,sep=hicsep)
    #bedpe.to_csv(tmp_path,index=False,header=False,single_file=True,sep=hicsep)
    ## Concatonate the bedpe files with no header 
    assert concatonation(input_paths,tmp_path,noheader=True), "ERROR: Unable to concatonate input files!"
    ## Call our sort command
    k = subprocess.run(f'sort -nb -k{chr1_ix},1 -k{chr2_ix},2 -k{pos1_ix},3 -k{pos2_ix},4 {tmp_path} >> {sorted_path}',shell=True)

    ## Set up if statements, if we are BOTH deduplicateing and soritng our inputs 
    if deduplicate:
        ## Appending?
        uniq_appending = False
        dups_appending = False
        last_chunk     = []
        ## Load in the sorted file for deduplication
        with pd.read_csv(sorted_path,sep=hicsep,chunksize=chunksize) as chunks:
            for chunk in chunks:
                ## Add to the counts
                Total_counts += chunk.Chrn1.count()
                ## Append the cunk to the last chunk 
                chunk = pd.concat([last_chunk,chunk],axis=0).reset_index(drop=True) if len(last_chunk) else chunk
                ## Get the last row 
                lastrow  = chunk.tail(1)
                
                ## Gather all rows within the last chunk 
                last_chunk = chunk[(chunk.Chrn2==lastrow.Chrn2.max()) & (chunk.Pos1==lastrow.Pos1.max()) & (chunk.Pos2==lastrow.Pos2.max())]
                ## Gather the rest of the chunk to dedup licat
                to_dedup   = chunk[(~chunk.index.isin(last_chunk.index.tolist()))]

                ## Set duplicate index 
                dups_ix    = to_dedup.duplicated(subset=drop_by)
                uniq_ix    = ~dups_ix

                ## Gather the duplicates and unique rows
                dups_rows  = to_dedup[(dups_ix)]
                uniq_rows  = to_dedup[(uniq_ix)]

                ## Count intra and inter dups
                chunk_dups = dups_rows.Chrn1.count()
                inter_dups = dups_rows.Inter.sum()
                intra_dups = chunk_dups - inter_dups

                ## Add to the totals 
                interdup_counts += inter_dups 
                intradup_counts += intra_dups

                ## Save out the non duplicates
                if sum(uniq_ix):
                    uniq_rows.to_csv(output_path,mode='a' if uniq_appending else 'w',index=False,header =not uniq_appending,sep=hicsep)
                    uniq_appending = True

                ## Save out the duplicates, if we are 
                if keep_dups and sum(dups_ix):
                    dups_rows.to_csv(dedupe_path,mode='a' if dups_appending else 'w',index=False,header =not dups_appending,sep=hicsep)
                    dups_appending = True

            ## Finish off the last row
            if len(last_chunk):
                ## Set the last row 
                to_dedup = last_chunk
                ## Set duplicate index 
                dups_ix    = to_dedup.duplicated(subset=drop_by)
                uniq_ix    = ~dups_ix

                ## Gather the duplicates and unique rows
                dups_rows  = to_dedup[(dups_ix)]
                uniq_rows  = to_dedup[(uniq_ix)]

                ## Count intra and inter dups
                chunk_dups = dups_rows.Chrn1.count()
                inter_dups = dups_rows.Inter.sum()
                intra_dups = chunk_dups - inter_dups

                ## Add to the totals 
                interdup_counts += inter_dups 
                intradup_counts += intra_dups

                ## Save out the non duplicates
                if sum(uniq_ix):
                    uniq_rows.to_csv(output_path,mode='a' if uniq_appending else 'w',index=False,header =not uniq_appending,sep=hicsep)

                ## Save out the duplicates, if we are 
                if keep_dups and sum(dups_ix):
                    dups_rows.to_csv(dedupe_path,mode='a' if dups_appending else 'w',index=False,header =not dups_appending,sep=hicsep)

        ## Remove the large temporary file
        remove(sorted_path)
        remove(tmp_path)

    ## Format new names and print counts
    new_names = ['InterDuplicates', 'IntraDuplicates']
    new_count = [interdup_counts, intradup_counts]
    ## Iterate thru and print the counts to log 
    [print('INFO: %s\t%s'%(a,b)) for a,b in zip(new_names,new_count)]
    ## Print to log
    print("Finished concatonation, sorting, and deduplicating bedpe files ending with %s"%filebackend)
## End of file 
