#!/usr/bin/env python
#################################################
##      Hi-C Chromosome Spliter (v 3.0.0)      ##
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

## Set the description
sort_desc = "Hi-C Chromosome spliter (v %s): Sorts an input txt file (space delimenated) representing Hi-C contacts."%version

## Set help messages 
i_help = "Path (or paths) to Hi-C text files for splitting."
o_help = "Path to output save file."
c_help = "Name of chromosome to parse."

## Bring in help messages from slurpy 
from defaults import mark_help, fileexists, hicsep, juicer_cols, juicer_types, dictzip, G_help, chunks

## Load in the short types form the juicer cols and types
short_columns = juicer_cols[:8]  + ['Chrn1','Chrn2','Qname1']
short_types   = juicer_types[:8] + [   int,   int,     str]

## Set the types of the columns 
short_type_dict = dictzip(short_columns,short_types)

## Load in tools from pysamtools
from pysamtools import writeset

## Bring in pandas 
import pandas as pd 

## Set the positional columns
pos_cols = ['Chrn1','Chrn2','Str1','Str2','Pos1','Pos2']

## ------------------------------------------------------ BODY of MAIN EXECUTABLE --------------------------------------------------------- ##
## If the script is called 
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse

    ## Make the parse
    parser = argparse.ArgumentParser(description = sort_desc)
    ## Add the required arguments
    parser.add_argument("-i", dest="I", type=str,  required=True, help=i_help,  metavar = "./path/to/input.txt")
    parser.add_argument("-o", dest="O", type=str,  required=True, help=o_help,  metavar = "./path/to/output.txt")
    parser.add_argument("-c", dest="C", type=str,  required=True, help=c_help,  metavar = 'chr1')
    parser.add_argument("-g", dest="G", nargs='+', required=True, help =G_help, metavar = 'chr1 chr2 ... chrX')
    ## Add boolean arguemnt
    parser.add_argument("--skip-dedup", dest="mark", help = mark_help, action='store_true')
    ## Set the paresed values as inputs
    inputs = parser.parse_args()

    ## Format input arguments
    inpath      = inputs.I
    savepath    = inputs.O
    chrom       = inputs.C
    chrlist     = inputs.G 
    notmarking  = inputs.mark 

    ## Check our work before we start 
    assert fileexists(inpath), "ERROR: A non-existant input path was passed!"

    ## Set output paths
    dupspath = savepath.split('.txt')[0] + '.duplicates.txt'

    ## Bring in remove 
    from os import remove

    ## remove previous attempts
    remove(dupspath) if fileexists(dupspath) else None 
    remove(savepath) if fileexists(savepath) else None 

    ## Use a with statement open 
    with pd.read_csv(inpath,sep=hicsep,usecols=short_columns,dtype=short_type_dict,chunksize=chunks) as hicloader:
        ## Iterate over the chunks
        cdf = pd.concat([df[(df.Chr1==chrom) & (df.Chr2.isin(chrlist))] for df in hicloader])
        ## Set fragmetn counts 
        cdf['Frag1'] = 0
        cdf['Frag2'] = 1

        ## If we are skipping duplciate marking 
        if not notmarking:
            ## Set the index as the Qname
            cdf.index = cdf.Qname1
            ## Gather the uniqe dataframe by position 
            uni_ix = cdf[pos_cols].drop_duplicates().index
            ## Calcualte the duplicates 
            duplicates  = set(cdf.Qname1.tolist()) - set(uni_ix.tolist())
            ## Write the duplicates out to file
            k = writeset(dupspath,duplicates) if len(duplicates) else None
            ## Gather the unique index 
            cdf = cdf[(cdf.Qname1.isin(uni_ix))]

        ## Left, right sort the chromosome dataframe by the genomic coordinates, take the juicer columns
        sdf = cdf.sort_values(['Chrn1','Chrn2','Pos1','Pos2'])[short_columns[:8]]
        ## Saveout to file
        sdf.to_csv(savepath,header=False,index=False,sep=hicsep)
## End of file 