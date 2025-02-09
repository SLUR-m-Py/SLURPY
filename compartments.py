#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
########################################
###    Hi-C Compartment Analysis     ###
########################################
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## ------------- Function Defining ------------------ ##
## Import fanc, sys, nupy and and padns 
import fanc, numpy as np, pandas as pd

## Ftn for chetting chormosome eigen vector 
def getchrev(evector:list,cdict:dict,coi:str) -> np.array:
    ## Return the eigen vector 
    return np.array(evector[cdict[coi][0]: cdict[coi][1]])

## Set column names 
chrom_cols = ['Chrom','Start','End','Eigenvalue','Compartment']
    
## Set the description
description = 'Calculates A and B compartmetns from Hi-C using FAN-C.'

## Set defaults
bin_size    = 100000
decplace    = 4

## Set help messages
i_help = 'Path to input .hic file.'
f_help = 'Optional path to an input .fasta (or .fa) file (used to align Hi-C data) for predicting GC content for compartments.'
b_help = 'The binsize or resolution (bp) to use in analysis (default: %s).'%bin_size
d_help = 'Decimal place to save out eigen values (default = %s).'%decplace
x_help = 'List of contigs or chromosomes to exclude from analysis'
P_help = 'Boolean flag to skip plotting compartment scores.'
D_help = 'Boolean flag to run script in verbose debugging mode.'


## ----------------------------------------------- MAIN EXECUTABLE --------------------------------------------------- ## 
## If the library is called as an executable
if __name__ == "__main__":
    ## ------------------------------------------- MODULE LOADING ---------------------------------------------------- ## 
    ## Load in pandas and arg parser
    import argparse 
    ## ------------------------------------------- PARSER SETTING ---------------------------------------------------- ## 
    ## Set parser
    parser = argparse.ArgumentParser(description=description)

    ## Set input arguments
    parser.add_argument("-i", "--inhic-file",   dest="i", required=True,  type=str,   help=i_help,  metavar= './path/to/in.hic') 
    parser.add_argument("-f", "--fasta-path",   dest="f", required=True, type=str,    help=f_help,  metavar='./path/to/in.fasta')

    ## Set default arguments 
    parser.add_argument("-b", "--binsize",      dest="b", required=False, type=int,   help=b_help,  default= bin_size)
    parser.add_argument("-d", "--decimals",     dest="d", required=False, type=int,   help=d_help,  default= decplace)
    parser.add_argument("-x", "--excludes",     dest="x", required=False, type=str,   help=x_help,  default= ['chrM'])

    ## Set boolean vars
    parser.add_argument("--debug",              dest="D",  help=D_help, action = 'store_true')

    ## Parse the arguments
    args = parser.parse_args()

    ## Print we are starting
    print('INFO: Starting A/B compartment correction and analysis.')

    ## Set needed input arguments 
    hic_path = args.i
    fas_path = args.f

    ## Set defaults 
    bin_size = args.b 
    decplace = args.d
    excludes = args.x

    ## Set input boolean vars
    debuging = args.D

    ## Reformat exclude list
    if  type(excludes) == str:
        print(excludes) if debuging else None 
        excludes = excludes.split(' ')
    else:
        print(excludes) if debuging else None 
        pass 

    ## Set the savepaths for the png and the 
    savepath = hic_path.split('/')[-1].split('.hic')[0] + '.compartment.scores.%s.csv'%bin_size
    saveppng = savepath.split('.csv')[0] + '.pdf'

    ## Print the paths 
    print('INFO: Saving results to dataframe: %s'%savepath) if debuging else None

    ## Set the bin str
    binstr = '@%skb'%int(bin_size/1000)

    ## Load in hic file 
    hic_obj = fanc.load(hic_path+binstr)

    ## Gather chromosome bins 
    chrlist = hic_obj.chromosomes()

    ## Set the bin dict of chromosomes
    chrbins = hic_obj.chromosome_bins

    ## Print we are starting fan-c
    print('INFO: Initilizing FAN-C:\n')
    ## Calcualte a and b compartment 
    ab = fanc.ABCompartmentMatrix.from_hic(hic_obj)

    ## Call ev ftn for gather eigen vectors
    ev = ab.eigenvector(genome=fas_path)

    ## Set chrom list
    chrom_dfs = []

    ## Print a line
    print(' ') if debuging else None 

    ## Iterate thru the chromosomes
    for cix,chrom in enumerate(chrlist): 

        ## Gather the bin bounds
        a,b = hic_obj.chromosome_bins[chrom]

        ## Gather bins 
        bins = [0] + list(np.cumsum(np.repeat(bin_size,b - a)))

        ## Calc chrom ev
        chrom_ev = getchrev(ev,chrbins,chrom)

        ## Correct the estimates of open or closed
        corrected = np.round(chrom_ev,decplace)

        ## Set the chromosome df
        chrom_df = pd.DataFrame([np.repeat(chrom,len(corrected)),
                                 bins[:-1],
                                 bins[1:],
                                 corrected,
                                 [(1 if ev >= 0 else -1) for ev in corrected]],
                                 index=chrom_cols).T

        ## Append to list 
        chrom_dfs.append(chrom_df)

    ## Merge the chrom dfs into genom
    genomic_df = pd.concat(chrom_dfs)

    ## Save-out the dataframe
    genomic_df[~(genomic_df.Chrom.isin(excludes))].to_csv(savepath,index=False)

    ## Print a line
    print(' ') if debuging else None 
## End of file 