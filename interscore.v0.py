#!/usr/bin/env python
#SBATCH --job-name=interscore           ## Name of job
#SBATCH --output=%x.%j.out              ## Name stdout
#SBATCH --error=%x.%j.err               ## Name stderr
#SBATCH --partition=tb                  ## Size of partitions
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=1               ## Number of tasks to be launched
##################################################
##      Inter-Chromosomal Compartment Score     ##
##################################################
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## Bring in defaults
from defaults import juicer_cols, juicer_types, hicsep, dictzip, G_help, fileexists, readtable, diagdir

## Bring in seaborn
import seaborn as sns

## Bring in mat plot lib
from matplotlib import pyplot as plt

## Set the short columns
short_columns = juicer_cols[:8]
short_types   = juicer_types[:8]

## Set the types of the columns 
short_type_dict = dictzip(short_columns,short_types)

## Ftn for calculating inter score 
def intscore(df,ci,cj,N):
    ## Interactions between chroms
    N_ij = df[(df.Chr1==ci) & (df.Chr2==cj)].Chr1.count().compute()
    ## All interactions envolving chrom 1
    N_i = df[(df.Chr1==ci) | (df.Chr2==ci)].Chr1.count().compute()
    ## All interaction envolving chrom 2
    N_j = df[(df.Chr1==cj) | (df.Chr2==cj)].Chr1.count().compute()
    ## Return the score
    return round(N_ij/(N *(N_i/N) * (N_j/N)),4)

## Set the description
inter_desc = "Inter-chromosomal score calculator."

## Set help messages 
i_help = "Path to Hi-C text files in short format."
O_help = "The output path to save .csv file with inter-chromosomal scores."

## ------------------------------------------------------ BODY of MAIN EXECUTABLE --------------------------------------------------------- ##
## If the script is called 
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse

    ## Make the parse
    parser = argparse.ArgumentParser(description = inter_desc)

    ## Add the required arguments
    parser.add_argument("-i", dest="I", type=str, required=True, help=i_help, metavar="./path/to/input.txt")

    ## Add optional arguments
    parser.add_argument("-G", "--genomelist", dest="G", type=str, required=False, help = G_help, metavar = './path/to/list.tsv',    default = None)
    parser.add_argument("-O", "--outputpath", dest="O", type=str, required=False, help = O_help, metavar = './path/to/saveout.csv', default = None)

    ## Set the paresed values as inputs
    inputs = parser.parse_args()

    ## Bring in needed mods 
    import dask.dataframe as dd, pandas as pd 

    ## Ftn for loading data
    def hicloader(toload):
        """Loads in a hic file (space delminated) in juicer long format with dask."""
        ## Loads in the dataframe with dask 
        return dd.read_csv(toload,sep=hicsep,names=short_columns,usecols=['Chr1','Chr2'],dtype=short_type_dict)

    ## Format input arguments
    inpath      = inputs.I  ## Input path
    pathtochrom = inputs.G  ## Path to chromosomes
    outputpath  = inputs.O  ## Output path 

    ## Reset output path
    outputpath =  outputpath if outputpath else f'./{diagdir}/'

    ## add the forward slash if needed 
    outputpath = outputpath if (outputpath[-1] == '/') else (outputpath + '/')

    ## Check our work before we start 
    assert fileexists(inpath), "ERROR: A non-existant input path was passed!"

    ## Set the output path
    headpath  = (outputpath + inpath.split('/')[-1].split('.txt')[0]) if outputpath else inpath.split('.txt')[0] 
    outcsv    = headpath + '.inter.score.csv' 
    outpng    = headpath + '.inter.score.png'

    ## Print out the csv and png paths 
    #print(outcsv)
    #print(outpng)
    
    ## load in the data
    df = hicloader(inpath)

    ## If a path was passed 
    if pathtochrom:
        ## Check that it exists 
        assert fileexists(pathtochrom), "ERROR: A non-existant input path for the list of chromosomes was passed!"
        ## Set the list of chromosomes 
        chrlist  = readtable(pathtochrom)[0].values
    else: ## Otherwise take the list of chroms from the dataframe, sort the df by the first chromosome, and Gather the chromlist
        chrlist = df['Chr1'].drop_duplicates().compute().values

    ## Initiate inter score df 
    inter_score_df = pd.DataFrame(index=chrlist,columns=chrlist,dtype=float)

    ## Calc total inter contacts    
    N = df.Chr1.count().compute()

    ## Iterate thru the chromosomes 
    for i,c1 in enumerate(chrlist):
        for j,c2 in enumerate(chrlist):
            ## If the chromosomes are not equal and we are going left to right. 
            if (c1 != c2) and (j>i):
                ## Calc the inter -chrom score 
                ics = intscore(df,c1,c2,N)

                ## Place score in df between c1 and c2 and their mirror pos
                inter_score_df.loc[c1,c2] = ics
                inter_score_df.loc[c2,c1] = ics 
    
    ## Save out the inter score df
    inter_score_df.to_csv(outcsv,index=True,header=True)

    ## Call a fig, set face color to w
    fig,ax = plt.subplots(1,1,figsize=(8,6.5))
    fig.set_facecolor('w')

    ## Call seaborn 
    k = sns.heatmap(inter_score_df,cmap='cividis',vmin=0,vmax=1)

    ## Saveout the figure
    plt.savefig(outpng,bbox_inches='tight',dpi=300)
## End of file 