#!/usr/bin/env python
#SBATCH --job-name=index.bedpe          ## Name of job
#SBATCH --output=%x.%j.out              ## Name stdout
#SBATCH --error=%x.%j.err               ## Name stderr
#SBATCH --partition=mpi                 ## Size of partitions
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=12              ## Number of tasks to be launched
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## Load in system
import sys 
## Load in current wd 
from os import getcwd
## append path path
sys.path.append(getcwd()) 
########################################################################################
    ##      Make an index dataframe of an input BEDPE file from SLUR(M)-py        ##
########################################################################################
short_cols = ['Seqrev1','Rname1','Pos1','Mapq1','Seqrev2','Rname2','Pos2','Mapq2']
to_use_col = ['Rname1','Rname2','Mapq2']
ST         = 'store_true'
desc       = "Generates a csv file representing the row indexes of chromosome contacts/sequenced read pairs."
B_help     = "Path to an input Hi-C text file from SLUR(M)-py in bedpe format or short file from Juicer. Must be unzipped."
G_help     = "Path to a genome map/dataframe/table with chromosome names in first column or path to a reference genome in fasta format."
S_help     = "Boolean flag to index a short file (from Juicer)."
D_help     = "Boolean flag to run in debug mode, printing a few variables and checks to screen."

## if the script is envoked
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse, dask.dataframe as dd, pandas as pd 

    ## Make the parse
    parser = argparse.ArgumentParser(description = desc)
    ## Add the required arguments
    parser.add_argument("-b", dest="B",  type=str, required=True,  help = B_help  ) 
    parser.add_argument("-g", dest="G",  type=str, required=True,  help = G_help  )
    parser.add_argument("-S", "--short", dest="S",  action = ST,   help = S_help  )
    parser.add_argument("-D", "--debug", dest="D",  action = ST,   help = D_help  )
    ## Set the parsed values as inputs
    inputs = parser.parse_args()
    
    ## Set input vars
    inpath = inputs.B
    gpath  = inputs.G 
    ishort = inputs.S 
    debug  = inputs.D

    ## Load in functions from defaults 
    from defaults import readtable
    from pysamtools import getchrlist

    ## If it is a fasta file 
    if gpath.split('.')[-1] in ['fasta','fa']:
        chrlist = getchrlist(gpath)
    else: 
        ## Load in genome tab
        genometab = readtable(gpath)
        ## Set chromosome list from table 
        chrlist = genometab[0].tolist()

    ## Print messages if debuging 
    print(len(chrlist)) if debug else None
    print(chrlist[:5]) if debug else None 

    ## Load in the dataframe 
    df = dd.read_csv(inpath,sep=' ',names=short_cols if ishort else None,usecols=to_use_col)

    ## Gather the row counts by chromosome combination
    testix = df.groupby(to_use_col[:-1],sort=False).count().compute().reset_index()
    ## Reorder the above dataframe by chromosome order 
    ## Initate ordered list
    ordered = []
    ## Itereat thru chormooems
    for i,c1 in enumerate(chrlist):
        for j,c2 in enumerate(chrlist):
            ## Only, left right sorted chromoomes, that are greater than the next are used
            if  i > j:
                pass
            else:
                ## Gather the chromosome info 
                tmp = testix[(testix.Rname1==c1) & (testix.Rname2==c2)]
                if  tmp.Mapq2.count() == 1:
                    ordered.append(tmp)
                elif tmp.Mapq2.count() > 1:
                    print("ERROR: The given files is not genomically left,right, top to bottom sorted!")
                else:
                    print(f'WARNING: No combination of {c1} and {c2}') if debug else None 
    ## Concat the ordered df               
    ordered = pd.concat(ordered)

    ## Calcualte cumulative count, stat and ends 
    ordered['Cumcount'] =  ordered.Mapq2.cumsum()
    ordered['Start'] = [0] + ordered.Cumcount.tolist()[:-1]
    ordered['End'] = ordered.Cumcount - 1

    ## saveout data
    ordered[to_use_col[:-1] + ['Start','End']].to_csv(inpath + '.index',header=True,index=False)
## Print to screen 
print("Finished indexing input file: %s"%inpath)
## EOF 