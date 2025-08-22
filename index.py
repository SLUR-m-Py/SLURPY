#!/usr/bin/env python
#SBATCH --job-name=index.bedpe          ## Name of job
#SBATCH --output=%x.%j.out              ## Name stdout
#SBATCH --error=%x.%j.err               ## Name stderr
#SBATCH --partition=mpi                 ## Size of partitions
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=12              ## Number of tasks to be launched
########################################################################################
    ##      Make an index dataframe of an input BEDPE file from SLUR(M)-py        ##
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
## Load in system
import sys 
## Load in current wd 
from os import getcwd
## append path path
sys.path.append(getcwd()) 

## Load in parameters
from parameters import hicsep,chunksize, ST, Z_help
## Set columns 
short_cols = ['Seqrev1','Rname1','Pos1','Mapq1','Seqrev2','Rname2','Pos2','Mapq2']
to_use_col = ['Rname1','Rname2','Mapq2']
out_cols   = to_use_col[:-1] + ['Start','End']

## Set messages 
desc       = "Generates a csv file representing the row indexes of chromosome contacts/sequenced read pairs."
B_help     = "Path to an input Hi-C text file from SLUR(M)-py in bedpe format or short file from Juicer. Can be gzipped."
G_help     = "Path to a genome map/dataframe/table with chromosome names in first column or path to a reference genome in fasta format."
S_help     = "Boolean flag to index a short file (from Juicer)."
D_help     = "Boolean flag to run in debug mode, printing a few variables and checks to screen."
K_help     = "Boolean flag to parse input bedpe file with dask dataframes."

## if the script is envoked
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse, dask.dataframe as dd, pandas as pd 

    ## Make the parse
    parser = argparse.ArgumentParser(description = desc)
    ## Add the required arguments
    parser.add_argument("-b", dest="B",  type=str,  required=True,   help = B_help  ) 
    parser.add_argument("-g", dest="G",  type=str,  required=False,  help = G_help, default=hicsep     )
    parser.add_argument("-Z", dest="Z",  type=int,  required=False,  help = Z_help, default=chunksize  )
    parser.add_argument("-S", "--short", dest="S",  action = ST,     help = S_help  )
    parser.add_argument("-D", "--debug", dest="D",  action = ST,     help = D_help  )
    parser.add_argument("--by-dask",     dest="K",  action = ST,     help = K_help  )
    ## Set the parsed values as inputs
    inputs = parser.parse_args()
    
    ## Set input vars
    inpath = inputs.B
    gpath  = inputs.G 
    ishort = inputs.S 
    debug  = inputs.D
    bydask = inputs.K

    ## Load in functions from defaults 
    from defaults import readtable, fileexists
    from pysamtools import getchrlist

    ## Check if the file is gzipped
    isgzipped = (inpath.split('.')[-1] == 'gz')
    ## Calcualte the mode of analysis 
    mode = bydask and (not isgzipped) and fileexists(gpath)
    ## Print the info 
    print("INFO: Using dask dataframes to parse input file: %s"%inpath) if (debug and mode) else None

    ## If the file is not gzipped, bydask has been promted, and the gpath is real
    if mode:
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
        df = dd.read_csv(inpath,sep=hicsep,names=short_cols if ishort else None,usecols=to_use_col)
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
        ordered['Start']    = [0] + ordered.Cumcount.tolist()[:-1]
        ordered['End']      = ordered.Cumcount - 1

    else: 
        ## initiate values 
        index_values = []
        ## Open with chunks 
        with pd.read_csv(inpath,sep=hicsep,names=short_cols if ishort else None,usecols=to_use_col,chunksize=chunksize) as chunks:
            for chunk in chunks:
                ## Group the chunk by chromosome pairs and Iterate thru 
                for (c1,c2),c in chunk.groupby(to_use_col[:-1],sort=False):
                    ## Append the index values 
                    index_values.append((c1,c2,c.index.min(),c.index.max()))
        ## set values into a datframe 
        index_values = pd.DataFrame(index_values,columns=out_cols)
        ## group the index values
        groups = index_values.groupby(to_use_col[:-1],sort=False)
        ## Gather stats and ends 
        starts = groups.Start.first().reset_index()
        ends   = groups.End.last().reset_index()
        ## merge into one dataframe 
        ordered = starts.merge(ends)

    ## saveout data
    ordered[out_cols].to_csv(inpath + '.index.csv',header=True,index=False)

## Print to screen 
print("Finished indexing input file: %s"%inpath)
## EOF 