#!/usr/bin/env python
#SBATCH --job-name=pairtosam            ## Name of job
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=12              ## Number of tasks to be launched
#SBATCH --partition=mpi                 ## Set the partition

## Load in mods 
import pandas as pd, sys 
"""
Qname1 Flag1 Rname1 Pos1 Mapq1 Cigar1 Seq1 Chimeric1 Seqrev1 Len1 End1 Chrn1 
Flag2 Rname2 Pos2 Mapq2 Cigar2 Seq2 Chimeric2 Seqrev2 Len2 End2 Chrn2 Distance Minmapq Orientation Inter 

"""
## Define read names 
read1 = ['Qname1','Flag1','Rname1','Pos1','Mapq1','Cigar1','Rname2','Pos2','Distance','Seq1']
read3 = [r.replace('2','3') for r in read1]
read2 = ['Qname1']+[r.replace('3','1') for r in [r.replace('1','2') for r in read3]][1:]

## Set other variables 
hicsep = ' '
samsep = '\t'
chunksize = 500000

## Define ftn
def formatsam(df):
    k1 = df[read1]
    k2 = df[read2]
    k2.columns = k1.columns 
    return pd.concat([k1,k2],axis=0).sort_values('Qname1').reset_index(drop=True)

def makesam(inpath,fend) -> str:
    return inpath.split(fend)[0]+ '.sam'

## Ftn for renaming files
def makeout(inpath) -> str:
    if '.txt' in inpath:
        outpath = makesam(inpath,'.txt')
    elif '.bedpe' in inpath:
        outpath = makesam(inpath,'.bedpe')
    else:
        outpath = inpath + '.sam'
    return outpath

## If the script is envoked 
if __name__ == "__main__":
    ## Set input path 
    input_path = sys.argv[1]
    ## Set output path 
    outpt_path = sys.argv[2] if (len(sys.argv) > 2) else makeout(input_path)

    ## Iniate i
    i = 0
    ## Open the pandas as chunks 
    with pd.read_csv(input_path,sep=hicsep,chunksize=chunksize) as chunks:
        ## Iterate thru
        for rows in chunks:
            ## Format sam 
            sam = formatsam(rows)
            ## Add qual field
            sam['Qual'] = '*'
            ## Check work
            sam.to_csv(outpt_path,sep=samsep,header=False,index=False,mode='a' if i else 'w')
            ## reset i 
            i += 1
## End file 