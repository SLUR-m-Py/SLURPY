#!/usr/bin/env python
#SBATCH --job-name=hictojuice           ## Name of job
#SBATCH --output=%x.%j.out              ## Name stdout
#SBATCH --error=%x.%j.err               ## Name stderr
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=12              ## Number of tasks to be launched
#SBATCH --partition=tb                  ## Set the partition
## bring in mods
import dask.dataframe as dd, pandas as pd 
## Brin in defaults
from gxgcounts import file_end, hicsep
## Load in params
from directories import macs3dir
## Lod in chunk size
from parameters import ST, chunksize, pairs_help, inter_help
"""
Juicer short format:

        'str1','chr1','pos1','frag1','str2','chr2','pos2','frag2'

str:  the strand of the contact, 0 for positive orientation anything else negtive
chr:  the chormosome of the contact
pos:  the position of the contact
frag: the fragment, if using dummy var must be different for the pair 
"""
## Set use columns for short format
short_cols = ['Seqrev1','Rname1','Pos1','Mapq1','Seqrev2','Rname2','Pos2','Mapq2','Inter']

"""
https://liz-fernandez.github.io/HiC-Langebio/04-matrix.html

Pairs format:

        readID chr1 pos1 chr2 pos2 strand1 strand2
"""
## Set use columns for pairs format
pairs_cols = ['Qname1','Rname1','Pos1','Rname2','Pos2','Seqrev1','Seqrev2']

## Generate a pre cursor to a .hic flie
desc = "Converts an input bedpe file (representing Hi-C contacts from SLURPY) to a short formated text file for juicer pre command."
## Set help message
I_help = "Input path to a bedpe file from SLURPY Hi-C pipeline."
B_help = "Format paired end bedpe into a single end bed file for ATAC-seq."
M_help = "Format input bedpe pairs file into macs3 compatible version."

## Def ftn for taking left size of the fragment 
def rowleft(row) -> int:
    return sorted(row)[1]

## Def ftn for taking right size of the fragment 
def rowright(row) -> int:
    return sorted(row)[2]

## Set position col
pos_cols = ['Rname1','Pos1','Pos2','End1','End2']

rpos1 = ['Rname1','Pos1','End1']
rpos2 = ['Rname2','Pos2','End2']
new_cols = ['Chrom','Left','Right','Strand']

## Fnt to make new chunk
def makechunk(df,poscols,strand='+') -> pd.DataFrame:
    ndf           = df[poscols].copy()
    ndf['Strand'] = strand
    ndf.columns   = new_cols
    return ndf 

## If the script is envoked 
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse
    ## Make the parse
    parser = argparse.ArgumentParser(description=desc)
    ## Add the required arguments
    parser.add_argument("-i",           dest="I", type=str, required=True,  help=I_help, metavar='./path/to/input.bedpe') 
    parser.add_argument('--bed',        dest="B", help=B_help,      action = ST)
    parser.add_argument('--bedpe',      dest="M", help=M_help,      action = ST)
    parser.add_argument('--pairs',      dest="P", help=pairs_help,  action = ST)
    parser.add_argument('--inter-only', dest="O", help=inter_help,  action = ST)

    ## Set the paresed values as inputs
    inputs = parser.parse_args()

    ## Set input
    input_path  = inputs.I 
    atacbed     = inputs.B 
    formacs3    = inputs.M
    makepairs   = inputs.P
    getinter    = inputs.O

    ## Check path
    assert file_end in input_path, "ERROR: We expected an input .bedpe file and didn't find that extension in: %s"%input_path

    ## IF this is an atac-seq sample 
    if atacbed:
        ## Forma the output path 
        output_path = f'{macs3dir}/{input_path.split('/')[-1]}' 
        ## Open with chunking 
        with pd.read_csv(input_path,sep=hicsep,usecols=pos_cols,chunksize=chunksize) as chunks:
            ## Iterate thru chunks 
            for i,chunk in enumerate(chunks):

                ## Assign the left and right chunk 
                chunk['Left']  = chunk[pos_cols[1:]].apply(rowleft, axis=1)
                chunk['Right'] = chunk[pos_cols[1:]].apply(rowright,axis=1)
                ## Save out the chunk
                chunk[['Rname1','Left','Right']].to_csv(output_path,header=False,index=False,mode='a' if i else 'w',sep='\t')

        """
        ## Open with chunking 
        with pd.read_csv(input_path,sep=hicsep,usecols=rpos1+rpos2,chunksize=chunksize) as chunks:
            ## Iterate thru chunks 
            for i,chunk in enumerate(chunks):
                ## Gather the positions
                new_chunk = pd.concat([makechunk(chunk,rpos1),makechunk(chunk,rpos2,strand='-')],axis=0).sort_values(by=new_cols)
                ## Save out the new chunk
                new_chunk.to_csv(output_path,header=False,index=False,mode='a' if i else 'w',sep='\t')
        """
    ## If in macs 3 mode 
    elif formacs3:
        ## Forma the output path 
        output_path = f'{macs3dir}/{input_path.split('/')[-1]}' 

        ## Open with chunking 
        with pd.read_csv(input_path,sep=hicsep,usecols=pos_cols,chunksize=chunksize) as chunks:
            ## Iterate thru chunks 
            for i,chunk in enumerate(chunks):
                ## ADd min and max of position columns 
                chunk['Min'] = chunk[pos_cols[1:]].min(axis=1)
                chunk['Max'] = chunk[pos_cols[1:]].max(axis=1)
               
                ## Save out the csv 
                chunk[['Rname1','Min','Max']].to_csv(output_path,header=False,index=False,mode='a' if i else 'w',sep='\t')
    ## Otherwise make pairs if make pairs is past
    elif makepairs:
        ## SEt output path
        output_path = input_path.split(file_end)[0] + '.pairs' 

        ## Load in data
        df = dd.read_csv(input_path,sep=hicsep,usecols=pairs_cols)
        ## Save out data in pairs format, should be tab seperated 
        df[pairs_cols].to_csv(output_path,single_file=True,header=False,index=False,sep='\t')

        ## Print to log
        print("Finished converting input bedpe file (%s) to pairs format."%output_path)
    else:
        ## SEt output path
        output_path = input_path.split(file_end)[0] + '.short' 

        ## Load in data
        df = dd.read_csv(input_path,sep=hicsep,usecols=short_cols)

        ## Remap mapq coulmns as our fragment dummy vars
        df['Mapq1'] = 0
        df['Mapq2'] = 1

        ## Parse inter chromosome
        df = df[(df.Inter>0)] if getinter else df

        ## Save out data
        df[short_cols[:-1]].to_csv(output_path,single_file=True,header=False,index=False,sep=hicsep)

        ## Print to log
        print("Finished converting input bedpe file (%s) to short format."%output_path)
## End of file 