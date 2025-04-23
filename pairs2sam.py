#!/usr/bin/env python
#SBATCH --job-name=pairstosam            ## Name of job
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=12              ## Number of tasks to be launched
#SBATCH --partition=mpi                 ## Set the partition

## Load in mods 
import pandas as pd, sys 

## Load in current wd 
from os import getcwd
## append path path
sys.path.append(getcwd()) 

## Bring in defaults
from defaults import submitter, fileexists, remove
## Load in parameters
from parameters import ST
## Load in directories
from directories import slurpydir, debugdir

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
description = 'pairs2sam.py converts input pairs file from SLUR(M)-py (.bedpe) to a SAM (tab seperated) file.'
threads = 8

## Define ftn
def formatsam(df:pd.DataFrame) -> pd.DataFrame:
    k1 = df[read1]
    k2 = df[read2]
    k2.columns = k1.columns 
    return pd.concat([k1,k2],axis=0).sort_values('Qname1').reset_index(drop=True)

def makesam(inpath:str,fend:str) -> str:
    return inpath.split(fend)[0]+ '.sam'

## Ftn for renaming files
def makeout(inpath:str) -> str:
    if '.txt' in inpath:
        outpath = makesam(inpath,'.txt')
    elif '.bedpe' in inpath:
        outpath = makesam(inpath,'.bedpe')
    else:
        outpath = inpath + '.sam'
    return outpath

## Ftn for correcting path
def correctpath(inpath:str,outpath) -> str:
    ## REturn outpath
    return outpath if outpath else makeout(inpath)

## Set load error message
loaderror = 'ERROR: The input genome map (%s) was not a tab or space deliminated map of the genome, with two columns representing the chromosome names and lengths!'

## Ftn for loading genome map
def loadmap(inpath:str) -> pd.DataFrame:
    ## Attempt one to load with tab sep 
    gmap = pd.read_csv(inpath,sep=samsep,header=None)
    ## IF the input genome map was not tab deliminated, try space delminiation 
    if gmap.shape[1] < 2:
        ## Load in again withs pace dlmin
        gmap = pd.read_csv(inpath,sep=hicsep,header=None)
        ## Check our work 
        assert gmap.shape[1] >= 2, loaderror%inpath
    ## Return gmap
    return gmap 

## Ftn for makign header
def samheader(df:pd.DataFrame) -> list:
    ## Gather a chromosome list form genome map (df), zip and format header lines 
    return [f'@SQ\tSN:{c}\tLN:{l}\n' for c,l in zip(df[0].tolist(),df[1].tolist())]

## Set help ftns
I_help = "Input path to a .bedpe file from SLURPY HI-C pipeline."
G_help = "Input paht to .txt file (tab or space deliminated) representing map of the genome (chormosome names and lengths)."
O_help = "Path of output sam file. If not given taken from input .bedpe file."
B_help = "Boolean flag to convert sam to bam."
T_help = "Thread count for sam to bam conversion."

## Ftn for calling this scrpt
def bedpetosam(inpath:str,genomepath:str,threads:int,tobam:bool,sname:str,pix=5) -> tuple:
    ## Return the formmatted commdn
    return [f'{slurpydir}/pairs2sam.py -i {inpath} -g {genomepath} -t {threads}' + (' --bam\n' if tobam else '\n')], f'{debugdir}/{pix}E.to.sam.{sname}.log'

## If the script is envoked 
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse
    ## Make the parse
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument("-i", dest="I", type=str,  required=True,   help=I_help, metavar='./path/to/input.bedpe'              ) 
    parser.add_argument("-g", dest="G", type=str,  required=True,   help=G_help, metavar='./path/to/chrom.map'                ) 
    parser.add_argument("-o", dest="O", type=str,  required=False,  help=O_help, metavar='./path/to/out.put.sam', default=None)
    parser.add_argument("-t", dest="T", type=int,  required=False,  help=T_help, metavar='n', default=threads)
    ## Add boolean
    parser.add_argument("--bam", dest="B",  help = B_help,    action = ST)

    ## Parse the inputs from parser obj 
    inputs = parser.parse_args()

    ## Set inputs
    input_path = inputs.I
    chrom_path = inputs.G
    outpt_path = inputs.O
    tobam      = inputs.B
    threads    = inputs.T 

    ## Reset output path if none was given  
    outpt_path = correctpath(input_path,outpt_path)

    ## Load in genome map
    genome_map = loadmap(chrom_path)

    ## Gather header lines
    header_lines = samheader(genome_map)

    ## Write out to file 
    with open(outpt_path,'w') as outfile:
        outfile.writelines(header_lines)
        outfile.close()

    ## Open the pandas as chunks 
    with pd.read_csv(input_path,sep=hicsep,chunksize=chunksize) as chunks:
        ## Iterate thru
        for rows in chunks:
            ## Format sam 
            sam = formatsam(rows)
            ## Add qual field
            sam['Qual'] = '*'
            ## Check work
            sam.to_csv(outpt_path,sep=samsep,header=False,index=False,mode='a')

    ## Converto bam
    if tobam:
        """
        samtools view -@ 8 -b 2401_054S.valid.sam | samtools sort - -@ 8 --write-index -o 2401_054S.valid.bam
        """
        ## Format path to bam file
        bam_path = outpt_path.replace('.sam','.bam')
        ## Format bam command 
        bamcommand = f'samtools view -@ {threads} -b {outpt_path} | samtools sort - -@ {threads} --write-index -o {bam_path}'
        ## Submit to os 
        submitter(bamcommand)

        ## If the bam path exisits, remove the sam file
        remove(outpt_path)  if fileexists(bam_path) else None 

    ## Format out file to print
    ins_file_name = input_path.split('/')[0]
    out_file_name = bam_path.split('/')[0] if tobam else outpt_path.split('/')[0]
    ## Print to log we are finished
    print(f'Finished bedpe pairs file conversion ({ins_file_name}) to {out_file_name}')
## End file 