#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
########################################
###    FROMATTING macs3 FUNCTIONS    ###
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
## ------------------------------------------------------------- macs3 Functions ----------------------------------------------------------------- ## 
## Load in pandas 
import pandas as pd 

## Load in scripts dir 
from directories import slurpydir, macs3dir 

## Ftn for formating length parameter
def formatlen(minlen):
    """Formats and returns the max length flag and parameter for a call to macs3."""
    ## Return the len
    return '--min-length %s'%minlen if minlen else ''

## Ftn for formating gap parameter
def formatgap(maxgap):
    """Formats and returns the max gap flag and parameter for a call to macs3."""
    ## Return the gap
    return '--max-gap %s'%maxgap if maxgap else ''

## Ftn for making input into a list of inputs
def makelist(input):
    """Makes an input variable into a list."""
    ## Return the listed input
    return input if (type(input) == list) else [input]

## Ftn for formating inputs
def sjoin(input):
    """Joins an input into a list, seperated by a space."""
    ## Make the input a list 
    inputs = makelist(input)
    ## Make the inputs into a joined list 
    return inputs[0] if (len(inputs) == 1) else ' '.join(inputs)

## Ftn for formatting control
def formatcontrol(incontrols):
    """Formats the input controls for a call to macs3 in chip mode."""
    ## Return the formated input control 
    return f'-c {sjoin(incontrols)}' if incontrols else ''

## Ftn for formatting inputs
def formatinput(inputs):
    """Formats the input bam files for a call to macs3."""    
    ## Return the joined inputs for macs3
    return f'-t {sjoin(inputs)}' if inputs else ''

## Set options for atac-seq opts
atacopts = '--keep-dup all --call-summits -B --SPMR --nolambda'

## Set options for chip-seq opts
narrow_chip = '--keep-dup all --call-summits -B --SPMR'
broad_chip  = '--keep-dup all --broad -B --SPMR'

## reformat the input bed names 
def threebedpe(inbedpe:list) -> list: 
    ## Inaite list, return the new file locations 
    return f'{macs3dir}/{inbedpe.split('/')[-1]}' 

## Write ftn for calling macs3 with atac seq data
def peakattack(bedpe:str,n,report,broad=False,gsize='hs',mg=None,ml=None,extraoptions=None,incontrols=None,outdir=f'./{macs3dir}',mode='BEDPE'):
    """Envokes macs3 callpeak function on for a run of the slurpy pipeline (n) on input bam file in bampe mode using the input genome size (g), maximum gap (ml), and minimum peak length (ml)."""
    ## Format the input options, if the options are set explicitly, and if the input control is set for chip experiment 
    if incontrols and broad:
        opts = broad_chip
    elif incontrols and not broad:
        opts = narrow_chip
    else: ## otherwise this is an atac exp 
        opts = atacopts
    ## Add the additional optsions 
    opts = opts + (' ' + extraoptions if extraoptions else '')
    
    ## Forma tthe conversion commands
    con_coms = [f'{slurpydir}/toshort.py --macs3 -i {bedpe}\n']
    ## Format the macs3 callpeak command
    macs_coms = [f'macs3 callpeak -t {threebedpe(bedpe)} {formatcontrol(incontrols)} -n {n} -g {gsize} -f {mode} --outdir {outdir} {opts} {formatgap(mg)} {formatlen(ml)} 2>> {report}\n', f'{slurpydir}/myecho.py Finished calling peaks in {bedpe} with macs3 {report}\n']
    return con_coms + macs_coms

## Set the narrow peak names
peaknames = ['Chrom','Start','End','Name','Score','Strand','Fold_change','-log10pvalue','-log10qvalue','Sumpos']

## Write ftn for loading in narrow peak file
def loadnarrowpeak(path,peakcols = peaknames,sep='\t') -> pd.DataFrame:
    """Loads in a narrow peak file from macs3."""
    ## Load in narrow peak file and return 
    return pd.read_csv(path,sep=sep,header=None,names=peakcols)

## Wrt ftn for slimming down narrow peak file
def slimdown(inbed,coi=['Chrom','Start','End']):
    """Returns a slimmer bed file from input, returning only columns of interst (COI) and droping duplicate rows."""
    ## Return only the columns of interest
    return inbed[coi].drop_duplicates()

## ---------------------------------------- Write ftns for peak analysis ------------------------------------------------ ## 
## Ftn for calc reads in peaks 
def readsinpeaks(inbam,peaks):
    """Summs the counts of pair-reads within an input peak file."""
    ## Return the counts 
    return sum([len([ r for r in inbam.fetch(contig=row.Chrom,start=row.Start,stop=row.End) if r.is_read1 and not r.is_supplementary]) for i,row in peaks.iterrows()])

## Ftn for calculating unique reads in a bam file 
def readsinbam(inbam):
    """Gathers read names from an input bam file and calculates the unique set of read names."""
    ## Return the length of the set of read names 
    return len(set([r.qname for r in inbam]))

## Ftn for calculating bam map from chromosome 
def makemap(inbam):
    """Generates a chromosome map from an input bam file."""
    ## Return the dataframe 
    return pd.DataFrame(zip(inbam.references,inbam.lengths),columns=['Chrom','Length'])

## ----------------------------------------------- MAIN EXECUTABLE --------------------------------------------------- ## 
## If the library is called as an executable
if __name__ == "__main__":
    ## ---------------------------------------------- DEFAULT VARIABLES ---------------------------------------------------- ## 
    ## Set description of this library and scirpt
    description = 'Calculates the fraction of reads within peaks from input bedpe and narrow peaks files.'

    ## Set defaults
    savename, dplace = './frip.stats.csv', 4

    ## Set help messages
    b_help = "Path to input bedpe file."
    p_help = "Path to input peak (BED) files from macs3."
    s_help = "Path and name of output diagnostic statistics."
    d_help = "Decimal place used to calcualte and save statistics (Default: %s)."%dplace

    ## ------------------------------------------- MODULE LOADING ---------------------------------------------------- ## 
    ## Load in pandas and arg parser
    import argparse, dask.dataframe as dd 
    ## Load in help messages from parameters 
    from parameters import g_help

    ## ------------------------------------------- PARSER SETTING ---------------------------------------------------- ## 
    ## Set parser
    parser = argparse.ArgumentParser(description=description)

    ## Add optional arguments
    parser.add_argument("-b",   dest="b",   required=True,    type=str,   help=b_help                       )
    parser.add_argument("-p",   dest="p",   required=True,    type=str,   help=p_help                       )
    parser.add_argument("-s",   dest="s",   required=False,   type=str,   help=s_help,    default=None      )  
    parser.add_argument("-g",   dest="g",   required=False,   type=int,   help=g_help,    default=False     )
    parser.add_argument("-d",   dest="d",   required=False,   type=int,   help=d_help,    default=dplace    )
    ## Parse the arguments
    args = parser.parse_args()

    ## ---------------------------------------- VARIABLE SETTING ---------------------------------------------------- ## 
    ## Gather inputs 
    bedpe_path, peak_path, save_path, genomesize,  dplace = args.b, args.p, args.s, args.g, args.d, 

    ## Load in dask 
    bedpe = dd.read_csv(bedpe_path,sep='\t',names=['Chrom','Left','Right'],header=None)

    ## Calc total
    total = bedpe.Chrom.count().compute()

    ## Load in narrow peak path, initate read count, and drop ducpliates to unique peaks 
    narrow          = loadnarrowpeak(peak_path)
    narrow['Reads'] = 0
    peaks           = narrow[['Chrom','Start','End','Reads']].drop_duplicates()

    ## Count the number of reads in peaks for each chromosome 
    for chrom,cdf in peaks.groupby('Chrom'):
        ## Set the tmporary creads df for this chromosome 
        creads = bedpe[(bedpe.Chrom==chrom)].compute()

        ## Iterate thru the cdf 
        for i,row in cdf.iterrows():
            ## Set the read count for each 
            peaks.loc[i,'Reads'] = creads[(creads.Left <= row.End) & (creads.Right >= row.Start)].Chrom.count()
    
    ## Calculate statists like the fript score, the number of peaks and sumits 
    fripscore = peaks.Reads.sum()/total
    bp        = (peaks.End - peaks.Start).sum()
    nsummits  = narrow.Chrom.count()
    npeaks    = peaks.Chrom.count()
    
    ## Iniate and fill in list for peak info 
    peak_info = [peak_path.split('/')[-1],total,nsummits,npeaks,fripscore,bp]
    ## Format into a df
    peak_info = pd.DataFrame(peak_info,index=['Peak File','Fragments','Summits','Peaks','FRiP','BP']).T
    ## IF genome size was givven
    if genomesize:
        ## Calculate the perecnt genome
        peak_info['Percent'] = 100*peak_info.BP/genomesize
    ## round the peak data
    peak_info = peak_info.round(dplace)
    ## Save the peak info
    peak_info.to_csv(save_path,index=False)
## End of file 