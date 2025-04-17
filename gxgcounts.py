#!/usr/bin/env python
#SBATCH --job-name=gene.count           ## Name of job
#SBATCH --output=%x.%j.out              ## Name stdout
#SBATCH --error=%x.%j.err               ## Name stderr
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=8               ## Number of tasks to be launched
#SBATCH --partition=tb                  ## Set the partition
## Load in mods 
import pandas as pd, sys, dask.dataframe as dd

## Load in current wd 
from os import getcwd
## append path path
sys.path.append(getcwd()) 

## Set variables used in loading and parsting 
poscols     = ['Rname1','Rname2','Pos1','Pos2']
chunksize   = 1000000
chrom_count = 24
file_end    = '.bedpe'
hicsep      = ' '

## Set gff column names
gffnames = ['Chrom','Source','Feature','Left','Right','Score','Strand','Frame','Attribute']

## Write ftn for loading in a gff file
def loadgff(inpath:str) -> pd.DataFrame:
    """Loads in an input gff file with column names using pandas read csv function."""
    ## Return the gff as a dataframe obj 
    return pd.read_csv(inpath,sep='\t',comment='#',names=gffnames,header=None)

## Write ftn if path is gzipped
def isgzip(inpath:str) -> bool:
    """Checks if the input file is g-zipped (i.e. ends with a gz)."""
    ## Return the check 
    return (inpath.split('.')[-1] == 'gz')

## Check if file is a gff file
def isgff(inpath:str) -> bool: 
    """Check if the file extension is a gff file."""
    ## Reutun boolean on iff this is a gff file 
    return (inpath.split('.')[-1] in ['gff','gtf'])

## Write ftn for laoding in gff file via dask df
def daskgff(inpath:str,coi:str,feature_space:list,gffsep='\t'):
    """Loads in a gff file with column names using dask dataframes."""
    ## Return the dask obj 
    gff = dd.read_csv(inpath,sep=gffsep,comment='#',names=gffnames,header=None)
    ## Return a subset of the gff for this chromosome 
    return gff[(gff.Chrom==coi) & (gff.Feature.isin(feature_space))].compute()

## Ftn for loading a three column bed file with dask df
def loadbed(inpath:str,coi:str) -> pd.DataFrame:
    """Loads in a bed file from path."""
    ## Return the tab seperated file
    tmp = dd.read_csv(inpath,sep='\t',names=['Chrom','Left','Right'],header=None)
    return tmp[(tmp.Chrom==coi)].compute()

## Ftn for loading gzipped hic data
def loadgzip(inpath:str,coi:str) -> pd.DataFrame:
    ## Iterate over chunks 
    tmp = []
    ## With the csv open in pandas 
    with pd.read_csv(inpath,sep=hicsep,chunksize=chunksize) as chunks:
        ## Iteraet over chuns 
        for df in chunks:
            ## Gather by chromosom 
            cdf = df[(df.Rname1==coi)][poscols]
            ## append if we have data 
            tmp.append(cdf) if cdf.shape[0] else None
    ## Concat the hic data and return 
    return pd.concat(tmp)

## Set the short columns
short_columns = ['Str1','Rname1','Pos1','Frag1','Str2','Rname2','Pos2','Frag2']

## Ftn for loading gzipped hic data
def shortload(inpath:str,coi:str) -> pd.DataFrame: 
    ## Iterate over chunks 
    tmp = []
    ## With the csv open in pandas 
    with pd.read_csv(inpath,sep=hicsep,header=None,chunksize=chunksize) as chunks:
        ## Iteraet over chuns 
        for df in chunks:
            df.columns = short_columns
            ## Gather by chromosom 
            cdf = df[(df.Rname1==coi)][poscols]
            ## append if we have data 
            tmp.append(cdf) if cdf.shape[0] else None
    ## Concat the hic data and return 
    return pd.concat(tmp)

## Ftn for load ing wht dask df 
def daskload(inpath:str,coi:str) -> pd.DataFrame:
    ## Make a dask df obj 
    df = dd.read_csv(inpath,sep=hicsep)
    ## Return only the chromosomes hits of interest, and onlty the positoin columsn 
    return df[(df.Rname1==coi)][poscols].compute()

## Set the feature list
feature_list = ['gene']

## Set dsecritpion
desc = "Processes and counts Hi-C contacts between features like genes from a GFF file or rows taken from a .bed file."

## Set help messages
I_help = "Input path to a bedpe file from SLURPY HI-C pipeline."
C_help = "The chromosome to limit analysis to. Used in parallel processing."
F_help = "Path to an input .gff or .bed file."
T_help = "Number of chromosomes expected and used to concatonate final results."
short_help = 'A boolean flag indicating that input Hi-C data is in juicer short format.'
merge_help = 'A boolean flag to activate the final merge protocal, collapsing across all genes.'

## If the script is envoked 
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse
    ## Make the parse
    parser = argparse.ArgumentParser(description =desc)
    ## Add the required arguments
    parser.add_argument("-i", dest="I", type=str,  required=True,   help=I_help, metavar='./path/to/input.bedpe') 
    parser.add_argument("-c", dest="C", type=str,  required=True,   help=C_help, metavar='chrX'                 )
    parser.add_argument("-f", dest="F", type=str,  required=True,   help=F_help, metavar='./path/to/t2t.gff'    )
    parser.add_argument("-t", dest="T", type=int,  required=False,  help=T_help, metavar=24, default=chrom_count)
    ## Add boolean for back version compatibilty
    parser.add_argument("--short", dest="short",  help = short_help,    action = 'store_true')
    parser.add_argument("--merge", dest="merge",  help = merge_help,    action = 'store_true')

    ## Set the paresed values as inputs
    inputs = parser.parse_args()

    ## Set the input text and gff path
    txt_path    = inputs.I     ## Set input text file path
    feat_path   = inputs.F     ## Set input feature file path 
    coi         = inputs.C     ## Set the chormosome of interest
    chrom_count = inputs.T     ## Set the number of chromosomes we are processing in total and other runs 
    short       = inputs.short ## Flag to work from short input file 
    merge       = inputs.merge ## Flag to merge outputs acorss genes from other runs 

    ## Check our ends
    if not short:
        assert file_end in txt_path, "ERROR: The file type was not recognized! Expected a .bedpe file representing hic counts."

    ## Gather head path from input txt file 
    head_path = txt_path.split(file_end)[0]
    ## Set the output path
    out_path = head_path + '.%s.gene.counts.csv'

    ## If this is a gff file
    if isgff(feat_path):
        ## Load in gff 
        cgenes = daskgff(feat_path,coi,feature_list)
        ## Set the name 
        cgenes['Name'] = [row.split('Name=')[-1].split(';')[0] for row in cgenes.Attribute.tolist()]
    ## Otherwise load a already made list of features 
    else: 
        ## Load in feature path 
        cgenes  = loadbed(feat_path,coi)
        ## Set the name for the fetures 
        cgenes['Name'] = ['feat%s'%i for i in cgenes.index.tolist()]

    ## Load in the hic data
    if short:
        chrhic = shortload(txt_path,coi)
    elif isgzip(txt_path) and (not short):
        chrhic = loadgzip(txt_path,coi)
    else:
        chrhic = daskload(txt_path,coi)

    ## Iniate the test 
    test = []
    for i, grow in cgenes.iterrows():
        ## Set gene name 
        gname = grow.Name

        ## Gather the contacts between within the gene
        ghic = chrhic[(chrhic.Pos1 <= grow.Right) & (chrhic.Pos1 >= grow.Left)]
        ## Gather those pairs not within the gene 
        whic = ghic[(ghic.Pos2<grow.Left) | (ghic.Pos2>grow.Right) & (ghic.Rname2==ghic.Rname1)]
        ## Or those between chromosomes 
        bhic =  ghic[(ghic.Rname2!=ghic.Rname1)]

        ## Append if we have data
        if whic.shape[0]:
            test.append((gname,whic,0))
        if bhic.shape[0]:
            test.append((gname,bhic,1))

    ## Initiate the gene by gene list 
    gbyg = []
    ## Iterate over test
    for i,df,j in test:
        ## Itereave over each df in tests 
        for ri,row in df.iterrows():
            ## Gather the position 
            pos  = row.Pos2
            chrm = row.Rname2
            ## Gather the gene they intersect with 
            k = cgenes[(cgenes.Left<=pos) & (cgenes.Right>=pos) & (cgenes.Chrom==chrm)]

            ## IF there are geens we intersex with 
            if k.shape[0]:
                for g in k.Name:
                    #gbyg.append((cgenes.loc[i,'Name'],g))
                    l = sorted([i+'_%s'%coi,g+'_%s'%chrm]) + [j]
                    ## Append to list 
                    gbyg.append(l)

    ## Make into a df
    gbyg = pd.DataFrame(gbyg,columns=['G1','G2','Inter'])
    ## initiate coutns 
    gbyg['Counts'] = 0

    ## Count the interactions between genes 
    gcounts = gbyg.groupby(['G1','G2','Inter']).count().reset_index()
    ## SAveout the data
    gcounts.to_csv(out_path%coi,index=False,header=True)

    ## If we are merging within this run 
    if merge:
        ## Free up memory 
        del gcounts
        del gbyg
        del cgenes 

         ## Bring in glob
        from glob import glob 
        import time

        ## Iniate list of all files 
        allfiles = []

        ## Check if all the chromosomes have been made 
        while len(allfiles) < chrom_count:
            ## Gather the files fro our output path 
            allfiles = sorted(glob(out_path%'*'))
            ## Sleep 5 secdons 
            time.sleep(5)

        ## set the outpath
        fin_out_path = head_path + '.gxg.%s.csv'%chrom_count
        ## Concat the files and saveout via dask dataframes 
        dd.read_csv(out_path%'*').to_csv(fin_out_path,index=False,single_file=True)
        
        ## Remove the previous files
        from os import remove
        from os.path import exists
        ## remove the fiels
        [remove(p) for p in allfiles] if exists(fin_out_path) else None 

    ## print to log
    print("Finished calculating featrue X feature counts for %s."%coi)
## EOF 