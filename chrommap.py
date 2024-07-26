
## Load in functions and variables from defaults
from defaults import aligndir, fileexists, readtable, readann

## Bring in chrom df from py-sam tools
from pysamtools import chromdf

## Set defalut messaging
chromgathering = 'INFO: Gathering chromosomes for processing.' 

## Set ftn for gathering chromosome map
def gathering(path_to_ref,path_to_chrom,list_exclude) -> tuple:
    """Formats a map of the refernece genome, which includes chromosome names and lengths"""
    ## Set fai path
    ann_path = path_to_ref + '.ann'
    ## Set new path to chrom
    newpath_to_chrom = aligndir + '/' + path_to_ref.split('/')[-1].split('.fa')[0] + '.txt'
    ## If a path of chrom file was passed and exists 
    if path_to_chrom and fileexists(path_to_chrom): 
        ## Set the list of chromosomes 
        chrbed = readtable(path_to_chrom)

    elif fileexists(ann_path):
        ## Patch path
        path_to_chrom = newpath_to_chrom
        ## Load in the fai file fromt he reference 
        chrbed = readann(ann_path)
        ## SAve out the bed file 
        chrbed[[0,1]].to_csv(path_to_chrom,sep=' ',index=False,header=False)

    ## Make a dataframe of the files 
    else: ## Bring in the list of chromosomes from the fasta file
        ## Patch chrom 
        path_to_chrom = newpath_to_chrom
        ## IF, the path to chrom form the reference file path already exists 
        if fileexists(path_to_chrom):
            chrbed = readtable(path_to_chrom)
        else:
            ## Gather tupes of chromosome ids and lengths 
            chrbed = chromdf(path_to_ref)
            ## Save out the chrbed as an .txt file for next time or further analysis 
            chrbed.to_csv(path_to_chrom,sep=' ',index=False,header=False)
      
    ## Gather the list of chromosomes from first column from the chr bed df
    chrlist = chrbed[0].tolist()

    ## Remove mito and other unwanted chromosomes
    chrlist = [c for c in chrlist if (c not in list_exclude)]

    ## Calcualte the genome size 
    gsize = chrbed[(chrbed[0].isin(chrlist))][1].sum()

    ## Return the chromosome list and genomesize
    return chrlist, gsize 