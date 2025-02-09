#!/usr/bin/env python
#####################################################################
##          Filter bedpe file representing Hi-C contacts           ##
#####################################################################
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## Set the descriptions
filter_desc = "Filters an input bedpe file (space deliminated) representing Hi-C contacts from a paired-end Hi-C experiment. Seperates file on chromosmes, unused contacts, unmapped reads, and errors."

## ----------------------------------- MODULE LOADING ------------------------------------ ##
## Bring in pandas
import dask.dataframe as dd
## Load in params
from defaults.parameters import Q_help, map_q_thres, hicsep, error_dist, L_help, E_help, r_help, X_help
## Bring in ftns from slurpy 
from defaults.defaults import fileexists
## Bring in seqIO
from Bio import SeqIO

## ----------------------------------- FUNCTIONS DEFINING --------------------------------- ##
## Ftn for setting output path
def makeoutpath(inpath,coi) -> str: 
    ## Gather head path 
    head_path = inpath.split('.bedpe')[0]
    ## Return the head path 
    return head_path + f'.{coi}.bedpe' ## NOTE: We didn't go with a zfill here b/c we don't know how many chromosomes / contigs a user will have

## Ftn for defninging restriciton sites 
def returnsite(enzyme) -> tuple: 
    """Returns the restriction site sequences and dangaling sequences based on input Hi-C library enzyme."""
    ## If we are working with the arima kit (which we are)
    if (enzyme.lower() == 'arima'):
        restsites, dangsites = ['GATC', 'GAATC', 'GAGTC', 'GATTC', 'GACTC'], ['GATC', 'AATC', 'AGTC','ATTC', 'ACTC']
    ## The MobI , DpnII, or Sau3AI enzymes
    elif enzyme.lower() in ['mboi', 'dpnii','sau3ai']:
        restsites, dangsites = ['GATC'], ['GATC']
    ## The HindIII enzyme
    elif (enzyme.lower() == 'hindiii'):
        restsites,dangsites = ['AAGCTT'],['AGCT']
    else: ## otherwise return none
        restsites, dangsites = False, False 
    ## Return the sites
    return restsites, dangsites

## Check the restriction sites 
def checkrests(gdf,restsites:list,refpath:str) -> tuple[list,list]:
    ## Iniate counts and index 
    rest_count, rest_index = [], []
    ## Gather the chroms 
    chroms = gdf.Rname1.unique()
    for ref in SeqIO.parse(refpath,format='fasta'):
        if (ref.id in chroms) or (ref.name in chroms):
        ## Group by the chromosome names
            df = gdf[(gdf.Rname1==ref.id) | (gdf.Rname1==ref.name)]
            ## Itreat thru the rows
            for rix,row in df.iterrows():
                ## Count the restriciton sites 
                k = 0
                for rsite in restsites:
                    ## calc start and end 
                    fragm_start = min([row['Pos1'],row['Pos2']]) + len(rsite)
                    fragm_end   = max([row['End1'],row['End2']]) - len(rsite)
                    k += str(ref.seq[fragm_start:fragm_end]).upper().count(rsite.upper())
                ## Append results for the row 
                rest_count.append(k)
                rest_index.append(rix)
    ## return nothing 
    return rest_index,rest_count

## Calc the max len fo resistance sites 
def maxrestlen(restsites) -> int: 
    return max([len(r) for r in restsites])

## Insure the rest sites are upper case 
def upperrest(restsites):
    ## Make sure they are uper
    return [r.upper() for r in restsites] ## NOTE: We don't need to search the space of invereted sequences # + [r[::-1] for r in restsites]

## Gather the return sequecnes, bed upper case ## NOTE: "find" ftn on strings is case sensitive!!!!!
def returnseq(chrseq,p1,p2,isleft):
    ## Return the reverse comp of seq if itis leftward of site
    return str(chrseq[p1:p2].reverse_complement() if isleft else chrseq[p1:p2]).upper()

## Ftn for finidnign sites
def findrest(seq,restsites):
    ## initate distandce aaway 
    dists_away = []
   
    ## Iterate thru restsites 
    for r in upperrest(restsites):
        k = seq.find(r)
        if k > -1:
            dists_away.append(k+len(r))
            
    ## Return lsit 
    return dists_away

## Ftn for formating distance as int
def findmindist(pos,indexes,isleft):
    ## Calc index and if it is left
    if len(indexes) and isleft:
        newpos = int(pos - min(indexes))
    ## Calc index and is not left 
    elif len(indexes) and not isleft:
        newpos = int(pos+min(indexes))
    ## empty indexes were passed
    else:
        newpos = -1 
    ## Return the index plus or minus pos
    return newpos

## Ftn for modyfying pos 
def modifypos(pos,rlen,isleft):
    ## Take the pos, adding the rest site len or minius if right side of read
    return int(pos+rlen) if isleft else int(pos-rlen)
    
## Clac left test, set pads 
def nearestrest(pos,restsites,chromseq,isleft,pads=[100,500,1000,2000,5000,10000]) -> int:
    ## GAther max len of restsits and modify the starting positoin based on max rest site
    modpos = modifypos(pos,maxrestlen(restsites),isleft)

    ## Iterate thru string sizes
    for p in pads:
        ## Gather the string of DNA to the left, by pad, of position, and rev comp
        theseq = returnseq(chromseq,modpos-p,modpos,isleft) if isleft else returnseq(chromseq,modpos,modpos+p,isleft)
        ## initate distandce aaway 
        dists_away = findrest(theseq,restsites)
        ## IF we found some the first try, or second, break 
        if len(dists_away):
            break
            
    ## Return the minimum index of distances if we have onces greater than -one pluse mpos 
    return findmindist(modpos,dists_away,isleft) 

## Check ends ftn 
def endcheck(read:str,strand:int,dangends:list) -> int:
    ## Checks the presesne of deanglign ends 
    return sum([read.upper().endswith(d) for d in dangends] if strand else [read.upper().startswith(d) for d in dangends])

## Ftn for checking if a ref obj is in list
def inref(refobj,inlist,outlist):
    chrom_id,chrom_name = refobj.id, refobj.name
    ## Return the id, name, and boolean
    return chrom_id,chrom_name,(chrom_id not in inlist) or (chrom_id in outlist) or (chrom_name not in inlist) or (chrom_name in outlist)

## ------------------------------------- Define Variables ------------------------------------------------- ##
B_help = "Path to an input bed paired-end file (bedpe)."
I_help = "List of chormosomes/contigs to only include in analysis"
dove_help = "Boolean flag to remove dovetailed paired-end reads (paired reads with overlapping mapped coordiantes) from analsyis (Default: is to remove these)."

## Set check names 
check_names = ['Rname1','Pos1','Pos2','End1','End2','Qname1']

## set the value tuples 
value_tuple = list(zip(['Left1','Right1','Left2','Right2'],['Pos1','End1','Pos2','End2'],[True,False,True,False]))


## -------------------------------------- MAIN EXECUTABLE -------------------------------------------------- ##
## if the script is envoked
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse
    ## Make the parse
    parser = argparse.ArgumentParser(description = filter_desc)

    ## Add the required arguments
    parser.add_argument("-b", dest="B", type=str,  required=True,   help=B_help, metavar='./path/to/input.bedpe'           )   
    parser.add_argument("-r", dest="R", type=str,  required=True,   help=r_help, metavar='./path/to/ref.fasta'             )
    
    ## Add optional args
    parser.add_argument("-e", dest="E", type=int,  required=False,  help=E_help, metavar='n',         default=error_dist   )
    parser.add_argument("-l", dest="L", type=str,  required=False,  help=L_help, metavar='Arima',     default='Arima'      )
    parser.add_argument("-q", dest="Q", type=int,  required=False,  help=Q_help, metavar='n',         default=map_q_thres  )
    parser.add_argument("-x", dest="X", nargs='+', required=False,  help=X_help, metavar='chrM',      default=['chrM']     )
    parser.add_argument("-i", dest="I", nargs='+', required=False,  help=I_help, metavar='chr1 chr2', default=[]           )
    
    ## Add boolean 
    parser.add_argument("--dovetails",  dest="D",  help = dove_help,    action = 'store_true')

    ## Set the paresed values as inputs
    inputs = parser.parse_args()

    ## Set inputs 
    bedpe_path = inputs.B  ## Path to bedpe file 
    refpath    = inputs.R  ## Set path to reference, needed to check restriction 
    excludes   = inputs.X  ## List of chromosomes/contigs to exclude 
    includos   = inputs.I  ## List of chormosomes to only includ
    minmapq    = inputs.Q  ## minimum mapping quality 
    errordist  = inputs.E  ## Set max distance to check for erros 
    library    = inputs.L  ## Set the library used 
    dovetail   = inputs.D  ## Flag to remove dovetail reads

    ## Check that we have a file
    assert fileexists(bedpe_path), "ERROR: No input file ( %s ) was found!"%bedpe_path

    ## Add the unmapped wild card to exclude list
    excludes.append('*')

    ## Load in the bedpe file via dask dataframes
    bedpe = dd.read_csv(bedpe_path,sep=hicsep)

    ## Gather the unmapped
    unmapped_qnames = bedpe[(bedpe.Rname1=='*') | (bedpe.Rname2=='*') | (bedpe.Cigar1=='*') | (bedpe.Cigar2=='*') | (bedpe.Pos1==0) | (bedpe.Pos2==0) | (bedpe.Inter<0)].Qname1.compute().tolist()

    ## Set dangling ends qnames, those that we know have dangling ends
    dang_end_qnames = bedpe[(bedpe.Dangend1>0) | (bedpe.Dangend2>0)].Qname1.compute().tolist()
    
    ## Calcualte the the reads that are dove tailed, NOTE if we do the above, this is pointless 
    dove_qnames = bedpe[(bedpe.End1>=bedpe.Pos2) & (bedpe.Inter==0)].Qname1.compute().tolist() if dovetail else []

    ## Gather the restriciton sites and dangling ends 
    restriciton_sites, dangling_ends = returnsite(library)

    ## If a library was passed, we need to check fragments for erros 
    if restriciton_sites: 
        ## Gather those intra chromosome pairs to check 
        to_check = bedpe[(bedpe.Distance<errordist) & (bedpe.Inter==0) & (~bedpe.Qname1.isin(unmapped_qnames+dang_end_qnames+dove_qnames))][check_names].compute().reset_index(drop=True)
        ## Initiate the fragment sites
        to_check['Left1'], to_check['Right1'], to_check['Left2'], to_check['Right2'] = -1, -1, -1, -1

        ## Gather the chromosome list
        chrlist = [c for c in to_check.Rname1.unique().tolist() if ((c in includos) and (not c in excludes))]

        ## Load in the reference 
        refs_parse = SeqIO.parse(refpath,format='fasta')

        ## Iterate thru the reference paths 
        for ref in refs_parse:
            ## Define the chromosome names and id, calc boolean
            chrom_id,chrom_name,chrom_cont = inref(ref,chrlist,excludes)
            ## If come across a chrom we don't recoginze 
            if chrom_cont:
                continue
           
            ## Set the chrom sequence
            chrseq = ref.seq
            ## Gather the contacts mapping for this chrom
            cdf = to_check[(to_check.Rname1 == chrom_id) | (to_check.Rname1==chrom_name)]

            ## Iterate thru the tupels 
            for (a,b,c) in value_tuple:
                to_check.loc[cdf.index,a] = cdf[b].apply(nearestrest, args=[restriciton_sites,chrseq,c])
        
        ## Gather the intra fragment read pairs, those with fragmetns / rest sites bounds that are equal or overlapping 5' to 3'
        intra_qnames = to_check[(to_check.Left1==to_check.Left2) | (to_check.Right1 == to_check.Right2) | (to_check.Right1>=to_check.Left2)].Qname1.tolist()
    else: ## Otherwisedo nothing 
        intra_qnames = []

    ## Calcualte the the reads that are dove tailed, NOTE if we do the above, this is pointless 
    dove_qnames = to_check[(to_check.End1>=to_check.Pos2)].Qname1.tolist() if dovetail else []

    ## Clac low quality qs
    lowqual_qnames = bedpe[(bedpe.Minmapq<minmapq)].Qname1.compute().tolist()

    ## combine all the errors
    error_qnames = dang_end_qnames + dove_qnames + intra_qnames + unmapped_qnames + lowqual_qnames

    ## Gather all other Qnames and Rnames
    keep_qnames = bedpe[(~bedpe.Qname1.isin(error_qnames)) & (bedpe.Rname1.isin(includos)) & (~bedpe.Rname1.isin(excludes))].Qname1.compute().tolist()

    ## Gather the not include list 
    noinclude_qnames = bedpe[(~bedpe.Qname1.isin(error_qnames)) & (~bedpe.Qname1.isin(keep_qnames))].Qname1.compute().tolist()

    ## Format boolean indexes 
    validix = bedpe.Qname1.isin(keep_qnames)
    notused = bedpe.Qname1.isin(noinclude_qnames)

    ## Compute a sorted list of valid chromosomes by reference indexed number
    chromosomes = sorted(bedpe[validix].Rname1.unique().compute())
    #chrom_n  = sorted(bedpe[validix].Chrn1.unique().compute())
    #print(chrom_n)
    ## Iterate thur and seperat on chromosomes, and saveout the contacts
    ## NOTE: We have to do this (rather than groupby) b/c dask dataframes can't precompute groups
    [bedpe[(validix) & (bedpe.Rname1==chrom)].to_csv(makeoutpath(bedpe_path,'valid.%s'%chrom),sep=hicsep,single_file=True,index=False) for chrom in chromosomes]

    ## Save out the unused reads 
    bedpe[notused].to_csv(makeoutpath(bedpe_path,'notused'),sep=hicsep,single_file=True,index=False) 

    ## Set error lables
    error_labels = ['Unmapped',       'Danglingend',  'Dovetailed', 'Intrafragment',   'Lowquality' ]
    error_groups = [unmapped_qnames, dang_end_qnames,  dove_qnames,   intra_qnames,   lowqual_qnames]
    ## Set empty list
    stored = []
    ## Save out the errors by group
    for i, (elabel,egroup) in enumerate(zip(error_labels,error_groups)):
        ## If the group has no length 
        if not len(egroup):
            continue
        ## Gather the rrors 
        tmp = bedpe[(bedpe.Qname1.isin(egroup)) & (~bedpe.Qname1.isin(stored))].compute()
        ## Add a label and save 
        tmp['Error'] = elabel
        tmp.to_csv(makeoutpath(bedpe_path,'errors'),mode='a' if i else 'w',sep=hicsep,index=False,header=not i)
        ## Update stored list 
        stored = stored + egroup
    
    ## print to log
    print("Finished filtering and splitting bedpe file: %s"%bedpe_path)
## End of file 