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
import pandas as pd, dask.dataframe as dd
## Bring in concat
from numpy import concatenate
## Load in params
from parameters import Q_help, map_q_thres, hicsep, error_dist, L_help, E_help, r_help, X_help, Z_help, chunksize
## Bring in ftns from slurpy 
from defaults import fileexists
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

## Set the use column names
use_columns = ['Qname1','Rname1','Rname2','Cigar1','Cigar2','End1','End2','Dangend1','Dangend2','Minmapq','Distance','Inter']

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
    parser.add_argument("-Z", dest="Z", type=int,  required=False,  help=Z_help, metavar='n',         default=chunksize    )
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
    error_dist = inputs.E  ## Set max distance to check for erros 
    library    = inputs.L  ## Set the library used 
    chunksize  = inputs.Z  ## set the chuck size 
    dovetail   = inputs.D  ## Flag to remove dovetail reads
    
    ## Check that we have a file
    assert fileexists(bedpe_path), "ERROR: No input file ( %s ) was found!"%bedpe_path

    ## Gather the restriciton sites and dangling ends 
    restriciton_sites, dangling_ends = returnsite(library)

    ## Iniate paths
    not_usede_paths = []
    too_check_paths = []
    hic_valid_paths = []

    ## Load in the bedpe file using chunks 
    with pd.read_csv(bedpe_path,sep=hicsep,chunksize=chunksize,usecols=use_columns) as chunks:
        ## Iteraet over the chunks 
        for i,bedpe in enumerate(chunks):
            ## Initate the lists of dataframes
            not_used = []
            ## Preset error 
            bedpe['Error'] = '' 

            ## Set output paths
            not_usede_path = makeoutpath(bedpe_path,f'erros.{i}')
            hic_valid_path = makeoutpath(bedpe_path,f'tohic.{i}')

            ## Append to path
            not_usede_paths.append(not_usede_path)
            hic_valid_paths.append(hic_valid_path)

            ## Gather the unmapped reads if any 
            unmapped = bedpe[(bedpe.Rname1=='*') | (bedpe.Rname2=='*') | (bedpe.Cigar1=='*') | (bedpe.Cigar2=='*') | (bedpe.Inter<0)].copy()
            unmapped['Error'] = 'unmapped'
            ## Append to not_used
            not_used.append(unmapped) if unmapped.shape[0] else None 
            ## Drop the unmapped 
            bedpe.drop(unmapped.index,axis=0,inplace=True)

            ## Set dangling ends, those that we know have dangling ends, and set error 
            dang_ends = bedpe[(bedpe.Dangend1>0) | (bedpe.Dangend2>0)].copy()
            dang_ends['Error'] = 'dangend'
            ## Append to not used
            not_used.append(dang_ends) if dang_ends.shape[0] else None 
            ## Drop the dang ends
            bedpe.drop(dang_ends.index,axis=0,inplace=True)

            ## Remove dovetailed reads, if doing so 
            if dovetail:
                ## Gather dovetailed reads, set error message 
                dovetailed = bedpe[(bedpe.End1>=bedpe.Pos2) & (bedpe.Inter==0)].copy()
                dovetailed['Error'] = 'dovetailed'
                ## append to the not used
                not_used.append(dovetailed) if dovetailed.shape[0] else None 
                ## Drop from bedpe chunk
                bedpe.drop(dovetailed.index,axis=0,inplace=True)
            else:
                pass 
            
            if len(excludes):
                ## Gather reads mapping to exclude list of chromosomes
                toexclude = bedpe[(bedpe.Qname1.isin(excludes) | bedpe.Qname2.isin(excludes))].copy()
                toexclude['Error'] = 'excluded'
                ## Appedn to the not used lsit
                not_used.append(toexclude) if toexclude.shape[0] else None 
                ## Drop out the unused apirs
                bedpe.drop(toexclude.index,axis=0,inplace=True)
            else:
                pass 

            if len(includos):
                ## Gather read pairs that we are not include
                toexclude = bedpe[~(bedpe.Qname1.isin(includos) & bedpe.Qname2.isin(includos))].copy()
                toexclude['Error'] = 'excluded'
                ## Appedn to the not used lsit
                not_used.append(toexclude) if toexclude.shape[0] else None 
                ## Drop out the unused apirs
                bedpe.drop(toexclude.index,axis=0,inplace=True)
            else:
                pass 

            ## Remove low quality mapping, set error to low qual
            lowqual = bedpe[(bedpe.Minmapq<minmapq)].copy()
            lowqual['Error'] = 'lowqual'
            ## Append to the not used list
            not_used.append(lowqual) if lowqual.shape[0] else None 
            ## Drop from bedpe chunk
            bedpe.drop(lowqual.index,axis=0,inplace=True)

            ## If restriction sites were passed 
            if restriciton_sites: 
                ## SEt the output path and append to list 
                too_check_path = makeoutpath(bedpe_path,f'tocheck.{i}')
                too_check_paths.append(too_check_paths)
                ## Gather the read pairs we plan to check for intra fragments
                tocheck = bedpe[(bedpe.Distance<error_dist) & (bedpe.Inter==0)].copy()
                ## SAve out the reads to check for intra fragments
                tocheck.to_csv(too_check_path,sep=hicsep,header=True,index=False) if tocheck.shape[0] else None 
                ## Drop these from bedpe
                bedpe.drop(tocheck.index,axis=0,inplace=True)
            else:
                pass 

            ## Concat the not used list so far and save out
            pd.concat(not_used,axis=0).to_csv(not_usede_path,sep=hicsep,header=True,index=False) if len(not_used) else None 
            ## Save out the not used reads and the to check names 
            bedpe.to_csv(hic_valid_path,sep=hicsep,header=True,index=False) if bedpe.shape[0] else None 

    ## Delet the last chunks, we don't need these
    del bedpe, not_used
    ## Set the next chunk
    last_chunk = i + 1
    ## Set the intra qnames
    intra_all_qnames = []

    ## If we are checking rest sites and we have dataframes to check 
    if restriciton_sites and len(too_check_paths):
        ## Bring in the reads to check 
        allcheck = dd.read_csv(too_check_paths,sep=hicsep)
        ## Gather the chromosome list
        chrlist = allcheck.Rname1.unique().compute()

        ## Load in the reference 
        refs_parse = SeqIO.parse(refpath,format='fasta')
        ## Iterate thru ref parser 
        for ref in refs_parse:
            if (ref.id not in chrlist) & (ref.name not in chrlist):
                continue
            ## set the dataframe to check 
            tocheck = allcheck[(allcheck.Rname1==ref.id) | (allcheck.Rname1==ref.name)][check_names].compute()
            ## Initiate the fragment sites
            tocheck['Left1'], tocheck['Right1'], tocheck['Left2'], tocheck['Right2'] = -1, -1, -1, -1

            ## Set the chrom sequence
            chrseq = ref.seq
 
            ## Iterate thru the tupels 
            for (a,b,c) in value_tuple:
                tocheck[a] = tocheck[b].apply(nearestrest, args=[restriciton_sites,chrseq,c])
        
            ## Gather the intra fragment read pairs, those with fragmetns / rest sites bounds that are equal or overlapping 5' to 3'
            intra_qnames = tocheck[(tocheck.Left1==tocheck.Left2) | (tocheck.Right1 == tocheck.Right2) | (tocheck.Right1>=tocheck.Left2)].Qname1.tolist()
            intra_all_qnames.append(intra_qnames)

        ## Concat the intra frag qnames
        intra_all_qnames = concatenate(intra_all_qnames)

        ## Set the qnames 
        errors = allcheck[(allcheck.Qname1.isin(intra_all_qnames))]
        errors['Error'] = 'intrafrag'
        ## Set path to save out the errors 
        out_error_path = makeoutpath(bedpe_path,f'errors.{last_chunk}')
        ## append to path
        not_usede_paths.append(out_error_path)
        ## SAve out the errors
        errors.to_csv(out_error_path,sep=hicsep,single_file=True,index=False)

        ## Gather the valid
        valid = allcheck[(~allcheck.Qname1.isin(intra_all_qnames))]
        ## SEt the output path to save the valid contacts
        out_valid_path = makeoutpath(bedpe_path,f'tohic.{last_chunk}')
        ## Append to list of paths of valid
        hic_valid_paths.append(out_valid_path)
        ## Save out the valid contacts
        valid.to_csv(out_valid_path,sep=hicsep,single_file=True,index=False)
    else:
        pass 

    ## Set output file path, gather the not used dataframes and save out as sicnle file
    not_usede_path = makeoutpath(bedpe_path,'notused')
    dd.read_csv(not_usede_paths,sep=hicsep).to_csv(not_usede_path,sep=hicsep,single_file=True,index=False)

    ## Gather the valid dataframes
    bedpe = dd.read_csv(hic_valid_paths,sep=hicsep)
    ## Compute a sorted list of valid chromosomes by reference indexed number
    chromosomes = sorted(bedpe.Rname1.unique().compute())
    ## Iterate thur and seperat on chromosomes, and saveout the contact. NOTE: We have to do this (rather than groupby) b/c dask dataframes can't precompute groups
    [bedpe[(bedpe.Rname1==chrom)].to_csv(makeoutpath(bedpe_path,'valid.%s'%chrom),sep=hicsep,single_file=True,index=False) for chrom in chromosomes]
    
    ## print to log
    print("Finished filtering and splitting bedpe file: %s"%bedpe_path)
## End of file 