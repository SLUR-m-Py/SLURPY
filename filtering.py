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
## Load in params
from parameters import chunksize, map_q_thres, hicsep, error_dist
## Bring in ftns from slurpy 
from defaults import fileexists
## Bring in seqIO
from Bio import SeqIO
from Bio.Seq import Seq

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

## ------------------------------------- Define Variables ------------------------------------------------- ##
B_help     = "Path to an input bed paired-end file (bedpe)."
I_help     = "List of chormosomes/contigs to only include in analysis"

## Import help messages
from parameters import ST, L_help, E_help, r_help, X_help, Z_help, intra_help, Q_help, m_help, hicex_help, dove_help

## Set check names 
check_names = ['Rname1','Pos1','Pos2','End1','End2','Qname1','Distance']

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
    parser.add_argument("-b", dest="b", type=str,  required=True,   help=B_help, metavar='./path/to/input.bedpe'           )   
    parser.add_argument("-r", dest="r", type=str,  required=True,   help=r_help, metavar='./path/to/ref.fasta'             )
    
    ## Add optional args
    parser.add_argument("-e", dest="e", type=int,  required=False,  help=E_help, metavar='n',         default=error_dist   )
    parser.add_argument("-l", dest="l", type=str,  required=False,  help=L_help, metavar='Arima',     default='Arima'      )
    parser.add_argument("-q", dest="q", type=int,  required=False,  help=Q_help, metavar='n',         default=map_q_thres  )
    parser.add_argument("-x", dest="x", nargs='+', required=False,  help=X_help, metavar='chrM',      default=['chrM']     )
    parser.add_argument("-i", dest="i", nargs='+', required=False,  help=I_help, metavar='chr1 chr2', default=[]           )
    parser.add_argument("-Z", dest="Z", type=int,  required=False,  help=Z_help, metavar='rown',      default=chunksize    )
    parser.add_argument("-M", dest="M", type=int,  required=False,  help=m_help, metavar='n (bp)',    default=0            )
    ## Add boolean 
    parser.add_argument("--dedovetail",     dest="D",  help = dove_help,    action = ST)
    parser.add_argument("--intra-only",     dest="I",  help = intra_help,   action = ST) 
    parser.add_argument("--hicexplorer",    dest="E",  help = hicex_help,   action = ST)

    ## Set the paresed values as inputs
    inputs = parser.parse_args()

    ## Set inputs 
    bedpe_path  = inputs.b   ## Path to bedpe file 
    refpath     = inputs.r   ## Set path to reference, needed to check restriction 
    error_dist  = inputs.e   ## Set max distance to check for erros 
    library     = inputs.l   ## Set the library used 
    minmapq     = inputs.q   ## minimum mapping quality 
    excludes    = inputs.x   ## List of chromosomes/contigs to exclude 
    includos    = inputs.i   ## List of chormosomes to only includ
    chunksize   = inputs.Z   ## set the chuck size 
    dovetail    = inputs.D   ## Flag to remove dovetail reads
    intraonly   = inputs.I   ## FLag to keep only intra
    max_dist    = inputs.M   ## Intiger to value to filter paired distances 
    hicexplorer = inputs.E   ## Boolean flag to run in hic explorer mode 
    
    ## Check that we have a file
    assert fileexists(bedpe_path), "ERROR: No input file ( %s ) was found!"%bedpe_path

    ## Gather the restriciton sites and dangling ends 
    restriciton_sites, dangling_ends = returnsite(library)

    ## Iniate paths
    not_usede_paths = []
    too_check_paths = []
    hic_valid_paths = []

    ## Set counts
    unmapped_counts = 0
    dovetail_counts = 0
    excluded_counts = 0
    lowquals_counts = 0
    interchr_counts = 0
    largefrg_counts = 0
    errorfrg_counts = 0
    validhic_counts = 0
    interhic_counts = 0
    intrahic_counts = 0

    ## Load in the bedpe file using chunks 
    with pd.read_csv(bedpe_path,sep=hicsep,chunksize=chunksize) as chunks:
        ## Iteraet over the chunks 
        for i,bedpe in enumerate(chunks):
            i = i + 1
            ## Initate the lists of dataframes
            not_used = []
            ## Preset error 
            bedpe['Error'] = '' 

            ## Set output paths of thoese not used, valid hic contacts, and those to check 
            not_usede_path = makeoutpath(bedpe_path,f'errors.{i}')
            hic_valid_path = makeoutpath(bedpe_path,f'tohic.{i}')
            too_check_path = makeoutpath(bedpe_path,f'tocheck.{i}')

            ## Gather the unmapped reads if any 
            unmapped = bedpe[(bedpe.Rname1=='*') | (bedpe.Rname2=='*') | (bedpe.Cigar1=='*') | (bedpe.Cigar2=='*') | (bedpe.Orientation=='Unmapped')].copy()
            ## Count the unmapped
            unmapped_count = unmapped.Rname1.count()
            ## If we have any unmapped
            if unmapped_count:
                ## Add to the overall count
                unmapped_counts += unmapped_count
                ## Set the error 
                unmapped['Error'] = 'unmapped'
                ## Append to not_used
                not_used.append(unmapped) 
                ## Drop the unmapped 
                bedpe.drop(unmapped.index,axis=0,inplace=True)

            ## Set dangling ends, those that we know have dangling ends, and set error 
            #dang_ends = bedpe[(bedpe.Dangend1>0) | (bedpe.Dangend2>0)].copy()
            #dang_ends['Error'] = 'dangend'
            ## Append to not used
            #not_used.append(dang_ends) if dang_ends.shape[0] else None 
            ### Drop the dang ends
            #bedpe.drop(dang_ends.index,axis=0,inplace=True)

            ## Remove dovetailed reads, if doing so 
            if dovetail:
                ## Gather dovetailed reads, set error message 
                dovetailed = bedpe[(bedpe.End1>=bedpe.Pos2) & (bedpe.Inter==0)].copy()
                ## Set the dovetail count
                dovetail_count = dovetailed.End1.count()
                ## If we have counts 
                if dovetail_count:
                    ## Add to the count of dove tail rads
                    dovetail_counts += dovetail_count
                    ## Set the error 
                    dovetailed['Error'] = 'dovetailed'
                    ## append to the not used
                    not_used.append(dovetailed)
                    ## Drop from bedpe chunk
                    bedpe.drop(dovetailed.index,axis=0,inplace=True)
            else:
                pass 
            
            ## Remove reads mapping to excluded chromosomes/contigs 
            if len(excludes):
                ## Gather reads mapping to exclude list of chromosomes
                toexclude = bedpe[(bedpe.Rname1.isin(excludes) | bedpe.Rname2.isin(excludes))].copy()
                ## Set the count 
                exclude_count = toexclude.Rname1.count()
                ## If we have counts
                if exclude_count:
                    ## Update the count
                    excluded_counts += exclude_count
                    ## Set the error message 
                    toexclude['Error'] = 'excluded'
                    ## Appedn to the not used lsit
                    not_used.append(toexclude) 
                    ## Drop out the unused apirs
                    bedpe.drop(toexclude.index,axis=0,inplace=True)
            else:
                pass 

            ## Remove reads not mapping to chromosomes/contigs of interest 
            if len(includos):
                ## Gather read pairs that we are not include
                toexclude = bedpe[~(bedpe.Rname1.isin(includos) & bedpe.Rname2.isin(includos))].copy()
                ## Calc rows 
                exclude_count = toexclude.Rname1.count()
                ## if we have chroms to not exclude
                if exclude_count:
                    ## update the count
                    excluded_counts += exclude_count
                    ## Add the error 
                    toexclude['Error'] = 'excluded'
                    ## Appedn to the not used lsit
                    not_used.append(toexclude) 
                    ## Drop out the unused apirs
                    bedpe.drop(toexclude.index,axis=0,inplace=True)
            else:
                pass 

            ## Remove low quality mapping, set error to low qual
            lowqual = bedpe[(bedpe.Minmapq<minmapq)].copy()
            ## Count the low qual
            lowqual_count = lowqual.Minmapq.count()
            ## If we have low quality reads
            if lowqual_count:
                ## Update counts
                lowquals_counts += lowqual_count
                ## Set error 
                lowqual['Error'] = 'lowqual'
                ## Append to the not used list
                not_used.append(lowqual)
                ## Drop from bedpe chunk
                bedpe.drop(lowqual.index,axis=0,inplace=True)

            ## Remove inter chromosome reads
            if intraonly:
                ## Gather inter
                interc = bedpe[(bedpe.Orientation=='Inter')].copy()
                ## Count 
                inter_count = interc.Rname1.count()
                ## If we have counts
                if inter_count:
                    ## Add to the count 
                    interchr_counts += inter_count
                    ## Add error 
                    interc['Error'] = 'inter'
                    ## Append to not used list
                    not_used.append(interc) 
                    ## Drop from bedpe chunk
                    bedpe.drop(interc.index,axis=0,inplace=True)
            else:
                pass 

            ## Filter those contacts larger distacen
            if max_dist:
                ## Gather those fragments to far away from eachother 
                tolarge = bedpe[(bedpe.Distance>max_dist) & (bedpe.Orientation!='Inter')].copy()
                ## Count the fragments
                tolarge_count = tolarge.Rname1.count()
                ## If we ahve counts 
                if tolarge_count:
                    ## Count the to large frags
                    largefrg_counts += tolarge_count
                    ## ADd the error 
                    tolarge['Error'] = 'largefrag'
                    ## Append to no not used list
                    not_used.append(tolarge) 
                    ## Drop from bedpe chunk
                    bedpe.drop(tolarge.index,axis=0,inplace=True)
            else:
                pass 

            ## Concat the not used list so far and save out
            if len(not_used):
                ## Append to path
                not_usede_paths.append(not_usede_path)
                ## save the output 
                pd.concat(not_used,axis=0).to_csv(not_usede_path,sep=hicsep,header=True,index=False) 
            
            ## If restriction sites were passed 
            if restriciton_sites: 
                ## Gather the read pairs we plan to check for intra fragments
                tocheck = bedpe[(bedpe.Distance<=error_dist) & (bedpe.Orientation.isin(['Outward','Inward'])) & (bedpe.Inter==0)].copy()
                ## Count the checks
                tocheck_count = tocheck.Rname1.count()
                ## If we have counts
                if tocheck_count:
                    ## Append to list 
                    too_check_paths.append(too_check_path)
                    ## Save out the reads to check for intra fragments
                    tocheck.drop('Error',axis=1).to_csv(too_check_path,sep=hicsep,header=True,index=False) 
                    ## Drop these from bedpe
                    bedpe.drop(tocheck.index,axis=0,inplace=True)
            else:
                pass 

            ## count the valid
            valid_count = bedpe.Rname1.count()
            ## If we have valid count
            if valid_count:
                ## Add to the total valid counts
                validhic_counts += valid_count
                ## append to path 
                hic_valid_paths.append(hic_valid_path)
                ## Save out the not used reads and the to check names 
                bedpe.drop('Error',axis=1).to_csv(hic_valid_path,sep=hicsep,header=True,index=False) 
                ## Count the inter
                interhic_count = bedpe.Inter.sum()
                interhic_counts += interhic_count
                ## Calc the intra hic 
                intrahic_count = valid_count - interhic_count
                intrahic_counts += intrahic_count

    ## Delet the last chunks, we don't need these
    del bedpe, not_used
    ## Set the intra qnames
    intra_all_qnames = set()
    #print("Above if statment")
    ## If we are checking rest sites and we have dataframes to check 
    if restriciton_sites and len(too_check_paths):
        ## Gather the complement of the restirction sites
        all_sites = (restriciton_sites + [str(Seq(r).complement()) for r in restriciton_sites]) if restriciton_sites else []

        #print("Within if statment")
        ## Bring in the reads to check 
        allcheck = dd.read_csv(too_check_paths,sep=hicsep)
        ## Gather the chromosome list
        chrlist = list(allcheck.Rname1.unique().compute()) ## NOTE: This must be of type list for boolean statemetns below to work
        #print(chrlist)
        #print(type(chrlist))
        ## Load in the reference 
        refs_parse = SeqIO.parse(refpath,format='fasta')
        ## Iterate thru ref parser 
        #print('Parsing per chromosome')
        for ref in refs_parse:
            #print('Within chromosome loop: %s %s'%(ref.id,ref.name))
            if (ref.id in chrlist) | (ref.name in chrlist):
                ## Calc the chr seq
                chrseq = ref.seq
                ## set the dataframe to check 
                tocheck = allcheck[(allcheck.Rname1==ref.id) | (allcheck.Rname1==ref.name)][check_names].compute()
                
                if hicexplorer: ## Perform our hic exploer filtering                 
                    ## Initiate the fragment sites
                    tocheck['Left1'], tocheck['Right1'], tocheck['Left2'], tocheck['Right2'] = -1, -1, -1, -1
        
                    ## Iterate thru the tupels 
                    for (a,b,c) in value_tuple:
                        tocheck[a] = tocheck[b].apply(nearestrest, args=[restriciton_sites,chrseq,c])
                
                    ## Gather the intra fragment read pairs, those with fragmetns / rest sites bounds that are equal or overlapping 5' to 3'
                    intra_frags = tocheck[(tocheck.Left1==tocheck.Left2) | (tocheck.Right1 == tocheck.Right2) | (tocheck.Right1>=tocheck.Left2)]
                    intra_count = intra_frags.Left1.count()

                    ## Gather qnames if intrafrags has hsape
                    if intra_count:
                        ## Gather qnames and update the list of qnames
                        intra_all_qnames.update(intra_frags.Qname1.tolist())
                else:
                    ## Format restrction site list and index 
                    n_restsites = []
                    rest_index = []

                    for i,row in tocheck.iterrows():
                        ## Gather the positions
                        a = row[['Pos1','Pos2','End1','End2']].min()
                        b = row[['Pos1','Pos2','End1','End2']].max()
                        n = sum([str(chrseq[a+len(r):b-len(r)]).upper().count(r) for r in all_sites])
                        
                        ## Append the count
                        n_restsites.append(n)
                        rest_index.append(i)

                    ## Format Rcount into df
                    rest_count = pd.DataFrame(n_restsites,index=rest_index,columns=['Rcount'])
                    ## SEt Rcount in to check 
                    tocheck['Rcount'] = rest_count['Rcount']

                    ## Gather the intra fragment read pairs, those with fragmetns / rest sites bounds that are equal or overlapping 5' to 3'
                    intra_frags = tocheck[(tocheck.Distance<0) | (tocheck.Rcount<1)]
                    intra_count = intra_frags.Pos1.count()

                    ## Gather qnames if intrafrags has hsape
                    if intra_count:
                        ## udpate intra qnames 
                        intra_all_qnames.update(intra_frags.Qname1.tolist())                    
            else: ## otherwise we didn't find the chromosmoe 
                pass 
            
        ## IF we have qnames that are errors
        if len(intra_all_qnames):
            #print('Yes we do!')
            ## Set the qnames 
            errors = allcheck[(allcheck.Qname1.isin(intra_all_qnames))]
            ## count errors
            error_count = errors.Inter.count().compute()
            ## add to the error count
            if error_count:
                errorfrg_counts += error_count
                errors['Error'] = 'intrafrag'
                ## Set path to save out the errors 
                out_error_path = makeoutpath(bedpe_path,f'errors.{0}')
                ## append to path
                not_usede_paths.append(out_error_path)
                ## SAve out the errors
                errors.to_csv(out_error_path,sep=hicsep,single_file=True,index=False)

        ## Gather the valid
        valid = allcheck[(~allcheck.Qname1.isin(intra_all_qnames))]
        ## Count the valid
        valid_count = valid.Inter.count().compute()
        ## If we have valid counts
        if valid_count:
            ## If we have valid counts 
            validhic_counts += valid_count
            ## SEt the output path to save the valid contacts
            out_valid_path = makeoutpath(bedpe_path,f'tohic.{0}')
            ## Append to list of paths of valid
            hic_valid_paths.append(out_valid_path)
            ## Save out the valid contacts
            valid.to_csv(out_valid_path,sep=hicsep,single_file=True,index=False)
            ## Count the inter chromosomes
            interhic_count = valid.Inter.sum()
            ## Aadd to the inter count 
            interhic_counts += interhic_count
            ## Calculate the intra
            intrahic_count = valid_count - interhic_count
            ## Add the intra
            intrahic_counts += intrahic_count
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
    
    ## Write out the counts
    new_names = [   'Unmapped',    'Dovetails',    'Excluded',    'Lowquality',   'Inter',     'Largefragment',  'Samefragment',    'Valid',       'Interhic',      'Intrahic']
    new_count = [unmapped_counts,dovetail_counts,excluded_counts,lowquals_counts,interchr_counts,largefrg_counts,errorfrg_counts,validhic_counts,interhic_counts,intrahic_counts]
    ## Iterate thru and print the counts to log 
    [print('INFO: %s\t%s'%(a,b)) for a,b in zip(new_names,new_count)]
    ## print to log
    print("Finished filtering and splitting bedpe file: %s"%bedpe_path)
## End of file 