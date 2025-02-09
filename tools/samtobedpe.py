#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
filterdesc = "Filters an input tsv file for properly mapped Hi-C contacts, filtering on mapping quality and intra-chromosomal distance."
## ----------------------------- v 1.0.0 ------------------------------ ## 
##   Sam to bedpe for Hi-C data processing 
## -------------------------------------------------------------------- ##
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545 
croth@lanl.gov
"""
## ----------------------------------- MODULE LOADING ------------------------------------ ##
## Bring in pandas and np 
import pandas as pd, numpy as np, re
## Load in params
from ..defaults.parameters import hicsep, chunksize, Z_help, r_help, L_help
## Bring in default ftns 
from ..defaults.defaults import pathexists, reset, basenosam
## Load in samnames 
from ..tools.pysamtools import getmatchsum, listzip, dictzip, loadref, loadsam
## Load in library restriction site 
from ..diethic.filterbedpe import returnsite, endcheck

## ------------------------ FUNCTION and VARIABLE FILTERING ------------------------------ ## 
## Make list of unique bits
ubits = [2**i for i in range(12)]

## List descripiotn of flags
descr = ['Multiseq',    # template having multiple segments in sequencing
         'Propalign',   # each segment properly aligned according to the aligner
         'Segmunmp',    # segment unmapped
         'Nextsegunmp', # next segment in the template unmapped
         'Seqrev',      # SEQ being reverse complemented
         'Nextseqrev',  # SEQ of the next segment in the template being reverse complemented
         'Read1',       # the first segment in the template
         'Read2',       # the last segment in the template
         'Secalign',    # secondary alignment
         'Nopass',      # not passing filters, such as platform/vendor quality controls
         'Pcrdup',      # PCR or optical duplicate
         'Supalign']    # supplementary alignment
         
## Zip and make into a dict 
bitflags = dict(zip(ubits,descr))

## Check our work
assert len(bitflags) == 12

## Write the cigar strings 
cigarlist = ['M','I','D','N','S','H','P','=','X']

## Format cigar dict
cigardict = dict(zip(cigarlist,range(len(cigarlist))))

## Ftn for spliting cigar
def splitcigar(s,strand):
    c = re.split('(\\d+)',s)[1:]
    d = np.array([*map(cigardict.get,c[1::2])],dtype=int)
    n = np.array(c[::2],dtype=int)
    l = listzip(d,n)
    return l if strand else l[::-1]

## Ftn for getting set of of flags from df
def getflags(df):
    ## Return list of flags 
    return list(set(df['Flag']))

## Ftn for adding alignment flags to sam df
def makeflagdf(df):
    ## Gather list of unique flags 
    flags = getflags(df)
    ## Calcualte a dataframe with bit inforamtion for each unqiue flag
    mbits = pd.DataFrame([[f & l for l in ubits] for f in flags],columns=descr)
    ## append the flags
    mbits['Flag'] = list(flags)
    ## Merge with origional df and return 
    return df.merge(mbits[['Flag','Seqrev','Read1','Read2']]).sort_values('Qname').reset_index(drop=True)

## Ftn for appending to df
def appenddf(df,outpath,thesep=hicsep):
    ## Append the df
    df.to_csv(outpath,mode='a' if pathexists(outpath) else 'w',index=False,header= not pathexists(outpath),sep=thesep)
    pass 

def checkread(df,c):
    assert np.sum(df['Read%s'%c]>0) == df.shape[0], "ERROR: Not all reads are the %s read in pair!"%c
    pass 

def renamecols(df,k):
    df.drop(['Rnext','Pnext','Tlen','Read1','Read2','Matchsum'],axis=1,inplace=True)
    df.columns = [c+str(k) for c in df.columns]
    return df    

## Ftn for makign dictionary of chromosome names and id numbers 
def makechromdict(inrefpath):
    if type(inrefpath) == str:
        ## If we the refeerenc has an fai index or is a refi index 
        if pathexists(inrefpath+'.fai') or (inrefpath.split('.') == 'fai'):
            ## Load in the chromosome dataframe
            chrdf = pd.read_csv(inrefpath+'.fai',sep='\t',header=None)[[0]]
            ## Format chromosome names into a list 
            chrdict = dictzip(chrdf.T.values[0],chrdf.index)
        else: ## otherwise, load in the ref and make a chromosome dictionary 
            ref = loadref(inrefpath)
            chrdict = dict([(r.name,i) for i,r in enumerate(ref)])
    ## Othereise it is a reference object, assuming
    else:
        assert type(inrefpath) == list, "ERROR: The expected input was a list for loadref ftn!"
        chrdict = dict([(r.name,i) for i,r in enumerate(inrefpath)])
    ## Assign unmapped reads 
    chrdict['*'] = -1
    return chrdict 

## Ftn for returning orientation of records 
def getoriented(s1,s2):
    """
    Returns the paired read orientation types, shown below:

    read 1                                  read 2
    -------------->      ...      <---------------
    inward


    <--------------      ...      --------------->
    outward


    <--------------      ...      <---------------
    left 


    -------------->      ...      --------------->
    right

    -----------------------------------------------
    |||||||||||||||||||||||||||||||||||||||||||||||
    --------------------genome---------------------

    """
    ## Gather the strand codes 
    if (s1 == 0) and (s2 ==0):
        orientation = 'Right'
    elif (s1 == 0) and (s2 >0):
        orientation = 'Inward'
    elif (s1 > 0) and (s2 == 0):
        orientation = 'Outward'
    elif (s1 > 0) and (s2 > 0):
        orientation = 'Left'
    else: ## Print an error 
        print("WARNING: We have encountered an orientation we did not recognize")
        orientation = 'Unknown'
    return orientation

def bycigar(seq,cigar):
    ## Check if cigar is unmapped
    if cigar == '*':
        retseq = seq
    else:
        ## split the cigar, with a positive strand
        cigs = splitcigar(cigar,1)
        ## Calc the front and back indexes 
        front_ = cigs[0][1] if cigs[0][0] >= 4 else 0
        backs_ = -cigs[-1][1] if cigs[-1][0] >= 4 else len(seq) + 1
        ## Set the return seq
        retseq = seq[front_:backs_]
    ## Return the sequences
    return retseq

## Ftn for formating long dataframe from juicer 
def formatlong(df:pd.DataFrame,r1_ix:list,r2_ix:list) -> pd.DataFrame:
    """
    Returns a list of mappign values from an input set of rows in sam format representing a Hi-C alignments. 
    The numeric value of the chormosome is also added.
    """
    ## Split into reads one and two 
    df1,df2 = df.loc[r1_ix].reset_index(drop=True), df.loc[r2_ix].reset_index(drop=True)
    ## Check reads
    checkread(df1,1),checkread(df2,2)
    ## Rename columns and concat the dfs
    df1,df2 = renamecols(df1,1), renamecols(df2,2)
    ## Concat the dfs 
    long = pd.concat([df1,df2],axis=1)
    ## Calculate distance
    long['Distance'] = long.Pos2 - long.Pos1
    ## Remap the distances
    long.loc[long.Distance<0] = long.loc[long.Distance<0,df2.columns.tolist()+df1.columns.tolist()+['Distance']].values
    ## Recalc distance metric
    long['Distance'] = long.Pos2 - long.Pos1
    ## Clac where the chromosomes are not left right sorted
    long['Test'] = long.Chrn1 - long.Chrn2
    ## Remap so chromosomes are left right sorted
    long.loc[long.Test>0] = long.loc[long.Test>0,df2.columns.tolist()+df1.columns.tolist()+['Test','Distance']].values
    ## Drop test
    long.drop('Test',axis=1,inplace=True)    
    ## Calculate min mapping quality
    long['Minmapq'] = long[['Mapq1','Mapq2']].min(axis=1)
    ## Add oriantation 
    long['Orientation'] = long.apply(lambda x: getoriented(x["Seqrev1"], x["Seqrev2"]), axis=1) if (long.shape[0] > 0) else None
    ## Return the longdf
    return long

## Ftn for post filtering recroeds 
def postfilter(inmapping:pd.DataFrame,outdfpath:str,chrdict:dict,danglingends,filename=None):  #,enzymelib:str,error_dist:int,frag_size:int,refseqio:str):
    ## add mapped flags
    mmapped = makeflagdf(inmapping)
    ## Check our work
    assert mmapped.shape[0] == inmapping.shape[0], "ERROR: There was a critical error in handeling read pairs."
    
    ## Gather the unique cigars from positivly orinted reads, calcualte and add the matcsum 
    poscigars = mmapped[(mmapped.Seqrev==0)][['Cigar']].drop_duplicates().copy().reset_index(drop=True)
    poscigars['Matchsum'] = poscigars.Cigar.apply(splitcigar,args=[1]).apply(getmatchsum)

    ## Gather the unique cigars form negativel orinted reads, calcualte and add the matcsum 
    negcigars = mmapped[(mmapped.Seqrev>0)][['Cigar']].drop_duplicates().copy().reset_index(drop=True)
    negcigars['Matchsum'] = negcigars.Cigar.apply(splitcigar,args=[0]).apply(getmatchsum)

    ## Calculate the number of positively and negativel oriented reads 
    pos_shape = mmapped[(mmapped.Seqrev==0)].shape[0]
    neg_shape = mmapped[(mmapped.Seqrev>0) ].shape[0]
    ## Merge the cigar values 
    pos = mmapped[(mmapped.Seqrev==0)].merge(poscigars)
    neg = mmapped[(mmapped.Seqrev>0) ].merge(negcigars)

    ## Check our shapes 
    assert (pos.shape[0] == pos_shape) and (neg.shape[0] == neg_shape)

    ## make a new map with the cigar values 
    newmap = pd.concat([pos,neg]).sort_values('Qname').reset_index(drop=True)
    ## Make the matsum column an int
    newmap['Matchsum'] = newmap['Matchsum'].apply(int)
    ## Gather and trim the seq
    newmap['Seq'] = newmap.apply(lambda row: bycigar(row['Seq'],row['Cigar']),axis=1)
    ## Add in seq len
    newmap['Len'] = newmap.Seq.apply(len)
    ## Add in end positoin
    newmap['End'] = newmap.Pos + newmap.Len 
    ## Check the danginling ends
    if danglingends:
        newmap['Dangend'] = newmap.apply(lambda row: endcheck(row['Seq'], row['Seqrev'], danglingends),axis=1) ## Removed, only check error fragments 
    else:
        newmap['Dangend'] = 0

    ## Check our shape and readnames
    assert newmap.shape[0] == mmapped.shape[0], "ERROR: Shape missmatch!"
    assert not len(set(newmap.Qname) - set(inmapping.Qname)), "ERROR: Missing reads!"
    
    ## Gather the minimized index for the match sum column and set to list 
    r1_ix = newmap[(newmap.Read1>0)].groupby('Qname').Matchsum.idxmin().tolist()
    r2_ix = newmap[(newmap.Read2>0)].groupby('Qname').Matchsum.idxmin().tolist()

    ## Add chromosomes integers 
    newmap['Chrn'] = newmap.Rname.map(chrdict)
    ## Make a long df
    long = formatlong(newmap,r1_ix,r2_ix)
    ## Check work 
    assert np.sum(long.Qname1 == long.Qname2), "ERROR: Read pair missmatch in chunk!"
    ## drop the qnames
    long.drop('Qname2',axis=1,inplace=True)

    ## calc intra vs intra contact, make inter -1 for unmapped reads
    long['Inter'] = 0                                               ## intra
    long.loc[(long.Chrn1!=long.Chrn2),'Inter'] = 1                  ## inter
    long.loc[(long.Chrn1<0) | (long.Chrn2<0),'Inter']  = -1         ## unmapped 

    ## Set distance of inter-chrom contacts to a dummy var
    long.loc[(long.Inter>0),'Distance'] = 10**8

    ## Set filename 
    if filename:
        long['File'] = filename

    ## Save out the pairs and pass back to main
    appenddf(long,outdfpath) if long.shape[0] else None
    pass 

## Set help strings
S_help = "Path to input .sam file from bwa. This can be a TSV file but must have a .sam extension."

## ---------------------------- MAIN SCRIPT & ARGUMENT PARSING -------------------------- ## 
## If the script is envoked
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse

    ## Make the parse
    parser = argparse.ArgumentParser(description = filterdesc)
    ## Add the required arguments
    parser.add_argument("-i", dest="I", type=str,  required=True,  help=S_help, metavar='./path/to/input.sam'    )
    parser.add_argument("-r", dest="R", type=str,  required=True,  help=r_help, metavar='./path/to/ref.fasta'    )
    parser.add_argument("-z", dest="Z", type=int,  required=False, help=Z_help, metavar='n',  default=chunksize  )
    parser.add_argument("-l", dest="L", type=str,  required=False, help=L_help, metavar='Arima', default='none'  )
    ## Set the paresed values as inputs
    inputs = parser.parse_args()

    ## Set inputs
    inpath       = inputs.I  ## Path to input sam file 
    refpath      = inputs.R  ## input path to reference 
    chunk_size   = inputs.Z  ## Chunk-size
    library      = inputs.L  ## Set the library size

    ## Gather file name 
    samname = basenosam(inpath)

    ## Load in file
    filehead = inpath.split('.sam')[0]
    outfile  = filehead + '.bedpe'
    ## Reset, removing previous runs
    reset([outfile])

    ## Format the chromosome dictionary from ref object
    chromdict = makechromdict(refpath)

    ## set dandling ends
    restriction_sites, dangling_ends = returnsite(library)

    ## Open the inpath as chunks 
    with loadsam(inpath,chunksize=chunk_size) as samchunks:
        lastc = []
        ## Iteraate thru the chunks 
        for inmap in samchunks:
            if inmap.shape[0] == 0:
                print("WARNING: This sam file ( %s ) was empty."%inpath)
                continue
            ## If this is the start of our second time or more thru this loop
            inmap = pd.concat([lastc,inmap],axis=0).reset_index(drop=True) if len(lastc) else inmap
            ## Gather the last Q name 
            lastq = inmap.loc[inmap.index.max(),'Qname']
            ## parse off the last Q
            lastc = inmap[(inmap.Qname==lastq)].copy()
            ## Cut off the end if this chunk is the smallest and their for the last
            tomap = inmap[(inmap.Qname!=lastq)]
            ## Call our filtering ftn
            postfilter(tomap,outfile,chromdict,dangling_ends,filename=samname) #,enzyme_lib,errorsize,min_frag_s,refpath)
        ## Make sure we get the last bit if it wasn't included, run the last chunk. I don't like this coding and will retrun to it later. 
        postfilter(lastc,outfile,chromdict,dangling_ends,filename=samname) if len(lastc) else None 

    ## Print to log
    print(f"Finished converting sam file {inpath} to bedpe format.")
## End of file 
