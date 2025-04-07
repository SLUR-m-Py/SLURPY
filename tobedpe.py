#!/usr/bin/env python
## Import mods 
import pandas as pd, numpy as np, re, sys 
## Loading defaults
from defaults import pathexists
## Load in variables and ftn from my other libs
from pysamtools import samnames, samtypes, dictzip, loadref, getmatchsum, listzip
## Bring in return site ftn form to bedp 
from filtering import returnsite
## Load in params
from parameters import hicsep

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

## Write the cigar strings 
cigarlist = ['M','I','D','N','S','H','P','=','X']

## Format cigar dict
cigardict = dict(zip(cigarlist,range(len(cigarlist))))

## Set the new types 
newtypes = samtypes + [int]

## Ftn for getting tab seperated fields from bwa line 
def getfields(inline,count=0) -> list:
    ## Strip the line 
    sept = inline.strip().split('\t')
    ## Return the seperated fields if a count was given, else all the fields 
    return [sept[i] for i in range(count)] if count else sept

## Ftn for returning bool if the alignment is chimeric
def ischimeric(line) -> bool: 
    ## REturn the boolean check 
    return 'SA' == getfields(line)[-1].split(':')[0]

## Ftn for making mapped list into a dataframe 
def makedf(mapped:list,totypes:list):
    tmp = pd.DataFrame(mapped,columns=samnames + ['Chimeric'])
    for i,c in enumerate(tmp.columns):
        tmp[c] = tmp[c].apply(totypes[i])
    return tmp 

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
        new_seq = []
        for i,c in cigs:
            ## If the the cigar is a match, insertion, deleation, or gap
            if i < 4:
                new_seq.append(seq[:c])
            ## trim the sequence 5' to 3' to match next cigar start 
            seq = seq[c:]
        ## Format the new seq
        retseq = ''.join(new_seq)
        ## Calc the front and back indexes 
        #front_ = cigs[0][1] if cigs[0][0] >= 4 else 0
        #backs_ = -cigs[-1][1] if cigs[-1][0] >= 4 else len(seq) + 1
        ## Set the return seq
        #retseq = seq[front_:backs_]
    ## Return the sequences
    return retseq

def lencigar(cigar):
    ## Check if cigar is unmapped
    if cigar == '*':
        retseq = 0
    else:
        ## split the cigar, with a positive strand
        cigs = splitcigar(cigar,1)
        new_len = []
        for i,c in cigs:
            ## If the the cigar is a match, insertion, deleation, or gap
            if i < 4:
                new_len.append(c)
        ## Calc new length 
        retseq = sum(new_len)
    ## Return the sequences
    return retseq

## Checks the read are in pair
def checkread(df,c):
    assert np.sum(df['Read%s'%c]>0) == df.shape[0], "ERROR: Not all reads are the %s read in pair!"%c
    pass 

## Ftn for renameing cols with int for read in pair 
def renamecols(df,k):
    df.drop(['Rnext','Pnext','Tlen','Read1','Read2','Matchsum'],axis=1,inplace=True)
    df.columns = [c+str(k) for c in df.columns]
    return df   

## Set orientation variable for easy use
orient = 'Orientation'

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
    ## Calculate distance beteween left of pairs, its negative if read 1 is to the right of read 2
    long['Distance'] = long.Pos2 - long.Pos1
    ## Remap the distances
    long.loc[long.Distance<0] = long.loc[long.Distance<0,df2.columns.tolist()+df1.columns.tolist()+['Distance']].values
    ## Recalc distance metric, take mins and maxs over the fragments 
    #long['Distance'] = long.Pos2 - long.End1
    dist_min = long[['Pos1','Pos2','End1','End2']].min(axis=1)
    dist_max = long[['Pos1','Pos2','End1','End2']].max(axis=1)
    long['Distance'] = dist_max - dist_min

    ## Clac where the chromosomes are not left right sorted
    long['Test'] = long.Chrn1 - long.Chrn2
    ## Remap so chromosomes are left right sorted
    long.loc[long.Test>0] = long.loc[long.Test>0,df2.columns.tolist()+df1.columns.tolist()+['Test','Distance']].values
    ## Drop test
    long.drop('Test',axis=1,inplace=True)    
    ## Calculate min mapping quality
    long['Minmapq'] = long[['Mapq1','Mapq2']].min(axis=1)
    ## Add oriantation 
    long[orient] = long.apply(lambda x: getoriented(x["Seqrev1"], x["Seqrev2"]), axis=1) if (long.shape[0] > 0) else None
    ## Return the longdf
    return long

## Check ends ftn 
def endcheck(read:str,strand:int,dangends:list) -> int:
    ## Checks the presesne of deanglign ends 
    return sum([read.upper().endswith(d) for d in dangends] if strand else [read.upper().startswith(d) for d in dangends])

## Ftn for appending to df
def appenddf(df:pd.DataFrame,outpath:str,toappend:bool,thesep=hicsep) -> None:
    ## Append the df
    df.to_csv(outpath,mode='a' if toappend else 'w',index=False,header= not toappend,sep=thesep)
    pass 

## Ftn for post filtering recroeds 
def postfilter(inmapping:pd.DataFrame,outdfpath:str,chrdict:dict,toappend:bool):  #,enzymelib:str,error_dist:int,frag_size:int,refseqio:str):
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
    #newmap['Seq'] = newmap.apply(lambda row: bycigar(row['Seq'],row['Cigar']),axis=1)
    ## Add in seq len
    newmap['Len'] = newmap['Cigar'].apply(lencigar)
    ## Add in end positoin
    newmap['End'] = newmap.Pos + newmap.Len 
    ## Check the danginling ends
    #if danglingends:
    #    newmap['Dangend'] = newmap.apply(lambda row: endcheck(row['Seq'], row['Seqrev'], danglingends),axis=1) ## Removed, only check error fragments 
    #else:
    #    newmap['Dangend'] = 0

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
    long['Inter'] = 0                                              
    ## Set inter bool
    inter_bool    = (long.Chrn1!=long.Chrn2)
    unmapped_bool = (long.Chrn1<0) | (long.Chrn2<0)

    ## Set inter and reset orientation 
    long.loc[inter_bool,'Inter']    = 1               
    long.loc[inter_bool,orient]     = 'Inter'
    long.loc[inter_bool,'Distance'] = 0

    ## Set unmpaped orientation 
    long.loc[unmapped_bool,'Inter'] = -1         
    long.loc[unmapped_bool,orient]  = 'Unmapped'

    ## Set distance of inter-chrom contacts to a dummy var
    ## Set filename 
    #if filename:
    #    long['File'] = filename

    ## Save out the pairs and pass back to main
    appenddf(long,outdfpath,toappend) if long.shape[0] else None
    pass 

## ------------------------ MAIN SCRIPT ------------------------------ ## 
## If the script is envoked by name 
if __name__ == "__main__":
    ## GAther input variabels 
    ref_path   = sys.argv[1]
    library    = sys.argv[2]
    outfile    = sys.argv[3]
    line_count = int(sys.argv[4])

    ## Set the sam name 
    samname    = outfile.split('/')[-1].split('.bedpe')[0]

    ## Initiate lists 
    filt_lines = []
    lastc      = []
    i          = 0
    c          = 0 

    ## Format the chromosome dictionary from ref object
    chromdict = makechromdict(ref_path)

    ## set dandling ends
    restriction_sites, dangling_ends = returnsite(library)

    ## Go thru stand input 
    for line in sys.stdin:
        if line.startswith('@'):
            continue
        i += 1
        ## Determine if read is cimeric 
        chimeric_bool = ischimeric(line)

        ## Append the first 7 or so fiels from bwa line 
        filt_lines.append(getfields(line,len(samnames)) + [int(chimeric_bool)])

        ## Trigger the processing into bedpe 
        if len(filt_lines) > line_count:
            ## Generate input map of filtered liens 
            inmap = makedf(filt_lines,newtypes)
            ## If this is the start of our second time or more thru this loop
            inmap = pd.concat([lastc,inmap],axis=0).reset_index(drop=True) if len(lastc) else inmap
            ## Gather the last Q name 
            lastq = inmap.loc[inmap.index.max(),'Qname']
            ## parse off the last Q
            lastc = inmap[(inmap.Qname==lastq)].copy()
            ## Cut off the end if this chunk is the smallest and their for the last
            tomap = inmap[(inmap.Qname!=lastq)]
            ## Call our filtering ftn
            postfilter(tomap,outfile,chromdict,toappend=bool(c)) #dangling_ends,,enzyme_lib,errorsize,min_frag_s,refpath)
            ## Rest filt lines 
            filt_lines = []
            c += 1

    ## If any are left from filt lines, append them and last chunk (if it has length too)
    if len(filt_lines):
        inmap = makedf(filt_lines,newtypes)
        ## If this is the start of our second time or more thru this loop
        tomap = pd.concat([lastc,inmap],axis=0).reset_index(drop=True) if len(lastc) else inmap
        ## Call our filtering ftn
        postfilter(tomap,outfile,chromdict,toappend=bool(c)) #,enzyme_lib,errorsize,min_frag_s,refpath)
        lastc = []

    ## If the last chunk still exsits and has alignments 
    if len(lastc):
        postfilter(tomap,outfile,chromdict,toappend=bool(c)) #,enzyme_lib,errorsize,min_frag_s,refpath)
    ## PRint to log 
    print("Finished parsing %s inputs from %s to %s"%(i,samname,outfile))
## EOF