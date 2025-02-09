#!/usr/bin/env python
## Set slurpy path 
#slurpy_path = '/panfs/biopan04/4DGENOMESEQ/EPICPIPELINE/SLURPY'
## filtering bwa standard output via python test
import sys
#sys.path.append(slurpy_path)
## Import subprocess 
import pandas as pd 
## Load in variables and ftn from my other libs
from ..tools.pysamtools import samnames, samtypes
from ..tools.samtobedpe import postfilter, makechromdict
from ..diethic.filterbedpe import returnsite

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
            postfilter(tomap,outfile,chromdict,dangling_ends,filename=samname) #,enzyme_lib,errorsize,min_frag_s,refpath)
            ## Rest filt lines 
            filt_lines = []

    ## If any are left from filt lines, append them and last chunk (if it has length too)
    if len(filt_lines):
        inmap = makedf(filt_lines,newtypes)
        ## If this is the start of our second time or more thru this loop
        tomap = pd.concat([lastc,inmap],axis=0).reset_index(drop=True) if len(lastc) else inmap
        ## Call our filtering ftn
        postfilter(tomap,outfile,chromdict,dangling_ends,filename=samname) #,enzyme_lib,errorsize,min_frag_s,refpath)
        lastc = []

    ## If the last chunk still exsits and has alignments 
    if len(lastc):
        postfilter(lastc,outfile,chromdict,dangling_ends,filename=samname) 
    ## PRint to log 
    print("Finished parsing %s inputs from %s to %s"%(i,samname,outfile))
## EOF