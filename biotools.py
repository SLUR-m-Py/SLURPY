#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## ----------------------------- v 8.0.0 ------------------------------ ## 
##   BIOTOOLS: Python functions for Working with biological data
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
## Load in defaults
from defaults import basename, sortglob, listzip, fileexists, basenoext, makelist, submitsbatch, ifprint, dictzip
## Load in parameters
from parameters import slurpydir, debugdir, runlocal, fakejobid, macs3dir, g_help
## Load in rm tree
from shutil import rmtree
## Load in sys
import sys, pysam, pandas as pd, argparse
## Brin in datetime
from datetime import datetime
## Load in os
from os import makedirs, remove
## Load in get size 
from os.path import getsize, exists
## Load in SeqIO
from Bio import SeqIO

## Write the sam error 
sam_error = "ERROR: The detected version of samtools is not greater than or equal to v 1.15.1!\nPlease update samtools and try again."
## Column names and types of sam file 
samnames = ['Qname','Flag','Rname','Pos','Mapq','Cigar','Rnext','Pnext','Tlen','Seq']
samtypes = [  str,   int,   str,    int,  int,    str,    str,    int,   int,   str]

## Write ftn for getting samtoosl version 
def getsamv(com='samtools --version') -> str:
    """Gathers the version of samtools as a string."""
    ## Returns the output of samools version
    return submitsbatch(com,returnid=False)[0].split(' ')[-1]

## Write ftn for spliting and making an int
def splitint(x,sep='.') -> list[int]:
    """Splits an input string x on SEP and converts values to integers."""
    ## Return the array of split ints
    return [int(z) for z in x.split(sep)]

## Check that we have the correct version of samtools
def checksam(version='1.15.1') -> bool:
    """Returns a boolean, checking if the version of samtools is correct."""
    ## Return the bool of the minimum difference
    return min([ a-b for a,b in zip(splitint(getsamv()),splitint(version))])>=0

## Ftn for loading in reference in fasta file format 
def loadref(inpath:str,format='fasta') -> list:
    """Returns a list of sequences from within input fasta file."""
    ## Return the records in the loaded reference 
    return [r for r in SeqIO.parse(inpath,format=format)]

## Reurn the match sumation of py records
def getmatchsum(cigartuples) -> int:
    """Returns the summed match score across cigar tuples from alignment record."""
    ## Set match sum to start at zero
    match_sum = 0
    ## Iterate over the cigar tuple
    for m, c in cigartuples[:]:
        ## If m is a match to the cmatch from pysam, break from loop
        if (m == pysam.CMATCH):
            break
        ## Add to match sum 
        match_sum += c
    ## Return the match sum
    return int(match_sum)

## Ftn for reseting clean bool
def patchclean(torerun:str,inclean:bool) -> bool:
    """Patches the clean boolean if the rerun step is clean or if the input clean boolean is true."""
    ## Patches the clean step
    return True if (torerun == 'clean') else inclean

## Ftn for setting run name
def setrunname(runname:str, incwd:str) -> str:
    """Sets the name of the pipeline run (if not given) from the current working directory."""
    ## Return the run name if it is defined, else use the basename of the cwd
    return runname if runname else basename(incwd)

## Ftn for getting fastqs
def getfastqs(fpath:str) -> list:
    """Lists fastq files in pair from input, wildcard path."""
    ## Gather all read fastq.gz files
    reads = sortglob(fpath)
    ## parse r1 and r2 
    r1,r2 = reads[::2], reads[1::2]
    ## Check r1 and r2, for proper naming convention
    for r in r1:
        assert '_R1_' in r, "ERROR: the first in pair needs to have the name _R1_ in the read name."
    for r in r2:
        assert '_R2_' in r, "ERROR: the second in pair needs to have the name _R2_ in the read name."
    ## Gather the pairs
    return listzip(r1,r2)

## Ftn for checking if a file is bwa indexed
def isbwaix(inref:str, indexends = ['amb', 'ann', 'bwt', 'pac', 'sa']) -> int:
    """Checks if the the input reference file is indexed."""
    ## Check to make sure the reference exists
    assert fileexists(inref), "ERROR: Unable to locate input reference file on path: %s"%inref
    ## Check if each of these files exist: amb, ann, bwt, pac and sa
    return sum([fileexists(inref+'.'+fe) for fe in indexends]) == len(indexends)

## Ftn for the begining of a call to fastp
def fastcut(i:str, I:str, o:str, O:str, r:str, n:int) -> tuple:
    """Formats the start of a call to fastp given inputs (i,I), outputs (o,O), and the report name (r)."""
    ## Format the fall and return
    return f'fastp -i {i} -I {I} -o {o} -O {O} -j {r}.{n}.json -h {r}.{n}.html', f'{r}.{n}.json', f'{r}.{n}.html'

## Write ftn for logging fastp has finished
def fastdry(r1:str, r2:str, report:str) -> str:
    """Formats an echo command for logging the completion of a fastp split."""
    ## Format and return the echo command
    return f'{slurpydir}/myecho.py Finished filtering and splitting on: {r1} {r2} {report}\n'

## Ftn for getting sample name
def getsamplename(read1:str, spliton='_R1_') -> str:
    """Get sample name from a fastq.gz file."""
    ## Return the name 
    return basename(read1).split(spliton)[0]

## Ftn for formating report
def reportname(infile:str, script:str, i=0) -> str:
    """Formats the report name given the input file name and type of script (str)."""
    ## Return the report name
    return f'{debugdir}/{i}.{script}.{basenoext(infile)}.log'

## Ftn for checking hard reset
def confirmreset(indirs:list) -> None:
    """Prompts input from the user to confirm the hard reset of pipeline."""
    ## While we are waiting for input 
    while True:
        ## Send a query to the user 
        query = input("INFO: A hard restart was detected!\n\tDo you wish to delete ALL previous files, directories, and runs?\n") 
        ## Gather the answere 
        anser = query[0].lower() 
        ## Make sure the queary and answer are defined 
        if (query == '') or not (anser in ['y','n']):
            ## Print info to user otherwise 
            print("INFO: Please answer with yes (y) or no (n)!") 
        else: ## If we have a defined answere break the loop 
            break 
    ## If we are cool with hard restart 
    if (anser == 'y'): 
        ## Print what we are doing
        print("INFO: Hard restart confirmed; Removing previous run.")
        ## Remove the dir and its contents 
        [rmtree(dir_path,ignore_errors=True) for dir_path in indirs]
    else: ## Otherwise, exit with the system
        sys.exit("WARNING: Aborting hard reset and run of the pipeline!\n\tExiting")
    pass 

## Write ftn for making a directory
def dirmaker(dirpath:str):
    """Makes a directory given an input path."""
    ## Return the os command
    return makedirs(dirpath,exist_ok=True)

## Write ftn for making folders used here
def makedirectories(folderpaths) -> list:
    """Generates the directories listed in folderpaths."""
    ## Iteratively make the dir in the folder paths
    return [dirmaker(f) for f in makelist(folderpaths)]

## Ftn for writing command output logs
def writeparams(script:str,runname:str,sstamp,inputs):
    ## Format the start time stamp 
    dt1 = datetime.fromtimestamp(sstamp)
    ## Set the output path
    outpath = f'{debugdir}/run.parameters.{script}.{runname}.{sstamp}.txt'
    ## Open and writeout 
    with open(outpath,'w') as fout:
        fout.write(f'{script} run on {dt1}\n')
        fout.write('\n'.join([str(i[0]) + "=" + str(i[1]) for i in (vars(inputs)).items()]))
    ## Close output path
    fout.close()
    pass 

## Ftn for returing chrom
def chromdf(inpath:str) -> list:
    """Returns a pandas dataframe made from the tuples of sequences ids and lenghts from within an input fasta file."""
    ## Return the records in the loaded reference 
    return pd.DataFrame([(r.id,len(r.seq)) for r in SeqIO.parse(inpath,format='fasta')])

## Ftn for sorting fastq files by size
def sortfastq(infastq:list,splitsizes:list):
    ## Format input fastq zipped list into df
    tmp = pd.DataFrame(infastq,columns=['Read1','Read2'])
    ## Gather sizes of fast
    tmp['Size1'] = [getsize(f1) for (f1,f2) in infastq]
    ## Sort values by size and reset index 
    tmp = tmp.sort_values('Size1').reset_index(drop=True)
    ## Modify split sizzes 
    splitsizes = [splitsizes[0] for s in range(tmp.shape[0])] if (len(splitsizes) < tmp.shape[0]) else splitsizes
    ## Add a sorted list to splitsize to match size of fastqs 
    tmp['Splitsize'] = sorted([int(s)*4 for s in splitsizes])
    ## Return tmp
    return tmp 

## Ftn for formating job ids
def formatids(cdf:pd.DataFrame, op:list, joinon=',') -> str:
    """Formats job ids within the command dataframe (cdf) given an operation."""
    ## Return the joined list of ides 
    return joinon.join([str(k) for k in cdf[(cdf.Operation.isin(op))].JobID.tolist() if len(str(k)) > 0])

## Ftn for fromating dependencies
def formatdepends(jobids:str, depend='afterok') -> str:
    """Formats job ids and dependencys for an sbatch submission."""
    ## Return the formated text
    return f' --dependency={depend}:{jobids} ' if (len(jobids) > 0) else ' '

## Ftn for formating job name 
def jobname(args:list) -> str:
    """Formats an input job name, joining inputs in args with a period."""
    ## Return the job name
    return '--job-name=%s'%'.'.join(map(str,args))

## Ftn for formating command control dataframe and restarting
def commandcontrol(commands:list, toreset:bool, pipelinesteps:list, rerunfrom:str, cols=['Jobfile','Sample','Experiment','Operation','Report','Torun','JobID']) -> tuple:
    """Generates a command and control dataframe for running pipeline."""
    ## Make into a dataframe
    commanddf = pd.DataFrame(commands,columns = cols)
    ## If a hard reset was called, remove the previous file reports 
    if_hard_reset = [(remove(fr) if fileexists(fr) else None) for fr in commanddf.Report.tolist()] if toreset else None 
    ## Check if the reports exist, iterate of the rows of command df
    for rix,row in commanddf.iterrows():
        commanddf.loc[rix,'Torun'] = 1 if fileexists(row.Report) else 0
    ## If we are re running any part of the pipeline we will re code from here
    if rerunfrom: ## Set the dict 
        rerun_dict = dictzip(pipelinesteps,[pipelinesteps[i:] for i in range(len(pipelinesteps))])
        #print(rerun_dict)
        ## Recode the command dataframe to run 
        for r in rerun_dict[rerunfrom]:
            commanddf.loc[(commanddf.Operation==r),'Torun'] = 0
    ## Reformat job id as a str 
    commanddf['JobID'] = commanddf.JobID.apply(str)
    ## Return the command dataframe and hard reset 
    return commanddf, if_hard_reset

## Ftn for sub mitting fastp 
def submitfastp(command_df:pd.DataFrame, subsbatchs:list, nodepartition:str, timestamp:str, debugmode=False) -> tuple:
    """Submits fastp command(s) to SLURM via sbatch."""
    ## Check theinput sub sbatchs is a list
    assert type(subsbatchs) == list, "ERROR: We expected a list for the input submission!"
    ## See if we have a master job id
    master = command_df[(command_df.Operation=='bwaix')].JobID.min() if ('bwaix' in command_df.Operation.tolist()) else 0
    ## Iterate thru the fastp opperations 
    for rix,row in command_df[(command_df.Operation=='fastp') & (command_df.Torun==0)].iterrows():
        ## Format the jobparts
        jobparts = jobname(['fastp',rix,nodepartition,timestamp])
        ## Format the sbatch  
        if master:
            intext = f'sbatch --dependency=afterok:{master} --partition={nodepartition} {row.Jobfile}'
        else:  
            intext = f'sbatch --dependency=singleton {jobparts} --partition={nodepartition} {row.Jobfile}'
        ## Append to sub list
        subsbatchs.append(intext)
        ## Append the fastp job id 
        command_df.loc[rix,'JobID'] = fakejobid + rix if (runlocal and debugmode) else submitsbatch(intext)
        ## Print if we are in debug mode
        ifprint(intext,debugmode)
    ## Return the command df and subsbatcsh 
    return command_df, subsbatchs

## Ftn for submission of jobs with dependencies 
def submitdependency(command_df:pd.DataFrame, operation:str, dependent:list, timestamp:str, clusterpart:str, bylast=False, group='Sample', debug=False) -> list:
    """Formats and submits sbatch jobs from command dataframe based on operations and dependents"""
    ## Initilzse cap and lists
    subsbatchs, dependent = [], makelist(dependent) + ['bwaix']
    ## group the 
    for sample_ix, sample_df in command_df[(command_df.Torun==0)].groupby(group):
        ## Gather the operaiton and dependent ix
        operations, dependents = sample_df[(sample_df.Operation==operation)], sample_df[(sample_df.Operation.isin(dependent))]
        ## Gather the operations indexes 
        op_ixs = operations.index.values
        ## Iterate thru the op indexs
        for m in op_ixs:
            ## gather the dependent ix
            dependent_ix = [i for i in dependents[(dependents.index < m)].index.tolist()]
            ## Take the last index if we are doing it by lst
            dependent_ix = dependent_ix[-1:] if bylast else dependent_ix
            ## Set the row
            row = command_df.loc[m,:]
            ## Format the bwa ids
            djobids, jobparts = formatids(command_df[(command_df.index.isin(dependent_ix)) | (command_df.Operation=='bwaix')],dependent), jobname([operation,m,clusterpart,timestamp])
            ## Set the sbatch text
            intext = f'sbatch{formatdepends(djobids)}--partition={clusterpart} {jobparts} {row.Jobfile}'
            ## print the comamnd
            ifprint(intext,debug)
            ## Append to sub list
            subsbatchs.append(intext)
            ## Submit the job and append the job id to the dataframe 
            command_df.loc[m,'JobID'] = fakejobid + m if (debug and runlocal) else submitsbatch(intext)
    ## Return submitted sbatchs 
    return subsbatchs

## reformat the input bed names 
def formatbymode(inbedpe:str,mode:str) -> list: 
    ## Inaite list, return the new file locations 
    return f'{macs3dir}/{basenoext(inbedpe)}.{mode.lower()}' 

## Ftn for formating values 
def formatval(valname,val) -> str:
    """Formats and returns the named parameter and its value for a call to macs3."""
    ## Format the named value
    return f' --{valname} {val} ' if val else ' '

## Write ftn for calling macs3 with atac seq data
def peakattack(bedpe:str,n:str,report:str,mode:str,gsize='hs',incontrols=[],shiftsize=0,extendsize=0,maxgap=0,minlen=0,keepdups='all',nolambda=False,broad=False,summits=False,outdir=f'./{macs3dir}') -> list[str]: 
    """Formats a call to the macs3 callpeak function for a run of the slurpy pipeline (n) on input bedpe file, using the input genome size (g), maximum gap (ml), and minimum peak length (ml)."""
    ## Format the no model paramater
    nomodel   = ' --nomodel '  if extendsize or shiftsize else ' '
    nolambda  = ' --nolambda ' if nolambda else ' '
    ## Format the borad option and call sumits opt
    isborad   = ' --broad ' if broad else ' '
    call_sums = ' --call-summits ' if summits else ' '
    ## Reformat controls
    if len(incontrols):
        ## make the input controls into bed format 
        short_controls = [f'{slurpydir}/toshort.py --{mode.lower()} -i {control} -s {shiftsize} -e {extendsize}\n' for control in incontrols]
        ## Repath the input controls
        controls = '-c ' + ' '.join([formatbymode(control,mode) for control in incontrols])
    else:
        short_controls = []
        controls = ''

    ## Format the conversion commands to the bedpe, the macs3 callpeak command, and the echo command 
    macs_coms = [f'{slurpydir}/toshort.py --{mode.lower()} -i {bedpe} -s {shiftsize} -e {extendsize}\n',
                 f'macs3 callpeak -t {formatbymode(bedpe,mode)} {controls}{formatval("keep-dup",keepdups)}-B --SPMR{nolambda}-n {n}{isborad}-g {gsize} -f {mode} --outdir {outdir}{formatval("max-gap",maxgap)}{formatval("min-length",minlen)}{call_sums}{nomodel}2>> {report}\n', 
                 f'{slurpydir}/myecho.py Finished calling peaks in {formatbymode(bedpe,mode)} with macs3 {report}\n']
    ## Return the macs coms 
    return short_controls + macs_coms

## Set the narrow peak names
peaknames = ['Chrom','Start','End','Name','Score','Strand','Fold_change','-log10pvalue','-log10qvalue','Sumpos']

## Write ftn for loading in narrow peak file
def loadnarrowpeak(path,peakcols = peaknames,sep='\t') -> pd.DataFrame:
    """Loads in a narrow peak file from macs3."""
    ## Load in narrow peak file and return 
    return pd.read_csv(path,sep=sep,header=None,names=peakcols)

## Set description of this library and scirpt
description = 'Calculates the fraction of reads within peaks from input bedpe and narrow peaks files.'
dplace = 4

## Set extension dict 
exten_dict = dict(zip(['csv','tsv','narrowPeak','bed','txt'],[',','\t','\t','\t',' ']))

## Set help messages
b_help = "Path to input bedpe file."
p_help = "Path to input peak (BED) files from macs3."
s_help = "Path and name of output diagnostic statistics."
d_help = "Decimal place used to calcualte and save statistics (Default: %s)."%dplace

## Define ftn for parsing variables
def parse_args():
    ## Set parser
    parser = argparse.ArgumentParser(description=description)
    ## Add optional arguments
    parser.add_argument("-b",   dest="b",   required=True,    type=str,   help=b_help                   )
    parser.add_argument("-p",   dest="p",   required=True,    type=str,   help=p_help                   )
    parser.add_argument("-s",   dest="s",   required=False,   type=str,   help=s_help,    default=None  )  
    parser.add_argument("-g",   dest="g",   required=False,   type=int,   help=g_help,    default=False )
    parser.add_argument("-d",   dest="d",   required=False,   type=int,   help=d_help,    default=dplace)
    ## Parse the arguments
    return parser.parse_args()

## ----------------------------------------------- MAIN EXECUTABLE --------------------------------------------------- ## 
## If the library is called as an executable
if __name__ == "__main__":
    ## ------------------------------------------- PARSER SETTING ---------------------------------------------------- ## 
    args = parse_args()
    ## Load in mods like dask dataframes 
    import dask.dataframe as dd 

    ## ---------------------------------------- VARIABLE SETTING ---------------------------------------------------- ## 
    ## Gather inputs 
    bedpe_path, peak_path, save_path, genomesize,  dplace = args.b, args.p, args.s, args.g, args.d, 

    ## Gather the file extension
    file_end = bedpe_path.split('.')[-1]
    ## Set column names 
    col_names = ['Chrom','Left','Right']
    col_names = col_names if (file_end == 'bedpe') else (col_names + ['Strand'])

    ## Load in dask 
    bedpe = dd.read_csv(bedpe_path,sep='\t',names=col_names,header=None)

    ## Calc total
    total = bedpe.Chrom.count().compute()

    ## Gather the extenstion of the input file
    extension = peak_path.split('.')[-1]

    ## Load in narrow peak path, initate read count, and drop ducpliates to unique peaks 
    narrow          = loadnarrowpeak(peak_path,sep=exten_dict[extension])
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
    peak_info.to_csv(save_path,index=False,float_format=f'%.{dplace}f')
## End of file 