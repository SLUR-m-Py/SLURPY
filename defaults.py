"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      SLURPY PIPELINE FUNCTIONS 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      MODULE LOADING 
## Load in directories
from directories import *
## Load in glob
from glob import glob 
## Bring in sub-process mod
import subprocess, sys 
## Bring in basename, get file size,  and pathexists 
from os.path import basename, exists as pathexists, getsize as getfilesize
## Bring in make dirs
from os import makedirs, remove 
## Bring in our costum ftns from pysamtools 
from pysamtools import ifprint, splitbam, splitsam, makelist, dictzip, getprimary
## Load in pandas
import pandas as pd 
## Load in date and time
from datetime import datetime
## Load in rm tree
from shutil import rmtree
## Bring in vars from prep run
from parameters import fakejobid, runlocal
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      SLURPY FUNCTIONS 
## Ftn for submitting que checks
def submitter(command:str) -> list:
    ## Sub process check output 
    return subprocess.check_output(command, shell=True).decode("utf-8").split('\n')   

## Ftn for cehcking queue 
#def checkqueue(step:str) -> int:
#    ## Sum the setps in l 
#    return sum([step in l for l in submitter('squeue')]) 

## Ftn for sorting fastq files by size
def sortfastq(infastq:list,splitsizes:list):
    ## Format input fastq zipped list into df
    tmp = pd.DataFrame(infastq,columns=['Read1','Read2'])
    ## Gather sizes of fast
    tmp['Size1'] = [getfilesize(f1) for (f1,f2) in infastq]
    ## Sort values by size and reset index 
    tmp = tmp.sort_values('Size1').reset_index(drop=True)
    ## Modify split sizzes 
    splitsizes = [splitsizes[0] for s in range(tmp.shape[0])] if (len(splitsizes) < tmp.shape[0]) else splitsizes
    ## Add a sorted list to splitsize to match size of fastqs 
    tmp['Splitsize'] = sorted([int(s)*4 for s in splitsizes])
    ## Return tmp
    return tmp 

## Ftn for returning a sorted glob
def sortglob(wildcard:str) -> list:
    """Retuns the sorted glob of input wild card."""
    ## Returns the sorted glob of input
    return sorted(glob(wildcard))

## Ftn for reseting file
def reset(infiles:list) -> list:
    """Removes the files given with in the intput infiles."""
    ## Return the none type
    return [(remove(f) if pathexists(f) else None) for f in infiles]

## Ftn for getting sample name
def getsamplename(read1:str, spliton='_R1_') -> str:
    """Get sample name from a fastq.gz file."""
    ## Return the name 
    return basename(read1).split(spliton)[0]

## Ftn for returning just the bam file name without the .bam
def basenobam(inbam:str) -> str:
    """Returns the basename of an input bam file without the .bam file extension."""
    ## Return the split basename
    return basename(splitbam(inbam))

## Ftn for returning just the bam file name without the .bam
def basenosam(inbam:str) -> str:
    """Returns the basename of an input bam file without the .bam file extension."""
    ## Return the split basename
    return basename(splitsam(inbam))

## Ftn for takeing basename 
def basenoext(inpath:str) -> str:
    """Returns the basename of an input file with no extension."""
    ## Returns the basename of in .sh file with no .sh
    return '.'.join(basename(inpath).split('.')[:-1])

## Write a ftn to check the file size meets our threshold
def checkfilesize(filepath:str, threshold=5) -> bool:
    """Returns a boolean on the input file path checking if a given file is larger (in bytes) than threshold."""
    ## Return the boolean of file size
    return getfilesize(filepath) > threshold

## Write ftn for checking file is a real deal!
def fileexists(filepath:str) -> bool:
    """Checks the existance and size of an input file (foodstuffs)."""
    ## Via os, check if the path to "food" exists and meets our size threshold
    return pathexists(filepath) and checkfilesize(filepath)

## Ftn to adds forward slash
def addfslash(inpath:str) -> str:
    """Corrects input inpath with a forward slash."""
    ## Returns a path if the last chracter is a forward path 
    return inpath if (inpath[-1] == '/') else inpath + '/'

## Write dict on extensions
extdict = dictzip(['tsv','fai','bed','txt','csv'],['\t','\t','\t',' ',','])

## Ftn set setp
def setsep(ext:str) -> str:
    """Given input extension returns the deliminator."""
    ## Return the extension 
    return extdict[ext]

## Ftn for reading in table 
def readtable(inpath:str,header=None,names=None) -> pd.DataFrame:
    """Reads in a table give its extension."""
    ## Return the table 
    return pd.read_csv(inpath,sep=setsep(inpath.split('.')[-1]),header=header,names=names)

## Ftn for checking fastp threads and splits
def checkfastp(fsplits:int,fthreads:int) -> tuple: 
    """Checks that the number of splits is divisable by the number of threads."""
    fthreads = fthreads if fthreads in [2,4,6,8] else 8
    ## Check our thread count 
    assert fthreads, "ERROR: The number of fastp threads must be an even integer between (and including) 2 and 8!"
    ## If needed reset that the fastp threads and splits such that they are a multiple of each
    if (fsplits%fthreads != 0):
        print(f'WARNING: The number of splits set for fastp ({fsplits}) is not an even multiple of the thread count ({fthreads}).')
        ## Reset the thread count to half the split number
        fsplits = fthreads
        ## Print we are resetting
        print(f'INFO: The number fastp threads to conduct {fsplits} splits was reset to {fthreads}.')
    ## Return the number of splits and threads for fastp 
    return fsplits,fthreads

## Ftn for reseting clean bool
def patchclean(torerun:str,inclean:bool) -> bool:
    """Patches the clean boolean if the rerun step is clean or if the input clean boolean is true."""
    ## Patches the clean step
    return True if (torerun == 'clean') else inclean

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

## Ftn for setting run name
def setrunname(runname:str, incwd:str) -> str:
    """Sets the name of the pipeline run (if not given) from the current working directory."""
    ## Return the run name if it is defined, else use the basename of the cwd
    return runname if runname else basename(incwd)

## Ftn for making a list from zipped items 
def listzip(a:list, b:list) -> list:
    """Returns a list of zipped items in a and b."""
    ## Return the list of zipped items 
    return list(zip(a,b))

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

## Ftn for calling sub-process
def submitsbatch(pathtoscript:str) -> int:
    """Checks output from sub-process submission of input path-to-script assuming it is an sbatch command. Returns the job id."""
    ## call sub process
    return subprocess.check_output(pathtoscript, shell=True).decode("utf-8").split('\n')[0].split(' ')[-1]

## Ftn for making items within a list into strings
def makestrs(inlist:list) -> list:
    """Formats all items within input list into a string."""
    ## Makes all items within input list into a string
    return [str(i) for i in inlist]

## Formats digit 
def formatn(n:float) -> str:
    """Formats input digit (float) to have 4 decimal places."""
    ## Format an integer to 4 deciimal places
    return "{:04d}".format(n)

## Ftn for checking if a file is bwa indexed
def isbwaix(inref:str, indexends = ['amb', 'ann', 'bwt', 'pac', 'sa']) -> int:
    """Checks if the the input reference file is indexed."""
    ## Check to make sure the reference exists
    assert fileexists(inref), "ERROR: Unable to locate input reference file on path: %s"%inref
    ## Check if each of these files exist: amb, ann, bwt, pac and sa
    return sum([pathexists(inref+'.'+fe) for fe in indexends]) == len(indexends)

## Ftn for readin an ann file from bwa index
def readann(inpath:str) -> list:
    """Parses an input reference .ann file from bwa index to generate chromosome dataframe of contig names and lengths."""
    ## open and read in the lines
    with open(inpath,'r') as inhandle:
        ## Take up to the third column of each line
        newlines = [l.split(' ')[:3] for l in inhandle]
        inhandle.close()
    ## Gather chromosomes
    contigs = [l[1] for l in newlines[1::2]]
    lengths = [int(l[1]) for l in newlines[2::2]]
    ## Format and return a dataframe
    return pd.DataFrame(list(zip(contigs,lengths)))
    
## Ftn for formating an sbatch text
def sbatch(nameojob:str, cpus:int, cwd:str, report:str, partition=None, nodes=1, tasks=1, runtime='200:00:00',nice=10**7, nodelist=None, memory=None) -> list:
    """Generates the sbatch settings for a script with a give input jobname, cpu count, working directory, and others."""
    ## Gather the extension
    jobext = nameojob.split('.')[-1] if nameojob else 'sh'
    ## Set the runner
    runner = 'python' if (jobext == 'py') else 'bash'
    ## return the formated sbatch txt
    settings = [f'#!/usr/bin/env {runner}\n',                                                 ##      The shebang
                 '#SBATCH --job-name=%s\n'%basenoext(nameojob) if nameojob else '\n',         ##      Name of job 
                 '#SBATCH -o %s\n'%report,                                                    ##      Set the debug dir and report anme 
                 '#SBATCH --nodes=%s\n'%str(nodes),                                           ##      Number of nodes
                 '#SBATCH --ntasks-per-node=%s\n'%str(tasks),                                 ##      Tasks per node
                 '#SBATCH --cpus-per-task=%s\n'%str(cpus),                                    ##      Number of cpus
                 '#SBATCH --time=%s\n'%runtime,                                               ##      Set nice parameter
                 '#SBATCH --nice=%s\n'%str(nice)]                                             ##      The allowed run time of the job 
    ## Add the current working dir                                                            ##  
    settings = settings + ['#SBATCH --chdir=%s\n'%cwd] if cwd else settings                   ##      Set the current workign dir
    ## Add the partition                                                                      ##      
    settings = settings + ['#SBATCH --partition=%s\n'%partition] if partition else settings   ##      The partitions
    ## Add nodes if they were passed
    settings = settings + ['#SBATCH --nodelist=%s\n'%','.join(nodelist)] if nodelist else settings ## List of nodes to run on
    ## Add the max memory
    settings = settings + ['#SBATCH --mem=%s'%memory] if memory else settings
    ## return settings 
    return settings

"""
 f'#SBATCH --output={errordir}/%x.%j.out' +'\n',                               ##      Output log 
 f'#SBATCH --error={errordir}/%x.%j.err' + '\n',                               ##      Error log
"""

## Ftn for the begining of a call to fastp
def fastcut(i:str, I:str, o:str, O:str, r:str, n:int) -> tuple:
    """Formats the start of a call to fastp given inputs (i,I), outputs (o,O), and the report name (r)."""
    ## Format the fall and return
    return f'fastp -i {i} -I {I} -o {o} -O {O} -j {r}.{n}.json -h {r}.{n}.html', f'{r}.{n}.json', f'{r}.{n}.html'

## Ftn for returning a boolean if in hic
def inhic(inexp:str) -> bool:
    """Boolean if the input experiment is hi-c."""
    ## Return the hic 
    return inexp == 'hic'

## Write ftn for logging fastp has finished
def fastdry(r1:str, r2:str, report:str) -> str:
    """Formats an echo command for logging the completion of a fastp split."""
    ## reformat the report 
    #report = report if report.split('.')[-1] == 'txt' else report + '.txt'
    ## Format and return the echo command
    return f'{slurpydir}/myecho.py Finished filtering and splitting on: {r1} {r2} {report}\n'

## Ftn for formating report
def reportname(inbam:str, script:str, i=0) -> str:
    """Formats the report name given the input bam file name and type of script (str)."""
    ## Return the report name
    return f'{debugdir}/{i}.{script}.{basenobam(inbam)}.log'

## Ftn for loading in script
def loadscript(inpath:str) -> list:
    """Loads in lines from a script."""
    ## Open the input path for reading 
    with open(inpath,'r') as infile:
        ## Read in lines 
        scriptlines = infile.readlines()
    ## Close the file 
    infile.close()
    ## Return all but the first line 
    return scriptlines[1:] + ['\n']

## Ftn for callign samblaster
def samblaster(inbam,outbam,report,threads) -> str:
    """Formats a call to samtools and samblaster given inputs."""
    ## Return the samtools and samblaster command 
    return f'samtools view -@ {threads} -h {inbam} | {slurpydir}/samblaster --ignoreUnmated -M 2>> {report} | samtools view -@ {threads} -Shb | samtools sort -@ {threads} - -o {outbam} -O BAM --write-index\n'

## Ftn for formating samblaster command
def markduplicates(inbam:str, threads:int, pix:int, script='mark') -> tuple:
    """Formats and submits a samtools and samblaster command to mark duplicates on input bam file (i) and saves to output (o)."""
    ## Format the output bams and the report name 
    outbam, report = f'{aligndir}/{basenobam(inbam)}.marked.bam', reportname(inbam,script,i=pix) 
    ## Format the sam-blaster and echo command
    blast_command, echo_command = samblaster(inbam,outbam,report,threads), f'{slurpydir}/myecho.py Finished marking duplicates in {outbam} {report}\n'
    ## Return the samblaster command and ecco chommand 
    return outbam, [blast_command, echo_command], report

## Ftn for formating echo split
def splitecho(inbam:str, report: str, script:str) -> str:
    """Formats an echo statment for the split command."""
    ## Return the formated command 
    return f'{slurpydir}/myecho.py Finished splitting {inbam} using {script} {report}\n'

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

## Write ftn for making head path
def headpath(inpath:str) -> str:
    """Returns the head path of a path of type string."""
    ## Split the path take up to the last item and reformat
    return '/'.join(inpath.split('/')[:-1]) 

## Ftn for formating bwa index
def bwaindex(refpath:str,script='index') -> tuple:
    """Formats command to index a reference via bwa."""
    ## return the index commands
    return [f'bwa index {refpath}\n', f'{slurpydir}/myecho.py Finished indexing reference on path {refpath} {reportname(refpath,script)}\n'], reportname(refpath,script)

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

## Ftn for filtering bam fie 
def filterbam(inbam:str, M:str, threads:int, chrlist:list, pix:int) -> tuple:
    """Formats a command to filter an input bam file seperating aligments on quality."""
    ## Format the report name and out bam name 
    report, outbam = reportname(inbam,'filter',i=pix), splitbam(inbam) + f'.primary.q{M}.bam'
    ## Format filter command and the echo command 
    bam_filter_command, echo_command = getprimary(inbam,M,threads,outbam,chroms=chrlist), f'{slurpydir}/myecho.py Finished filtering {inbam} at mapping quality of {M} {report}\n'
    ## Format and return commands
    return outbam, [bam_filter_command, echo_command], report

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

## Set the sub list error message
sublist_err = "ERROR: We expected a list for the input submission!"

## Ftn for sub mitting fastp 
def submitfastp(command_df:pd.DataFrame, subsbatchs:list, nodepartition:str, timestamp:str, debugmode=False) -> tuple:
    """Submits fastp command(s) to SLURM via sbatch."""
    ## Check theinput sub sbatchs is a list
    assert type(subsbatchs) == list, sublist_err
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

## Ftn for sub mitting bwa commadns 
def submitbwa(command_df:pd.DataFrame, subsbatchs:list, nodepartition:str, timestamp:str, nruns:int, debugmode=False) -> tuple:
    """Submits bwa commands to SLURM via sbatch."""
    ## Check theinput sub sbatchs is a list
    assert type(subsbatchs) == list, sublist_err
    ## Iterate thru the rows of the command file 
    for s,rows in command_df[(command_df.Torun==0) & (command_df.Operation=='bwa')].groupby('Sample'):
        ## Gather the bwa indexs for this sampe 
        bwaindexs = rows.index
        ## Group the bwa runs for this sample
        groupedbwa = [bwaindexs[i::nruns] for i in range(nruns)]
        ## Iterate over the grouped bwa smaples 
        for gi, gbwa in enumerate(groupedbwa):
            for i, bwai in enumerate(gbwa):
                ## Set the row from the command files 
                row = command_df.loc[bwai,:]
                ## Format the job parts
                jobparts = jobname(['bwa',gi,nodepartition,timestamp])
                ## Gather the fastpids 
                fastpids = formatids(command_df[(command_df.Sample==s) | (command_df.Operation=='bwaix')],['fastp','bwaix']) 
                ## If this is the first in the bwa submission, add the fastp dependency, otherwise make the bwa a named singleton 
                intext = f'sbatch{formatdepends(fastpids)}{jobparts} --partition={nodepartition} {row.Jobfile}' if (i==0) else f'sbatch --dependency=singleton {jobparts} --partition={nodepartition} {row.Jobfile}'
                ## Append to sub list
                subsbatchs.append(intext)
                ## Submit the job and append the job id to the dataframe 
                command_df.loc[bwai,'JobID'] = fakejobid + bwai if (runlocal and debugmode) else submitsbatch(intext)
                ## Print if we are in debug mode
                ifprint(intext,debugmode)
    ## Return command df and subsbatchs 
    return command_df, subsbatchs 

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

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## 
## JUICER PARAMETERS
## Set Juicer columns and data types 
juicer_cols  = [ 'Str1','Chr1','Pos1','Frag1','Str2','Chr2','Pos2','Frag2','Mapq1','Cigar1','Seq1','Mapq2','Cigar2','Seq2','Qname1']
juicer_types = [   int,  str,   int,    int,   int,    str,  int,   int,    int,     str,    str,    int,    str,    str,    str  ]

## Format datatypes into dict
juicer_type_dict = dictzip(juicer_cols,juicer_types)
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## End of file 