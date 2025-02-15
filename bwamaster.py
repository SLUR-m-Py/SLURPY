#!/usr/bin/env python
#SBATCH --job-name=bwa.master           ## Name of job
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=1               ## Number of tasks to be launched
#SBATCH --nice=2147483645               ## Nice parameter, sets job to lowest priority 
## Bring in ftns and variables from defaluts 
from defaults import sortglob, sbatch, submitsbatch, fileexists, getfilesize
from directories import splitsdir, comsdir, debugdir, slurpydir, bedtmpdir
## Load in write to file from pysam tools 
from pysamtools import writetofile, listzip
## load in sleep
from time import sleep
## Load input vars from params
from parameters import refmetavar, bwathreads, lib_default, nice, hic_options, waittime

## Set opttions 
line_count  = 5000
pix         = 1     

## Ftn for checkign if single 
def issingle(r:str)-> bool:
    return '.singletons.' in r.lower()

## Ftn for checking if failed in 
def isfailed(r:str) -> bool:
    return '.failed.' in r.lower() 

## Dftn for getting read1
def getread1(wc:str) -> list:
    return [r for r in sortglob(wc) if not (issingle(r) or isfailed(r))]

## Ftn for defining read two
def formatread2(firstreads:list) -> list:
    return ['_R2_'.join(r.split('_R1_')) for r in firstreads]

## Def for checkign report 
def reportcheck(reportpaths) -> bool:
    finreports = [r for r in reportpaths if fileexists(r)]
    ## If the report list lengths match 
    if len(finreports) == len(reportpaths):
        finished = []
        for r in finreports:
            with open(r,'r') as infile:
                finished.append(infile.readlines()[-1].startswith('Finished'))
                infile.close()
        ## Sum the kicker 
        kicker = (sum(finished) == len(finreports))
    else: ## other wise keep kicken 
        kicker = False
    return kicker

## Ftn for formating the bwa master 
def bwamaster(sname:str,refpath:str,threads:int,cwd:str,partition:str,debug:bool,nice:int,inhic=False,pix=pix,linecount=line_count,library=None,forced=False,nodelist=None):
    ## Format command 
    command = f'{slurpydir}/bwamaster.py -s {sname} -r {refpath} -b {threads} -c {cwd} -P {partition} -N {nice} -l {linecount}' + (f' -L {library}' if library else '') + (' --debug' if debug else '') + (' --hic' if inhic else '') + (' --force' if forced else '') + (' --nodelist %s'%' '.join(nodelist) if nodelist else '')
    ## Format report 
    report  = f'{debugdir}/{pix}.bwa.master.{sname}.log'
    return [command], report 

## Ftn for echo bwa finish to a log
def bwaecho(o,l=None) -> str:
    """Formats an echo command to print finishing statment of bwa alignment to log."""
    ## Format the message given the output (o) file name and log (l)
    return  f'{slurpydir}/myecho.py Finished alignment of split: {o} {l}' if l else ''
    
## Write ftn for formating bwa mem command with paired reads
def bwamem_paired(r1,r2,ref,outbam,log,vmode=1,threads=4,opts='-M') -> list[str]:
    """Formats a bwa command and requires as input paired-reads, a reference name, an output bam file, a log file, experiment type."""
    ## Set the bwa mem option based on experiment type and format the bwa mem call, submit to shell via submit command
    return [f'bwa mem {opts} -v {vmode} -t {threads} {ref} {r1} {r2} 2>> {log} | samtools view -hb -@ {threads} -o {outbam} -O BAM\n', bwaecho(outbam,log)]

## Write ftn to check filesize fo fastq
def sizecheck(read1,read2) -> list: 
    ## Iniate list 
    read_pairs = []
    ## Iterat thru the read pairs 
    for r1,r2 in listzip(read1,read2):
        ## If either is zero 
        if (getfilesize(r1) == 0) | (getfilesize(r2) == 0):
            ## Print a warning 
            print("WARNING: the given number of read pairs from %s or %s was zero."%(r1,r2))
        else:
            read_pairs.append((r1,r2))
    ## Return the read pairs 
    return read_pairs

## Set description
bwadescr = 'A submission script that formats bwa/bedpe commands for paired fastq file from fastp splits of a given sample.'

## Load in help messages from parameters 
from parameters import s_help, r_help, b_help, P_help, L_help, N_help, debug_help, force_help, node_help

## Set help messages 
c_help     = 'The current working directory'
l_help     = 'The number of lines from bwa to buffer in list. Default is: %s'%line_count
hic_flag   = 'Flag to run in Hi-C mode.'

##      MAIN SCRIPT & ARGUMENT PARSING 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## If the script is envoked by name 
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse
    ## Make the parse
    parser = argparse.ArgumentParser(description = bwadescr)

    ## Add the required argument
    parser.add_argument("-s", "--sample",         dest="s",     type=str,  required=True,  help = s_help, metavar = 'sample name'                        )
    parser.add_argument("-r", "--refix",          dest="r",     type=str,  required=True,  help = r_help, metavar = refmetavar                           ) 
    parser.add_argument("-c", "--cwd",            dest="c",     type=str,  required=True,  help = c_help, metavar = './the/cwd'                          )
    parser.add_argument("-b", "--bwa-threads",    dest="b",     type=int,  required=False, help = b_help, metavar = bwathreads,   default = bwathreads   )
    parser.add_argument("-P", "--partition",      dest="P",     type=str,  required=False, help = P_help, metavar = 'tb',         default = 'tb'         ) 
    parser.add_argument("-L", "--library",        dest="L",     type=str,  required=False, help = L_help, metavar = 'MboI',       default = lib_default  )
    parser.add_argument("-l", "--line-count",     dest="l",     type=int,  required=False, help = l_help, metavar = 'n',          default = line_count   )
    parser.add_argument("-N", "--nice",           dest="N",     type=int,  required=False, help = N_help, metavar = 'n',          default = nice         )
    parser.add_argument("--nodelist",             dest="nodes", nargs='+', required=False, help = node_help,                      default = None         )
    parser.add_argument("--debug",                dest="debug",     help = debug_help,    action = 'store_true'                                          )
    parser.add_argument("--hic",                  dest="hic",       help = hic_flag,      action = 'store_true'                                          )
    parser.add_argument("--force",                dest="force",     help = force_help,    action = 'store_true'                                          )
    ## Set the paresed values as inputs
    inputs = parser.parse_args() 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    
    ## Set inputs 
    sample_name  = inputs.s
    ref_path     = inputs.r
    thread_count = inputs.b
    the_cwd      = inputs.c 
    partitions   = inputs.P
    library      = inputs.L
    line_count   = inputs.l
    nice         = inputs.N
    nodes        = inputs.nodes
    debug        = inputs.debug 
    ishic        = inputs.hic 
    forced       = inputs.force 

    ## Gather the first reads 
    read_ones  = getread1(f'{splitsdir}/*.{sample_name}_R1_*.fastq.gz')
    read_twos  = formatread2(read_ones)
    
    ## Print to file 
    print('INFO: Spawning %s calls to bwa.'%len(read_twos))
    ## Check the size of the read pairs 
    read_pairs = sizecheck(read_ones,read_twos)

    ## Iniate list 
    bwa_reports = []
    ## if we are formating hic run
    options = hic_options if ishic else '-M'
        
    ## Iterate thru the read pairs 
    for i, (r1,r2) in enumerate(read_pairs):
        ## SEt report and file name 
        bwa_repo = f'{debugdir}/{pix}.bwa.{i}.{sample_name}.log'
        bwa_file = f'{comsdir}/{pix}.bwa.{i}.{sample_name}.sh' 
        outfile  = f'{bedtmpdir}/{i}.{sample_name}.bedpe'

        ## format the command 
        bwa_coms = [f'bwa mem -v 1 -t {thread_count-1} {options} {ref_path} {r1} {r2} | {slurpydir}/tobedpe.py {ref_path} {library} {outfile} {line_count}\n## EOF']

        ## If the report exists and has alredy been run, just skip
        prekick = reportcheck([bwa_repo])
        if (prekick and fileexists(bwa_file) and fileexists(outfile)) and (not forced):
            print(f'WARNING: Detected a finished run ({outfile}) from {bwa_file} in {bwa_repo}.\nINFO: Skipping.\n')
            continue

        ## Write the bwa command to file 
        writetofile(bwa_file, sbatch(None,thread_count,the_cwd,bwa_repo,nice=nice,nodelist=nodes) + bwa_coms, debug)
        ## Sleep here to wait
        sleep(waittime)
        ## append the report
        bwa_reports.append(bwa_repo)

        ## Submit the command to SLURM
        submitsbatch(f'sbatch --partition={partitions} {bwa_file}')

    ## Check the reports
    kicker = True 

    ## While thekicker is true 
    while kicker and len(bwa_reports):
        kicker = not reportcheck(bwa_reports)
        ## Wait a minitue 
        sleep(2*waittime)

    ## Print to log 
    print("Finished %s bwa submissions for sample: %s"%(len(bwa_reports),sample_name))
## EOF 