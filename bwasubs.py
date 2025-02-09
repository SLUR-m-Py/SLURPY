#!/usr/bin/env python
#SBATCH --job-name=bwa.master           ## Name of job
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=1               ## Number of tasks to be launched
#SBATCH --nice=2147483645               ## Nice parameter, sets job to lowest priority 
## Load in sys
import sys 
## Append pipelien
sys.path.append('/SLURPY/pipeline')
## Bring in ftns and variables from defaluts 
from defaults.defaults import sortglob, sbatch, submitsbatch, fileexists, splitsdir, comsdir, debugdir, bamtmpdir, pipelinedir
## Load in write to file from pysam tools 
from defaults.tools.pysamtools import writetofile
## load in sleep
from time import sleep

## Set opts
bwa_options = '-5SMP'   
line_count  = 5000
pix         = 1     

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
def bwamaster(sname:str,refpath:str,library:str,threads:int,cwd:str,partition:str,debug:bool,pix=pix,linecount=line_count):
    command = f'{pipelinedir}/bwasubs.py -s {sname} -r {refpath} -b {threads} -c {cwd} -P {partition} -L {library} -l {linecount} ' +  ('--debug' if debug else '')
    report  = f'{debugdir}/{pix}.bwa.master.{sname}.log'
    return [command], report 

## Set description
bwadescr = 'A submission script that formats bwa/bedpe commands for paired fastq file from fastp splits of a given sample.'

## Load inputs
from .defaults.parameters import s_help,r_help,b_help,P_help,L_help, debug_help, refmetavar,bwathreads,part, lib_default

c_help = 'The current working directory'
l_help = 'The number of lines from bwa to buffer in list. Default is: %s'%line_count

##      MAIN SCRIPT & ARGUMENT PARSING 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## If the script is envoked by name 
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse
    ## Make the parse
    parser = argparse.ArgumentParser(description = bwadescr)

    ## Add the required argument
    parser.add_argument("-s", "--sample",         dest="s",     type=str,  required=True,  help = s_help, metavar = 'sample name'                                  )
    parser.add_argument("-r", "--refix",          dest="r",     type=str,  required=True,  help = r_help, metavar = refmetavar                                     ) 
    parser.add_argument("-c", "--cwd",            dest="c",     type=str,  required=True,  help = c_help, metavar = './the/cwd'                                    )
    parser.add_argument("-b", "--bwa-threads",    dest="b",     type=int,  required=False, help = b_help, metavar = bwathreads,             default = bwathreads   )
    parser.add_argument("-P", "--partition",      dest="P",     type=str,  required=False, help = P_help, metavar = part,                   default = part         ) 
    parser.add_argument("-L", "--library",        dest="L",     type=str,  required=False, help = L_help, metavar = 'MboI',                 default = lib_default  )
    parser.add_argument("-l", "--line-count",     dest="l",     type=int,  required=False, help = l_help, metavar = 'n',                    default = line_count   )
    parser.add_argument("--debug",                dest="debug",     help = debug_help,    action = 'store_true'                                                    )
    
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
    debug        = inputs.debug 

    ## Gather the first reads 
    read_ones  = sortglob(f'{splitsdir}/*{sample_name}*_R1_*.fastq.gz')
    read_twos  = ['_R2_'.join(r.split('_R1_')) for r in read_ones]

    ## List the pairs 
    read_pairs = list(zip(read_ones,read_twos))

    ## Iniate list 
    bwa_reports = []

    ## Iterate thru the read pairs 
    for i, (r1,r2) in enumerate(read_pairs):
        ## Set the out file
        outfile   = f'{bamtmpdir}/{i}.{sample_name}.bedpe'
        ## format the command 
        bwa_com   = f'bwa mem -v 1 -t {thread_count-1} {bwa_options} {ref_path} {r1} {r2} | {pipelinedir}/tobedpe.py {ref_path} {library} {outfile} {line_count}\n## EOF'
        bwa_repo  = f'{debugdir}/{pix}.bwa.{i}.{sample_name}.log'
        bwa_file  = f'{comsdir}/{pix}.bwa.{i}.{sample_name}.sh' 

        ## Write the bwa command to file 
        writetofile(bwa_file, sbatch(None,thread_count,the_cwd,bwa_repo) + [bwa_com], debug)

        ## Submit the command to SLURM
        submitsbatch(f'sbatch --partition={partitions} {bwa_file}')

        ## append the report
        bwa_reports.append(bwa_repo)
        ## Sleep
        sleep(1)

    ## Check the reports
    kicker = True 

    ## While thekicker is true 
    while kicker:
        kicker = not reportcheck(bwa_reports)
        ## Wait a minitue 
        sleep(10)

    ## Print to log 
    print("Finished %s bwa submissions for sample: %s"%(len(bwa_reports),sample_name))
## EOF 