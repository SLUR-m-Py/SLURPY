#!/usr/bin/env python
#SBATCH --job-name=bwa.master           ## Name of job
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=1               ## Number of tasks to be launched
#SBATCH --nice=2147483645               ## Nice parameter, sets job to lowest priority 
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## Bring in ftns and variables from defaluts 
from defaults import sortglob, fileexists, submitsbatch, sbatch, listzip, ifprint
## Load input vars from params
from parameters import refmetavar, bwathreads, nice, hic_options, waittime, nparallel
## Load in directiory
from parameters import splitsdir, comsdir, debugdir, slurpydir, bedtmpdir
## Load in help messages from parameters 
from parameters import ST, s_help, r_help, b_help, P_help, L_help, N_help, debug_help, force_help, node_help, j_help, B_help, slurmem_help
## load in sleep
from time import sleep
## Bring in tile and arange 
from numpy import tile, arange
## Load in file size from os path
from os.path import getsize
## load in remove
from os import remove
## Bring in myecho
from myecho import writetofile
## Bring in unfinished
from checkwork import unfinished
## Bring in argparse and set parser
import argparse

## Set opttions for line count and step in pipelien 
line_count  = 5000
pix         = 1     

## Set description and help messages for this script
bwadescr   = 'A submission script that formats bwa/bedpe commands for paired fastq file from fastp splits of a given sample.'
c_help     = 'The current working directory'
l_help     = 'The number of lines from bwa to buffer in list. Default is: %s'%line_count
hic_flag   = 'Flag to run in Hi-C mode.'

## Ftn for checkign if single 
def issingle(r:str)-> bool: return '.singletons.' in r.lower()

## Ftn for checking if failed in 
def isfailed(r:str) -> bool: return '.failed.' in r.lower() 

## Dftn for getting read1
def getread1(wc:str) -> list: return [r for r in sortglob(wc) if not (issingle(r) or isfailed(r))]

## Ftn for defining read two
def formatread2(firstreads:list) -> list: return ['_R2_'.join(r.split('_R1_')) for r in firstreads]

## Write ftn to check filesize fo fastq
def sizecheck(read1,read2) -> list: 
    ## Iniate list 
    read_pairs = []
    ## Iterat thru the read pairs 
    for r1,r2 in listzip(read1,read2):
        ## If either is zero 
        if (getsize(r1) == 0) | (getsize(r2) == 0):
            ## Print a warning 
            print("WARNING: the given number of read pairs from %s or %s was zero."%(r1,r2))
        else:
            read_pairs.append((r1,r2))
    ## Return the read pairs 
    return read_pairs

## Ftn for tiling jobs
def vectortile(k,n): return tile(arange(k)+1,n)

## ftn for making sure reports are finished
def unfinishedreports(reports) -> bool: return bool(sum([unfinished(f) for f in reports]))

## Ftn for formating the bwa master 
def bwamaster(sname:str,refpath:str,threads:int,cwd:str,partition:str,debug:bool,nice:int,njobs:int,pix=pix,linecount=line_count,forced=False,nodelist=None,bwaopts='',memory=None) -> tuple[list[str], str]:
    ## Format command to the bwatobedpe.py 
    command = f'{slurpydir}/bwatobedpe.py -s {sname} -r {refpath} -b {threads} -c {cwd} -P {partition} -N {nice} -l {linecount} -j {njobs}' \
            + (' --debug' if debug else '') + (' --force' if forced else '') + (' --nodelist %s'%' '.join(nodelist) if nodelist else '') \
            + (' -B %s'%bwaopts if len(bwaopts) else '') + (' --memory %s'%memory if memory else '')
    ## Format report 
    report  = f'{debugdir}/{pix}.bwa.to.bedpe.{sname}.log'
    ## Return the command and report 
    return [command+'\n'], report 

def parse_args():
    ## Make the parse
    parser = argparse.ArgumentParser(description = bwadescr)

    ## Add the required argument
    parser.add_argument("-s", "--sample",         dest="s",     type=str,  required=True,  help = s_help, metavar = 'sample name'                        )
    parser.add_argument("-r", "--refix",          dest="r",     type=str,  required=True,  help = r_help, metavar = refmetavar                           ) 
    parser.add_argument("-c", "--cwd",            dest="c",     type=str,  required=True,  help = c_help, metavar = './the/cwd'                          )
    parser.add_argument("-b", "--bwa-threads",    dest="b",     type=int,  required=False, help = b_help, metavar = bwathreads,   default = bwathreads   )
    parser.add_argument("-P", "--partition",      dest="P",     type=str,  required=False, help = P_help, metavar = 'tb',         default = 'tb'         ) 
    parser.add_argument("-l", "--line-count",     dest="l",     type=int,  required=False, help = l_help, metavar = 'n',          default = line_count   )
    parser.add_argument("-N", "--nice",           dest="N",     type=int,  required=False, help = N_help, metavar = 'n',          default = nice         )
    parser.add_argument("-j", "--njobs",          dest="j",     type=int,  required=False, help = j_help, metavar = 'n',          default = nparallel    )
    parser.add_argument("-B", "--bwa-opts",       dest="bwa",   type=str,  required=False, help = B_help, metavar = hic_options,  default = hic_options  )
    parser.add_argument("--nodelist",             dest="nodes", nargs='+', required=False, help = node_help,                      default = None         )
    parser.add_argument("--memory",               dest="mem",   type=str,  required=False, help = slurmem_help, metavar = '40G',  default = None         )

    ## Add boolean vars 
    parser.add_argument("--debug",                dest="debug",     help = debug_help,    action = ST                                          )
    parser.add_argument("--force",                dest="force",     help = force_help,    action = ST                                          )
    ## Set the paresed values as inputs
    return parser.parse_args() 

## MAIN SCRIPT & ARGUMENT PARSING 
def main():
    ## Set the paresed values as inputs
    inputs = parse_args() 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Set inputs 
    sample_name  = inputs.s         ## Set sample name
    ref_path     = inputs.r         ## Path to reference file
    thread_count = inputs.b         ## Number of threads
    the_cwd      = inputs.c         ## The cwd
    partitions   = inputs.P         ## Paritions to run nodes on 
    line_count   = inputs.l         ## Number of lines to parse from file
    nice         = inputs.N         ## Nice parameter, "Its nice to be nice"
    nparallel    = inputs.j         ## Number of parallel bwa jobs to run 
    bwa_opts     = inputs.bwa       ## String of options to also feed into bwa 
    nodes        = inputs.nodes     ## Node list 
    memory       = inputs.mem       ## Job memory in slurm batchs 
    debug        = inputs.debug     ## To debug, or not debug
    forced       = inputs.force     ## Use the force Luke, let go Luke. Luke trust me.
    memory       = inputs.mem       ## Set SLURM mem limit 

    ## Format the options for bwa mem
    options = f'-v 1 -t {thread_count} ' + ' '.join(bwa_opts.split(","))

    ## Gather the first reads 
    read_ones  = getread1(f'{splitsdir}/*.{sample_name}_R1_*.fastq.gz')
    ## Format the second read in pair and calculate len of read two
    read_twos  = formatread2(read_ones)
    nreadtwo = len(read_twos)

    ## Check the size of the read pairs 
    read_pairs = sizecheck(read_ones,read_twos)
    nreads  = len(read_pairs)
    assert nreads, "ERROR: No read pairs were found to have reads!"

    ## Calc difference in size 
    sizedif = nreadtwo-nreads 
    sizewarn = 'WARNING: %s read pairs (in .fastq.gz) were detected to have no reads.'%(sizedif)
    ## PRint if the size dif is nonzero 
    ifprint(sizewarn,sizedif)

    ## Iniate list of bwa files and reports 
    bwa_files   = []
    bwa_repos   = []

     ## gather job numbers / names
    job_numbers = vectortile(nparallel,nreads)
    job_names   = [f'{job_numbers[i]}.bwa.sh' for i in range(nreads)]

    ## Iterate thru the read pairs 
    for i, (r1,r2) in enumerate(read_pairs):
        out_file_prefix = f'{i}.{sample_name}'
        ## SEt report and file name 
        bwa_repo  = f'{debugdir}/{pix}.bwa.{out_file_prefix}.log'
        bwa_file  = f'{comsdir}/{pix}.bwa.{out_file_prefix}.sh' 
        outfile   = f'{bedtmpdir}/{out_file_prefix}.bedpe'

        ## format the command to bwa mem 
        bwa_coms = [f'refpath={ref_path}\n',
                    f'fread={r1}\n',
                    f'rread={r2}\n',
                    f'bwa mem {options} $refpath $fread $rread | {slurpydir}/tobedpe.py $refpath {outfile} {line_count}\n']
        
        ## If we are not forcing the run, then check if it exists, and the report exists and has alredy been run, just skip
        if (not forced) and (not unfinished(bwa_repo)) and fileexists(outfile):
            print(f'INFO: Detected a finished run ({outfile}) from {bwa_file} in {bwa_repo}, skipping.\n')
        else:
            ## Write the bwa command to file 
            writetofile(bwa_file, sbatch(job_names[i],thread_count,the_cwd,bwa_repo,nice=nice,nodelist=nodes,memory=memory) + bwa_coms, debug)
            ## append the report and files
            bwa_files.append(bwa_file)
            bwa_repos.append(bwa_repo)
            ## Remove previous filtering log if it exists
            filter_report = f'{debugdir}/{pix}.filter.bedpe.{out_file_prefix}.log'
            ## Remove the filter log
            remove(filter_report) if not unfinished(filter_report) else None 
                
    ## Recount the nubmer submitted 
    to_submit = len(bwa_files)
    ## Iniate submitted 
    submitted = 0
    ## While sub miiting 
    while submitted < to_submit:
        ## Try to Submit the command to SLURM
        try:
            job_id = submitsbatch(f'sbatch --partition={partitions} --dependency=singleton {bwa_files[submitted]}')
            ## Add to sub
            submitted += 1
            ## Print job id
            print('Job ID: %s'%job_id)
            ## Wiat a few seconds
            sleep(2)
        except Exception as error:
            print(error)

    ## Check the reports are finished 
    while unfinishedreports(bwa_repos):
        ## Wait a minitue 
        sleep(waittime)

    ## Print to log 
    print("Finished %s bwa submissions for sample: %s"%(to_submit,sample_name))

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## If the script is envoked by name 
if __name__ == "__main__":
    main()
## EOF 