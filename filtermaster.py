#!/usr/bin/env python
#SBATCH --job-name=filtermaster         ## Name of job
#SBATCH --output=%x.%j.out              ## Name stdout
#SBATCH --error=%x.%j.err               ## Name stderr
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=1               ## Number of tasks to be launched
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## Bring in ftns and variables from defaluts 
from defaults import sortglob, sbatch, submitsbatch
## Load in variables from parameters
from parameters import map_q_thres, error_dist, daskthreads, nice, chunksize, waittime, nparallel, lib_default
## Load in directories from parameters
from parameters import comsdir, debugdir, bedtmpdir, slurpydir
## Load in help messages
from parameters import ST, Q_help, L_help, E_help, r_help, X_help, t_help, N_help, Z_help, m_help, P_help, force_help, node_help, dove_help, intra_help, j_help, hicex_help, slurmem_help
## load in sleep
from time import sleep
## Load in report check 
from bwatobedpe import vectortile, unfinishedreports, hic_flag
## Bring in write to file
from myecho import writetofile
## Bring in finished check
from checkwork import unfinished
## Load in argparse
import argparse

## Set stage in piepline
pix = 2

## Ftn for formating / joinign a list 
def formatinput(inlist) -> str: return ' '.join([str(x) for x in inlist])

## Ftn for formating commands to this script 
def filtermaster(sname:str,refpath:str,cwd:str,xcludes:list,includes:list,mapq:int,errordistance:int,threads:int,library:list,partitions:str,debug:bool,nice:int,njobs:int,pix=pix,forced=False,chunksize=chunksize,maxdist=0,nodelist=None,dovetail=False,removeinter=False,hicexplorer=False,memory=None) -> tuple[list[str], str]:
    command = f'{slurpydir}/filtermaster.py -s {sname} -r {refpath} -c {cwd} -q {mapq} -e {errordistance} -t {threads} -N {nice} -Z {chunksize} -M {maxdist} -x {formatinput(xcludes)} -i {formatinput(includes)} -l {formatinput(library)} -P {partitions} -j {njobs}'\
            + (' --dedovetail' if dovetail else '') + (' --debug' if debug else '') + (' --force' if forced else '') \
            + (' --intra-only' if removeinter else '')  + (' --nodelist %s'%' '.join(nodelist) if nodelist else '') \
            + (' --hicexplorer' if hicexplorer else '') + (' --memory %s'%memory if memory else '') 
    ## Format the report 
    report  = f'{debugdir}/{pix}.filter.master.{sname}.log'
    ## Return the 
    return [command+'\n'], report 

## Set description of sub script and help messages 
filtdescr  = 'The submission script of filterbedpe across sample paritions'
s_help     = 'Sample starting name to form wild card extraction of paths'
I_help     = "List of chormosomes/contigs to only include in analysis"
c_help     = 'The current working directory'

def parse_args():
    ## Make the parse
    parser = argparse.ArgumentParser(description = filtdescr)
    ## Add the required arguments
    parser.add_argument("-s", dest="S", type=str,  required=True,   help=s_help, metavar='./path/to/input.bedpe'           )   
    parser.add_argument("-r", dest="R", type=str,  required=True,   help=r_help, metavar='./path/to/ref.fasta'             )
    parser.add_argument("-c", dest="C", type=str,  required=True,   help=c_help, metavar = './the/cwd'                     )
    ## Add optional args
    parser.add_argument("-x", dest="X", nargs='+', required=False,  help=X_help, metavar='chrM',      default=['chrM']     )
    parser.add_argument("-i", dest="i", nargs='+', required=False,  help=I_help, metavar='chr1 chr2', default=[]           )
    parser.add_argument("-q", dest="Q", type=int,  required=False,  help=Q_help, metavar='n',         default=map_q_thres  )
    parser.add_argument("-e", dest="E", type=int,  required=False,  help=E_help, metavar='n',         default=error_dist   )
    parser.add_argument("-l", dest="L", nargs='+', required=False,  help=L_help, metavar='arima',     default=lib_default  )
    parser.add_argument("-P", dest="P", type=str,  required=False,  help=P_help, metavar='tb',        default = 'tb'       ) 
    parser.add_argument("-t", dest="T", type=int,  required=False,  help=t_help, metavar=daskthreads, default=daskthreads  )
    parser.add_argument("-N", dest="N", type=int,  required=False,  help=N_help, metavar = 'n',       default = nice       )
    parser.add_argument("-Z", dest="Z", type=int,  required=False,  help=Z_help, metavar = 'n',       default=chunksize    )
    parser.add_argument("-M", dest="M", type=int,  required=False,  help=m_help, metavar = 'n',       default=0            )
    parser.add_argument("-j", dest="j", type=int,  required=False,  help=j_help, metavar = 'n',       default=nparallel    )
    parser.add_argument("--nodelist",   dest="nodes", nargs='+', required=False, help = node_help,    default = None       )
    parser.add_argument("--memory",     dest="mem",   type=str,  required=False, help = slurmem_help, metavar = '40G',  default = None)
    ## Add boolean 
    parser.add_argument("--dedovetail",     dest="tails",  help = dove_help,    action = ST )
    parser.add_argument("--debug",          dest="debug",  help = dove_help,    action = ST )
    parser.add_argument("--force",          dest="force",  help = force_help,   action = ST )
    parser.add_argument("--hic",            dest="hic",    help = hic_flag,     action = ST )
    parser.add_argument("--intra-only",     dest="Intra",  help = intra_help,   action = ST )
    parser.add_argument("--hicexplorer",    dest="hicexp", help = hicex_help,   action = ST )
    ## Set the paresed values as inputs
    return parser.parse_args()

## define main ftn 
def main():
    ## Set the paresed values as inputs
    inputs = parse_args()

    ## Set iinputs
    sample_name = inputs.S          ## Sample name 
    ref_path    = inputs.R          ## Path to reference genome fasta file
    the_cwd     = inputs.C          ## The current working directory
    xcludos     = inputs.X          ## Chromosomes to exclude from analysis 
    includos    = inputs.i          ## Chromosomes to include in analysis
    map_q_thres = inputs.Q          ## Mapping quality threshold
    error_dist  = inputs.E          ## Error distance of Hi-C contacts
    elibrary    = inputs.L          ## The enzamatic library to run on 
    partitions  = inputs.P          ## Cluster paritions to run jobs on 
    threads     = inputs.T          ## The number of threads, thread them all!
    nice        = inputs.N          ## How nice are we?
    chunksize   = inputs.Z          ## Chunksize 
    max_dist    = inputs.M          ## Max linear distance of read pairs
    nparallel   = inputs.j          ## Set number of parallel jobs 
    nodes       = inputs.nodes      ## Set the nodes to run on 
    dovetail    = inputs.tails      ## Flag to remove dovetail reads
    debug       = inputs.debug      ## Flag to debug 
    forced      = inputs.force      ## Flag to force 
    intra_only  = inputs.Intra      ## Flag for intra chromosomal contacts only 
    hicexplorer = inputs.hicexp     ## Flag to run filtering in hicexplorer mode
    memory      = inputs.mem        ## SLURM memory limit 
    
    ## Bring in bedpe paths, calc len
    bedpe_paths = sortglob(f'{bedtmpdir}/*.{sample_name}.bedpe')
    nbedpe = len(bedpe_paths)
    ## Check work
    assert nbedpe, "ERROR: Unable to find bedpe files associated with sample: %s"%sample_name

    ## Initiate filter files (.sh) and reports (.txt)
    filter_files    = []
    filter_repos    = []

    ## gather job numbers / names
    job_numbers = vectortile(nparallel,nbedpe)
    job_names   = [f'{job_numbers[i]}.filter.sh' for i in range(nbedpe)]

    ## Iterate thru the paths
    for i,bedpe in enumerate(bedpe_paths):
        ## Format the report, file and check for this filt of bedpe split 
        filter_repo  = f'{debugdir}/{pix}.filter.bedpe.{i}.{sample_name}.log'
        filter_file  = f'{comsdir}/{pix}.filter.bedpe.{i}.{sample_name}.sh' 

        ## Format commands 
        filter_coms   = [f'{slurpydir}/filtering.py -b {bedpe} -e {error_dist} -l {" ".join(elibrary)} -q {map_q_thres} -r {ref_path} -x {formatinput(xcludos)} -i {formatinput(includos)} -Z {chunksize} -M {max_dist}' + (' --dedovetail' if dovetail else ' ') + (' --intra-only' if intra_only else '') +  (' --hicexplorer' if hicexplorer else '') + '\n']

        ## If we are not forcing the run, then check if it exists
        if (not forced) and (not unfinished(filter_repo)):
            print(f'INFO: Detected a finished run from {filter_file} in {filter_repo}, skipping.\n',flush=True)
        else:
            ## Write the bwa command to file 
            writetofile(filter_file, sbatch(job_names[i],threads,the_cwd,filter_repo,nice=nice,nodelist=nodes,memory=memory) + filter_coms, debug)
            ## append the report and files 
            filter_files.append(filter_file)
            filter_repos.append(filter_repo)

    ## Recount the nubmer submitted 
    to_submit = len(filter_files)
    ## Iniate submitted 
    submitted = 0
    ## While sub miiting 
    while submitted < to_submit:
        ## Try to Submit the command to SLURM
        try:
            job_id = submitsbatch(f'sbatch --partition={partitions} --dependency=singleton {filter_files[submitted]}')
            ## Add to sub
            submitted += 1
            ## Pring job id
            print('Job ID: %s'%job_id,flush=True)
            ## Wiat a few seconds
            sleep(3)
        except Exception as error:
            print(error,flush=True)

    ## Check the reports are finished 
    while unfinishedreports(filter_repos):
        ## Wait a minitue 
        sleep(waittime)

    ## Print to log 
    print("Finished %s bwa submissions for sample: %s"%(to_submit,sample_name),flush=True)
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## if the script is envoked
if __name__ == "__main__":
    main()
## EOF 