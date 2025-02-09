#!/usr/bin/env python
#SBATCH --job-name=bwa.master           ## Name of job
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=1               ## Number of tasks to be launched
#SBATCH --nice=2147483645               ## Nice parameter, sets job to lowest priority 

## Bring in ftns and variables from defaluts 
from defaults import sortglob, sbatch, submitsbatch, comsdir, debugdir, bedtmpdir, slurpydir
## Load in params
from parameters import Q_help, map_q_thres, error_dist, L_help, E_help, r_help, X_help, t_help, N_help, daskthreads, part, P_help, nice
## Load in write to file from pysam tools 
from pysamtools import writetofile
## load in sleep
from time import sleep
## Load in report check 
from bwamaster import reportcheck

## Set stage in piepline
pix = 2

## Ftn for formating / joinign a list 
def formatinput(inlist):
    return ' '.join([str(x) for x in inlist])

## Ftn for formating commands to this script 
def filtermaster(sname:str,refpath:str,cwd:str,xcludes:list,includes:list,mapq:int,errordistance:int,threads:int,library:str,partitions:str,todovetail:bool,debug:bool,pix=2):
    command = f'{slurpydir}/filtersubs.py -s {sname} -r {refpath} -c {cwd} -q {mapq} -e {errordistance} -t {threads} -x {formatinput(xcludes)} -i {formatinput(includes)} -l {library} -P {partitions}' + (' --dovetails' if todovetail else ' ') + (' --debug' if debug else ' ')
    report  = f'{debugdir}/{pix}.filter.master.{sname}.log'
    return [command+'\n'], report 

## Set description of sub script and help messages 
filtdescr  = 'The submission script of filterbedpe across sample paritions'
s_help     = 'Sample starting name to form wild card extraction of paths'
I_help     = "List of chormosomes/contigs to only include in analysis"
c_help     = 'The current working directory'
dove_help  = "Boolean flag to remove dovetailed paired-end reads (paired reads with overlapping mapped coordiantes) from analsyis (Default: is to remove these)."

## -------------------------------------- MAIN EXECUTABLE -------------------------------------------------- ##
## if the script is envoked
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse
    ## Make the parse
    parser = argparse.ArgumentParser(description = filtdescr)

    ## Add the required arguments
    parser.add_argument("-s", dest="S", type=str,  required=True,   help=s_help, metavar='./path/to/input.bedpe'           )   
    parser.add_argument("-r", dest="R", type=str,  required=True,   help=r_help, metavar='./path/to/ref.fasta'             )
    parser.add_argument("-c", dest="C", type=str,  required=True,   help=c_help, metavar = './the/cwd'                     )

    ## Add optional args
    parser.add_argument("-x", dest="X", nargs='+', required=False,  help=X_help, metavar='chrM',      default=['chrM']     )
    parser.add_argument("-i", dest="I", nargs='+', required=False,  help=I_help, metavar='chr1 chr2', default=[]           )
    parser.add_argument("-q", dest="Q", type=int,  required=False,  help=Q_help, metavar='n',         default=map_q_thres  )
    parser.add_argument("-e", dest="E", type=int,  required=False,  help=E_help, metavar='n',         default=error_dist   )
    parser.add_argument("-l", dest="L", type=str,  required=False,  help=L_help, metavar='Arima',     default='Arima'      )
    parser.add_argument("-P", dest="P", type=str,  required=False,  help=P_help, metavar=part,        default = part       ) 
    parser.add_argument("-t", dest="T", type=int,  required=False,  help=t_help, metavar=daskthreads, default=daskthreads  )
    parser.add_argument("-N", dest="N", type=int,  required=False,  help=N_help, metavar = 'n',       default = nice       )

    ## Add boolean 
    parser.add_argument("--dovetails",  dest="tails",  help = dove_help,    action = 'store_true')
    parser.add_argument("--debug",      dest="debug",  help = dove_help,    action = 'store_true')

    ## Set the paresed values as inputs
    inputs = parser.parse_args()

    ## Set iinputs
    sample_name = inputs.S
    ref_path    = inputs.R
    the_cwd     = inputs.C
    xcludos     = inputs.X
    includos    = inputs.I
    map_q_thres = inputs.Q
    error_dist  = inputs.E
    elibrary    = inputs.L
    partitions  = inputs.P
    threads     = inputs.T
    nice        = inputs.N
    dovetail    = inputs.tails  ## Flag to remove dovetail reads
    debug       = inputs.debug  ## Flag to debug 

    ## Bring in bedpe paths
    bedpe_paths = sortglob(f'{bedtmpdir}/*.{sample_name}.bedpe')
    ## Check work
    assert len(bedpe_paths), "ERROR: Unable to find bedpe files associated with sample: %s"%sample_name

    ## INiate filter repos
    filter_reports = []
    ## Iterate thru the paths
    for i,bedpe in enumerate(bedpe_paths):
        ## format the command 
        filter_com   = f'{slurpydir}/filterbedpe.py -b {bedpe} -e {error_dist} -l {elibrary} -q {map_q_thres} -r {ref_path} -x {formatinput(xcludos)} -i {formatinput(includos)}' + (' --dovetails' if dovetail else ' ')
        filter_repo  = f'{debugdir}/{pix}.filter.bedpe.{i}.{sample_name}.log'
        filter_file  = f'{comsdir}/{pix}.filter.bedpe.{i}.{sample_name}.sh' 

        ## Write the bwa command to file 
        writetofile(filter_file, sbatch(None,threads,the_cwd,filter_repo,nice=nice) + [filter_com+'\n'], debug)

        ## Submit the command to SLURM
        submitsbatch(f'sbatch --partition={partitions} {filter_file}')

        ## append the report
        filter_reports.append(filter_repo)
        ## Sleep a second
        sleep(1)

    ## Check the reports
    kicker = True 

    ## While thekicker is true 
    while kicker:
        kicker = not reportcheck(filter_reports)
        ## Wait a minitue 
        sleep(10)

    ## Print to log 
    print("Finished %s bwa submissions for sample: %s"%(len(filter_reports),sample_name))

