#!/usr/bin/env python
#SBATCH --job-name=filter.master        ## Name of job
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=1               ## Number of tasks to be launched
#SBATCH --nice=2147483645               ## Nice parameter, sets job to lowest priority 

## Bring in ftns and variables from defaluts 
from defaults import sortglob, sbatch, submitsbatch, fileexists, remove
## Load in directories
from directories import comsdir, debugdir, bedtmpdir, slurpydir, checkerdir
## Load in write to file from pysam tools 
from pysamtools import writetofile
## Load in params
from parameters import map_q_thres, error_dist, daskthreads, nice, chunksize, waittime, nparallel
## load in sleep
from time import sleep
## Load in report check 
from bwatobedpe import reportcheck, vectortile, hic_flag

## Set stage in piepline
pix = 2

## Ftn for formating / joinign a list 
def formatinput(inlist) -> str:
    return ' '.join([str(x) for x in inlist])

## Ftn for formating commands to this script 
def filtermaster(sname:str,refpath:str,cwd:str,xcludes:list,includes:list,mapq:int,errordistance:int,threads:int,library:str,partitions:str,debug:bool,nice:int,njobs:int,pix=pix,forced=False,chunksize=chunksize,maxdist=0,nodelist=None,dovetail=False,removeinter=False,hicexplorer=False) -> tuple[list[str], str]:
    command = f'{slurpydir}/filtermaster.py -s {sname} -r {refpath} -c {cwd} -q {mapq} -e {errordistance} -t {threads} -N {nice} -Z {chunksize} -M {maxdist} -x {formatinput(xcludes)} -i {formatinput(includes)} -l {library} -P {partitions} -j {njobs}' + (' --dedovetail' if dovetail else '') + (' --debug' if debug else '') + (' --force' if forced else '') + (' --intra-only' if removeinter else '')  + (' --nodelist %s'%' '.join(nodelist) if nodelist else '') + (' --hicexplorer' if hicexplorer else '')
    report  = f'{debugdir}/{pix}.filter.master.{sname}.log'
    return [command+'\n'], report 

## Load in help messages
from parameters import ST, Q_help, L_help, E_help, r_help, X_help, t_help, N_help, Z_help, m_help, P_help, force_help, node_help, dove_help, intra_help, j_help, hicex_help

## Set description of sub script and help messages 
filtdescr  = 'The submission script of filterbedpe across sample paritions'
s_help     = 'Sample starting name to form wild card extraction of paths'
I_help     = "List of chormosomes/contigs to only include in analysis"
c_help     = 'The current working directory'

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
    parser.add_argument("-i", dest="i", nargs='+', required=False,  help=I_help, metavar='chr1 chr2', default=[]           )
    parser.add_argument("-q", dest="Q", type=int,  required=False,  help=Q_help, metavar='n',         default=map_q_thres  )
    parser.add_argument("-e", dest="E", type=int,  required=False,  help=E_help, metavar='n',         default=error_dist   )
    parser.add_argument("-l", dest="L", type=str,  required=False,  help=L_help, metavar='Arima',     default='Arima'      )
    parser.add_argument("-P", dest="P", type=str,  required=False,  help=P_help, metavar='tb',        default = 'tb'       ) 
    parser.add_argument("-t", dest="T", type=int,  required=False,  help=t_help, metavar=daskthreads, default=daskthreads  )
    parser.add_argument("-N", dest="N", type=int,  required=False,  help=N_help, metavar = 'n',       default = nice       )
    parser.add_argument("-Z", dest="Z", type=int,  required=False,  help=Z_help, metavar = 'n',       default=chunksize    )
    parser.add_argument("-M", dest="M", type=int,  required=False,  help=m_help, metavar = 'n',       default=0            )
    parser.add_argument("-j", dest="j", type=int,  required=False,  help=j_help, metavar = 'n',       default=nparallel    )
    parser.add_argument("--nodelist",   dest="nodes", nargs='+', required=False, help = node_help,    default = None       )

    ## Add boolean 
    parser.add_argument("--dedovetail",     dest="tails",  help = dove_help,    action = ST )
    parser.add_argument("--debug",          dest="debug",  help = dove_help,    action = ST )
    parser.add_argument("--force",          dest="force",  help = force_help,   action = ST )
    parser.add_argument("--hic",            dest="hic",    help = hic_flag,     action = ST )
    parser.add_argument("--intra-only",     dest="Intra",  help = intra_help,   action = ST )
    parser.add_argument("--hicexplorer",    dest="hicexp", help = hicex_help,   action = ST )

    ## Set the paresed values as inputs
    inputs = parser.parse_args()

    ## Set iinputs
    sample_name = inputs.S
    ref_path    = inputs.R
    the_cwd     = inputs.C
    xcludos     = inputs.X
    includos    = inputs.i
    map_q_thres = inputs.Q
    error_dist  = inputs.E
    elibrary    = inputs.L
    partitions  = inputs.P
    threads     = inputs.T
    nice        = inputs.N
    chunksize   = inputs.Z 
    max_dist    = inputs.M 
    nparallel   = inputs.j 
    nodes       = inputs.nodes
    dovetail    = inputs.tails      ## Flag to remove dovetail reads
    debug       = inputs.debug      ## Flag to debug 
    forced      = inputs.force      ## Flag to force 
    intra_only  = inputs.Intra      ## Flag for intra chromosomal contacts only 
    hicexplorer = inputs.hicexp     ## Flag to run filtering in hicexplorer mode

    ## Bring in bedpe paths, calc len
    bedpe_paths = sortglob(f'{bedtmpdir}/*.{sample_name}.bedpe')
    nbedpe = len(bedpe_paths)
    ## Check work
    assert nbedpe, "ERROR: Unable to find bedpe files associated with sample: %s"%sample_name

    ## Format bed pe checks
    bed_checkers = [f'{checkerdir}/{i}.{sample_name}.bedpe.log' for i in range(nbedpe)]
    ## Remove the previous checks if any
    [remove(bed_check) for bed_check in bed_checkers if fileexists(bed_check)]

    ## Initiate filter files (.sh)
    filter_files    = []

    ## gather job numbers / names
    job_numbers = vectortile(nparallel,nbedpe)
    job_names   = [f'{job_numbers[i]}.filter.sh' for i in range(nbedpe)]

    ## Iterate thru the paths
    for i,bedpe in enumerate(bedpe_paths):
        ## Format the report, file and check for this filt of bedpe split 
        filter_repo  = f'{debugdir}/{pix}.filter.bedpe.{i}.{sample_name}.log'
        filter_file  = f'{comsdir}/{pix}.filter.bedpe.{i}.{sample_name}.sh' 
        bed_check    = bed_checkers[i]

        ## Format commands 
        filter_coms   = [f'{slurpydir}/filtering.py -b {bedpe} -e {error_dist} -l {elibrary} -q {map_q_thres} -r {ref_path} -x {formatinput(xcludos)} -i {formatinput(includos)} -Z {chunksize} -M {max_dist}' + (' --dedovetail' if dovetail else ' ') + (' --intra-only' if intra_only else '') +  (' --hicexplorer' if hicexplorer else '') + '\n',
                         f'{slurpydir}/myecho.py Finished bedpe filtering of split {i} {bed_check}\n## EOF']

        ## If we are not forcing the run, then check if it exists
        if not forced:
            ## If the report exists and has alredy been run, just skip
            if fileexists(filter_repo) and fileexists(filter_file) and reportcheck(filter_repo):
                print(f'WARNING: Detected a finished run from {filter_file} in {filter_repo}.\nINFO: Skipping.\n')
                ## Reformat the commands
                filter_coms = [filter_coms[-1]]

        ## Write the bwa command to file 
        writetofile(filter_file, sbatch(job_names[i],threads,the_cwd,filter_repo,nice=nice,nodelist=nodes) + filter_coms, debug)
        ## append the report and files 
        filter_files.append(filter_file)

    ## Iniate submitted 
    submitted = 0
    ## While sub miiting 
    while submitted < nbedpe:
        ## Try to Submit the command to SLURM
        try:
            submitsbatch(f'sbatch --partition={partitions} --dependency=singleton {filter_files[submitted]}')
            ## Add to sub
            submitted += 1
            ## Wiat a few seconds
            sleep(waittime)
        except Exception as error:
            print(error)

    ## While thekicker is true 
    while not (sum([fileexists(f) for f in bed_checkers]) == nbedpe):
        ## Wait a minitue 
        sleep(2*waittime)

    ## Print to log 
    print("Finished %s bwa submissions for sample: %s"%(nbedpe,sample_name))
## EOF 