#!/usr/bin/env python
#SBATCH --job-name=bwa.master           ## Name of job
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=1               ## Number of tasks to be launched
#SBATCH --nice=2147483645               ## Nice parameter, sets job to lowest priority 

## Bring in ftns and variables from defaluts 
from defaults import sortglob, sbatch, submitsbatch, fileexists, getfilesize, remove
## Load in directories
from directories import splitsdir, comsdir, debugdir, slurpydir, bedtmpdir, checkerdir
## Load in write to file from pysam tools 
from pysamtools import writetofile, listzip, ifprint
## Load input vars from params
from parameters import refmetavar, bwathreads, lib_default, nice, hic_options, waittime, nparallel
## load in sleep
from time import sleep
## Bring in tile and arange 
from numpy import tile, arange

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
def reportcheck(reportpath) -> bool:
    finished = False 
    if fileexists(reportpath):
        with open(reportpath,'r') as infile:
            for l in infile:
                if 'finished' in l.lower():
                    finished = True
                    break
            infile.close()
    return finished

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

## Ftn for tiling jobs
def vectortile(k,n):
    return tile(arange(k)+1,n)

## Ftn for formating the bwa master 
def bwamaster(sname:str,refpath:str,threads:int,cwd:str,partition:str,debug:bool,nice:int,njobs:int,pix=pix,linecount=line_count,library=None,forced=False,nodelist=None,bwaopts='') -> tuple[list[str], str]:
    ## Format command 
    command = f'{slurpydir}/bwatobedpe.py -s {sname} -r {refpath} -b {threads} -c {cwd} -P {partition} -N {nice} -l {linecount} -j {njobs}' + (f' -L {library}' if library else '') + (' --debug' if debug else '') + (' --force' if forced else '') + (' --nodelist %s'%' '.join(nodelist) if nodelist else '') + (' -B %s'%bwaopts if len(bwaopts) else '')
    ## Format report 
    report  = f'{debugdir}/{pix}.bwa.to.bedpe.{sname}.log'
    ## Return the command and report 
    return [command], report 

## Load in help messages from parameters 
from parameters import ST, s_help, r_help, b_help, P_help, L_help, N_help, debug_help, force_help, node_help, j_help, B_help

## Set description and help messages 
bwadescr = 'A submission script that formats bwa/bedpe commands for paired fastq file from fastp splits of a given sample.'
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
    parser.add_argument("-j", "--njobs",          dest="j",     type=int,  required=False, help = j_help, metavar = 'n',          default = nparallel    )
    parser.add_argument("-B", "--bwa-opts",       dest="bwa",   type=str,  required=False, help = B_help, metavar = hic_options,  default = hic_options  )
    parser.add_argument("--nodelist",             dest="nodes", nargs='+', required=False, help = node_help,                      default = None         )

    ## Add boolean vars 
    parser.add_argument("--debug",                dest="debug",     help = debug_help,    action = ST                                          )
    parser.add_argument("--force",                dest="force",     help = force_help,    action = ST                                          )
    ## Set the paresed values as inputs
    inputs = parser.parse_args() 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    
    ## Set inputs 
    sample_name  = inputs.s         ## Set sample name
    ref_path     = inputs.r         ## Path to reference file
    thread_count = inputs.b         ## Number of threads
    the_cwd      = inputs.c         ## The cwd
    partitions   = inputs.P         ## Paritions to run nodes on 
    library      = inputs.L         ## Library used in hic construcition
    line_count   = inputs.l         ## Number of lines to parse from file
    nice         = inputs.N         ## Nice parameter, "Its nice to be nice"
    nparallel    = inputs.j         ## Number of parallel bwa jobs to run 
    nodes        = inputs.nodes     ## Node list 
    debug        = inputs.debug     ## To debug, or not debug
    forced       = inputs.force     ## Use the force Luke, let go Luke. Luke trust me.
    bwa_opts     = inputs.bwa       ## String of options to also feed into bwa 

    ## Format the options for bwa mem
    options = f'-v 1 -t {thread_count} ' + ' '.join(bwa_opts.split(','))

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
    
    ## format list of checks
    bwa_checkers = [f'{checkerdir}/{i}.{sample_name}.bwa.log' for i in range(nreads)]
    ## Remove the previous checks if any
    [remove(bwa_check) for bwa_check in bwa_checkers if fileexists(bwa_check)]

    ## Iniate list of bwa files
    bwa_files   = []

     ## gather job numbers / names
    job_numbers = vectortile(nparallel,nreads)
    job_names   = [f'{job_numbers[i]}.bwa.sh' for i in range(nreads)]

    ## Iterate thru the read pairs 
    for i, (r1,r2) in enumerate(read_pairs):
        ## SEt report and file name 
        bwa_repo  = f'{debugdir}/{pix}.bwa.{i}.{sample_name}.log'
        bwa_file  = f'{comsdir}/{pix}.bwa.{i}.{sample_name}.sh' 
        outfile   = f'{bedtmpdir}/{i}.{sample_name}.bedpe'
        bwa_check = bwa_checkers[i]

        ## format the command 
        bwa_coms = [f'refpath={ref_path}\n',
                    f'bwa mem {options} $refpath {r1} {r2} | {slurpydir}/tobedpe.py $refpath {library} {outfile} {line_count}\n',
                    f'{slurpydir}/myecho.py Finished bwa alignment of split {i} {bwa_check}\n## EOF']

        ## If we are not forcing the run, then check if it exists
        if not forced:
            ## If the report exists and has alredy been run, just skip
            if fileexists(outfile) and fileexists(bwa_repo) and reportcheck(bwa_repo):
                print(f'WARNING: Detected a finished run ({outfile}) from {bwa_file} in {bwa_repo}.\nINFO: Skipping.\n')
                ## REformat coms so its just the check 
                bwa_coms = [bwa_coms[-1]]
        
        ## Write the bwa command to file 
        writetofile(bwa_file, sbatch(job_names[i],thread_count,the_cwd,bwa_repo,nice=nice,nodelist=nodes) + bwa_coms, debug)
        ## append the report and files
        bwa_files.append(bwa_file)

    ## Iniate submitted 
    submitted = 0
    ## While sub miiting 
    while submitted < nreads:
        ## Try to Submit the command to SLURM
        try:
            submitsbatch(f'sbatch --partition={partitions} --dependency=singleton {bwa_files[submitted]}')
            ## Add to sub
            submitted += 1
            ## Wiat a few seconds
            sleep(waittime)
        except Exception as error:
            print(error)

    ## While thekicker is true 
    while not (sum([fileexists(f) for f in bwa_checkers]) == nreads):
        ## Wait a minitue 
        sleep(2*waittime)

    ## Print to log 
    print("Finished %s bwa submissions for sample: %s"%(nreads,sample_name))
## EOF 