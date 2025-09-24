## Bring inmods
import argparse
## Set the run local var
runlocal = False 
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      SLUR(M)-py: A SLURM Powered Pythonic Pipeline for Parallel Processing of 3D (Epi)genomic Profiles
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
## List command for canceling all
"""
squeue -u croth | grep 'croth' | awk '{print $1}' | xargs -n 1 scancel
squeue -u croth | grep 'croth' | grep tb | awk '{print $1}' | xargs -n 1 scancel
squeue -u croth | grep 'croth' | grep gpu | grep "(DependencyNeverSatisfied)" | awk '{print $1}' | xargs -n 1 scancel
"""
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      VARIABLE SETTING
## Fake job id used for debugging only 
fakejobid = 666666         
## Path to human t2t ref on canopus and my local machine 
if runlocal:
    t2t_refpath = '/Users/croth/Desktop/T2T/GCA_009914755.4_T2T-CHM13v2.0_genomic.named.fasta'
    t2t_gffpath = '/Users/croth/Desktop/T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.named.gff'
    t2t_gtfpath = '/Users/croth/Desktop/T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.named.gtf'
else:
    t2t_refpath = '/panfs/biopan04/4DGENOMESEQ/REFERENCES/T2T/GCA_009914755.4_T2T-CHM13v2.0_genomic.named.fasta' 
    t2t_gffpath = '/panfs/biopan04/4DGENOMESEQ/REFERENCES/T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.named.gff'
    t2t_gtfpath = '/panfs/biopan04/4DGENOMESEQ/REFERENCES/T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.named.gtf'

## Set local path to vero gm ref
vero_refpath = '/panfs/biopan04/4DGENOMESEQ/REFERENCES/GreenMVA/Chlorocebus_sabeus_mva.fasta'

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## SET thread cound 
threads = 8
##      Hi-C and ATAC-seq DEFAULT VARIABLE SETTING  
splitsize    = 10**7         ##     The number of splits made by fastp 
bwathreads   = threads       ##     Number of threads used by calls to bwa 
samthreads   = threads       ##     Number of threads used by calls to samtools 
daskthreads  = threads       ##     Number of threads used by calls to dask df 
fastpthreads = threads       ##     Number of threads in fastp 
parts        = ['tb']        ##     Defalut partition 
map_q_thres  = 30            ##     Minimum mapping quality threhosld 
waittime     = 1             ##     Seconds of buffer time 
error_dist   = 25000         ##     The idstance to check for erros 
lib_default  = 'Arima'       ##     Defalut library used to make Hi-C experimetns 
hic_options  = ',-5SMP'      ##     Hi-C options passed to bwa 
chunksize    = 950000        ##     Chunks size for parsing with pandas
hicsep       = ' '           ##     Text deliminator 
line_count   = 10**7         ##     Number of lines 
mito         = 'chrM'        ##     The name of the mitocondrial contig (in humns)
xmemory      = 49152         ##     Sets the memory used by juicer pre command 
nice         = 10**8         ##     Set the nice parameter 
nparallel    = 24            ##     Number of jobs to run in parallele 
max_dist     = 0             ##     Maximum distance allowed for paired end mapping
max_nchrom   = 200           ##     Maximum number of chromosomes allowed to parse
binsizes     = [2500000,     ##     Set the binsizes of resolution for Hi-C analysis 
                2000000,
                1000000,
                 750000,
                 500000,
                 250000,
                 200000,
                 100000,
                  75000,
                  50000,
                  25000,
                  10000,
                  5000]

## Set store tru var
ST = 'store_true'
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      METAVARS
## Set metavars
c_metavar  = './path/to/control.bam'
g_metavar  = 'bp'
refmetavar = './path/to/reference.fasta'  
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

## Set pipeline of hic analysis ## NOTE defined after defaults import 
hic_pipeline  = ['fastp', 'bwa', 'filter','dedup','concat','gxg','toshort','hic' ,'macs3','sam','count','clean']
##                  0        1      2        3       4      5a      5b      5c      5d      5e     6B      6C

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      HELP MESSAGES 
## Set help messages for input and default variables 
r_help = "Path to input reference referecne (in .fasta or .fa format) with an assoicated (.fai) bwa index to use for alignment." 
fqhelp = "Path to the directory holding zipped, paired-end seqeuncing data in fastq.gz format."
F_help = "The approximate number of reads per split made by fastp on input fastq files. Default is: %s."%splitsize
T_help = "The number of threads used across all applications of the run (fastp, bwa, dask.dataframes). Default is: %s."%threads 
f_help = "The number of threads used in fastp to split input fastq files. Default is: %s. Note: must be an even multiple of the number of splits."%fastpthreads
b_help = "The number of threads used per bwa alignment on split input fastq files. Default is: %s."%bwathreads
B_help = "A comma starting and seperated list (no spaces) of options (and their values) for the bwa mem algorithm (for example ,-t,2,-k,10,-y,5,-S). See bwa mem for help and a list of options."
n_help = "Run name used to name output files. Default behavior is to take the common name from the input read pairs."
M_help = "Name of the mitochondrial contig. Default is: %s."%mito
X_help = "List of chromosomes/contigs to exclude from analysis. Default behavior is to process all within the passed .fasta or .fa file."
P_help = "The type of partition jobs formatted by slurpy run on. Default is %s."%parts[0]
Q_help = "Mapping quality threshold to filter alignments. Default is: %s."%map_q_thres
c_help = "Path to control or input files used in ChIP-seq experiments. Must be a .bam file or a .bedpe file or processed via SLUR(M)-py."
G_help = "Path to list of chromosomes (by name) to include in final analysis. Default behavior expects a tab seperated tsv or bed, comma seperated csv, or space seperated txt file with no header."
g_help = "Size of the genome being analyzed, used as parameter for macs3. Inputs can be integers in bp or two letter short hand, for e.g. hs for homo sapiens. Default behavior is to calculate this value from the reference file."
E_help = "Minimum fragment size of read pairs scanned for an intersecting restriction fragment site (if passed thru library parameter). Default is %s. These pairs are also marked for dangling ends and self-circles."%error_dist
L_help = "The name of the restriction site enzyme (or library prep) used in Hi-C sample creation. Default is %s. Options include Arima, MboI, DpnII, Sau3AI, and HindIII. Note: passing none (i.e. Dovetail) is also allowed, but checks for restriction sites and dangling ends will be skipped."%lib_default
Z_help = "Number of rows loaded into pandas at a time. Default is: %s. WARNING: while increasing could speed up pipeline it could also cause memeory issues."%chunksize
t_help = "The number of threads used in calls to functions and calculations with pandas and dask dataframes. Default is: %s."%daskthreads
s_help = "The number of threads used in calls to samtools. Default is: %s"%samthreads
J_help = "Path to a juicer jar file with the juicer pre command. Required for .hic file creation."
S_help = "Space seperated list of chromosome resolutions (i.e. bin sizes) for .hic files. Default: %s"%' '.join(map(str,binsizes))
x_help = "Amount of Xmx and Xms memory passed to juicer\'s pre command. Default is: %s."%xmemory
A_help = "The path to a gff or bed file with a feature space (i.e. genes) to count gene x gene interactions. Must have columns named Chrom, Left, and Right specifying genomic coordiantes."
a_help = "A SLURM job ID, used as a dependency, specifying all jobs in this run to start after succssful termination."
N_help = "The SLURM nice parameter, an integer value lowering the job priority of submissions. Default is: %s"%nice
R_help = "Step within the pipeline to re-run from. Options include: %s"%', '.join(hic_pipeline) 
m_help = "Maximum allowed distance between intra-chromosomal pairs. Default is zero, setting will activate filter."
j_help = "Number of bwa and filtering jobs allowed to run in parallel. Default: %s"%nparallel

## Specific help messages
node_help  = "Space seperated list of nodes to run jobs on."
dove_help  = "Boolean flag to remove dovetailed paired-end reads (paired reads with overlapping mapped coordiantes) from analsyis (Default: is not to remove these)."
atac_help  = "Preset mode to run in ATAC-seq analysis mode."
rnas_help  = "Preset mode to run in RNA-seq analysis mode."
intra_help = "Boolean flag to remove read pairs spanning multiple chromosomes."
mcool_help = "Flag to make an mcool file with cooler."
pairs_help = "Convert final output to pairs format defined by the 4DNucleome consortium."
inter_help = "Flag to return only inter-chromosomal contacts"
hicex_help = "Flag to run stricter intra-fragment filtering."
save_help  = 'Flag to not save out duplicate read pairs to file.'
maxnc_help = "Total number of allowed chromosomes to parse and process from fasta file. Must be lower than %s to avoid over submitting jobs to SLURM."%max_nchrom

## MACS3 specific variables
shift_size   = 75
extendsize   = 150
shift_help   = "Size (in bp) to shift (3' to 5') end of a read for peak analysis in MACS3 (example 75 bp). Only used when --nomodel is passed. Default: %s"%shift_size
extend_help  = "Size (in bp) to extend (5' to 3') read position for peak analysis in MACS3 (exmale 150 bp). Only used when --nomodel is passed. Default: %s"%extendsize
macs_help    = "Mode and file type for peak calling with MACS3. Options include BED or BEDPE. Default: BEDPE"
summits_help = "Boolean flag to all summits with MACS3. Default: False"
nomodel_help = "Boolean flag to turn off shifting model in MACS3. Default: False"
lambda_help  = "Boolean flag to turn off local lambda calculation in MACS3. Default: False"
mgap_help    = "Max gap between peaks called in MACS3."
mlen_help    = "Minimum length of peaks called in MACS3."
slurmem_help = "The max amount of memory (for e.g. 4OG) passed to SLURM sbatch processes. If set, this value is applied acorss all subprocesses. Default: None"
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      BOOLEAN HELP MESSAGES
## Set help messages for bollean vars
restart_help  = "Flag to force the pipeline to reset and run from start."
runlocal_help = "Disables sbatch submission and submits the script via bash to a local os."
force_help    = "Flag to force the overwrite of output files generated from bwa."
debug_help    = "A flag to run in verbose mode, printing sbatch commands. Default behavior is false."
mark_help     = "Pass this flag to skip marking and removing duplicates. Default behavior is false (conduct duplicate marking)."
broad_help    = "Flag to call broad peaks using the --broad-cutoff=0.1 setting in macs3. See macs3 callpeak --help for more details."
clean_help    = "If included will run clean up script at end of run. The default behavior is false, can be run after pipeline."
count_help    = "Boolean flag to turn off final diagnostics on Hi-C and ATAC-seq samples. Default is to perform the diagnostics"
skipq_help    = "Flag to skip initial quality control and filtering with fastp (i.e. only split reads)."
merge_help    = "Passing this flag will keep replicates/samples seperate acorss (n) pairs of input fastqs, generating (n) outputs rather than one final output."
peaks_help    = "A boolean flag to skip peak calling via macs3."
short_help    = "A boolean flag to convert output bedpe file to short format for hic creation with juicer pre. Defaults to True if the juicer jarpath (-J) is specificed."
tosam_help    = "Flag to convert output .bedpe file from SLUR(M)-py to .sam format via samtools."
tobam_help    = "Flag to convert output .bedpe file from SLUR(M)-py to .bam format via samtools."
wgs_help      = "Passing this flag will run SLUR(M)-py in whole-genome sequencing (wgs) mode, parsing paired-end reads as if from wgs experiments."
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      DIRECTORY NAMES 
## Set directory names 
debugdir    = 'logs'                 ##       Hold logs for debug 
fastqdir    = 'fastqs'               ##       The directory holding fastq files 
aligndir    = 'aligned'              ##       Holds final aligments 
splitsdir   = 'splits'               ##       Temporary dir for split fastq 
comsdir     = 'commands'             ##       Folder for holding all command files 
macs3dir    = 'peaks'                ##       Has results from macs3 
hicdir      = 'merged'               ##       Has hic resluts 
diagdir     = 'diagnostics'          ##       Plots for diagnostics are held here 
bedtmpdir   = 'bedpe'                ##       Tempeory bedpe dir
checkerdir  = 'checks'               ##       Temporary log to check for script outputs
slurpydir   = './SLURPY'             ##       The script directory holding this file
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

## Set the descripton 
slurpy_description = 'A SLURM Powered, Pythonic Pipeline, Performing Parallel Processing of Piared-end Sequenced Reads Prepaired from 3D Epigenomic Profiles.'

## Write ftn for parsing args 
def slurpy_args(descr=slurpy_description):
    ## Make the parse
    parser = argparse.ArgumentParser(description = descr)

    ## Add the required argument
    parser.add_argument("-r", "--refix",          dest="r",       type=str,  required=True,  help = r_help, metavar = refmetavar                                     ) 
    ## Add default arguments
    parser.add_argument("--fastq","--fq",         dest="fq",      type=str,  required=False, help = fqhelp, metavar = './path/to/fastqs',     default = fastqdir     )
    parser.add_argument("-F", "--fastp-splits",   dest="F",       nargs='+', required=False, help = F_help, metavar = splitsize,              default = [splitsize]  )
    parser.add_argument("-T", "--threads",        dest="T",       type=int,  required=False, help = T_help, metavar = 'n',                    default = 0            )
    parser.add_argument("-P", "--partition",      dest="P",       nargs='+', required=False, help = P_help, metavar = 'tb gpu fast',          default = parts        ) 
    parser.add_argument("-M", "--mtDNA",          dest="M",       type=str,  required=False, help = M_help, metavar = mito,                   default = mito         )
    parser.add_argument("-X", "--exclude",        dest="X",       nargs='+', required=False, help = X_help, metavar = 'chrX, chrY ...',       default = []           )
    parser.add_argument("-Q", "--map-threshold",  dest="Q",       type=int,  required=False, help = Q_help, metavar = map_q_thres,            default = map_q_thres  )
    parser.add_argument("-R", "--rerun-from",     dest="R",       type=str,  required=False, help = R_help, metavar = 'step',                 default = None         )
    parser.add_argument("-a", "--afterok",        dest="a",       type=int,  required=False, help = a_help, metavar = fakejobid,              default = 0            )
    parser.add_argument("-N", "--nice",           dest="N",       type=int,  required=False, help = N_help, metavar = 'n',                    default = nice         )

    ## Set number of threads across software 
    parser.add_argument("-j", "--n-parallel",     dest="j",       type=int,  required=False, help = j_help, metavar = 'n',                    default = nparallel    )
    parser.add_argument("-f", "--fastp-threads",  dest="f",       type=int,  required=False, help = f_help, metavar = fastpthreads,           default = fastpthreads )
    parser.add_argument("-t", "--dask-threads",   dest="t",       type=int,  required=False, help = t_help, metavar = daskthreads,            default = daskthreads  )
    parser.add_argument("-b", "--bwa-threads",    dest="b",       type=int,  required=False, help = b_help, metavar = bwathreads,             default = bwathreads   )
    parser.add_argument("-B", "--bwa-options",    dest="B",       type=str,  required=False, help = B_help, metavar = hic_options,            default = hic_options  )

    ## Set values for Hi-C analysis 
    parser.add_argument("-n", "--run-name",       dest="n",       type=str,  required=False, help = n_help, metavar = 'name',                 default = None         )
    parser.add_argument("-E", "--error-distance", dest="E",       type=int,  required=False, help = E_help, metavar = 'bp',                   default = error_dist   )
    parser.add_argument("-L", "--library",        dest="L",       type=str,  required=False, help = L_help, metavar = 'MboI',                 default = lib_default  )
    parser.add_argument("-Z", "--chunksize",      dest="Z",       type=int,  required=False, help = Z_help, metavar = 'n',                    default = chunksize    )
    parser.add_argument("-G", "--genomelist",     dest="G",       type=str,  required=False, help = G_help, metavar = './path/to/list.tsv',   default = False        )
    parser.add_argument("-J", "--jar-path",       dest="J",       type=str,  required=False, help = J_help, metavar = './path/to/juicer.jar', default = None         )
    parser.add_argument("-xmx", "--Xmemory",      dest="xmx",     type=int,  required=False, help = x_help, metavar = xmemory,                default = xmemory      )
    parser.add_argument("-S", "--bin-sizes",      dest="S",       nargs='+', required=False, help = S_help, metavar = '25000, 10000, ...',    default = binsizes     )
    parser.add_argument("-m", "--max-dist",       dest="m",       type=int,  required=False, help = m_help, metavar = '1000',                 default = max_dist     )
    #parser.add_argument("-gtf","--features",      dest="gxg",     type=str,  required=False, help = A_help, metavar= './path/to/my.gff',      default = 'none'       )
    parser.add_argument("--nodelist",             dest="nodes",   nargs='+', required=False, help = node_help,                                default = None         )
    parser.add_argument("--memory",               dest="slurmem", type=str,  required=False, help = slurmem_help, metavar='40G',              default = None         )

    ## Set Pipeline boolean blags 
    parser.add_argument("--restart",              dest="restart",   help = restart_help,   action = ST)
    parser.add_argument("--nomerge",              dest="merge",     help = merge_help,     action = ST)
    parser.add_argument("--force",                dest="force",     help = force_help,     action = ST)
    parser.add_argument("--debug",                dest="debug",     help = debug_help,     action = ST)
    parser.add_argument("--clean",                dest="clean",     help = clean_help,     action = ST)
    parser.add_argument("--nocount",              dest="count",     help = count_help,     action = ST)

    ## Set Hi-C related boolean flags 
    parser.add_argument("--toshort",              dest="toshort",   help = short_help,     action = ST)
    parser.add_argument("--pairs",                dest="makepairs", help = pairs_help,     action = ST)
    parser.add_argument("--mcool",                dest="mcool",     help = mcool_help,     action = ST)
    parser.add_argument("--hicexplorer",          dest="hicexp",    help = hicex_help,     action = ST)
    parser.add_argument("--inter-only",           dest="inter",     help = inter_help,     action = ST)
    
    ## Set ATAC-seq specifict
    parser.add_argument("--atac-seq",             dest="atac",        help = atac_help,     action = ST)
    parser.add_argument("--broad",                dest="broad",       help = broad_help,    action = ST)
    parser.add_argument("--skipmacs3",            dest="peaks",       help = peaks_help,    action = ST)
    parser.add_argument("--nolambda",             dest="nolambda",    help = lambda_help,   action = ST)
    parser.add_argument("--nomodel",              dest="nomodel",     help = nomodel_help,  action = ST)
    parser.add_argument("--call-summits",         dest="summits",     help = summits_help,  action = ST)
    parser.add_argument("--shift-size",           dest="shiftsize",   help = shift_help,    default=shift_size, required=False, metavar = 'bp')
    parser.add_argument("--extend-size",          dest="extendsize",  help = extend_help,   default=extendsize, required=False, metavar = 'bp')
    parser.add_argument("--macs-mode",            dest="macmode",     help = macs_help,     default='BEDPE',    required=False, metavar = 'BEDPE')
    parser.add_argument("-c", "--controls",       dest="c",           help = c_help,        default=[],         required=False, metavar = c_metavar, nargs='+')
    parser.add_argument("--max-gap",              dest="maxgap",      help = mgap_help,     default=0,          required=False, type=int)
    parser.add_argument("--max-length",           dest="minlen",      help = mlen_help,     default=0,          required=False, type=int)
    parser.add_argument("--max-number-chroms",    dest="maxnc",       help = maxnc_help,    default=max_nchrom, required=False, type=int)

    ## RNA-seq and other vairables 
    parser.add_argument("--rna-seq",              dest="rnas",       help = rnas_help,     action = ST)
    parser.add_argument("--skipfastp",            dest="sfast",      help = skipq_help,    action = ST)
    parser.add_argument("--skipdedup",            dest="skipdedup",  help = mark_help,     action = ST)
    parser.add_argument("--dont-save-dups",       dest="save",       help = save_help,     action = ST)
    parser.add_argument("--dedovetail",           dest="tails",      help = dove_help,     action = ST)
    parser.add_argument("--sam",                  dest="tosam",      help = tosam_help,    action = ST)
    parser.add_argument("--bam",                  dest="tobam",      help = tobam_help,    action = ST)
    parser.add_argument("--wgs",                  dest="wgs",        help = wgs_help,      action = ST)

    ## Set the paresed values as inputs
    return parser.parse_args() 
## EOF 