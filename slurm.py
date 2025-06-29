#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      SLUR(M).py ; A Hi-C, Multiomic, Python Powered Processing Pipeline 
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
./SLURPY/slurm.py -r vero -G ../Chlorocebus_sabeus_mva.genome.sizes.autosome.filtered.bed -F 2000000 -P fast gpu -j 20
squeue -u croth | grep mpi | awk '{print $1}' | xargs -n 1 scancel
squeue -u croth | grep 'croth' | grep gpu | grep "(DependencyNeverSatisfied)" | awk '{print $1}' | xargs -n 1 scancel
""" 
## List slurpy command for VERO analysis
"""
./SLURPY/slurm.py -r vero -P tb fast gpu -G ../Chlorocebus_sabeus_mva.genome.sizes.autosome.filtered.bed -M NC_008066.1 -F 150000 15000000 --merge --restart 

 --nodelist c1003 c1004 c1005 c1006 c0823 c0825 c0825
"""
##      SET PIPELINES
## 
## HIC
## Set pipeline of hic analysis ## NOTE defined after defaults import 
hic_pipeline  = ['fastp', 'bwa', 'filter','dedup','concat','gxg','toshort','hic' ,'macs3','sam','count','clean']
##                  0        1      2        3       4      5a      5b      5c      5d      5e     6B      6C
## Join pipeline names by commas
h_pipe = ', '.join(hic_pipeline) 
inhic = True 
count_mod = ''
## 
##      MODULE LOADING and VARIABLE SETTING 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## Load in vars from parameter space
from parameters import * 
## Load in default directories from defaults 
from directories import * 
## Load in ftns from other libraries
from pysamtools import checksam, writetofile
## Load in panda cat ftn
from pandacat import pandacat
## Load in bwa master
from bwatobedpe import bwamaster
## Load in filter master
from filtermaster import filtermaster
## Load in bedpe to sam
from pairs2sam import bedpetosam

## Set the ftn descritption
hiclite_descr = "Processing and analysis pipeline for paired-end sequencing data from Hi-C experiments."
## Define help messages
R_help = R_help%h_pipe
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


##      FUNCTION DEFINING
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## Load in ftns from defaluts
from defaults import basename, getsamplename, reportname, fastcut, fastdry, ifprint

## Ftn for formating fastp command to filter and split reads
def fastpeel(r1:str, r2:str, w:int, s:int, ishic=True, pix=0, options=fastp_opts, toremoveend='.toberemoved.fastq.gz', singleend='.singletons.fastq.gz', failend='.failed.fastq.gz',z=4) -> tuple:
    """Formats calls to fastp given input fastq files."""
    ## Set variable names 
    r1_bn    = basename(r1).split('.fastq')[0]               ##      Generate basename of first read 
    r2_bn    = basename(r2).split('.fastq')[0]               ##      Generate basename of second read 
    temp1    = f'{splitsdir}/{r1_bn}{toremoveend}'           ##      Format the first temp file for splits to be remove
    temp2    = f'{splitsdir}/{r2_bn}{toremoveend}'           ##      Format the second temp file for splits to be remove 
    single1  = f'{splitsdir}/{r1_bn}{singleend}'             ##      Make singleton file from first split 
    single2  = f'{splitsdir}/{r2_bn}{singleend}'             ##      Make singleton file from second split 
    split1   = f'{splitsdir}/{basename(r1)}'                 ##      Set the first split name
    split2   = f'{splitsdir}/{basename(r2)}'                 ##      Create the second split name 
    report   = reportname(getsamplename(r1),'fastp',i=pix)   ##      Format the report name
    failed   = f'{splitsdir}/{getsamplename(r1)}{failend}'   ##      Generate the list of the failed reads
    ## Format commands if filtering, then splitting is being performed
    ## Call the first run of fastp, json, and json are also returend
    fastp_header0, fastp_json0, fastp_html0 = fastcut(r1,r2,temp1,temp2,report,0)
    ## Set initial filerting
    initial_filtering = fastp_header0 + f' --unpaired1 {single1} --unpaired2 {single2} --failed_out {failed} --thread {w} ' + ' '.join(options[:2]) + '\n'
    ## Call the second run of fastp to split reads into splits, return aslo json, and html 
    fastp_header1, fastp_json1, fastp_html1 = fastcut(temp1,temp2,split1,split2,report,1)
    ## Call fastp the second time to split reads
    split_filtered = fastp_header1 + f' -S {s} -z {z} --thread {w} ' + ' '.join(options) + '\n'

    ## Format command if just splitting is begin performed
    fastp_header, fastp_json, fastp_html = fastcut(r1,r2,split1,split2,report,0)
    ## Call fastp the second time to split reads
    just_split = fastp_header + f' -S {s} -z {z} --thread {w} ' + ' '.join(options) + '\n'

    ## Call the fast dry command and format the remove command 
    dry, throwout = fastdry(r1,r2,report), f'rm {temp1} {temp2}\n'
    ## Format the mv commands for the html and json
    mvh, mvj = f'mv {fastp_json0} {diagdir}/\n', f'mv {fastp_html0} {diagdir}/\n'

    ## Return the comands depenent aupon if we are spliting or filtering, and the reprot 
    return [just_split,dry,mvh,mvj] if ishic else [initial_filtering,split_filtered,throwout,dry,mvh,mvj], report

## Ftn for calling juicer's pre command
def juicerpre(intxt:str, outhic:str, Xmemory:int, jarfile:str, threadcount:int, bins:list, genomepath:str) -> tuple:
    """
    java -Xmx49152m -Xms49152m -jar $jarpath pre -j 5 -r 500000,250000,200000,150000,100000,50000,25000,10000,5000,1000 $1 ${1}.hic $2
    """
    ## Format the report name 
    report = reportname(outhic+'.bam','hic.pre',i='6C')
    ## Set the java command for the passed juicer jar file 
    prestr = ['java -Xmx%sm -Xms%sm -jar %s pre -j %s -r %s %s %s %s\n'%(Xmemory,Xmemory,jarfile,threadcount,','.join(map(str,bins)),intxt,outhic,genomepath),
              f'{slurpydir}/myecho.py Finished formating Hi-C contacts into a .hic file on path: {outhic} {report}\n']

    ## Return the pre and report
    return prestr, report
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


##      MAIN SCRIPT & ARGUMENT PARSING 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## If the script is envoked by name 
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse
    ## Make the parse
    parser = argparse.ArgumentParser(description = slurpy_descr)

    ## Add the required argument
    parser.add_argument("-r", "--refix",          dest="r",       type=str,  required=True,  help = r_help, metavar = refmetavar                                     ) 
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
    parser.add_argument("-x", "--Xmemory",        dest="x",       type=int,  required=False, help = x_help, metavar = xmemory,                default = xmemory      )
    parser.add_argument("-S", "--bin-sizes",      dest="S",       nargs='+', required=False, help = S_help, metavar = '25000, 10000, ...',    default = binsizes     )
    parser.add_argument("-m", "--max-dist",       dest="m",       type=int,  required=False, help = m_help, metavar = '1000',                 default = max_dist     )
    parser.add_argument("-gtf","--features",      dest="gxg",     type=str,  required=False, help = A_help, metavar= './path/to/my.gff',      default = 'none'       )
    parser.add_argument("--nodelist",             dest="nodes",   nargs='+', required=False, help = node_help,                                default = None         )

    ## Set Pipeline boolean blags 
    parser.add_argument("--restart",              dest="restart",   help = restart_help,   action = ST)
    parser.add_argument("--nomerge",              dest="merge",     help = merge_help,     action = ST)
    parser.add_argument("--force",                dest="force",     help = force_help,     action = ST)
    parser.add_argument("--debug",                dest="debug",     help = debug_help,     action = ST)
    parser.add_argument("--clean",                dest="clean",     help = clean_help,     action = ST)
    parser.add_argument("--count",                dest="count",     help = count_help,     action = ST)

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
    parser.add_argument("-c", "--controls",       dest="c",           help = c_help,        default=None,       required=False, metavar = c_metavar, nargs='+')
    parser.add_argument("--max-gap",              dest="maxgap",      help = mgap_help,     default=0,          required=False)
    parser.add_argument("--max-length",           dest="minlen",      help = mlen_help,     default=0,          required=False)

    ## RNA-seq and other vairables 
    parser.add_argument("--rna-seq",              dest="rnas",       help = rnas_help,     action = ST)
    parser.add_argument("--skipfastp",            dest="sfast",      help = skipq_help,    action = ST)
    parser.add_argument("--skipdedup",            dest="skipdedup",  help = mark_help,     action = ST)
    parser.add_argument("--dedovetail",           dest="tails",      help = dove_help,     action = ST)
    parser.add_argument("--sam",                  dest="tosam",      help = tosam_help,    action = ST)
    parser.add_argument("--bam",                  dest="tobam",      help = tobam_help,    action = ST)

    ## Set the paresed values as inputs
    inputs = parser.parse_args() 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
   

    ##      PASS ARGUMENTS & SET VARIABLES 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Set required variables
    reference_path  = inputs.r            ##     Set path to the reference genome

    ## Set default vairables              ##
    splitsize       = inputs.F            ##     Number of splits in fastp, this is a line count, multiple by four to get actual read counts
    threadn         = inputs.T            ##     Set the number of parallel runs of bwa 
    partitions      = inputs.P            ##     Set the partition 
    mito            = inputs.M            ##     Set the mito contig name 
    excludes        = inputs.X            ##     List of chromosomes to exclude from analysis 
    mapq            = inputs.Q            ##     Set the mapping quality threshold 
    rerun           = inputs.R            ##     Setp to rerun pipeline from 
    bwaix_jobid     = inputs.a            ##     The job id to have all submissions wait on   
    nice            = inputs.N            ##     Sets the nice parameter 

    ## Set threads across softwares                         
    nparallel       = inputs.j            ##     Number of parallel jobs to run at once acorss bwa and filter 
    fastp_threads   = inputs.f            ##     Number of fastp threads
    daskthreads     = inputs.t            ##     Number of threads in dask
    bwa_threads     = inputs.b            ##     Number of threads in bwa alignments
    bwa_opts        = inputs.B            ##     Set options for bwa 

    ## Set variables for Hi-C                
    run_name        = inputs.n            ##     The name of the samples 
    error_dist      = inputs.E            ##     Distance in bp to check hi-C contacts for errors 
    enzymelib       = inputs.L            ##     Restriction-enzymelib used in Hi-C prep 
    chunksize       = inputs.Z            ##     Chunk size (row number) to load in with pandas 
    pathtochrom     = inputs.G            ##     Path to list of chromosomes to use 
    jarpath         = inputs.J            ##     Path to juicer jar file 
    xmemory         = inputs.x            ##     Amount of memory passed to juicer pre command 
    binsizes        = inputs.S            ##     Bins / resolutions used in hi-c analysis 
    max_dist        = inputs.m            ##     Maximum distance of paired end reads
    feature_space   = inputs.gxg          ##     Path to a gff or bed file used in g x g interaction matrix / df 
    nodes           = inputs.nodes        ##     List of nodes 
    
    ## Set pipeline boolean vars                    
    hardreset       = inputs.restart      ##     Resetart the slurpy run, removing previous
    postmerging     = not inputs.merge    ##     Forces premerge of outputs 
    force           = inputs.force        ##     Force overwrite of output alignment files 
    debug           = inputs.debug        ##     Run in debug mode 
    ifclean         = inputs.clean        ##     Flag to run clean up script 
    counting        = inputs.count        ##     Flag to count the read pairs

    ## Set Hi-C related boolean flasgs                  
    toshort         = inputs.toshort      ##     Flag the make short file, kicks if jarpath was given 
    makepairs       = inputs.makepairs    ##     Flag to make pairs file 
    make_mcool      = inputs.mcool        ##     Flag to make mcool file  
    hicexplorer     = inputs.hicexp       ##     Flag to run stricter filtering like Hi-C explorer
    get_inter       = inputs.inter        ##     Flag to return inter chromosome counts in bedpe file
    
    ## MACS3 related vars
    atac_seq        = inputs.atac         ##     Boolean flag to run in atac-seq mode 
    ifbroad         = inputs.broad        ##     Boolean to activate broader peak calling in macs3 
    skippeaks       = inputs.peaks        ##     Skips peak calling with macs3 
    nolambda        = inputs.nolambda     ##     Flag to skip local lambda control in macs3
    nomodel         = inputs.nomodel      ##     Flag to skip shifting model
    callsummits     = inputs.summits      ##     Flag to call summits in macs3
    shift_size      = inputs.shiftsize    ##     Value of shift size (bp)
    extendsize      = inputs.extendsize   ##     Value of read extension 
    macs3mode       = inputs.macmode      ##     Format in MACS3, eg. bedpe.
    chip_control    = inputs.c            ##     Set the input control for chip experimetn
    max_gap         = inputs.maxgap       ##     Max gap in used in MACS3    
    min_len         = inputs.minlen       ##     Min lenght of peaks in MACS3

    ## Set RNA-seq like vars 
    rna_seq         = inputs.rnas         ##     Boolean flag to run in rna-seq mode 
    sfastp          = inputs.sfast        ##     Flag to skip fastp filtering 
    skipduplicates  = inputs.skipdedup    ##     Boolean to mark duplicates       
    dedovetail      = inputs.tails        ##     Boolean for dove tailing 
    tosam           = inputs.tosam        ##     Boolean flag to convert .bedpe file to sam format
    tobam           = inputs.tobam        ##     Boolean flag to convert .bedpe file to bam format

    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    
    ##      CORRECT SPLIT / Chunksize
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    max_splitsize = max([int(s) for s in splitsize])
    if chunksize > max_splitsize:
        ## Reset chunksize
        chunksize = int(round(max_splitsize/3,0))
        print("INFO: Reseting chunksize (Z) to: %s"%chunksize)

    ##      PRESET Thread counts
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Check if we are in atac seq mode    
    if threadn:
        ## reset threads                           
        fastp_threads   = threadn           ##     Number of fastp threads
        bwa_threads     = threadn           ##     Number of threads in bwa alignments
        daskthreads     = threadn           ##     Number of threads in dask 

    ## Reset paortions into a comma joined list
    partition = ','.join(partitions)

    ##      ROTH SETTINGS
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Reset reference with presets if we are running as Cullen Roth (These can be edited if you are not me :-))
    ist2t  = (reference_path.lower() == 't2t')
    isvero = (reference_path.lower() == 'vero')

    ## If human reference 
    reference_path = t2t_refpath if ist2t else reference_path
    feature_space  = t2t_gffpath if (feature_space.lower()  == 't2t') else feature_space

    ## If vero was pass
    reference_path = vero_refpath if isvero else reference_path

    ## Set bwa master partition to mpi (FOR ROTH only)
    bwa_partition   = 'mpi' if ist2t or isvero else partition
    filt_partition  = 'mpi' if ist2t or isvero else partition
    clean_partition = 'mpi' if ist2t or isvero else partition
    time_partition  = 'mpi' if ist2t or isvero else partition
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      INITILIZATION 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Load in ftns from defalut
    from defaults import pathexists, patchclean, fileexists
    ## Check that the fastq and reference path exists
    assert pathexists(fastqdir), fastqserror
    assert pathexists(reference_path), noref_path%reference_path 
    ## Check the versions of samtools, the user email is an email and the experiment mode is one we know
    assert checksam(), not_sam_err 
    ## Reformat clean boolean if clean was passed from restart
    ifclean = patchclean(rerun,ifclean)
    ## Initilizse list for sbatch
    sub_sbatchs = []

    ## Reset bwa options
    bwa_opts = bwa_opts if (bwa_opts[0] == ',') else (','+bwa_opts)

    ## Set the lib error
    lib_error = "ERROR: The passed library name of enzyme(s) used in Hi-C prep was not recognized."
    ## Check library if it was pased
    assert enzymelib.lower() in ['mboi', 'dpnii','sau3ai','hindiii','arima','none'], lib_error

    ## Check if we have a jarpath
    if jarpath:
        assert fileexists(jarpath), "ERROR: The given jarpath %s could not be found! Please check path the given path and try again."%jarpath
        ## Set othe booleans 
        toshort    = True 
        make_mcool = False
        makepairs  = False 

    ## If we are making a cool file 
    if make_mcool:
        makepairs = True 
        jarpath   = False
        toshort   = False

    ## Set feature space to a boolean var, needed to make my settings work 
    if feature_space == 'none':
        feature_space = False
    ## Check is a path the give feature space 
    if feature_space:
        assert fileexists(feature_space), "ERROR: The given features data %s could not be found! Please check the given path and try agian."%feature_space
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      PRESET ATAC and ChIP-seq mode
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Check if we are mapping RNA-seq reads
    if rna_seq:
        keep_dovetail = True
        inhic         = False
        enzymelib     = 'none'
        max_dist      = 0
        count_mod     = '--atac-seq'
        bwa_opts      = ',-M' if (bwa_opts == hic_options) else bwa_opts
        skippeaks     = True 
    
    ## Check if we are in atac seq mode    
    if atac_seq or chip_control:
        keep_dovetail = True
        inhic         = False
        enzymelib     = 'none'
        max_dist      = max_dist if max_dist else 1000
        count_mod     = '--atac-seq'
        bwa_opts      = ',-M' if (bwa_opts == hic_options) else bwa_opts
    
        ## Load in macs3 ftns
        from pymacs3 import peakattack
        ## Set the broad pkeack
        broadpeak = '--broad' if ifbroad else ''

    ## Check our control files if they were passed 
    if chip_control:
        assert (type(chip_control) == list) and (type(chip_control[0]) == str), 'ERROR: Inputs for chip control are not a type we recognize (should be a list of strings)!'
        ## Iterate over the controls inputs 
        for chip_con in chip_control:
            ## Check they exist 
            assert fileexists(chip_con), f'ERROR: Bad path to bam file! Unable to locate input control for chip experiment: {chip_con}'
    
    ## Reset experiment mode 
    if atac_seq:
        experi_mode = 'atac'
    if rna_seq:
        experi_mode = 'rnas'
    elif chip_control:
        experi_mode = 'chip'
    else: ## Set the hardset hic mode 
        experi_mode = 'hic'

    ## Set peakcalling boolean 
    peakcalling = atac_seq and (not skippeaks)

    ## Moved to trigger if statment correctly
    print("INFO: Running analysis with the %s Hi-C library."%enzymelib) if inhic else None 

    ## Reset shift and extention size 
    shift_size,extendsize = (shift_size,extendsize) if nomodel else (0,0)
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      CONFIRM THE HARD RESTART 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Load in defluats 
    from defaults import confirmreset
    ## Format group dirs 
    grouped_dirs = [debugdir,aligndir,splitsdir,comsdir,diagdir,bedtmpdir,hicdir,checkerdir] + ([macs3dir] if ((atac_seq or chip_control) and not skippeaks) else [])
    ## If there is a hard reset passed 
    confirmreset(grouped_dirs) if hardreset else None
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      MORE MODULE LOADING (and time stamp setting)
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Load in the list dir ftn 
    from os import getcwd as gcwd
    ## Load in time mod
    import time 
    ## Set the unique time stamp
    stamp = time.time()
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      DIRECTORY MAKING & TIME STAMP SUBMISSION 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Load in ftns for directory and run name 
    from defaults import setrunname, makedirectories, sortglob, remove
    ## Let the user know we are making directories 
    print(directormaking)
    ## Get the current working dir
    the_cwd = gcwd()
    ## If the run name is none
    run_name = setrunname(run_name,the_cwd)
    ## Make the directories 
    makedirectories(grouped_dirs)
    ## Gather and remove the old command files, error, and output logs if they are there 
    [remove(f) for f in sortglob(f'{comsdir}/*.sh') + sortglob(f'{debugdir}/*.err') + sortglob(f'{debugdir}/*.out')]
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      WRITING OUT PARAMS
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Load in ftn from defaluts
    from defaults import writeparams
    ## Write out params 
    writeparams(f'{experi_mode}.py',run_name,stamp,inputs)
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      CHROMOSOME GATHERING 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Load in more mods
    from chrommap import gathering, chromgathering
    ## Gather chromosomes and inform the user we are gathering chromosomes
    print(chromgathering)
    ## Expand exlcude list to include mitochondria contig
    excludes.append(mito)
    ## Gather a list of chromosomes 
    chrlist,genome_size,pathtochrom = gathering(reference_path,pathtochrom,excludes,not make_mcool)
    ## Calculate chrlist
    nchrom = len(chrlist)
    ## Print if there are fewer than ten
    print('INFO: Running on -> %s'%' '.join(chrlist)) if (nchrom < 10) else None
    print('INFO: Processing %s chromosomes from list.'%nchrom)
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      CHECK FOR A BWA INDEXING  
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Set index error
    index_error = index_error%reference_path
    ## Load in bwaix check
    from defaults import isbwaix
    ## Assert our truth
    assert isbwaix(reference_path),index_error
    
    ## Format the command file and initilize sample names
    command_files, samplenames = [],[]

    ## Append to command
    command_files.append(('jobfile','sample',experi_mode,'bwaix','report',0,bwaix_jobid)) if bwaix_jobid else None 
    #print(command_files)
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      FASTQ GATHERING
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Bring in get fastqs fromd efaluts
    from defaults import getfastqs, sortfastq
    ## Inform user we are formating jobs
    print(formatingfastq)
    ## Gather the fastqs 
    in_fastqs = getfastqs(fastqdir+'/*.fastq.gz')
    ## Assert we have fastq files
    assert len(in_fastqs), missingfqs
    ## Sort by fastq size
    in_fastqs = sortfastq(in_fastqs,splitsize)
    ## Saveout the sizes to debug dir 
    in_fastqs.to_csv(f'{debugdir}/fastq.sizes.csv',header=True,index=False)
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      SAMPLE NAME MAKING AND HANDELING
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Load in sbatch from defaults
    from defaults import sbatch
    ## iterate thru the pairs 
    for fastqix, fastrow in in_fastqs.iterrows():
        ## Split off the row values from our sorted dataframe 
        r1,r2,sizeosplit = fastrow.Read1,fastrow.Read2,fastrow.Splitsize
        ## Gather the sample name, exp mode, and fastp command file 
        sample_name = getsamplename(r1) 
        ## Append the sample names
        samplenames.append(sample_name)
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


        ##      FASTP SPLITTING 
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## Iniate pix at 0, set fastp command interface
        pix = 0
        ##  Set the the fastp command file name 
        fastp_command_file =  f'{comsdir}/{pix}.fastp.{sample_name}.sh'  
        ##  Gather the fastp command and report 
        fastp_coms, fastp_repo = fastpeel(r1,r2,fastp_threads,sizeosplit,ishic=inhic or sfastp)
        ##  Write the command to file
        writetofile(fastp_command_file, sbatch(None,fastp_threads,the_cwd,fastp_repo,nice=1,nodelist=nodes) + fastp_coms, debug)
        ##  Append command to file
        command_files.append((fastp_command_file,sample_name,experi_mode,hic_pipeline[pix],fastp_repo,0,''))
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

        
        ##      BWA MASTER JOB SUBMISSION
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## Set the point of pipeline
        pix = 1
        ## Call the bwa master command
        bwa_master_file = f'{comsdir}/{pix}.bwa.master.{sample_name}.sh'
        ## Gahter the bwa master command and report
        bwa_master_commands, bwa_master_repo = bwamaster(sample_name,reference_path,bwa_threads,the_cwd,partition,debug,nice,nparallel,library=enzymelib,forced=force,nodelist=nodes,bwaopts=bwa_opts)
        ## Write command to file
        writetofile(bwa_master_file, sbatch('bwa.master',1,the_cwd,bwa_master_repo,nice=nice,nodelist=nodes) + bwa_master_commands, debug)
        ## Append to command fil
        command_files.append((bwa_master_file,sample_name,experi_mode,hic_pipeline[pix],bwa_master_repo,0,''))
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


        ##      MASTER FILTER BEDPE / Hi-C CONTACTS JOB SUBMISSION 
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------ ##
        ## 2. Set the setp in pipeline and the new bedbe file name
        pix = 2
        ## Call the master filter command
        filter_master_file = f'{comsdir}/{pix}.filter.bedpe.master.{sample_name}.sh'
        ## Gather the filter master commadn and report
        filter_master_commands, filter_master_repo = filtermaster(sample_name,reference_path,the_cwd,excludes,chrlist,mapq,error_dist,daskthreads,enzymelib,partition,debug,nice,nparallel,forced=force,maxdist=max_dist,chunksize=chunksize,nodelist=nodes,dovetail=dedovetail,removeinter= not inhic,hicexplorer=hicexplorer)
        ## Write command to file
        writetofile(filter_master_file, sbatch('filter.bedpe.master',1,the_cwd,filter_master_repo,nice=nice,nodelist=nodes) + filter_master_commands, debug)
        ## Append to command file
        command_files.append((filter_master_file,sample_name,experi_mode,hic_pipeline[pix],filter_master_repo,0,''))
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

  
    ##      MERGING, CONCATONATION, DEDUPLICATION, SORTING 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Iniate new cat files
    new_catfiles = []
    ## Iterate over the sample names
    for si,sname in enumerate(samplenames):
        ## Set the smaple name 
        sample_name = run_name if postmerging else sname 
        ## If the sample is the first and we are merging all smaples 
        if (si > 0) and postmerging:
            continue
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## 3. Deduplicating and sorintg 
        ## Initiate list of out put hic files for each chromosome
        hiccats_outs = []
        hicdups_outs = []
        ## Set pipeline step 
        pix = 3
        ## Set the new concat file name by chromosome index 
        for chrom in chrlist:
            ## Set start name for wildcard use to bring in inputs  to ftn 
            sample_start = f'valid.{chrom}.bedpe' if postmerging else f'{sample_name}.valid.{chrom}.bedpe'
            hiccat_out = f'{hicdir}/{sample_name}.valid.{chrom}.bedpe'
            hicdup_out = f'{hicdir}/{sample_name}.duplicates.{chrom}.bedpe'
            ## Append to hic cat outs to preserve order 
            hiccats_outs.append(hiccat_out)
            hicdups_outs.append(hicdup_out)
            ## SEt the report name
            hiccat_repo = reportname(sample_start,hic_pipeline[pix],i=pix)
            ## Set the command
            hiccat_coms = [f'{slurpydir}/deduphic.py -b {sample_start} -o {hiccat_out} -d {hicdup_out} --sort' +  ('\n' if skipduplicates else ' --dedup\n')]

            ## make concat file name
            hiccat_file = f'{comsdir}/{pix}.dedup.{sample_name}.valid.{chrom}.sh'
            ## Write the concat command to file
            writetofile(hiccat_file, sbatch(hiccat_file,daskthreads,the_cwd,hiccat_repo,nice=nice,nodelist=nodes) + hiccat_coms, debug)
            ## Append the concat command
            command_files.append((hiccat_file,sample_name,experi_mode,hic_pipeline[pix],hiccat_repo,0,''))    
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

        
        ## Hi-C CONCATONATION
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## 4. Merging unused contacts, and concatonation of duplicates and valid sets 
        ## Set pipeline step 
        pix = 4
        ## Set start name for wildcard use to bring in inputs to ftn , the output file and the ouput report, set ouptu name  
        unused_start  = [bedtmpdir + ('/*.notused.bedpe' if postmerging else f'/*.{sample_name}.notused.bedpe')]
        
        ## Set file names, concat the inputs, set output names, repor tnames, and files 
        hic_fileends   = ['notused','duplicates','valid']
        concat_inputs  = [unused_start,hicdups_outs,hiccats_outs]
        concat_outputs = [f'{aligndir}/{sample_name}.{f}.bedpe' for f in hic_fileends]
        concat_reports = [reportname(f,hic_pipeline[pix],i=f'{pix}{a}') for a,f in zip(['A','B','C'],concat_outputs)]
        concat_files   = [f'{comsdir}/{pix}.concat.{sample_name}.{f}.sh' for f in hic_fileends]

        ## Iterate thru the zipped concat files 
        for (con_ins,con_out,con_rep,con_fil) in zip(concat_inputs,concat_outputs,concat_reports,concat_files):
            ## Format the commands 
            concat_coms = pandacat(con_ins,con_out,rmheader=True)
            ## Write commands to file 
            writetofile(con_fil, sbatch(con_fil,1,the_cwd,con_rep,nice=nice,nodelist=nodes) + concat_coms, debug)
            ## append command to command file list
            command_files.append((con_fil,sample_name,experi_mode,hic_pipeline[pix],con_rep,0,'')) 

        ## Patch new cat file to next analysis step, must be the valid contacts/alignments 
        newcatfile = concat_outputs[-1]
        new_catfiles.append(newcatfile)
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


        ## OPTIONAL ANALYZES
        ## Set pipeline step 
        pix = 5
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        """
        Optional analyzes
        gene x gene counts
        bedpe ---> short format
        jucier pre command for hic creation 
        """
        ## GENE X GENE from HIC
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## 5A. Create a feature space
        if inhic:  
            if feature_space:
                ## Iterate over chromosome list 
                for i,coi in enumerate(chrlist):
                    ## Set the report, commands, and gxg script file name 
                    gxg_repo   = reportname(sample_start,'gxg%s'%i,i=f'{pix}A')
                    gxg_commands = [f'{slurpydir}/gxgcounts.py -i {newcatfile} -c {coi} -f {feature_space} -t {nchrom}' + (' --merge\n' if not i else '\n')]
                    gxg_file     = f'{comsdir}/{pix}A.gxg.{i}.{sample_name}.sh'
                    ## Write to file for the gxg script, pasing dask thread count, the cwd, commands and debug mode 
                    writetofile(gxg_file,sbatch(gxg_file,daskthreads,the_cwd,gxg_repo,nice=nice,nodelist=nodes) + gxg_commands, debug)
                    ## Append to command list for calling/submission later 
                    command_files.append((gxg_file,sample_name,experi_mode,'gxg',gxg_repo,0,''))
            ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

            ## BEDPE TO SHORT or PAIRS
            ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
            ## 5B. Make the bedpe file into a short file for juicer pre command 
            if toshort: 
                ## make a report
                short_repo = reportname(sample_name,'toshort',i=f'{pix}B')
                ## Format the command
                short_commands = [f'{slurpydir}/toshort.py -i {newcatfile}\n'] 
                ## Set command file 
                short_file = f'{comsdir}/{pix}B.toshort.{sample_name}.sh'
                ## Wriet the short command to file
                writetofile(short_file,sbatch(short_file,daskthreads,the_cwd,short_repo,nice=nice,nodelist=nodes) + short_commands, debug)
                ## append the short command
                command_files.append((short_file,sample_name,experi_mode,'toshort',short_repo,0,''))

            if get_inter:
                ## make a report
                short_repo = reportname(sample_name,'inter.short',i=f'{pix}B')
                ## Format the command
                short_commands = [f'{slurpydir}/toshort.py --inter-only -i {newcatfile}\n'] 
                ## Set command file 
                short_file = f'{comsdir}/{pix}B.inter.short.{sample_name}.sh'
                ## Wriet the short command to file
                writetofile(short_file,sbatch(short_file,daskthreads,the_cwd,short_repo,nice=nice,nodelist=nodes) + short_commands, debug)
                ## append the short command
                command_files.append((short_file,sample_name,experi_mode,'toshort',short_repo,0,''))

            ## 5B. Make the bedpe file into a short file for juicer pre command 
            if makepairs:
                ## make a report
                short_repo = reportname(sample_name,'pairs',i=f'{pix}B')
                ## Format the command
                short_commands = [f'{slurpydir}/toshort.py --pairs -i {newcatfile}\n'] 
                ## Set command file 
                short_file = f'{comsdir}/{pix}B.pairs.{sample_name}.sh'
                ## Wriet the short command to file
                writetofile(short_file,sbatch(short_file,daskthreads,the_cwd,short_repo,nice=nice,nodelist=nodes) + short_commands, debug)
                ## append the short command
                command_files.append((short_file,sample_name,experi_mode,'toshort',short_repo,0,''))
                ## Format the pairs file name
                newpairsfile = newcatfile.split('.bedpe')[0] + '.pairs'
            ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

            ##  HIC FILE CREATION / JUICER PRE 
            ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
            ## 5C. IF we have a valid jar path
            if jarpath:
                ## Set the output hic path
                shortfile, outhicpath = newcatfile.split('.bedpe')[0] + '.short', newcatfile.split('.bedpe')[0] +'.hic'
                ## Sort and saveout the the concatonated hic file 
                jpre_coms, jpre_repo = juicerpre(shortfile,outhicpath,xmemory,jarpath,fastp_threads,binsizes,pathtochrom)
                ## make concat file name
                jpre_file = f'{comsdir}/{pix}C.juicerpre.{sample_name}.sh'
                ## Write the concat command to file
                writetofile(jpre_file, sbatch(jpre_file,fastp_threads,the_cwd,jpre_repo,nice=nice,nodelist=nodes) + jpre_coms, debug)
                ## Append the concat command
                command_files.append((jpre_file,sample_name,experi_mode,'hic',jpre_repo,0,''))

            ## 5C. Otherwise, make a cooler and mcool file
            if make_mcool: 
                ## format call to cooler
                """
                From: https://liz-fernandez.github.io/HiC-Langebio/04-matrix.html
                cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 chr_file.txt:10000 ZmEn_HiC_2_1_2.hicup.bsorted.pairs ZmEn_2_10k.cool
                """
                ## Set out cool and mcool files
                outcool  = newpairsfile + '.cool'
                outmcool = newpairsfile + '.mcool'
                ## Set cooler report and file name 
                coolrepo = reportname(sample_name,'mcool',i=f'{pix}C')
                coolfile = f'{comsdir}/{pix}C.mcool.{sample_name}.sh'
                ## Set cooler commands 
                cooler_coms = [f'cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 {pathtochrom}:{min(binsizes)} {newpairsfile} {outcool}\n',
                            f'cooler zoomify {outcool} -r {",".join(map(str,binsizes))} -o {outmcool}\n'
                            f'rm {outcool}\n',
                            f'{slurpydir}/myecho.py Finished formating mcool file! {coolrepo}\n']
                ## Write the concat command to file
                writetofile(coolfile, sbatch(coolfile,fastp_threads,the_cwd,coolrepo,nice=nice,nodelist=nodes) + cooler_coms, debug)
                ## Append the concat command
                command_files.append((coolfile,sample_name,experi_mode,'hic',coolrepo,0,''))
            ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

        ##      PEAK CALLING WITH MACS3
        ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## 5D. If we are running analysis on atac-seq experiments and the peak calling is taking place 
        elif peakcalling:
            ## Format the macs3 call report name
            macs3_report, macs3_filename = reportname(sample_name,'macs3',i=f'{pix}D'), f'{comsdir}/{pix}D.macs3.{sample_name}.sh'
            ## Format the command to macs3
            macs3_commands = peakattack(newcatfile,sample_name,macs3_report,macs3mode.upper(),gsize=genome_size,incontrols=chip_control,shiftsize=shift_size,extendsize=extendsize,maxgap=max_gap,minlen=min_len,nolambda=nolambda,broad=broadpeak,summits=callsummits) + [f'{slurpydir}/pymacs3.py -b {macs3dir}/{sample_name}*.valid.bed* -p {macs3dir}/{sample_name}_peaks.narrowPeak -m {macs3mode} -s {diagdir}/{sample_name}.frip.stats.csv -g {genome_size}\n',f'{slurpydir}/myecho.py Finished calculating FrIP from macs3 {macs3_report}\n']
            ## Write the macs3 commands to file
            writetofile(macs3_filename, sbatch(macs3_filename,1,the_cwd,macs3_report) + macs3_commands, debug)
            ## Append the macs3 command 
            command_files.append((macs3_filename,sample_name,experi_mode,'macs3',macs3_report,0,''))

        ##      SAM or BAM CONVERSION
        ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## 5E. If we are converting bedpe to sam or bam
        if tosam or tobam:
            ## Format the sam file name
            sam_filename = f'{comsdir}/{pix}E.sam.{sample_name}.sh'
            ## format the commadn to pairs2sam
            sam_commands, sam_report = bedpetosam(newcatfile,pathtochrom,samthreads,tobam,sample_name)
            ## Wriet the mac
            writetofile(sam_filename, sbatch(sam_filename,samthreads,the_cwd,sam_report) + sam_commands, debug)
            ## Append the sam command to command file list
            command_files.append((sam_filename,sample_name,experi_mode,'sam',sam_report,0,''))

    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## <- 
    
    ##      SUBMITTING TIME-STAMP   
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## 6A. Final timestamp, out clean up
    pix = 6
    ## Format count commands 
    ## Set the timesampe assocaited file names 
    timestamp_file   = f'{diagdir}/{run_name}.timestamp.{stamp}.txt'             ##     Name of the output file 
    timestampsh      = f'{comsdir}/{pix}A.time.stamp.sh'                         ##     Name of the .sh bash file 
    timestamp_repo = reportname(run_name,f'timestamp.{stamp}',i=f'{pix}A')       ##     Name of the log to report to 
    ## Formath time stamp and echo commands 
    times_commands = [ f'{slurpydir}/totalcount.py {run_name} {diagdir}\n', f'{slurpydir}/endstamp.py {timestamp_file} {stamp}\n']
    ## Format the command file name and write to sbatch, we will always ask the timestamp to run even in debug mode 
    writetofile(timestampsh, sbatch(timestampsh,1,the_cwd,timestamp_repo,nice=1,nodelist=nodes) + times_commands, False)
    ## Append the timestamp command to file
    command_files.append((timestampsh,run_name,experi_mode,'timestamp',timestamp_repo,0,''))
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

    ##      COUNTING 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## 6B. If the clean boolean or rerun vars were set 
    if counting:
        ## Format command to remove uneedeed files 
        counting_sh   = f'{comsdir}/{pix}B.thecount.sh'             ##   Set the bash file name 
        counting_repo = reportname(run_name,'thecount',i=f'{pix}B')   ##   Set the report 
        counting_coms = [f'{slurpydir}/totalcount.py {run_name} {diagdir}\n'] + [f'{slurpydir}/pairs2count.py -i {newf} -c {chunksize} {count_mod}\n' for newf in new_catfiles ]
        ## Format the command to clean up          
        writetofile(counting_sh, sbatch(counting_sh,4,the_cwd,counting_repo,nice=nice,nodelist=nodes)+ counting_coms, debug)
        ## Append the clean up command to file
        command_files.append((counting_sh,run_name,experi_mode,'count',counting_repo,0,''))
    else: ## Otherwise do nothing
        pass

    ##      CLEANING UP FILES   
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## 6C. If the clean boolean or rerun vars were set 
    if ifclean or (rerun == 'clean'): 
        ## Format command to remove uneedeed files 
        remove_sh   = f'{comsdir}/{pix}C.cleanup.sh'             ##   Set the bash file name 
        remove_repo = reportname(run_name,'clean',i=f'{pix}C')   ##   Set the report 
        remove_comm = [f'{slurpydir}/remove.py {bedtmpdir} {splitsdir} {hicdir} {checkerdir}\n', f'{slurpydir}/gzipy.py ./{aligndir}/*.bedpe\n', f'{slurpydir}/gzipy.py ./{aligndir}/*.short\n', f'{slurpydir}/gzipy.py ./{aligndir}/*.valid.pairs\n']
        ## Format the command to clean up          
        writetofile(remove_sh, sbatch(remove_sh,1,the_cwd,remove_repo,nice=nice,nodelist=nodes)+ remove_comm, debug)
        ## Append the clean up command to file
        command_files.append((remove_sh,run_name,experi_mode,'clean',remove_repo,0,''))
    else: ## Otherwise do nothing
        pass
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      PIPELINE COMMAND SUBMISSION TO SLURM    
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Bring in commmand control, fastp submitter, and submit depends
    from defaults import commandcontrol, submitfastp, submitdependency
    ##      A) RESTARTING                          
    ## Call the command dataframe, remove previous logs if hard reset was passed 
    command_files, was_hard_reset = commandcontrol(command_files,hardreset,hic_pipeline,rerun)
    #print(command_files)
    ## 
    ##      0) SUBMITTING FASTP SPLITTING          
    ## Call the submission ftn for fastp commands 
    command_files, sub_sbatchs = submitfastp(command_files,sub_sbatchs,partition,stamp,debugmode=debug)
    ## 
    ##      1) SUBMITTING BWA ALIGNMENTS           
    ## Call the submit bwa ftn
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,hic_pipeline[1],hic_pipeline[0],stamp,bwa_partition,debug=debug)
    ##
    ##      2) SUBMITTING BEDPE FILTERING and SPLITTING
    ## Submit the bedpe filtering scirpt 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,hic_pipeline[2],hic_pipeline[1],stamp,filt_partition,debug=debug)
    ##
    ##      3) SUBMITTING DEDUPLICATING AND SORTING 
    ## Submit the dedup and sort script 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,hic_pipeline[3],hic_pipeline[2],stamp,partition,debug=debug,group='Experiment' if postmerging else 'Sample')
    ##
    ##      4) SUBMITTING CONCATONATION 
    ## Submit the concatonation sciprt 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,hic_pipeline[4],hic_pipeline[3],stamp,partition,debug=debug,group='Experiment' if postmerging else 'Sample')
    ##
    ##      5A) SUBMITTING Gene X Gene interaction script 
    ## Submit the g x g script
    sub_sbatchs = sub_sbatchs + (submitdependency(command_files,'gxg',hic_pipeline[4],stamp,partition,debug=debug,group='Experiment' if postmerging else 'Sample') if feature_space else []) 
    ##
    ##      5B) SUBMITTING CONVERSION FROM BEDPE to JUICER SHORT 
    ## Submit the toshort.py command 
    sub_sbatchs = sub_sbatchs + (submitdependency(command_files,'toshort',hic_pipeline[4],stamp,partition,debug=debug,group='Experiment' if postmerging else 'Sample') if (toshort or makepairs or get_inter) else []) 
    ##
    ##      5C) Hi-C FILE CREATION 
    ## Call the juicer pre command for hic file creation if jarpath was passed 
    sub_sbatchs = sub_sbatchs + (submitdependency(command_files,'hic','toshort',stamp,partition,debug=debug,group='Experiment' if postmerging else 'Sample') if (jarpath or make_mcool) else [])
    ##
    ##      5D) PEAK CALLING w/ MACS3
    ## Call the peak calling command 
    sub_sbatchs = sub_sbatchs + (submitdependency(command_files,'macs3',hic_pipeline[4],stamp,partition,debug=debug,group='Experiment' if postmerging else 'Sample') if peakcalling else [])
    ##
    ##      5E) SAM or BAM CONVERSION
    ## Call the sam conversion script 
    sub_sbatchs = sub_sbatchs + (submitdependency(command_files,'sam',hic_pipeline[4],stamp,partition,debug=debug,group='Experiment' if postmerging else 'Sample') if (tosam or tobam) else [])
    ##
    ## Set the last step
    #if feature_space:
    #    last_step = 'gxg'
    #elif jarpath or toshort:
    #    last_step = 'toshort'
    #elif peakcalling:
    #    last_step = 'macs3'
    #elif make_mcool:
    #    last_step = 'hic'
    #else:
    #    last_step = hic_pipeline[4]   
    last_step = hic_pipeline[4]   
    ##
    ##      6) SUBMITTING TIME STOP COMMANDS 
    ## Submit time stamp 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'timestamp',last_step,stamp,time_partition,group='Experiment',debug=debug) 
    ##
    ##      7) COUNTING COMMANDS
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'count',last_step,stamp,clean_partition,group='Experiment',debug=debug) 

    ##      8) CLEAN UP COMMANDS 
    ## Submit the clean up command if the flag was passed 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'clean','count' if counting else 'timestamp',stamp,clean_partition,group='Experiment',debug=debug) 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        

    ##      SAVING OUT COMMAND FILE 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Write the sbatch commands
    writetofile(f'{debugdir}/sbatch.log.txt',[l+'\n' for l in sub_sbatchs], False)
    ## Calc number of jobs submitted
    njobstosub, njobssubbed = command_files[(command_files.Torun==0)].shape[0], len(sub_sbatchs)
    ## Check our work
    ifprint(f'WARNING: The number of expected jobs to run ({njobstosub}) and number of jobs submitted {njobssubbed} do not match!' ,not (njobstosub == njobssubbed))
    ## Print the number of commands being submitted
    print(f'INFO: A total of {njobssubbed} jobs were submitted with this run' + (f' and will start after {bwaix_jobid}.' if bwaix_jobid else '.')  )
    ## If zero jobs were submitted
    ifprint('WARNING: Zero jobs were submitted; If this was unexpected try running slurpy again, including --restart flag.',(njobssubbed==0))
    ## Save out the command files 
    command_files.to_csv(f'{debugdir}/command.file.csv', index=False)
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## End of file 