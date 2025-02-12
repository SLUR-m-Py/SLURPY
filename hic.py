#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      hic.py; A Hi-C Processing Pipeline 
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
## List slurpy command for VERO analysis
"""
./SLURPY/hic.py -r ../GreenMVA/Chlorocebus_sabeus_mva.fasta -P tb,fast,gpu -G ../Chlorocebus_sabeus_mva.genome.sizes.autosome.filtered.bed -M NC_008066.1 --merge -F 200000 10000000 --restart -a 873337
"""
##      SET PIPELINES
## 
## HIC
## Set pipeline of hic analysis ## NOTE defined after defaults import 
hic_pipeline  = ['fastp', 'bwa', 'filter','dedup','concat','gxg','toshort','hic','clean']
##                  0        1        2       3      4      5a      5b      5c     6
## Join pipeline names by commas
h_pipe = ', '.join(hic_pipeline) 
## 
##      MODULE LOADING and VARIABLE SETTING 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## Load in vars from parameter space
from parameters import * 
## Load in default directories from defaults 
from directories import * 
## Load in ftns from other libraries
from pysamtools import checksam, writetofile
## Bring in bwa mem ftn for hic
from pybwatools import bwamem_hic
## Load in panda cat ftn
from pandacat import pandacat
## Load in bwa master
from bwamaster import bwamaster
## Load in filter master
from filtermaster import filtermaster

## Set the ftn descritption
hiclite_descr = "Processing and analysis pipeline for paired-end sequencing data from Hi-C experiments."
## Set the hardset hic mode 
experi_mode = 'hic'
## Define help messages
R_help = pipe_help%h_pipe
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


##      FUNCTION DEFINING
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## Load in ftns from defaluts
from defaults import basename, getsamplename, reportname, fastcut, fastdry, ifprint

## Ftn for formating fastp command to filter and split reads
def hicpeel(r1:str, r2:str, w:int, s:int, z=4,  options=fastp_opts, script = 'fastp') -> tuple:
    """Formats calls to fastp given input fastq files."""
    ## Format the splits
    split1, split2   = f'{splitsdir}/{basename(r1)}', f'{splitsdir}/{basename(r2)}'  
    ## Set the report name
    report = reportname(getsamplename(r1),script,i=0)
    ## Gather the fastp header, json file, and html 
    fastp_header,fastp_json,fastp_html = fastcut(r1,r2,split1,split2,report,1)
    ## Call fastp again to split the files
    slice = fastp_header + f' -S {s} -z {z} --thread {w} ' + ' '.join(options) + '\n'
    ## Call the fast dry command and format the remove command 
    dry = fastdry(r1,r2,report)
    ## Format and return the command based on the experiment type (ie if it is hi-c or not)
    return [slice,dry,f'mv {fastp_html} {diagdir}/\n',f'mv {fastp_json} {diagdir}/\n'], report 

## Ftn for formating fastq files for bwa mem in hic mode
def prepbwamem(r1:str, r2:str, refix:str, t:int,) -> tuple:
    """Formats a bwa mem command given the paired reads, reference index, thread count, and mode."""
    ## Format the sam file name and report name
    sam, report = f'{bedtmpdir}/{getsamplename(r1)}.sam', f'{debugdir}/1.bwa.{getsamplename(r1)}.log'
    ## Format the bwa mem command 
    makebams = bwamem_hic(r1,r2,refix,sam,report,threads=t)
    ## Return the sam, make bam command and reprot 
    return sam, makebams, report  

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
    parser = argparse.ArgumentParser(description = hiclite_descr)

    ## Add the required argument
    parser.add_argument("-r", "--refix",          dest="r",     type=str,  required=True,  help = r_help, metavar = refmetavar                                     ) 
    parser.add_argument("-F", "--fastp-splits",   dest="F",     nargs='+', required=False, help = F_help, metavar = splitsize,              default = [splitsize]  )
    parser.add_argument("-B", "--parallel-bwa",   dest="B",     type=int,  required=False, help = B_help, metavar = parallelbwa,            default = parallelbwa  )
    parser.add_argument("-P", "--partition",      dest="P",     type=str,  required=False, help = P_help, metavar = part,                   default = part         ) 
    parser.add_argument("-M", "--mtDNA",          dest="M",     type=str,  required=False, help = M_help, metavar = mito,                   default = mito         )
    parser.add_argument("-X", "--exclude",        dest="X",     nargs='+', required=False, help = X_help, metavar = 'chrX, chrY ...',       default = []           )
    parser.add_argument("-Q", "--map-threshold",  dest="Q",     type=int,  required=False, help = Q_help, metavar = map_q_thres,            default = map_q_thres  )
    parser.add_argument("-R", "--rerun-from",     dest="R",     type=str,  required=False, help = R_help, metavar = 'step',                 default = None         )
    parser.add_argument("-q", "--fastq",          dest="q",     type=str,  required=False, help = q_help, metavar = '.fastq.gz',            default = fends        )
    parser.add_argument("-a", "--afterok",        dest="a",     type=int,  required=False, help = a_help, metavar = fakejobid,              default = 0            )
    parser.add_argument("-N", "--nice",           dest="N",     type=int,  required=False, help = N_help, metavar = 'n',                    default = nice         )

    ## Set number of threads 
    parser.add_argument("-f", "--fastp-threads",  dest="f",     type=int,  required=False, help = f_help, metavar = fastpthreads,           default = fastpthreads )
    parser.add_argument("-b", "--bwa-threads",    dest="b",     type=int,  required=False, help = b_help, metavar = bwathreads,             default = bwathreads   )
    parser.add_argument("-t", "--dask-threads",   dest="t",     type=int,  required=False, help = t_help, metavar = daskthreads,            default = daskthreads  )

    ## Set values for Hi-C analysis 
    parser.add_argument("-n", "--run-name",       dest="n",     type=str,  required=False, help = n_help, metavar = 'name',                 default = None         )
    parser.add_argument("-E", "--error-distance", dest="E",     type=int,  required=False, help = E_help, metavar = 'bp',                   default = error_dist   )
    parser.add_argument("-L", "--library",        dest="L",     type=str,  required=False, help = L_help, metavar = 'MboI',                 default = lib_default  )
    parser.add_argument("-Z", "--chunksize",      dest="Z",     type=int,  required=False, help = Z_help, metavar = 'n',                    default = chunksize    )
    parser.add_argument("-G", "--genomelist",     dest="G",     type=str,  required=False, help = G_help, metavar = './path/to/list.tsv',   default = False        )
    parser.add_argument("-J", "--jar-path",       dest="J",     type=str,  required=False, help = J_help, metavar = './path/to/juicer.jar', default = None         )
    parser.add_argument("-x", "--Xmemory",        dest="x",     type=int,  required=False, help = x_help, metavar = xmemory,                default = xmemory      )
    parser.add_argument("-S", "--bin-sizes",      dest="S",     nargs='+', required=False, help = S_help, metavar = '25000, 10000, ...',    default = binsizes     )
    parser.add_argument("-gxg","--features",      dest="gxg",   type=str,  required=False, help = A_help, metavar= './path/to/my.gff',      default = 'none'       )

    ## Set boolean flags 
    parser.add_argument("--toshort",              dest="toshort",   help = short_help,    action = 'store_true')
    parser.add_argument("--restart",              dest="restart",   help = restart_help,  action = 'store_true')
    parser.add_argument("--force",                dest="force",     help = force_help,    action = 'store_true')
    parser.add_argument("--debug",                dest="debug",     help = debug_help,    action = 'store_true')
    parser.add_argument("--skipdedup",            dest="skipdedup", help = mark_help,     action = 'store_true')
    parser.add_argument("--clean",                dest="clean",     help = clean_help,    action = 'store_true')
    parser.add_argument("--merge",                dest="merge",     help = merge_help,    action = 'store_true')

    ## Set the paresed values as inputs
    inputs = parser.parse_args() 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
   

    ##      PASS ARGUMENTS & SET VARIABLES 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Set required variables
    reference_path  = inputs.r            ##     Set path to the reference genome

    ## Set default vairables              ##
    splitsize       = inputs.F            ##     Number of splits in fastp 
    bwa_runs        = inputs.B            ##     Set the number of parallel runs of bwa 
    partition       = inputs.P            ##     Set the partition 
    mito            = inputs.M            ##     Set the mito contig name 
    excludes        = inputs.X            ##     List of chromosomes to exclude from analysis 
    mapq            = inputs.Q            ##     Set the mapping quality threshold 
    rerun           = inputs.R            ##     Setp to rerun pipeline from 
    fend            = inputs.q            ##     End of the input fastq files
    bwaix_jobid     = inputs.a            ##     The job id to have all submissions wait on   
    nice            = inputs.N            ##     Sets the nice parameter 

    ## Set threads                           
    fastp_threads   = inputs.f            ##     Number of fastp threads
    bwa_threads     = inputs.b            ##     Number of threads in bwa alignments
    daskthreads     = inputs.t            ##     Number of threads in dask

    ## Set variables for Hi-C                
    run_name        = inputs.n            ##     The name of the samples 
    error_dist      = inputs.E            ##     Distance in bp to check hi-C contacts for errors 
    enzymelib       = inputs.L            ##     Restriction-enzymelib used in Hi-C prep 
    chunksize       = inputs.Z            ##     Chunk size (row number) to load in with pandas 
    pathtochrom     = inputs.G            ##     Path to list of chromosomes to use 
    jarpath         = inputs.J            ##     Path to juicer jar file 
    xmemory         = inputs.x            ##     Amount of memory passed to juicer pre command 
    binsizes        = inputs.S            ##     Bins / resolutions used in hi-c analysis 
    feature_space   = inputs.gxg          ##     Path to a gff or bed file used in g x g interaction matrix / df 
                                            
    ## Set boolean vars                    
    toshort         = inputs.toshort      ##     Flag the make short file, kicks if jarpath was given 
    hardreset       = inputs.restart      ##     Resetart the slurpy run, removing previous
    force           = inputs.force        ##     Force overwrite of output alignment files 
    debug           = inputs.debug        ##     Run in debug mode 
    skipduplicates  = inputs.skipdedup    ##     Boolean to mark duplicates 
    ifclean         = inputs.clean        ##     Flag to run clean up script 
    postmerging     = inputs.merge        ##     Forces premerge of outputs 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    
    
    ##      ROTH SETTINGS
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Reset reference with presets if we are running as Cullen Roth (These can be edited if you are not me :-))
    reference_path = t2t_refpath if (reference_path.lower() == 't2t') else reference_path
    feature_space  = t2t_gffpath if (feature_space.lower()  == 't2t') else feature_space
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
    ## Set the hic file ends
    hicfileends = [mito] + hicfileends_tmp
    ## Initilizse list for sbatch
    sub_sbatchs = []

    ## Set the lib error
    lib_error = "ERROR: The passed library name of enzyme(s) used in Hi-C prep was not recognized."
    ## Check library if it was pased
    assert enzymelib.lower() in ['mboi', 'dpnii','sau3ai','hindiii','arima','none'], lib_error
    print("INFO: Running analysis with the %s Hi-C library."%enzymelib)

    ## Check if we have a jarpath
    if jarpath:
        assert fileexists(jarpath), "ERROR: The given jarpath %s could not be found! Please check path the given path and try again."%jarpath
        toshort = True 

    ## Set feature space to a boolean var, needed to make my settings work 
    if feature_space == 'none':
        feature_space = False
    ## Check is a path the give feature space 
    if feature_space:
        assert fileexists(feature_space), "ERROR: The given features data %s could not be found! Please check the given path and try agian."%feature_space
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      CONFIRM THE HARD RESTART 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Load in defluats 
    from defaults import confirmreset
    ## Format group dirs 
    grouped_dirs = [debugdir,aligndir,splitsdir,comsdir,diagdir,bedtmpdir,hicdir]
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
    chrlist,genome_size,pathtochrom = gathering(reference_path,pathtochrom,excludes)
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
    command_files, samplenames, filteredbams = [], [], []

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
    in_fastqs = getfastqs(fastqdir+'/*.gz')
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
        fastp_coms, fastp_repo = hicpeel(r1,r2,fastp_threads,sizeosplit)
        ##  Write the command to file
        writetofile(fastp_command_file, sbatch(None,fastp_threads,the_cwd,fastp_repo,nice=1) + fastp_coms, debug)
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
        bwa_master_commands, bwa_master_repo = bwamaster(sample_name,reference_path,bwa_threads,the_cwd,partition,debug,nice,library=enzymelib,inhic=True,forced=force)
        ## Write command to file
        writetofile(bwa_master_file, sbatch('bwa.master',1,the_cwd,bwa_master_repo,nice=nice) + bwa_master_commands, debug)
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
        filter_master_commands, filter_master_repo = filtermaster(sample_name,reference_path,the_cwd,excludes,chrlist,mapq,error_dist,daskthreads,enzymelib,partition,True,debug,nice,forced=force)
        ## Write command to file
        writetofile(filter_master_file, sbatch('filter.bedpe.master',1,the_cwd,filter_master_repo,nice=nice) + filter_master_commands, debug)
        ## Append to command file
        command_files.append((filter_master_file,sample_name,experi_mode,hic_pipeline[pix],filter_master_repo,0,''))
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

  
    ##      FILE and TXT MERGING 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## 3. merging,  deduplicating and sorintg
    pix = 3
    ## Iterate over the sample names
    for si,sname in enumerate(samplenames):
        ## Set the smaple name 
        sample_name = run_name if postmerging else sname 
        ## If the sample is the first and we are merging all smaples 
        if (si > 0) and postmerging:
            break 
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## Initiate list of out put hic files for each chromosome
        hiccats_outs = []
        hicdups_outs = []
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
            hiccat_file = f'{comsdir}/{pix}.hiccat.{sample_name}.valid.{chrom}.sh'
            ## Write the concat command to file
            writetofile(hiccat_file, sbatch(hiccat_file,daskthreads,the_cwd,hiccat_repo,nice=nice) + hiccat_coms, debug)
            ## Append the concat command
            command_files.append((hiccat_file,sample_name,experi_mode,hic_pipeline[pix],hiccat_repo,0,''))    

        ## Set start name for wildcard use to bring in inputs  to ftn 
        sample_start = 'notused.bedpe' if postmerging else f'{sample_name}.notused.bedpe'
        hiccat_out = f'{aligndir}/{sample_name}.notused.bedpe'
        ## SEt the report name
        hiccat_repo = reportname(sample_start,hic_pipeline[pix],i=pix)
        ## Set the command
        hiccat_coms = [f'{slurpydir}/deduphic.py -b {sample_start} -o {hiccat_out}\n']
        
        ## make concat file name
        hiccat_file = f'{comsdir}/{pix}.hiccat.{sample_name}.notused.sh'
        ## Write the concat command to file
        writetofile(hiccat_file, sbatch(hiccat_file,daskthreads,the_cwd,hiccat_repo,nice=nice) + hiccat_coms, debug)
        ## Append the concat command
        command_files.append((hiccat_file,sample_name,experi_mode,hic_pipeline[pix],hiccat_repo,0,''))
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

        
        ## Hi-C CONCATONATION
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## 4B. Concat the duplicate Hi-C contacts from valid files, set pipeline. NOTE we reversed the order there 
        pix = 4
        ## Format the infilename, we need this last filename in the loop for the hic
        newcatfile = f'{aligndir}/{sample_name}.duplicates.bedpe'
        ## Format report for pandacat
        concat_repo = reportname(newcatfile,hic_pipeline[pix],i=f'{pix}B')
        ## Format command to merge the chromosome bedpe files to each, the chromosome order is important 
        concat_coms = pandacat(hicdups_outs,newcatfile,report=concat_repo,rmheader=True)
        ## make concat file name
        concat_file = f'{comsdir}/{pix}.concat.{sample_name}.duplicates.sh'
        ## Write the concat command to file
        writetofile(concat_file, sbatch(concat_file,1,the_cwd,concat_repo,nice=nice) + concat_coms, debug)
        ## Append the concat command
        command_files.append((concat_file,sample_name,experi_mode,hic_pipeline[pix],concat_repo,0,''))  
        
        ## 4A. Concat the input Hi-C contacts from valid files, set pipeline 
        ## Format the infilename, we need this last filename in the loop for the hic
        newcatfile = f'{aligndir}/{sample_name}.valid.bedpe'
        ## Format report for pandacat
        concat_repo = reportname(newcatfile,hic_pipeline[pix],i=f'{pix}A')
        ## Format command to merge the chromosome bedpe files to each, the chromosome order is important 
        concat_coms = pandacat(hiccats_outs,newcatfile,report=concat_repo,rmheader=True)
        ## make concat file name
        concat_file = f'{comsdir}/{pix}.concat.{sample_name}.valid.sh'
        ## Write the concat command to file
        writetofile(concat_file, sbatch(concat_file,1,the_cwd,concat_repo,nice=nice) + concat_coms, debug)
        ## Append the concat command
        command_files.append((concat_file,sample_name,experi_mode,hic_pipeline[pix],concat_repo,0,'')) 
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


        ## OPTIONAL ANALYSES
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        """
        Optional analyzes
        gene x gene counts
        bedpe ---> short format
        jucier pre command for hic creation 
        """
        ## GENE X GENE interaction 
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        pix = 5
        ## 5A. Create a feature space 
        if feature_space:
            ## Iterate over chromosome list 
            for i,coi in enumerate(chrlist):
                ## Set the report, commands, and gxg script file name 
                gxg_repo   = reportname(sample_start,'gxg%s'%i,i=f'{pix}A')
                gxg_commands = [f'{slurpydir}/gxgcounts.py -i {newcatfile} -c {coi} -f {feature_space} -t {nchrom}' + (' --merge\n' if not i else '\n')]
                gxg_file     = f'{comsdir}/{pix}A.gxg.{i}.{sample_name}.sh'
                ## Write to file for the gxg script, pasing dask thread count, the cwd, commands and debug mode 
                writetofile(gxg_file,sbatch(gxg_file,daskthreads,the_cwd,gxg_repo,nice=nice) + gxg_commands, debug)
                ## Append to command list for calling/submission later 
                command_files.append((gxg_file,sample_name,experi_mode,'gxg',gxg_repo,0,''))
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


        ## BEDPE TO SHORT
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
            writetofile(short_file,sbatch(short_file,daskthreads,the_cwd,short_repo,nice=nice) + short_commands, debug)
            ## append the short command
            command_files.append((short_file,sample_name,experi_mode,'toshort',short_repo,0,''))
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
            writetofile(jpre_file, sbatch(jpre_file,fastp_threads,the_cwd,jpre_repo,nice=nice) + jpre_coms, debug)
            ## Append the concat command
            command_files.append((jpre_file,sample_name,experi_mode,'hic',jpre_repo,0,''))
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      SUBMITTING TIME-STAMP   
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## 6A. Final timestamp, out clean up
    pix = 6
    ## Set the timesampe assocaited file names 
    timestamp_file   = f'{diagdir}/{run_name}.timestamp.{stamp}.txt'             ##     Name of the output file 
    timestampsh      = f'{comsdir}/{pix}A.time.stamp.sh'                         ##     Name of the .sh bash file 
    timestamp_repo = reportname(run_name,f'timestamp.{stamp}',i=f'{pix}A')       ##     Name of the log to report to 
    ## Formath time stamp and echo commands 
    times_commands = [f'{slurpydir}/endstamp.py {timestamp_file} {stamp}\n']
    ## Format the command file name and write to sbatch, we will always ask the timestamp to run even in debug mode 
    writetofile(timestampsh, sbatch(timestampsh,1,the_cwd,timestamp_repo,nice=1) + times_commands, False)
    ## Append the timestamp command to file
    command_files.append((timestampsh,run_name,experi_mode,'timestamp',timestamp_repo,0,''))
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      CLEANING UP FILES   
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## 6B. If the clean boolean or rerun vars were set 
    if ifclean or (rerun == 'clean'): 
        ## Format command to remove uneedeed files 
        remove_sh   = f'{comsdir}/{pix}B.cleanup.sh'             ##   Set the bash file name 
        remove_repo = reportname(run_name,'clean',i=f'{pix}B')   ##   Set the report 
        remove_comm = [f'{slurpydir}/remove.py {bedtmpdir} {splitsdir} {hicdir}\n', f'{slurpydir}/gzipy.py ./{aligndir}/*.bedpe\n', f'{slurpydir}/gzipy.py ./{aligndir}/*.short\n']
        ## Format the command to clean up          
        writetofile(remove_sh, sbatch(remove_sh,1,the_cwd,remove_repo,nice=nice) + remove_comm, debug)
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
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,hic_pipeline[1],hic_pipeline[0],stamp,partition,debug=debug)

    ##      2) SUBMITTING BEDPE FILTERING and SPLITTING
    ## Submit the bedpe filtering scirpt 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,hic_pipeline[2],hic_pipeline[1],stamp,partition,debug=debug)

    ##      3) SUBMITTING DEDUPLICATING AND SORTING 
    ## Submit the dedup and sort script 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,hic_pipeline[3],hic_pipeline[2],stamp,partition,debug=debug,group='Experiment' if postmerging else 'Sample')

    ##      4) SUBMITTING CONCATONATION 
    ## Submit the concatonation sciprt 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,hic_pipeline[4],hic_pipeline[3],stamp,partition,debug=debug,group='Experiment' if postmerging else 'Sample')
    
    ##      5A) SUBMITTING Gene X Gene interaction script 
    ## Submit the g x g script
    sub_sbatchs = sub_sbatchs + (submitdependency(command_files,'gxg',hic_pipeline[4],stamp,partition,debug=debug,group='Experiment' if postmerging else 'Sample') if feature_space else []) 
    ## Set the next step in pipeline 
    last_step = 'gxg' if jarpath else hic_pipeline[4]

    ##      5B) SUBMITTING CONVERSION FROM BEDPE to JUICER SHORT 
    ## Submit the toshort.py command 
    sub_sbatchs = sub_sbatchs + (submitdependency(command_files,'toshort',hic_pipeline[4],stamp,partition,debug=debug,group='Experiment' if postmerging else 'Sample') if toshort else []) 
    ## Set the next step in pipeline 
    last_step = 'toshort' if jarpath else hic_pipeline[4]

    ##      5C) Hi-C FILE CREATION 
    ## Call the juicer pre command for hic file creation if jarpath was passed 
    sub_sbatchs = sub_sbatchs + (submitdependency(command_files,'hic','toshort',stamp,partition,debug=debug,group='Experiment' if postmerging else 'Sample') if jarpath else [])
    ## Set the next step in pipeline 
    last_step = 'hic' if jarpath else hic_pipeline[4]

    ##      6) SUBMITTING TIME STOP COMMANDS 
    ## Submit time stamp 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'timestamp',last_step,stamp,partition,group='Experiment',debug=debug) 

    ##      6) CLEAN UP COMMANDS 
    ## Submit the clean up command if the flag was passed 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'clean','timestamp',stamp,partition,group='Experiment',debug=debug) 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        

    ##      SAVING OUT COMMAND FILE 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Write the sbatch commands
    writetofile(f'{debugdir}/sbatch.log.txt',[l+'\n' for l in sub_sbatchs], False)
    ## Calc number of jobs submitted
    njobstosub, njobssubbed = command_files[(command_files.Torun==0)].shape[0], len(sub_sbatchs)
    ## Check our work
    ifprint(f'WARNING: The number of expected jobs to run ({njobstosub}) and number of jobs submitted {njobssubbed} does not match!' ,not (njobstosub == njobssubbed))
    ## Print the number of commands being submitted
    print(f'INFO: A total of {njobssubbed} jobs were submitted with this run' + (f' and will start after {bwaix_jobid}.' if bwaix_jobid else '.')  )
    ## If zero jobs were submitted
    ifprint('WARNING: Zero jobs were submitted; If this was unexpected try running slurpy again, including --restart flag.',(njobssubbed==0))
    ## Save out the command files 
    command_files.to_csv(f'{debugdir}/command.file.csv', index=False)
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## End of file 