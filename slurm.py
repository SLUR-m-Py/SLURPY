#!/usr/bin/env python
## Load in parase args ftn, optional variables, and reference paths 
from parameters import slurpy_args, hic_options, hic_pipeline, samthreads, t2t_refpath, t2t_gtfpath, vero_refpath
## Load in directoris
from parameters import fastqdir, splitsdir, diagdir, slurpydir, debugdir, aligndir, comsdir, bedtmpdir, hicdir, macs3dir
## Load in ftns from bio tools for formating calls to fastp
from biotools import fastcut, reportname, fastdry, getsamplename, sortfastq
## Load in ftn from bio tools for writing files and submittting jobs 
from biotools import commandcontrol, submitdependency, submitfastp, writeparams
## Load in other ftn from bio tools 
from biotools import peakattack, isbwaix, sortglob, setrunname, gathering
## Load in ftns from biotools for formating paths, directories, and reseting
from biotools import makedirectories, checksam, patchclean, confirmreset
## Load in panda cat 
from pandacat import pandacat
## Load in bwa master
from bwatobedpe import bwamaster
## Load in filter master
from filtermaster import filtermaster
## Load in bedpe to sam
from pairs2sam import bedpetosam
## Load in ftn from os path
from os.path import exists, basename, isdir
## laod in os
from os import remove, getcwd
## Bring in defautls
from defaults import sbatch, ifprint, basenoext
## Load in write to file
from myecho import writetofile
## Load in time
from time import time
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
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

##      FUNCTION DEFINING
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## Define options for fastpeel ftn
fastp_opts = ['--dont_eval_duplication','--disable_length_filtering','--disable_adapter_trimming','--disable_quality_filtering','--disable_trim_poly_g']

## Ftn for formating fastp command to filter and split reads
def fastpeel(r1:str, r2:str, w:int, s:int, ishic=True, pix=0, options:str=fastp_opts, toremoveend='.toberemoved.fastq.gz', singleend='.singletons.fastq.gz', failend='.failed.fastq.gz',z=4, fastp_splits_dir=splitsdir,diagnos=diagdir) -> tuple:
    """Formats calls to fastp given input fastq files."""
    ## Set variable names 
    r1_bn    = basenoext(r1)                                 ##      Generate basename of first read 
    r2_bn    = basenoext(r2)                                 ##      Generate basename of second read 
    temp1    = f'{fastp_splits_dir}/{r1_bn}{toremoveend}'           ##      Format the first temp file for splits to be remove
    temp2    = f'{fastp_splits_dir}/{r2_bn}{toremoveend}'           ##      Format the second temp file for splits to be remove 
    single1  = f'{fastp_splits_dir}/{r1_bn}{singleend}'             ##      Make singleton file from first split 
    single2  = f'{fastp_splits_dir}/{r2_bn}{singleend}'             ##      Make singleton file from second split 
    split1   = f'{fastp_splits_dir}/{basename(r1)}'                 ##      Set the first split name
    split2   = f'{fastp_splits_dir}/{basename(r2)}'                 ##      Create the second split name 
    report   = reportname(getsamplename(r1),'fastp',i=pix)   ##      Format the report name
    failed   = f'{fastp_splits_dir}/{getsamplename(r1)}{failend}'   ##      Generate the list of the failed reads
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
    mv_fastp = [f'mkdir -p {diagnos}/fastp\n',f'mv {fastp_json0} {diagnos}/fastp/\n',f'mv {fastp_html0} {diagnos}/fastp/\n']

    ## Return the comands depenent aupon if we are spliting or filtering, and the reprot 
    if ishic:
        return [just_split,dry] + mv_fastp, report
    else:
        return [initial_filtering,split_filtered,throwout,dry] + mv_fastp, report

## Ftn for calling juicer's pre command
def juicerpre(intxt:str, outhic:str, Xmemory:int, jarfile:str, threadcount:int, bins:list, genomepath:str,pix:str, exdir=slurpydir) -> tuple:
    """
    java -Xmx49152m -Xms49152m -jar $jarpath pre -j 5 -r 500000,250000,200000,150000,100000,50000,25000,10000,5000,1000 $1 ${1}.hic $2
    """
    ## Format the report name 
    report = reportname(outhic+'.bam','hic.pre',i=pix)
    ## Set the java command for the passed juicer jar file 
    prestr = ['java -Xmx%sm -Xms%sm -jar %s pre -j %s -r %s %s %s %s\n'%(Xmemory,Xmemory,jarfile,threadcount,','.join(map(str,bins)),intxt,outhic,genomepath),
              f'{exdir}/myecho.py Finished formating Hi-C contacts into a .hic file on path: {outhic} {report}\n']
    ## Return the pre and report
    return prestr, report
    
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      MAIN SCRIPT & ARGUMENT PARSING 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
def main(executive_dir:str=slurpydir,
         fastqs_dir:str=fastqdir,
         fastpsplits_dir:str=splitsdir,
         alignment_dir:str=aligndir, 
         commands_dir:str=comsdir,
         logs_dir:str=debugdir,
         diagnostics_dir:str=diagdir,
         bedpe_dir:str=bedtmpdir,
         hic_dir:str=hicdir,
         peaks_dir:str=macs3dir,
         bwa_options:str=hic_options,
         samtools_ncpu:int=samthreads,
         human_reference:str=t2t_refpath,
         human_gtf:str=t2t_gtfpath,
         vero_reference:str=vero_refpath,
         pipeline_steps:list=hic_pipeline,
         inhic:bool=False,
         peakcalling:bool=False,
         count_mode:str='--atac-seq'):

    ## Set inputs from parser
    inputs = slurpy_args()
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
   
    ##      PASS ARGUMENTS & SET VARIABLES 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Set required variables
    reference_path  = inputs.r              ##     Set path to the reference genome

    ## Set default vairables              
    path_to_fastqs  = inputs.fq             ##     Set the path to the input fastqs, if non was passed, assumes they are at the level of the call to slurpy
    splitsize       = inputs.F              ##     Number of splits in fastp, this is a line count, multiple by four to get actual read counts
    threadn         = inputs.T              ##     Set the number of parallel runs of bwa 
    partitions      = inputs.P              ##     Set the partition 
    mito            = inputs.M              ##     Set the mito contig name 
    excludes        = inputs.X              ##     List of chromosomes to exclude from analysis 
    mapq            = inputs.Q              ##     Set the mapping quality threshold 
    rerun           = inputs.R              ##     Setp to rerun pipeline from 
    bwaix_jobid     = inputs.a              ##     The job id to have all submissions wait on   
    nice            = inputs.N              ##     Sets the nice parameter 

    ## Set threads across softwares                         
    nparallel       = inputs.j              ##     Number of parallel jobs to run at once acorss bwa and filter 
    fastp_threads   = inputs.f              ##     Number of fastp threads
    daskthreads     = inputs.t              ##     Number of threads in dask
    bwa_threads     = inputs.b              ##     Number of threads in bwa alignments
    bwa_opts        = inputs.B              ##     Set options for bwa 

    ## Set variables for Hi-C                
    run_name        = inputs.n              ##     The name of the samples 
    error_dist      = inputs.E              ##     Distance in bp to check hi-C contacts for errors 
    enzymelib       = inputs.L              ##     Restriction-enzyme library used in Hi-C prep 
    chunksize       = inputs.Z              ##     Chunk size (row number) to load in with pandas 
    pathtochrom     = inputs.G              ##     Path to list of chromosomes to use 
    jarpath         = inputs.J              ##     Path to juicer jar file 
    xmemory         = inputs.xmx            ##     Amount of memory passed to juicer pre command 
    binsizes        = inputs.S              ##     Bins / resolutions used in hi-c analysis 
    max_dist        = inputs.m              ##     Maximum distance of paired end reads
    feature_space   = 'none'                ##     Path to a gff or bed file used in g x g interaction matrix / df 
    nodes           = inputs.nodes          ##     List of nodes 
    slurm_mem       = inputs.slurmem        ##     Amount of memory setting in SLURM
    
    ## Set pipeline boolean vars                    
    hardreset       = inputs.restart        ##     Resetart the slurpy run, removing previous
    postmerging     = not inputs.merge      ##     Forces premerge of outputs 
    force           = inputs.force          ##     Force overwrite of output alignment files 
    debug           = inputs.debug          ##     Run in debug mode 
    ifclean         = inputs.clean          ##     Flag to run clean up script 
    counting        = not inputs.count      ##     Flag to count the read pairs

    ## Set Hi-C related boolean flasgs                  
    toshort         = inputs.toshort        ##     Flag the make short file, kicks if jarpath was given 
    makepairs       = inputs.makepairs      ##     Flag to make pairs file 
    make_mcool      = inputs.mcool          ##     Flag to make mcool file  
    hicexplorer     = inputs.hicexp         ##     Flag to run stricter filtering like Hi-C explorer
    get_inter       = inputs.inter          ##     Flag to return inter chromosome counts in bedpe file
    
    ## MACS3 related vars
    atac_seq        = inputs.atac           ##     Boolean flag to run in atac-seq mode 
    ifbroad         = inputs.broad          ##     Boolean to activate broader peak calling in macs3 
    skippeaks       = inputs.peaks          ##     Skips peak calling with macs3 for atac-seq experiments
    nolambda        = inputs.nolambda       ##     Flag to skip local lambda control in macs3
    nomodel         = inputs.nomodel        ##     Flag to skip shifting model
    callsummits     = inputs.summits        ##     Flag to call summits in macs3
    shift_size      = inputs.shiftsize      ##     Value of shift size (bp)
    extendsize      = inputs.extendsize     ##     Value of read extension 
    macs3mode       = inputs.macmode        ##     Format in MACS3, eg. bedpe.
    chip_control    = inputs.c              ##     Set the input controls for chip experiments
    max_gap         = inputs.maxgap         ##     Max gap in used in MACS3    
    min_len         = inputs.minlen         ##     Min lenght of peaks in MACS3
    max_nchrom      = inputs.maxnc          ##     Max number of chromosomes 

    ## Set RNA-seq like vars 
    rna_seq         = inputs.rnas           ##     Boolean flag to run in rna-seq mode 
    wgs_seq         = inputs.wgs            ##     Boolean flag to run in wgs seq mode
    sfastp          = inputs.sfast          ##     Flag to skip fastp filtering 
    skipduplicates  = inputs.skipdedup      ##     Boolean to mark duplicates       
    keep_dups       = inputs.save           ##     Boolean flag to save out duplicates once identified
    dedovetail      = inputs.tails          ##     Boolean for dove tailing 
    tosam           = inputs.tosam          ##     Boolean flag to convert .bedpe file to sam format
    tobam           = inputs.tobam          ##     Boolean flag to convert .bedpe file to bam format

    ## Set the time stamp
    stamp = time()
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

    ##      CORRECT SPLIT / Chunksize
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    max_splitsize = max([int(s) for s in splitsize])
    chunk_check   = (chunksize > max_splitsize)
    chunksize = int(round(max_splitsize/3,0)) if chunk_check else chunksize
    ifprint("INFO: Reseting chunksize (Z) to: %s"%chunksize,chunk_check)

    ##      PRESET Thread counts
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Check if we are in atac seq mode    
    if threadn:
        ## reset threads                           
        fastp_threads   = threadn           ##     Number of fastp threads
        bwa_threads     = threadn           ##     Number of threads in bwa alignments
        daskthreads     = threadn           ##     Number of threads in dask 

    ## Reset paortions into a comma joined list
    new_parts = []
    for part in partitions:
        for p in part.split():
            for c in p.split(','):
                if len(c):
                    new_parts.append(c)

    ## Make a single commo sepereated string
    partition = ','.join(new_parts)


    ##      ROTH SETTINGS
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Reset reference with presets if we are running as Cullen Roth (These can be edited if you are not me :-))
    ist2t  = (reference_path.lower() == 't2t')
    isvero = (reference_path.lower() == 'vero')

    ## If human reference 
    reference_path = human_reference if ist2t else reference_path
    feature_space  = human_gtf if (feature_space.lower()  == 't2t') else feature_space
    ## If vero was pass
    reference_path = vero_reference if isvero else reference_path
    #print(reference_path)

    ## Set bwa master partition to mpi (FOR ROTH only)
    bwa_partition   = 'gpu,fast,mpi,tb' if ist2t or isvero else partition
    filt_partition  = bwa_partition     if ist2t or isvero else partition
    clean_partition = bwa_partition     if ist2t or isvero else partition
    time_partition  = bwa_partition     if ist2t or isvero else partition
    dedup_partition = bwa_partition     if ist2t or isvero else partition

    ## Set nodelists to run bwa and filter master
    bwanodes = ['c0826', 'c1002', 'c0701', 'c0702']
    filternodes = bwanodes
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

    ##      INITILIZATION 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Modify the fastq dir
    if isdir(path_to_fastqs):
        ## If it isn't a relaive path and the last character is a forward slash
        if path_to_fastqs != './' and (path_to_fastqs[-1] == '/'):
            ## Remove the extra forward slash
            path_to_fastqs = path_to_fastqs[:-1]
        ## Set the path to the fastqs
        fastqs_dir = path_to_fastqs
    else:
        print("WANRING: Unable to find given fastq directory: %s\n\tDefaulting to %s"%(path_to_fastqs,fastqs_dir))
    
    ## Check that the fastq and reference path exists
    assert isdir(fastqs_dir), "ERROR: Unable to detect a fastqs directory holding fastq files!"
    assert exists(reference_path), "ERROR: We could not find the provided input path -- %s -- of the reference file!"%reference_path 
    ## Check the versions of samtools, the user email is an email and the experiment mode is one we know
    assert checksam(), "ERROR: The detected version of samtools is not greater than or equal to v 1.15.1!\nPlease update samtools and try again."
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
        assert exists(jarpath), "ERROR: The given jarpath %s could not be found! Please check path the given path and try again."%jarpath
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
        assert exists(feature_space), "ERROR: The given features data %s could not be found! Please check the given path and try agian."%feature_space
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      PRESET ATAC and ChIP-seq mode
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Check if we are mapping RNA-seq reads
    if rna_seq:
        experi_mode   = 'rna-seq'
        bwa_opts      = ',-M' if (bwa_opts == bwa_options) else bwa_opts
        print("INFO: Processing RNA-seq sample.")

    ## Check if we are in atac seq mode    
    elif atac_seq or len(chip_control) or ifbroad:
        ## Check our control files if they were passed 
        if len(chip_control):
            experi_mode = 'chrp-seq'
            assert (type(chip_control) == list) and (type(chip_control[0]) == str), 'ERROR: Inputs for chip control are not a type we recognize (should be a list of strings)!'
            ## Iterate over the controls inputs 
            for chip_con in chip_control:
                ## Check they exist 
                assert exists(chip_con), f'ERROR: Bad path to bam file! Unable to locate input control for chip experiment: {chip_con}'
            ## Otherwise it is an atac-seq 
        else:
            experi_mode = 'atac-seq'
        ## SEt hic vars to false and none 
        max_dist      = max_dist if max_dist else 1000
        bwa_opts      = ',-M' if (bwa_opts == bwa_options) else bwa_opts
        peakcalling   = not skippeaks
        print("INFO: Processing ATAC-seq sample, calling peaks with MACS3.")
        
    ## set wgs mode 
    elif wgs_seq:
        bwa_opts      = ',-M' if (bwa_opts == bwa_options) else bwa_opts
        print("INFO: Processing sample from whole-genome sequencing.")

    else:
        ## It is a Hi-C experiment
        experi_mode   = 'hic'
        inhic         = True 
        count_mode    = ''
        ## Moved to trigger if statment correctly
        print("INFO: Running analysis with the %s Hi-C library."%enzymelib)

    ## Redefine enzymatic lib
    enzymelib = 'none' if not inhic else enzymelib

    ## Reset shift and extention size 
    shift_size,extendsize = (shift_size,extendsize) if nomodel else (0,0)
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

    ##      CONFIRM THE HARD RESTART 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Format group dirs 
    grouped_dirs = [logs_dir,alignment_dir,fastpsplits_dir,commands_dir,diagnostics_dir,bedpe_dir,hic_dir] + ([peaks_dir] if peakcalling else [])
    ## If there is a hard reset passed 
    confirmreset(grouped_dirs) if hardreset else None
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

    ##      DIRECTORY MAKING & TIME STAMP SUBMISSION 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Let the user know we are making directories 
    print('INFO: Making local directories.')
    ## Get the current working dir
    the_cwd = getcwd()
    ## If the run name is none
    run_name = setrunname(run_name,the_cwd)
    ## Make the directories 
    makedirectories(grouped_dirs)
    ## Gather and remove the old command files, error, and output logs if they are there 
    [remove(f) for f in sortglob(f'{commands_dir}/*.sh') + sortglob(f'{logs_dir}/*.err') + sortglob(f'{logs_dir}/*.out')]
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      WRITING OUT PARAMS
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Write out params 
    writeparams(f'{experi_mode}.py',run_name,stamp,inputs)
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      CHROMOSOME GATHERING 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Expand exlcude list to include mitochondria contig
    excludes.append(mito)
    ## Gather a list of chromosomes 
    chrlist,genome_size,pathtochrom = gathering(reference_path,pathtochrom,excludes,not make_mcool)
    ## Calculate chrlist
    nchrom = len(chrlist)
    ## Set error messave
    nchrom_err = "ERROR: The number of chromosmes (%s) in the passed reference file are larger than the allowed number (%s)! Consider passing the -G argument or using a more complete reference. You can also increase the --max-number-chroms argument."
    ## Set chromosome cut off
    assert nchrom < max_nchrom, nchrom_err%(nchrom,max_nchrom)
    ## Print if there are fewer than ten
    print('INFO: Running on -> %s'%' '.join(chrlist)) if (nchrom < 10) else None
    print('INFO: Processing %s chromosomes from list.'%nchrom)
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      CHECK FOR A BWA INDEXING  
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Set index error
    index_error = "ERROR: We could not detect an index associated with the input reference path: %s\nINFO: Index the reference (bwa index) and try again."
    ## Assert our truth
    assert isbwaix(reference_path),index_error%reference_path
    
    ## Format the command file and initilize sample names
    command_files, samplenames = [],[]

    ## Append to command
    command_files.append(('jobfile','sample',experi_mode,'bwaix','report',0,bwaix_jobid)) if bwaix_jobid else None 
    #print(command_files)
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      FASTQ GATHERING
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Sort by fastq size
    in_fastqs = sortfastq(fastqs_dir+'/*.fastq.gz',splitsize)
    ## Saveout the sizes to debug dir 
    in_fastqs.to_csv(f'{logs_dir}/fastq.sizes.csv',header=True,index=False)
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      SAMPLE NAME MAKING AND HANDELING
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
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
        fastp_command_file =  f'{commands_dir}/{pix}.fastp.{sample_name}.sh'  
        ##  Gather the fastp command and report 
        fastp_coms, fastp_repo = fastpeel(r1,r2,fastp_threads,sizeosplit,ishic=inhic or sfastp)
        ##  Write the command to file
        writetofile(fastp_command_file, sbatch(None,fastp_threads,the_cwd,fastp_repo,nice=1,nodelist=nodes,memory=slurm_mem) + fastp_coms, debug)
        ##  Append command to file
        command_files.append((fastp_command_file,sample_name,experi_mode,pipeline_steps[pix],fastp_repo,0,''))
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

        
        ##      BWA MASTER JOB SUBMISSION
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## Set the point of pipeline
        pix = 1
        ## Call the bwa master command
        bwa_master_file = f'{commands_dir}/{pix}.bwa.master.{sample_name}.sh'
        ## Gahter the bwa master command and report
        bwa_master_commands, bwa_master_repo = bwamaster(sample_name,reference_path,bwa_threads,the_cwd,partition,debug,nice,nparallel,forced=force,nodelist=nodes,bwaopts=bwa_opts,memory=slurm_mem)
        ## Write command to file
        writetofile(bwa_master_file, sbatch('bwa.master',1,the_cwd,bwa_master_repo,nice=nice,nodelist=bwanodes) + bwa_master_commands, debug)
        ## Append to command fil
        command_files.append((bwa_master_file,sample_name,experi_mode,pipeline_steps[pix],bwa_master_repo,0,''))
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


        ##      MASTER FILTER BEDPE / Hi-C CONTACTS JOB SUBMISSION 
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------ ##
        ## 2. Set the setp in pipeline and the new bedbe file name
        pix = 2
        ## Call the master filter command
        filter_master_file = f'{commands_dir}/{pix}.filter.bedpe.master.{sample_name}.sh'
        ## Gather the filter master commadn and report
        filter_master_commands, filter_master_repo = filtermaster(sample_name,reference_path,the_cwd,excludes,chrlist,mapq,error_dist,daskthreads,enzymelib,partition,debug,nice,nparallel,forced=force,maxdist=max_dist,chunksize=chunksize,nodelist=nodes,dovetail=dedovetail,removeinter= not inhic,hicexplorer=hicexplorer,memory=slurm_mem)
        ## Write command to file
        writetofile(filter_master_file, sbatch('filter.bedpe.master',1,the_cwd,filter_master_repo,nice=nice,nodelist=filternodes) + filter_master_commands, debug)
        ## Append to command file
        command_files.append((filter_master_file,sample_name,experi_mode,pipeline_steps[pix],filter_master_repo,0,''))
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
        ## Set the new concat file name by chromosome name, in reverse order, this is to free up resrouces as smaller chromosomes finish first
        for chrom in chrlist[::-1]:
            ## Set valid ead
            valid_end = f'_valid_{chrom}.bedpe'
            ## Set start name for wildcard use to bring in inputs  to ftn 
            sample_start = valid_end if postmerging else f'{sample_name}' + valid_end
            hiccat_out = f'{hic_dir}/{sample_name}.valid.{chrom}.bedpe'
            hicdup_out = f'{hic_dir}/{sample_name}.duplicates.{chrom}.bedpe'
            ## Append to hic cat outs to preserve order 
            hiccats_outs.append(hiccat_out)
            hicdups_outs.append(hicdup_out)
            ## SEt the report name
            hiccat_repo = reportname(sample_start,pipeline_steps[pix],i=pix)
            ## Set the command
            hiccat_coms = [f'{executive_dir}/deduphic.py -b {sample_start} -o {hiccat_out} -z {chunksize}' + (' ' if skipduplicates else f' -d {hicdup_out} ') + ('--save-dups' if not keep_dups else '')]

            ## make deduplication file name
            hiccat_file = f'{commands_dir}/{pix}.dedup.{sample_name}.valid.{chrom}.sh'
            ## Write the deduplicat command to file
            writetofile(hiccat_file, sbatch(hiccat_file,daskthreads,the_cwd,hiccat_repo,nice=nice,nodelist=nodes,memory=slurm_mem) + hiccat_coms, debug)
            ## Append the deduplication command
            command_files.append((hiccat_file,sample_name,experi_mode,pipeline_steps[pix],hiccat_repo,0,''))    
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

        
        ## Hi-C CONCATONATION
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## 4. Merging unused contacts, and concatonation of duplicates and valid sets 
        ## Set pipeline step 
        pix = 4
        ## Set start name for wildcard use to bring in inputs to ftn , the output file and the ouput report, set ouptu name  
        unused_start  = [bedpe_dir + ('/*.notused.bedpe' if postmerging else f'/*.{sample_name}.notused.bedpe')]
        
        ## Set file names, concat the inputs, set output names, repor tnames, and files 
        hic_fileends   = ['notused','duplicates','valid']
        concat_inputs  = [unused_start,hicdups_outs,hiccats_outs[::-1]] ## Reverse this to match reverse above 
        concat_outputs = [f'{alignment_dir}/{sample_name}.{f}.bedpe' for f in hic_fileends]
        concat_reports = [reportname(f,pipeline_steps[pix],i=f'{pix}{a}') for a,f in zip(['A','B','C'],concat_outputs)]
        concat_files   = [f'{commands_dir}/{pix}.concat.{sample_name}.{f}.sh' for f in hic_fileends]

        ## Iterate thru the zipped concat files 
        for (con_ins,con_out,con_rep,con_fil) in zip(concat_inputs,concat_outputs,concat_reports,concat_files):
            ## Format the commands 
            concat_coms = pandacat(con_ins,con_out)
            ## Write commands to file 
            writetofile(con_fil, sbatch(con_fil,1,the_cwd,con_rep,nice=nice,nodelist=nodes,memory=slurm_mem) + concat_coms, debug)
            ## append command to command file list
            command_files.append((con_fil,sample_name,experi_mode,pipeline_steps[pix],con_rep,0,'')) 

        ## Patch new cat file to next analysis step, must be the valid contacts/alignments 
        newcatfile = concat_outputs[-1]
        new_catfiles.append(newcatfile)
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


        ## OPTIONAL File conversions and analyses
        ## Set pipeline step 
        pix = 5
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        """
        Optional conversions / analyzes for Hi-C
        5A: BEDPE to .short format
        5B: BEDPE to pairs file format
        5C: Short to Hi-C
        5D: Pairs to mcool
        5E: Inter-chromosomal parsing
        5F: gene x gene counts
        5G: Peak calling with MACS3 (for atac-seq)
        5H: BEDPE to SAM/BAM
        """
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        if inhic:  
            ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
            ## BEDPE TO SHORT or PAIRS
            ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
            ## 5A. Make the bedpe file into a short file for juicer pre command 
            if toshort: 
                ## make a report
                short_repo = reportname(sample_name,'toshort',i=f'{pix}A')
                ## Format the command
                short_commands = [f'{executive_dir}/toshort.py -i {newcatfile}\n'] 
                ## Set command file 
                short_file = f'{commands_dir}/{pix}A.toshort.{sample_name}.sh'
                ## Wriet the short command to file
                writetofile(short_file,sbatch(short_file,daskthreads,the_cwd,short_repo,nice=nice,nodelist=nodes,memory=slurm_mem) + short_commands, debug)
                ## append the short command
                command_files.append((short_file,sample_name,experi_mode,'toshort',short_repo,0,''))
            
            ## 5B. Make the bedpe file into a short file for juicer pre command 
            if makepairs:
                ## make a report
                short_repo = reportname(sample_name,'pairs',i=f'{pix}B')
                ## Format the command
                short_commands = [f'{executive_dir}/toshort.py --pairs -i {newcatfile}\n'] 
                ## Set command file 
                short_file = f'{commands_dir}/{pix}B.pairs.{sample_name}.sh'
                ## Wriet the short command to file
                writetofile(short_file,sbatch(short_file,daskthreads,the_cwd,short_repo,nice=nice,nodelist=nodes,memory=slurm_mem) + short_commands, debug)
                ## append the short command
                command_files.append((short_file,sample_name,experi_mode,'toshort',short_repo,0,''))
                ## Format the pairs file name
                newpairsfile = newcatfile.split('.bedpe')[0] + '.pairs'

            ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
            ##  HIC/COOL FILE CREATION / JUICER PRE 
            ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
            ## 5C. IF we have a valid jar path
            if jarpath:
                ## Set the output hic path
                shortfile, outhicpath = newcatfile.split('.bedpe')[0] + '.short', newcatfile.split('.bedpe')[0] +'.hic'
                ## Sort and saveout the the concatonated hic file 
                jpre_coms, jpre_repo = juicerpre(shortfile,outhicpath,xmemory,jarpath,fastp_threads,binsizes,pathtochrom,f'{pix}C')
                ## make concat file name
                jpre_file = f'{commands_dir}/{pix}C.juicerpre.{sample_name}.sh'
                ## Write the concat command to file
                writetofile(jpre_file, sbatch(jpre_file,fastp_threads,the_cwd,jpre_repo,nice=nice,nodelist=nodes,memory=slurm_mem) + jpre_coms, debug)
                ## Append the concat command
                command_files.append((jpre_file,sample_name,experi_mode,'hic',jpre_repo,0,''))

            ## 5D. Otherwise, make a cooler and mcool file
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
                coolrepo = reportname(sample_name,'mcool',i=f'{pix}D')
                coolfile = f'{commands_dir}/{pix}D.mcool.{sample_name}.sh'
                ## Set cooler commands 
                cooler_coms = [f'cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 {pathtochrom}:{min(binsizes)} {newpairsfile} {outcool}\n',
                            f'cooler zoomify {outcool} -r {",".join(map(str,binsizes))} -o {outmcool}\n'
                            f'rm {outcool}\n',
                            f'{executive_dir}/myecho.py Finished formating mcool file! {coolrepo}\n']
                ## Write the concat command to file
                writetofile(coolfile, sbatch(coolfile,fastp_threads,the_cwd,coolrepo,nice=nice,nodelist=nodes,memory=slurm_mem) + cooler_coms, debug)
                ## Append the concat command
                command_files.append((coolfile,sample_name,experi_mode,'hic',coolrepo,0,''))

            ## 5E. Inter-chromosomsomal analysis 
            if get_inter:
                ## make a report
                short_repo = reportname(sample_name,'inter.short',i=f'{pix}E')
                ## Format the command
                short_commands = [f'{executive_dir}/toshort.py --inter-only -i {newcatfile}\n'] 
                ## Set command file 
                short_file = f'{commands_dir}/{pix}E.inter.short.{sample_name}.sh'
                ## Wriet the short command to file
                writetofile(short_file,sbatch(short_file,daskthreads,the_cwd,short_repo,nice=nice,nodelist=nodes,memory=slurm_mem) + short_commands, debug)
                ## append the short command
                command_files.append((short_file,sample_name,experi_mode,'toshort',short_repo,0,''))


            ## 5F. Feature space analysis, currently disabled
            if feature_space:
                ## Iterate over chromosome list 
                for i,coi in enumerate(chrlist):
                    ## Set the report, commands, and gxg script file name 
                    gxg_repo   = reportname(sample_start,'gxg%s'%i,i=f'{pix}F')
                    gxg_commands = [f'{executive_dir}/gxgcounts.py -i {newcatfile} -c {coi} -f {feature_space} -t {nchrom}' + (' --merge\n' if not i else '\n')]
                    gxg_file     = f'{commands_dir}/{pix}F.gxg.{i}.{sample_name}.sh'
                    ## Write to file for the gxg script, pasing dask thread count, the cwd, commands and debug mode 
                    writetofile(gxg_file,sbatch(gxg_file,daskthreads,the_cwd,gxg_repo,nice=nice,nodelist=nodes,memory=slurm_mem) + gxg_commands, debug)
                    ## Append to command list for calling/submission later 
                    command_files.append((gxg_file,sample_name,experi_mode,'gxg',gxg_repo,0,''))
                
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ##      PEAK CALLING WITH MACS3
        ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## 5G. If we are running analysis on atac-seq experiments and the peak calling is taking place 
        elif peakcalling:
            ## Format the macs3 call report name
            macs3_report, macs3_filename = reportname(sample_name,'macs3',i=f'{pix}G'), f'{commands_dir}/{pix}G.macs3.{sample_name}.sh'
            ## Format the name of the output peaks,  output bed or bedpe file
            peak_path   = f'{peaks_dir}/{sample_name}_peaks.broadPeak' if ifbroad else f'{peaks_dir}/{sample_name}_peaks.narrowPeak'
            outbed_file = f'{peaks_dir}/{sample_name}*.valid.{macs3mode.lower()}'
            ## Format the command to macs3
            macs3_commands = peakattack(newcatfile,sample_name,macs3_report,macs3mode.upper(),gsize=genome_size,incontrols=chip_control,
                                        shiftsize=shift_size,extendsize=extendsize,maxgap=max_gap,minlen=min_len,nolambda=nolambda,
                                        broad=ifbroad,summits=callsummits) +\
                                        [f'{executive_dir}/biotools.py -b {outbed_file} -p {peak_path} -s {diagnostics_dir}/{sample_name}.frip.stats.csv -g {genome_size}\n',
                                         f'{executive_dir}/myecho.py Finished calculating FRiP from macs3 {macs3_report}\n'] + \
                                        ([f'{executive_dir}/annotator.py -i {peak_path} -g {feature_space} --plot\n',
                                         f'{executive_dir}/myecho.py Finished annotating peaks {macs3_report}'] if feature_space else [])
            ## Write the macs3 commands to file
            writetofile(macs3_filename, sbatch(macs3_filename,1,the_cwd,macs3_report,memory=slurm_mem) + macs3_commands, debug)
            ## Append the macs3 command 
            command_files.append((macs3_filename,sample_name,experi_mode,'macs3',macs3_report,0,''))

        ##      SAM or BAM CONVERSION
        ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## 5H. If we are converting bedpe to sam or bam
        if tosam or tobam:
            ## Format the sam file name
            sam_filename = f'{commands_dir}/{pix}H.sam.{sample_name}.sh'
            ## format the commadn to pairs2sam
            sam_commands, sam_report = bedpetosam(newcatfile,pathtochrom,samtools_ncpu,tobam,sample_name)
            ## Wriet the mac
            writetofile(sam_filename, sbatch(sam_filename,samtools_ncpu,the_cwd,sam_report,memory=slurm_mem) + sam_commands, debug)
            ## Append the sam command to command file list
            command_files.append((sam_filename,sample_name,experi_mode,'sam',sam_report,0,''))

    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## <- 
    pix = 6
    ##      COUNTING 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## 6A. If the clean boolean or rerun vars were set 
    if counting:
        ## Format command to remove uneedeed files 
        counting_sh   = f'{commands_dir}/{pix}A.thecount.sh'             ##   Set the bash file name 
        counting_repo = reportname(run_name,'thecount',i=f'{pix}A')   ##   Set the report 
        counting_coms = [f'{executive_dir}/pairs2count.py -i {newf} -c {chunksize} {count_mode}\n' for newf in new_catfiles ]
        ## Format the command to clean up          
        writetofile(counting_sh, sbatch(counting_sh,4,the_cwd,counting_repo,nice=nice,nodelist=nodes) + counting_coms, debug)
        ## Append the clean up command to file
        command_files.append((counting_sh,run_name,experi_mode,'count',counting_repo,0,''))
    else: ## Otherwise do nothing
        pass
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ##
    ##      SUBMITTING TIME-STAMP   
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## 6B. Final timestamp, out clean up
    ## Format count commands 
    ## Set the timesampe assocaited file names 
    timestamp_file   = f'{diagnostics_dir}/{run_name}.timestamp.{stamp}.txt'     ##     Name of the output file 
    timestampsh      = f'{commands_dir}/{pix}B.time.stamp.sh'                    ##     Name of the .sh bash file 
    timestamp_repo = reportname(run_name,f'timestamp.{stamp}',i=f'{pix}B')       ##     Name of the log to report to 
    ## Formath time stamp and echo commands 
    times_commands = [f'{executive_dir}/endstamp.py {timestamp_file} {stamp}\n',
                      f'{executive_dir}/memoryprofile.py\n']
    ## Format the command file name and write to sbatch, we will always ask the timestamp to run even in debug mode 
    writetofile(timestampsh, sbatch(timestampsh,1,the_cwd,timestamp_repo,nice=1,nodelist=nodes) + times_commands, False)
    ## Append the timestamp command to file
    command_files.append((timestampsh,run_name,experi_mode,'timestamp',timestamp_repo,0,''))
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ##
    ##      CLEANING UP FILES   
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## 6C. If the clean boolean or rerun vars were set 
    if ifclean or (rerun == 'clean'): 
        ## Format command to remove uneedeed files 
        remove_sh   = f'{commands_dir}/{pix}C.cleanup.sh'             ##   Set the bash file name 
        remove_repo = reportname(run_name,'clean',i=f'{pix}C')   ##   Set the report 
        remove_comm = [f'{executive_dir}/checkwork.py clean\n', f'{executive_dir}/checkwork.py gzip ./{alignment_dir}/*.bedpe\n', f'{executive_dir}/checkwork.py gzip ./{alignment_dir}/*.short\n', f'{executive_dir}/checkwork.py gzip ./{alignment_dir}/*.valid.pairs\n']
        ## Format the command to clean up          
        writetofile(remove_sh, sbatch(remove_sh,1,the_cwd,remove_repo,nice=nice,nodelist=nodes) + remove_comm, debug)
        ## Append the clean up command to file
        command_files.append((remove_sh,run_name,experi_mode,'clean',remove_repo,0,''))
    else: ## Otherwise do nothing
        pass
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## 
    ##      PIPELINE COMMAND SUBMISSION TO SLURM    
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ##      A) RESTARTING                          
    ## Call the command dataframe, remove previous logs if hard reset was passed 
    command_files, was_hard_reset = commandcontrol(command_files,hardreset,pipeline_steps,rerun)
    #print(command_files)
    ## 
    ##      0) SUBMITTING FASTP SPLITTING          
    ## Call the submission ftn for fastp commands 
    command_files, sub_sbatchs = submitfastp(command_files,sub_sbatchs,partition,stamp,debugmode=debug)
    ## 
    ##      1) SUBMITTING BWA ALIGNMENTS           
    ## Call the submit bwa ftn
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,pipeline_steps[1],pipeline_steps[0],stamp,bwa_partition,debug=debug)
    ##
    ##      2) SUBMITTING BEDPE FILTERING and SPLITTING
    ## Submit the bedpe filtering scirpt 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,pipeline_steps[2],pipeline_steps[1],stamp,filt_partition,debug=debug)
    ##
    ##      3) SUBMITTING DEDUPLICATING AND SORTING 
    ## Submit the dedup and sort script 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,pipeline_steps[3],pipeline_steps[2],stamp,dedup_partition,debug=debug,group='Experiment' if postmerging else 'Sample')
    ##
    ##      4) SUBMITTING CONCATONATION 
    ## Submit the concatonation sciprt 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,pipeline_steps[4],pipeline_steps[3],stamp,partition,debug=debug,group='Experiment' if postmerging else 'Sample')
    ##
    ##      5A) SUBMITTING Gene X Gene interaction script 
    ## Submit the g x g script
    sub_sbatchs = sub_sbatchs + (submitdependency(command_files,'gxg',pipeline_steps[4],stamp,partition,debug=debug,group='Experiment' if postmerging else 'Sample') if feature_space else []) 
    ##
    ##      5B) SUBMITTING CONVERSION FROM BEDPE to JUICER SHORT 
    ## Submit the toshort command 
    sub_sbatchs = sub_sbatchs + (submitdependency(command_files,'toshort',pipeline_steps[4],stamp,partition,debug=debug,group='Experiment' if postmerging else 'Sample') if (toshort or makepairs or get_inter) else []) 
    ##
    ##      5C) Hi-C FILE CREATION 
    ## Call the juicer pre command for hic file creation if jarpath was passed 
    sub_sbatchs = sub_sbatchs + (submitdependency(command_files,'hic','toshort',stamp,partition,debug=debug,group='Experiment' if postmerging else 'Sample') if (jarpath or make_mcool) else [])
    ##
    ##      5D) PEAK CALLING w/ MACS3
    ## Call the peak calling command 
    sub_sbatchs = sub_sbatchs + (submitdependency(command_files,'macs3',pipeline_steps[4],stamp,partition,debug=debug,group='Experiment' if postmerging else 'Sample') if peakcalling else [])
    ##
    ##      5E) SAM or BAM CONVERSION
    ## Call the sam conversion script 
    sub_sbatchs = sub_sbatchs + (submitdependency(command_files,'sam',pipeline_steps[4],stamp,partition,debug=debug,group='Experiment' if postmerging else 'Sample') if (tosam or tobam) else [])
    ##
    ##      6A) COUNTING COMMANDS
    ## Sumbit the count command 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'count',pipeline_steps[4],stamp,clean_partition,group='Experiment',debug=debug)
    ##   
    ##      6B) SUBMITTING TIME STOP COMMANDS 
    ## Submit time stamp 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'timestamp',pipeline_steps[4:-1],stamp,time_partition,group='Experiment',debug=debug) 
    ## 
    ##      6C) CLEAN UP COMMANDS 
    ## Submit the clean up command if the flag was passed 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'clean','count' if counting else 'timestamp',stamp,clean_partition,group='Experiment',debug=debug) 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        

    ##      SAVING OUT COMMAND FILE 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Write the sbatch commands
    writetofile(f'{logs_dir}/sbatch.log.txt',[l+'\n' for l in sub_sbatchs], False)
    ## Calc number of jobs submitted
    njobstosub, njobssubbed = command_files[(command_files.Torun==0)].shape[0], len(sub_sbatchs)
    ## Check our work
    ifprint(f'WARNING: The number of expected jobs to run ({njobstosub}) and number of jobs submitted {njobssubbed} do not match!' , (njobstosub < njobssubbed))
    ## Print the number of commands being submitted
    print(f'INFO: A total of {njobssubbed} jobs were submitted with this run' + (f' and will start after {bwaix_jobid}.' if bwaix_jobid else '.')  )
    ## If zero jobs were submitted
    ifprint('WARNING: Zero jobs were submitted; If this was unexpected try running slurpy again, including --restart flag.',(njobssubbed==0))
    ## Save out the command files, once for patching to other command, and a time stamped one
    command_files.to_csv(f'{logs_dir}/command.file.csv', index=False)
    command_files.to_csv(f'{logs_dir}/command.file.{stamp}.csv', index=False)
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

## If the script is envoked by name 
if __name__ == "__main__":
    main()
## End of file 