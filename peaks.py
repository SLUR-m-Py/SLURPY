#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      atac.py; ATAC- & ChIP-seq Processing & Analysis Pipeline 
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
squeue -u croth | grep 'mpi' | grep 'croth' | awk '{print $1}' | xargs -n 1 scancel



"""
## List slurpy command for VERO analysis
"""
./SLURPY/peaks.py -r ../Cryo_ATAC/GreenMVA_Ref/Chlorocebus_sabeus_mva.fasta -M NC_008066.1 -g 2744115311 -G ../../Chlorocebus_sabeus_mva.genome.sizes.autosome.filtered.bed

"""

##      MODULE LOADING and VARIABLE SETTING
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## Load in vars from prep run
from defaults.parameters import *
## Load in ftns and variables from defaults
from defaults.defaults import *
## Load in ftns from other libraries
from tools.pysamtools import checksam, writetofile, txttobam, outnames, bambyreadname
## Bring in bwa mem ftn for atac-seq
from tools.pybwatools import bwamem_paired
## Load in panda cat
from pipeline.pandacat import pandacat


## Set the ftn descritption
peakatac_descr = "Processing and analysis pipeline for paired-end sequencing data from ATAC-seq (and ChIP-seq) experiments."
## Set the hardset hic mode 
experi_mode = 'peaks'
## Define help messages
R_help = pipe_help%a_pipe
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


##      FUNCTION DEFINING
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## Ftn for formating fastp command to filter and split reads
def fastpeel(r1:str, r2:str, exptype:str, w:int, s:int, pix:int, z=4, options=fastp_opts, toremoveend='.toberemoved.fastq.gz', singleend='.singletons.fastq.gz', failend='.failed.fastq.gz') -> tuple:
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
    split_filtered = fastp_header1 + f' -s {s} -z {z} --thread {w} ' + ' '.join(options) + '\n'

    ## Format command if just splitting is begin performed
    fastp_header, fastp_json, fastp_html = fastcut(r1,r2,split1,split2,report,0)
    ## Call fastp the second time to split reads
    just_split = fastp_header + f' -s {s} -z {z} --thread {w} ' + ' '.join(options) + '\n'

    ## Call the fast dry command and format the remove command 
    dry, throwout = fastdry(r1,r2,report), f'rm {temp1} {temp2}\n'
    ## Format the mv commands for the html and json
    mvh, mvj = f'mv {fastp_json0} {diagdir}/\n', f'mv {fastp_html0} {diagdir}/\n'

    ## Return the comands depenent aupon if we are spliting or filtering, and the reprot 
    return [just_split,dry,mvh,mvj] if inhic(exptype) else [initial_filtering,split_filtered,throwout,dry,mvh,mvj], report

## Define ftn for formaiting fastq files for bwa mem in paried end mode 
def prepbwamem(r1:str, r2:str, refix:str, t:int, pix:int) -> tuple:
    """Formats a bwa mem command given the paired reads, reference index, thread count, and mode."""
    ## Format the bam file name and report name
    bam, report = f'{bamtmpdir}/{getsamplename(r1)}.bam', reportname(getsamplename(r1),'bwa',i=pix) # f'{debugdir}/{getsamplename(r1)}.bwa.log.txt'
    ## Format the bwa mem command, we are in paired end mode 
    makebams = bwamem_paired(r1,r2,refix,bam,report,threads=t)
    ## Return the list of bams to make and the count of those made
    return bam, makebams, report  

## Ftn for formating call to our split command
def splitcommand(bam:str, mito:str, threads:int, pix:int) -> tuple:
    """Formats and returns a command to our protocol for splitting an input bam file."""
    ## Call the report and echo command
    report = reportname(bam,'split',i=pix)
    ## Set the output names fromt he bam file 
    mapped_txt, placed_txt, mito_txt, unmap_txt, bedpe_txt = outnames(bam,mito)
    ## Format the command for splitting the input bam file, index the output bam file, and samtools commands to form bam files from read names
    commands = [f'{scriptsdir}/readspliter.py -b {bam} -M {mito} >> {report}\n'] + [bambyreadname(bam,r,threads)[0] for r in [mapped_txt, placed_txt, mito_txt]] + [splitecho(bam,report,'readspliter.py')]
    ## Return the command
    return commands, report 





## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


##      MAIN SCRIPT & ARGUMENT PARSING 
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## If the script is envoked by name 
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse
    ## Make the parse
    parser = argparse.ArgumentParser(description = peakatac_descr)

    ## Add the required argument
    parser.add_argument("-r", "--refix",         dest="r", type=str,  required=True, help = r_help,  metavar = refmetavar                                 )
    
    ## Add the default arguments
    parser.add_argument("-F", "--fastp-splits",  dest="F", type=int,  required=False, help = F_help, metavar = splitsize,            default = splitsize  )
    parser.add_argument("-B", "--parallel-bwa",  dest="B", type=int,  required=False, help = B_help, metavar = parallelbwa,          default = parallelbwa)
    parser.add_argument("-P", "--partition",     dest="P", type=str,  required=False, help = P_help, metavar = part,                 default = part       ) 
    parser.add_argument("-M", "--mtDNA",         dest="M", type=str,  required=False, help = M_help, metavar = mito,                 default = mito       )
    parser.add_argument("-X", "--exclude",       dest="X", nargs='+', required=False, help = X_help, metavar = 'chrX, chrY ...',     default = []         )
    parser.add_argument("-Q", "--map-threshold", dest="Q", type=int,  required=False, help = Q_help, metavar = map_q_thres,          default = map_q_thres)
    parser.add_argument("-R", "--rerun-from",    dest="R", type=str,  required=False, help = R_help, metavar = 'step',               default = None       )
    parser.add_argument("-q", "--fastq",         dest="q", type=str,  required=False, help = q_help, metavar = '.fastq.gz',          default = fends      )
    parser.add_argument("-a", "--afterok",       dest="a", type=int,  required=False, help = a_help, metavar = fakejobid,            default = 0          )

    ## Set number of threads 
    parser.add_argument("-f", "--fastp-threads", dest="f", type=int,  required=False, help = f_help, metavar = splitsize,            default = splitsize  )
    parser.add_argument("-b", "--bwa-threads",   dest="b", type=int,  required=False, help = b_help, metavar = bwathreads,           default = bwathreads )
    parser.add_argument("-t", "--sam-threads",   dest="t", type=int,  required=False, help = s_help, metavar = samthreads,           default = samthreads )

    ## Set values for macs3
    parser.add_argument("-n", "--run-name",      dest="n", type=str,  required=False, help = n_help, metavar = 'name',               default = None       )
    parser.add_argument("-g", "--genome-size",   dest="g", type=str,  required=False, help = g_help, metavar = g_metavar,            default = None       )
    parser.add_argument("-c", "--controls",      dest="c", nargs='+', required=False, help = c_help, metavar = c_metavar,            default = None       )
    parser.add_argument("-G", "--genomelist",    dest="G", type=str,  required=False, help = G_help, metavar = './path/to/list.tsv', default = None       )

    ## Set boolean flags 
    parser.add_argument("--restart",             dest="start",  help = restart_help,  action = 'store_true')
    parser.add_argument("--debug",               dest="debug",  help = debug_help,    action = 'store_true')
    parser.add_argument("--skipdedup",           dest="mark",   help = mark_help,     action = 'store_true')
    parser.add_argument("--clean",               dest="clean",  help = clean_help,    action = 'store_true')
    parser.add_argument("--merge",               dest="merge",  help = merge_help,    action = 'store_true')

    ## Set ATAC-seq specifict
    parser.add_argument("--skipfastp",           dest="sfast",  help = skipq_help,    action = 'store_true')
    parser.add_argument("--broad",               dest="broad",  help = broad_help,    action = 'store_true')
    parser.add_argument("--skipmacs3",           dest="peaks",  help = peaks_help,    action = 'store_true')

    ## Set the paresed values as inputs
    inputs = parser.parse_args() 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      PASS ARGUMENTS & SET VARIABLES 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Set required variables
    reference_path  = inputs.r       ##     Set path to the reference genome 
                                     ##
    ## Set default vairables         ##
    fastp_splits    = inputs.F       ##     Number of splits in fastp 
    bwa_runs        = inputs.B       ##     Set the number of parallel runs of bwa 
    partition       = inputs.P       ##     Set the partition 
    mito            = inputs.M       ##     Set the mito contig name 
    mapq            = inputs.Q       ##     Set the mapping quality threshold 
    excludes        = inputs.X       ##     List of chromosomes to exclude from analysis 
    rerun           = inputs.R       ##     Setp to rerun pipeline from 
    fend            = inputs.q       ##     End of the input fastq files 
    bwaix_jobid     = inputs.a       ##     The job id to have all submissions wait on   
                                     ##
    ## Set threads                   ##
    fastp_threads   = inputs.f       ##     Number of fastp threads
    bwa_threads     = inputs.b       ##     Number of threads in bwa alignments
    sam_threads     = inputs.t       ##     Number of samtools threads 
                                     ##
    ## Set values for macs3          ## 
    run_name        = inputs.n       ##     The name of the samples 
    g_size          = inputs.g       ##     Size of the genome   
    chip_control    = inputs.c       ##     Set the input control for chip experimetn
    pathtochrom     = inputs.G       ##     Path to list of chromosomes to use 
                                     ##
    ## Set boolean vars              ## 
    hardreset       = inputs.start   ##     Resetart the slurpy run, removing previous
    debug           = inputs.debug   ##     Run in debug mode 
    skipduplicates  = inputs.mark    ##     Boolean to mark duplicates
    sfastp          = inputs.sfast   ##     Flag to skip fastp filtering 
    ifbroad         = inputs.broad   ##     Boolean to activate broader peak calling in macs3 
    ifclean         = inputs.clean   ##     Flag to run clean up script 
    postmerging     = inputs.merge   ##     Forces premerge of outputs 
    skippeaks       = inputs.peaks   ##     Skips peak calling with macs3 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

    
    ##      ROTH SETTINGS
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Reset reference with presets if we are running as Cullen Roth (These can be edited if you are not me :-))
    reference_path = t2t_refpath if (reference_path.lower() == 't2t') else reference_path
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      INITILIZATION 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Check that the fastq and reference paths exists
    assert pathexists(fastqdir), fastqserror
    assert pathexists(reference_path), noref_path%reference_path 
    ## Check the versions of samtools, the user email is an email and the experiment mode is one we know
    assert checksam(), not_sam_err 
    ## If needed reset that the fastp threads and splits such that they are a multiple of each
    fastp_splits,fastp_threads = checkfastp(fastp_splits,fastp_threads)
    ## Reformat clean boolean if clean was passed from restart
    ifclean = patchclean(rerun,ifclean)
    ## Initilizse list for sbatch
    sub_sbatchs = []

    ## Load in macs3 ftns
    from scripts.pymacs3 import peakattack
    ## Set the broad pkeack
    broadpeak = '--broad' if ifbroad else ''

    ## Check our control files if they were passed 
    if chip_control:
        assert (type(chip_control) == list) and (type(chip_control[0]) == str), f'ERROR: Inputs for chip control are not a type we recognize (should be a list of strings)!'
        ## Iterate over the controls inputs 
        for chip_con in chip_control:
            ## Check they exist 
            assert fileexists(chip_con), f'ERROR: Bad path to bam file! Unable to locate input control for chip experiment: {chip_con}'
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      CONFIRM THE HARD RESTART 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## apppend the macs3 dir to the group
    grouped_dirs = [debugdir,aligndir,splitsdir,comsdir,diagdir,bamtmpdir,macs3dir]
    ## If there is a hard reset passed 
    confirmreset(grouped_dirs) if hardreset else None 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      MORE MODULE LOADING 
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
    ## Let the user know we are making directories 
    print(directormaking)
    ## Get the current working dir
    the_cwd = gcwd()
    ## If the run name is none
    run_name = setrunname(run_name,the_cwd) 
    ## Make the directories 
    makedirectories(grouped_dirs)
    ## Gather and remove the old command files, error, and output logs
    [remove(f) for f in sortglob(f'{comsdir}/*.sh') + sortglob(f'{debugdir}/*.err') + sortglob(f'{debugdir}/*.out')]
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      WRITING OUT PARAMS
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    writeparams(f'{experi_mode}.py',run_name,stamp,inputs)

    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      CHROMOSOME GATHERING 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Load in more mods
    ## Gather chromoosoes
    from tools.chrommap import gathering, chromgathering
    ## Inform the user we are gathering chromosomes
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
    index_error = "ERROR: We could not detect an index associated with the input reference path: %s.\nINFO: Index the reference (bwa index) and try again."%reference_path
    ## Assert our truth
    assert isbwaix(reference_path),index_error
    
    ## Format the command file and initilize sample names
    command_files, samplenames, filteredbams = [], [], []

    ## Append to command
    command_files.append(('jobfile','sample',experi_mode,'bwaix','report',0,bwaix_jobid)) if bwaix_jobid else None 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      FASTQ GATHERING
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Inform user we are formating jobs
    print(formatingfastq)
    ## Gather the fastqs 
    in_fastqs = getfastqs(fastqdir+'/*.gz')
    ## Assert we have fastq files
    assert len(in_fastqs), missingfqs
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      SAMPLE NAME HANDELING
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## iterate thru the pairs 
    for fastqix, (r1,r2) in enumerate(in_fastqs):
        ## Gather the sample name, exp mode, and fastp command file 
        sample_name = getsamplename(r1) 
        ## Append the sample names
        samplenames.append(sample_name)
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


        ##      FASTP SPLITTING 
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## Iniate pix at 0, set fastp command interface
        pix = 0
        ## Set the the fastp command file name 
        fastp_command_file =  f'{comsdir}/fastp.{sample_name}.sh'  
        ## Gather the fastp command and report 
        fastp_coms, fastp_report = fastpeel(r1,r2,'hic' if sfastp else experi_mode,fastp_threads,fastp_splits,pix)
        ## Write the command to file
        writetofile(fastp_command_file, sbatch(None,fastp_threads,the_cwd,fastp_report) + fastp_coms, debug)
        ## Append command to file
        command_files.append((fastp_command_file,sample_name,experi_mode,peak_pipeline[pix],fastp_report,0,0))
        ## Format the split fastq files
        fastq_splits = [(f'{splitsdir}/{formatn(i+1)}.{basename(r1)}', f'{splitsdir}/{formatn(i+1)}.{basename(r2)}') for i in range(fastp_splits) ]
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


        ##      PARALLIZATION of ALIGNMENT 
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## Format the bwa commands
        for (s1,s2) in fastq_splits:
            ## Set the step in pipeline 
            pix = 1
            ## Call the bwa meme prep ftn, return the split bam file, the commands, and their report
            sbam, bwa_coms, bwa_report = prepbwamem(s1,s2,reference_path,bwa_threads,pix)
            ## Format the bwa commands
            bwa_command_file =  f'{comsdir}/bwa.{basenobam(sbam)}.sh' 
            ## Write the bwa command to file 
            writetofile(bwa_command_file, sbatch(None,bwa_threads,the_cwd,bwa_report) + bwa_coms, debug)
            ## Append bwa command to list of commands 
            command_files.append((bwa_command_file,sample_name,experi_mode,peak_pipeline[pix],bwa_report,0,0))
            ## --------------------------------------------------------------------------------------------------------------------------------------------------------- ##


            ##      BAM SPLITTING (post bwa) 
            ## --------------------------------------------------------------------------------------------------------------------------------------------------------- ##
            ## 2. Set the setp in pipelien 
            pix = 2
            ## Format report name, commands, and file name 
            split_filename = f'{comsdir}/split.{basenobam(sbam)}.sh'
            ## Call the split ftn
            split_coms, split_report = splitcommand(sbam,mito,sam_threads,pix)
            ## Write the split command
            writetofile(split_filename, sbatch(split_filename,sam_threads,the_cwd,split_report) + split_coms, debug)
            ## Apppend the seperating command
            command_files.append((split_filename,sample_name,experi_mode,peak_pipeline[pix],split_report,0,0))

            ## Gather output bam name and txt names
            name_mapd_txt, name_plac_txt, name_mito_txt, name_unmp_txt, bedpe_txt = outnames(sbam,mito)    
            ## Remap the names to a bam file name 
            mapd_bam, plac_bam, mito_bam = txttobam(name_mapd_txt), txttobam(name_plac_txt), txttobam(name_mito_txt)
            ## --------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      FILE and TXT MERGING 
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Iterate over the sample names
    for si,sname in enumerate(samplenames):
        ## Set the smaple name 
        sample_name = run_name if postmerging else sname 
        ## If the sample is the first and we are merging all smaples 
        if (si > 0) and postmerging:
            break 
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


        ##      BAM MERGING 
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## 3. Set the setp in pipeline and the new bedbe file name
        pix = 3
        ## Format commands to merge on mapped, placed, and mitochondrial 
        for bamend in ['mapped','placed',mito]:
            ## format the merged bam name and wild card
            tmp_inwild_bam = f'*.{bamend}.bam' if postmerging else f'*.{sample_name}.{bamend}.bam'
            ## Format the merge command 
            merge_coms, merge_report = mergebam(f'{sample_name}.{bamend}.bam', tmp_inwild_bam, sam_threads,pix)
            ## Format the name of the merge command 
            merge_file = f'{comsdir}/concat.{sample_name}.{bamend}.sh'
            ## Write the merge command to file 
            writetofile(merge_file, sbatch(merge_file,sam_threads,the_cwd,merge_report) + merge_coms, debug)
            ## Append the merge command
            command_files.append((merge_file,sample_name,experi_mode,peak_pipeline[pix],merge_report,0,0))

        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


        ##      MERGE UNMAPPED TXT 
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## Set the bam end file name
        bamend = 'unmapped'
        ## Format the infilename, we need this last bedpe for the hic
        newcatfile =  f'{aligndir}/{sample_name}.{bamend}.txt'
        ## Format the outfile 
        wildcard_filename = f'{bamtmpdir}/*.{bamend}.txt' if postmerging else f'{bamtmpdir}/*.{sample_name}*.{bamend}.txt'
        ## Set report for panda cat with new coding
        concat_report = reportname(newcatfile,peak_pipeline[pix],i=pix)
        ## Format command to merge the unmapped txt file
        concat_coms = pandacat([wildcard_filename],newcatfile,concat_report)
        ## make concat file name
        concat_file = f'{comsdir}/concat.{sample_name}.{bamend}.sh'
        ## Write the concat command to file
        writetofile(concat_file, sbatch(concat_file,1,the_cwd,concat_report) + concat_coms, debug)
        ## Append the concat command
        command_files.append((concat_file,sample_name,experi_mode,peak_pipeline[pix],concat_report,0,0))
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


        ##      SET MERGED BAM FILE
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## Format the merged bam file
        merged_map = f'{aligndir}/{sample_name}.mapped.bam'
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


        ##      DUPLICATE MARKING
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## Format the merged bam file
        pix = 4
        ## If we are skipping duplicates 
        if skipduplicates: 
            ## Skip duplicate marking 
            ifprint('INFO: Skipping duplicate marking.',skipduplicates)
            ## Patch the merged to the marked bam 
            to_filt_bam = merged_map

        else: ## Otherwise call samblaster, make the file name 
            markdup_filename = f'{comsdir}/mark.{sample_name}.sh' 
            ## Call the mark duplicates ftn, returning the output name of the makred bam, the commands for makring, and the assoicated report 
            to_filt_bam, mark_coms, mark_report = markduplicates(merged_map,sam_threads,pix)
            ## Write to file the mark dups commands
            writetofile(markdup_filename, sbatch(markdup_filename,sam_threads,the_cwd,mark_report) + mark_coms, debug)
            ## Append the marking duplicates command 
            command_files.append((markdup_filename,sample_name,experi_mode,peak_pipeline[pix],mark_report,0,0))
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


        ##      FILTER BAM FILE 
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
        ## Filtering step
        pix = 5
        ## Set the filter command file name 
        filter_filename = f'{comsdir}/filter.{sample_name}.sh'
        ## Call the filter command 
        filt_bam, filt_coms, filt_report = filterbam(to_filt_bam,mapq,sam_threads,chrlist,pix)
        ## Write to file the filt commands and sbatch
        writetofile(filter_filename, sbatch(filter_filename,sam_threads,the_cwd,filt_report) + filt_coms, debug)
        ## Append the filtering command
        command_files.append((filter_filename,sample_name,experi_mode,peak_pipeline[pix],filt_report,0,0))

        ## Append the output filtered bam
        filteredbams.append(filt_bam)
        ## ------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## <- This indent marks all the following commands run post aligment and filtering of individual samples 


    ##      COUNTING BAM and TXT FILES   
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## 6A. Start counting
    pix = 6
    ## Format counting command, the count report, and the command
    counting_filename, count_report = f'{comsdir}/{pix}A.countbams.{run_name}.sh', reportname('fragments','count',i=f'{pix}A')  #reportname(f'{aligndir}/bam.counts.bam','count')
    ## Format the save file name
    save_dist_name = f'-S ./{diagdir}/{run_name}.bam.fragment.dist.png\n'
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      FORMATING FRAGMENT DIST CALC FOR CHIP EXPER
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## If the chip control was given 
    if chip_control:
        ## Format the contorl samples
        joined_controls = ' '.join(chip_control)
        ## Format the fragment commands
        frag_calc_commands = f'{scriptsdir}/fragmentdist.py -b ./{aligndir}/*.primary.*.bam {joined_controls} ' + save_dist_name
    else:
        ## Format the fragment histogram distribution for atac seq samples 
        frag_calc_commands = f'{scriptsdir}/fragmentdist.py -b ./{aligndir}/*.primary.*.bam ' + save_dist_name
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      COUNTING COMMANDS
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## 6A cont. Start analysis 
    ## List the count commands 
    count_commands = [f'{scriptsdir}/countbams.py {run_name}\n', frag_calc_commands, f'{scriptsdir}/myecho.py Finished counting bam files in {aligndir} dir. {count_report}\n']
    ## Wriet the coutn command to file
    writetofile(counting_filename,sbatch(counting_filename,1,the_cwd,count_report) + count_commands, debug)
    ## Append the counting command
    command_files.append((counting_filename,run_name,experi_mode,peak_pipeline[pix],count_report,0,0))
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      PEAK CALLING WITH MAC2
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## 6B. If we are running analysis on atac-seq experiments and the peak calling is taking place 
    if skippeaks: ## Print we are skipping peak calling 
        print("INFO: Skipping peak calling.")
    else: 
        ## Gather the genome size 
        #gsize = genomesize(g_size,reference_path,mito)
        gsize = g_size if g_size else genome_size
        ## Format the macs3 call report name
        macs3_report, macs3_filename = reportname(run_name,'macs3',i=f'{pix}B'), f'{comsdir}/{pix}B.macs3.{run_name}.sh'
        ## Format the command to macs3
        macs3_commands = peakattack(filteredbams,run_name,macs3_report,gsize=gsize,broad=broadpeak,incontrols=chip_control) + [f'{scriptsdir}/pymacs3.py -s {diagdir}/{run_name}.frip.stats.csv\n',f'{scriptsdir}/myecho.py Finished calculating FrIP from macs3 {macs3_report}\n']
        ## Write the macs3 commands to file
        writetofile(macs3_filename, sbatch(macs3_filename,1,the_cwd,macs3_report) + macs3_commands, debug)
        ## Append the macs3 command 
        command_files.append((macs3_filename,run_name,experi_mode,'macs3',macs3_report,0,0))
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      SUBMITTING TIME-STAMP   
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    pix = 7
    ## 7A. Final timestamp
    ## Set the timesampe assocaited file names 
    timestamp_file   = f'{diagdir}/{run_name}.timestamp.{stamp}.txt'             ##     Name of the output file 
    timestampsh      = f'{comsdir}/{pix}A.time.stamp.sh'                                ##     Name of the .sh bash file 
    timestamp_report = reportname(run_name,f'timestamp.{stamp}',i=f'{pix}A')     ##     Name of the log to report to 
    ## Formath time stamp and echo commands 
    times_commands = [f'{scriptsdir}/endstamp.py {timestamp_file} {stamp}\n', f'{scriptsdir}/myecho.py Finished SLURPY run of sample: {run_name}. {timestamp_report}\n']
    ## Format the command file name and write to sbatch, we will always ask the timestamp to run even in debug mode 
    writetofile(timestampsh, sbatch(timestampsh,1,the_cwd,timestamp_report) + times_commands, False)
    ## Append the timestamp command to file
    command_files.append((timestampsh,run_name,experi_mode,'timestamp',timestamp_report,0,0))
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      CLEANING UP FILES   
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## 7B. Cleaning
    ## If the clean boolean or rerun vars were set 
    if ifclean or (rerun == 'clean'): 
        ## Format command to remove uneedeed files 
        remove_sh     = f'{comsdir}/{pix}B.cleanup.sh'                      ##   Set the bash file name 
        remove_report = reportname(run_name,'clean',i = f'{pix}B')   ##   Set the report 
        ## set remove coms list 
        remove_coms = [f'{scriptsdir}/remove.py {bamtmpdir} {splitsdir}\n', f'{scriptsdir}/gzipy.py ./{aligndir}/*.txt\n'] + [f'rm {aligndir}/*.mapped.bam\n' if not skipduplicates else '']
        ## Format the command to clean up          
        writetofile(remove_sh, sbatch(remove_sh,1,the_cwd,remove_report) + remove_coms, debug)
        ## Append the clean up command to file
        command_files.append((remove_sh,run_name,experi_mode,'clean',remove_report,0,0))
    else: ## Otherwise do nothing
        pass
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


    ##      PIPELINE COMMAND SUBMISSION TO SLURM    
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ##      1) RESTARTING                          
    ## Call the command dataframe, remove previous logs if hard reset was passed 
    command_files, was_hard_reset = commandcontrol(command_files,hardreset,peak_pipeline,rerun)
    ## 
    ##      2) SUBMITTING FASTP SPLITTING          
    ## Call the submission ftn for fastp commands 
    command_files, sub_sbatchs = submitfastp(command_files,sub_sbatchs,partition,stamp,debugmode=debug)
    ## 
    ##      3) SUBMITTING BWA ALIGNMENTS           
    ## Call the submit bwa ftn
    command_files, sub_sbatchs = submitbwa(command_files,sub_sbatchs,partition,stamp,bwa_runs,debugmode=debug)
    ## 
    ##      4) SUBMITTING SPLITTING COMMANDS
    ## Submit the bam split commands and append the submitted jobs
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'split','bwa',stamp,partition,bylast=True,debug=debug)
    ## 
    ##      5) SUMITTING THE CONCATONATION COMMANDS   
    ## Submit the concat commands and append the submitted jobs
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'concat','split',stamp,partition,debug=debug,group='Experiment' if postmerging else 'Sample')
    ##
    ##      6) SUBMITTING MARKING COMMANDS      
    ## Submit the duplicate marking commands
    sub_mark = [] if skipduplicates else submitdependency(command_files,'mark','concat',stamp,partition,debug=debug,group='Experiment' if postmerging else 'Sample')
    ## Add the sbatchs 
    sub_sbatchs = sub_sbatchs + sub_mark 
    ##
    ##      7) SUBMITTING BAM FILTERING COMMANDS 
    ## Submit the filtering commands and submit to file
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'filter','concat' if skipduplicates else 'mark',stamp,partition,debug=debug,group='Experiment' if postmerging else 'Sample')
    ##
    ##      8) SUBMITTING COUNTING COMMANDS      
    ## Submit the counting commands and append to the list of submitted jobs 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'count','filter',stamp,partition,group='Experiment',debug=debug)
    ##
    ##      9) SUBMITTING COUNTING COMMANDS      
    ## Submit the call to macs3 if calling peaks 
    sub_macs3 = [] if skippeaks else submitdependency(command_files,'macs3','filter',stamp,partition,group='Experiment',debug=debug)
    ## Append macs3 calls
    sub_sbatchs = sub_sbatchs + sub_macs3
    ## Set next step
    next_step = ['count'] if skippeaks else ['macs3','count']
    ##
    ##      10) SUBMITTING TIME COMMANDS 
    ## Submit time stamp, map the above steps/commands the time stamp needs to wait on 
    sub_sbatchs = sub_sbatchs + submitdependency(command_files,'timestamp',next_step,stamp,partition,group='Experiment',debug=debug)
    ## 
    ##      11) CLEAN UP COMMANDS 
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
    print(f'INFO: A total number of {njobssubbed} jobs were submitted with this run.')
    ## If zero jobs were submitted
    ifprint('WARNING: Zero jobs were submitted; If this was unexpected try running slurpy again, including --restart flag.',(njobssubbed==0))
    ## Save out the command files 
    command_files.to_csv(f'{debugdir}/command.file.csv', index=False)
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
## End of file 