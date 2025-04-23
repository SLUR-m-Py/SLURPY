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
else:
    t2t_refpath = '/panfs/biopan04/4DGENOMESEQ/REFERENCES/T2T/GCA_009914755.4_T2T-CHM13v2.0_genomic.named.fasta' 
    t2t_gffpath = '/panfs/biopan04/4DGENOMESEQ/REFERENCES/T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.named.gff'

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
nparallel    = 33            ##     Number of jobs to run in parallele 
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

## Define options for fastpeel ftn
fastp_opts = ['--dont_eval_duplication','--disable_length_filtering','--disable_adapter_trimming','--disable_quality_filtering','--disable_trim_poly_g']

## Set list of experiments 
explist = ['wgs','atac','chip','hic']
ST = 'store_true'
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      METAVARS
## Set metavars
c_metavar  = './path/to/control.bam'
g_metavar  = 'bp'
refmetavar = './path/to/reference.fasta'  
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      HELP MESSAGES 
## Set help messages for input and default variables 
e_help = "The type of epigenomic expeirment (for exmale ATAC or Hi-C) used to generate sequenced pair-end reads (options include wgs, atac, chip, or hic)."
r_help = "Path to input reference referecne (in .fasta or .fa format) with an assoicated (.fai) bwa index to use for alignment." 
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
c_help = "Path to control or input bam/bedpe files used in ChIP-seq experiments."
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
R_help = "Step within the pipeline to re-run from. Options include: %s"
m_help = "Maximum allowed distance between intra-chromosomal pairs. Default is zero, setting will activate filter."
j_help = "Number of bwa and filtering jobs allowed to run in parallel. Default: %s"%nparallel

## Specific help messages
node_help = "Space seperated list of nodes to run jobs on."
dove_help  = "Boolean flag to remove dovetailed paired-end reads (paired reads with overlapping mapped coordiantes) from analsyis (Default: is not to remove these)."
atac_help  = "Preset mode to run in ATAC-seq analysis mode."
intra_help = "Boolean flag to remove read pairs spanning multiple chromosomes."
mcool_help = "Flag to make an mcool file with cooler."
pairs_help = "Convert final output to pairs format defined by the 4DNucleome consortium."
inter_help = "Flag to return only inter-chromosomal contacts"
hicex_help = "Flag to run stricter intra-fragment filtering."

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
count_help    = "Boolean flag to performe diagnostics on Hi-C and ATAC-seq samples."
skipq_help    = "Flag to skip initial quality control and filtering with fastp (i.e. only split reads)."
merge_help    = "Passing this flag will keep replicates/samples seperate acorss (n) pairs of input fastqs, generating (n) outputs rather than one final output."
peaks_help    = "A boolean flag to skip peak calling via macs3."
short_help    = "A boolean flag to convert output bedpe file to short format for hic creation with juicer pre. Defaults to True if the juicer jarpath (-J) is specificed."
tosam_help    = "Flag to convert output .bedpe file from SLUR(M)-py to .sam format via samtools."
tobam_help    = "Flag to convert output .bedpe file from SLUR(M)-py to .bam format via samtools."

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      INFO MESSAGES
## Set messages printed to user 
directormaking  = 'INFO: Making local directories.'
formatingfastq  = 'INFO: Formatting jobs.'

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      ERROR MESSAGES
## Write error messages
fastqserror = "ERROR: Unable to detect a fastqs directory holding fastq files!"
missingfqs  = "ERROR: No fastq.gz files were detected!"
not_sam_err = "ERROR: The detected version of samtools is not greater than or equal to v 1.15.1!\nPlease update samtools and try again."
noref_path  = "ERROR: We could not find the provided input path -- %s -- of the reference file!"
index_error = "ERROR: We could not detect an index associated with the input reference path: %s.\nINFO: Index the reference (bwa index) and try again."

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
slurpy_descr = 'A SLURM Powered, Pythonic Pipeline, Performing Parallel Processing of Piared-end Sequenced Reads Prepaired from 3D Epigenomic Profiles.'
## EOF 