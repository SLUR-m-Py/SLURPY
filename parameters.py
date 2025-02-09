#!/usr/bin/env python
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
## Set the run local var
runlocal = False 
## Fake job id used for debugging only 
fakejobid = 666666         
## Path to human t2t ref on canopus and my local machine 
if runlocal:
    t2t_refpath = '/Users/croth/Desktop/T2T/GCA_009914755.4_T2T-CHM13v2.0_genomic.named.fasta'
    t2t_gffpath = '/Users/croth/Desktop/T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.named.gff'
else:
    t2t_refpath = '/panfs/biopan04/4DGENOMESEQ/REFERENCES/T2T/GCA_009914755.4_T2T-CHM13v2.0_genomic.named.fasta' 
    t2t_gffpath = '/panfs/biopan04/4DGENOMESEQ/REFERENCES/T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.named.gff'

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      Hi-C and ATAC-seq DEFAULT VARIABLE SETTING  
splitsize    = 10000000      ##     The number of splits made by fastp 
bwathreads   = 4             ##     Number of threads used by calls to bwa 
samthreads   = 4             ##     Number of threads used by calls to samtools 
daskthreads  = 4             ##     Number of threads used by calls to dask df 
parallelbwa  = splitsize     ##     Number of parallel runs of bwa 
fastpthreads = 12            ##     Number of threads in fastp 
part         = 'tb'          ##     Defalut partition 
map_q_thres  = 30            ##     Minimum mapping quality threhosld 
error_dist   = 10000         ##     The idstance to check for erros 
circle_dist  = 10000         ##     The distance to check for self circles 
lib_default  = 'Arima'       ##     Defalut library used to make Hi-C experimetns 
chunksize    = 250000        ##     Chunks size for parsing with pandas
hicsep       = ' '           ##     Text deliminator 
line_count   = 10**7         ##     Number of lines 
fends        = '.fastq.gz'   ##     End of fastq fiels 
mito         = 'chrM'        ##     The name of the mitocondrial contig (in humns)
xmemory      = 49152         ##     Sets the memory used by juicer pre command 
nice         = 10**7         ##     Set the nice parameter 
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
F_help = "The number of splits to make for each pair of input fastq files. Default is: %s. Controls the total number of splits across the run."%splitsize
f_help = "The number of threads used in fastp to split input fastq files. Default is: %s. Note: must be an even multiple of the number of splits."%fastpthreads
b_help = "The number of threads used per bwa alignment on split input fastq files. Default is: %s."%bwathreads
n_help = "Run name used to name output files. Default behavior is to take the common name from the input read pairs."
M_help = "Name of the mitochondrial contig. Default is: %s."%mito
X_help = "List of chromosomes/contigs to exclude from analysis. Default behavior is to process all within the passed .fasta or .fa file."
B_help = "Number of parallel bwa alignments to run. Defaults to %s. Controls the number of bwa jobs submitted at once to slurm."%parallelbwa
P_help = "The type of partition jobs formatted by slurpy run on. Default is %s."%part
Q_help = "Mapping quality threshold to filter alignments. Default is: %s."%map_q_thres
c_help = "Path to control or input bam files used in ChIP-seq experiments."
G_help = "Path to list of chromosomes (by name) to include in final analysis. Default behavior expects a tab seperated tsv or bed, comma seperated csv, or space seperated txt file with no header."
g_help = "Size of the genome being analyzed, used as parameter for macs3. Inputs can be integers in bp or two letter short hand, for e.g. hs for homo sapiens. Default behavior is to calculate this value from the reference file."
C_help = "Linear genomic distance to check outward facing, intra-chromosomal Hi-C contacts for self-circle artifacts. Default is %s. Passing zero (0) will skip this check."%circle_dist 
E_help = "Minimum fragment size of read pairs scanned for an intersecting restriction fragment site (if passed thru library parameter). Default is %s. These pairs are also marked for dangling ends and self-circles."%error_dist
L_help = "The name of the restriction site enzyme (or library prep) used in Hi-C sample creation. Default is %s. Options include Arima, MboI, DpnII, Sau3AI, and HindIII. Note: passing none (i.e. Dovetail) is also allowed, but checks for restriction sites and dangling ends will be skipped."%lib_default
Z_help = "Number of rows loaded into pandas at a time. Default is: %s. WARNING: while increasing could speed up pipeline it could also cause memeory issues."%chunksize
q_help = "The expected file extension of input fastq files. Default is: %s"%fends
t_help = "The number of threads used in calls to functions and calculations with pandas and dask dataframes. Default is: %s."%daskthreads
s_help = "The number of threads used in calls to samtools. Default is: %s"%samthreads
J_help = "Path to juicer jar file for juicer pre command. Required for .hic file creation."
S_help = "Space seperated list of chromosome resolutions (i.e. bin sizes) for .hic files. Default: %s"%' '.join(map(str,binsizes))
x_help = "Amount of Xmx and Xms memory passed to juicer\'s pre command. Default is: %s."%xmemory
A_help = "The path to a gff or bed file with a feature space (i.e. genes) to count gene x gene interactions. Must have columns named Chrom, Left, and Right specifying genomic coordiantes."
a_help = "A SLURM job ID, used as a dependency, specifying all jobs in this run to start after succssful termination."
N_help = "The SLURM nice parameter, an integer value lowering the job priority of submissions. Default is: %s"%nice
#R_help = "Step within the pipeline to re-run from. Default behavior is to recheck logs to determin position of run. Options for Hi-C analysis include: %s. Options for ATAC-seq include: %s"%(h_pipe,a_pipe)
pipe_help = "Step within the pipeline to re-run from. Options include: %s"
R_help   = pipe_help
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
##      BOOLEAN HELP MESSAGES
## Set help messages for bollean vars
restart_help  = "Flag to force the pipeline to reset and run from start."
runlocal_help = "Disables sbatch submission and submits the script via bash to a local os."
debug_help    = "A flag to run in verbose mode, printing sbatch commands. Default behavior is false."
mark_help     = "Pass this flag to skip marking and removing duplicates. Default behavior is false (conduct duplicate marking)."
broad_help    = "Flag to call broad peaks using the --broad-cutoff=0.1 setting in macs3. See macs3 callpeak --help for more details."
clean_help    = "If included will run clean up script at end of run. The default behavior is false, can be run after pipeline."
skipq_help    = "Flag to skip initial quality control and filtering with fastp (i.e. only split reads)."
merge_help    = "Passing this flag will merge across all pairs of fastqs for final output."
peaks_help    = "A boolean flag to skip peak calling via macs3."
short_help    = "A boolean flag to convert output bedpe file to short format for hic creation with juicer pre. Defaults to True if the juicer jarpath (-J) is specificed."
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
slurpy_descr = 'Prepares a run of the SLUR(M)-py pipeline given user inputs (such as path to genome reference and experiment type) and generates an ingredients.tsv file for the slurpy.py function'

## Set file ends used in this script and other filtering stages 
hicfileends_tmp = ['unmapped','oddling','lowqual','distance','dangling','errors','selfcircle','tohic'] 


##      FUNCTIONS
##
## Write ftns for parsing
## Ftn to format inputs and help messages 
def formatdf(inkeys,inputs,helps):
    return pd.DataFrame([(k,inputs[k],helps[k]) for k in inkeys],columns=['Parameter','Input','Help'])

## Set output file name
params_file = 'ingredients.tsv'

## Set definitions
## Ftn for converting list into string
def spacelist(inrow):
    return ''.join([a for a in inrow if a not in ['[','\'',']']]).split(', ') #' '.join(inrow[1:-1].split(','))

## Ftn for detecting if input row is a weird list 
def isstrlist(inrow):
    return (inrow[0] == '[') and (inrow[-1] == ']')

## Bring in pandas 
import pandas as pd 

## Def for resetting inputs 
def resetinputs(indict,expmod):
    idf = pd.read_csv(params_file,sep='\t')
    ## Check the experiment
    #assert indict['e'] == expmod, 'ERROR: The parameters that were passed do not match the input parameters from file'
    ## Iterate thru the input tsv file 
    for i,row in idf.iterrows():
        if i: ## Skip first row 
            newinput = row.Input

            ## If none was passed from the tsv file of parameters, just pass 
            if newinput == 'none':
                continue

            ## Reassign booleans 
            if newinput == 'True':
                newinput = True
            elif newinput == 'False':
                newinput = False 
            elif isstrlist(newinput):
                newinput = spacelist(newinput)

            ## Assign new input from row
            indict[row.Parameter] = newinput
    ## Return modifyied dictionary 
    return indict

## If the script is envoked by name 
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse
    ## Make the parse
    parser = argparse.ArgumentParser(description = slurpy_descr)

    ## Add the required argument
    parser.add_argument("-e","--experiment",      dest="e", type=str,  required=True, help = e_help,  metavar = 'hic, atac, chip, or wgs'                      )
    parser.add_argument("-r", "--reference-path", dest="r", type=str,  required=True, help = r_help,  metavar = refmetavar                                     ) 
    
    ## Add the default arguments
    parser.add_argument("-F", "--fastp-splits",   dest="F", type=int,  required=False, help = F_help, metavar = splitsize,              default = splitsize    )
    parser.add_argument("-B", "--parallel-bwa",   dest="B", type=int,  required=False, help = B_help, metavar = parallelbwa,            default = parallelbwa  )
    parser.add_argument("-P", "--partition",      dest="P", type=str,  required=False, help = P_help, metavar = part,                   default = part         ) 
    parser.add_argument("-M", "--mtDNA",          dest="M", type=str,  required=False, help = M_help, metavar = mito,                   default = mito         )
    parser.add_argument("-X", "--exclude",        dest="X", nargs='+', required=False, help = X_help, metavar = 'chrX, chrY ...',       default = 'none'       )
    parser.add_argument("-Q", "--map-threshold",  dest="Q", type=int,  required=False, help = Q_help, metavar = map_q_thres,            default = map_q_thres  )
    parser.add_argument("-R", "--rerun-from",     dest="R", type=str,  required=False, help = R_help, metavar = 'step',                 default = 'none'       )
    parser.add_argument("-q", "--fastq",          dest="q", type=str,  required=False, help = q_help, metavar = '.fastq.gz',            default = fends        )
    parser.add_argument("-a", "--afterok",        dest="a", type=int,  required=False, help = a_help, metavar = fakejobid,              default = 0            )

    ## Set number of threads 
    parser.add_argument("-f", "--fastp-threads",  dest="f", type=int,  required=False, help = f_help, metavar = fastpthreads,           default = fastpthreads )
    parser.add_argument("-b", "--bwa-threads",    dest="b", type=int,  required=False, help = b_help, metavar = bwathreads,             default = bwathreads   )
    parser.add_argument("-t", "--dask-threads",   dest="t", type=int,  required=False, help = t_help, metavar = daskthreads,            default = daskthreads  )

    ## Set values for Hi-C analysis 
    parser.add_argument("-n", "--run-name",       dest="n", type=str,  required=False, help = n_help, metavar = 'name',                 default = 'none'       )
    parser.add_argument("-E", "--error-distance", dest="E", type=int,  required=False, help = E_help, metavar = 'bp',                   default = error_dist   )
    parser.add_argument("-L", "--library",        dest="L", type=str,  required=False, help = L_help, metavar = 'MboI',                 default = lib_default  )
    #parser.add_argument("-D", "--mindist",        dest="D", type=int,  required=False, help = D_help, metavar = 'n',                    default = same_fragm   )
    parser.add_argument("-Z", "--chunksize",      dest="Z", type=int,  required=False, help = Z_help, metavar = 'n',                    default = chunksize    )
    parser.add_argument("-G", "--genomelist",     dest="G", type=str,  required=False, help = G_help, metavar = './path/to/list.tsv',   default = 'none'       )
    parser.add_argument("-J", "--jar-path",       dest="J", type=str,  required=False, help = J_help, metavar = './path/to/juicer.jar', default = 'none'       )
    parser.add_argument("-x", "--Xmemory",        dest="x", type=int,  required=False, help = x_help, metavar = xmemory,                default = xmemory      )
    parser.add_argument("-S", "--bin-sizes",      dest="S",  nargs='+', required=False, help = S_help, metavar = '25000, 10000, ...',    default = binsizes     )
    parser.add_argument("-gxg","--features",      dest="gxg", type=str,  required=False, help = A_help, metavar= './path/to/my.gff',      default = 'none'       )

    ## Set value for ATAC-seq / peaks analysis
    parser.add_argument("-g", "--genome-size",   dest="g", type=str,   required=False, help = g_help, metavar = g_metavar,              default = None         )
    parser.add_argument("-c", "--controls",      dest="c", nargs='+',  required=False, help = c_help, metavar = c_metavar,              default = None         )

    ## Set boolean flags for Hi-C 
    parser.add_argument("--toshort",              dest="toshort",     help = short_help,    action = 'store_true')

    ## Boolean flags for all analsyis 
    parser.add_argument("--restart",              dest="restart",     help = restart_help,  action = 'store_true')
    parser.add_argument("--debug",                dest="debug",       help = debug_help,    action = 'store_true')
    parser.add_argument("--skipdedup",            dest="skipdedup",   help = mark_help,     action = 'store_true')
    parser.add_argument("--clean",                dest="clean",       help = clean_help,    action = 'store_true')
    parser.add_argument("--merge",                dest="merge",       help = merge_help,    action = 'store_true')

    ## Set additional boolean flags for atac
    parser.add_argument("--skipfastp",            dest="skipfastp",   help = skipq_help,    action = 'store_true')
    parser.add_argument("--broad",                dest="broad",       help = broad_help,    action = 'store_true')
    parser.add_argument("--skipmacs3",            dest="skipmacs3",   help = peaks_help,    action = 'store_true')

    ## Format inputs
    inputs = parser.parse_args()
    ## Call a dctionary 
    vars_dict = vars(inputs)

    from os.path import exists

    ## Set experiment mode
    expmod = inputs.e.lower()

    ## Check the expmod
    assert expmod in explist, "ERROR: The type of experiment ( %s ) was not recognized. Options include %s"%(expmod,', '.join(explist))

    ## Reset reference spaces for Cullen 
    reference_path = inputs.r
    feature_space  = inputs.gxg

    ## reassign expmod
    vars_dict['e'] = './SLURPY/hic.py' if (expmod == 'hic') else './SLURPY/peaks.py'

    ## Format dictionar of inputs
    all_help = [e_help,r_help,
                F_help,B_help,P_help,M_help,X_help,Q_help,R_help,q_help,a_help,
                f_help,b_help,t_help,
                n_help,E_help,L_help,Z_help,G_help,J_help,x_help,S_help,A_help,
                g_help,c_help,
                restart_help,debug_help,mark_help,clean_help,merge_help,short_help,
                skipq_help,broad_help,peaks_help]

    ##      ROTH SETTINGS
    ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------- ##
    ## Reset reference with presets if we are running as Cullen Roth (These can be edited if you are not me :-))
    vars_dict['r']   = t2t_refpath if (reference_path.lower() == 't2t') else reference_path
    vars_dict['gxg'] = t2t_gffpath if (feature_space.lower()  == 't2t')  else feature_space

    ## Check refernece
    assert exists(vars_dict['r']), "ERROR: The given path to the reference file ( %s ) does not exist!"%vars_dict['r']
    if vars_dict['gxg'] != 'none':
        assert exists(vars_dict['gxg']), "ERROR: The given path to the features file ( %s ) does not exist!"%vars_dict['gxg']

    ## Format a help dict
    help_dict = dict(zip(vars_dict.keys(),[a.split('Default')[0] for a in all_help]))

    ## Set parameters we don't need for hi-C analysis or atac-seq by name
    hi_c_novar = ['skipmacs3','broad','skipfastp','c','g']
    peak_novar = ['short','A','S','x','J','Z','D','L','E']

    ## Set hic and atacseq variable list 
    hi_c_vars = [i for i in vars_dict.keys() if i not in hi_c_novar]
    peak_vars = [i for i in vars_dict.keys() if i not in peak_novar]

    ## SEt the parameters df
    idf = formatdf(hi_c_vars,vars_dict,help_dict) if (expmod == 'hic') else formatdf(peak_vars,vars_dict,help_dict)

    ## Save out the df
    idf.to_csv(params_file,sep='\t',index=False)
## End of file