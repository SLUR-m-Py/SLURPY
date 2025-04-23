# SLURPY
[SLUR(M)-py: A SLURM Powered Pythonic Pipeline for Parallel Processing of 3D (Epi)genomic Profiles](https://www.biorxiv.org/content/10.1101/2024.05.18.594827v2)

## Setting up the computing environment
SLUR(M)-py (pronounced slurpy) was developed using anaconda (python v 3.10.13) for a linex OS on a high-performance computing cluster. 
We recommend using conda to manage the python environment needed by SLUR(M)-py. 
Below are commands needed to set up the "bioenv" for running slurpy. 

```
## Make a python environment named bioenv 
conda create -n bioenv 

## Activate the environment
conda activate bioenv 

## After making a new conda environment install needed packages.
## Use conda to install needed mods, and bring in mods from bioconda
conda install numpy pandas matplotlib seaborn dask biopython pysam samtools bwa -c bioconda 
```

For calling peaks in ATAC-seq and narrow peak ChIP-seq data sets we use macs3 (the latest version of macs2).
We recommend installing macs3 via pip (shown below). 

```
pip install macs3
```

A full list of the python libraries and their versions used to develop slurpy are listed within [python.dependencies.txt](https://github.com/SLUR-m-Py/SLURPY/blob/main/python.dependencies.txt).

## Installation
Downloading the repository as a .zip archive is easiest. For developers a simple clone command with git works too:

```
git clone https://github.com/SLUR-m-Py/SLURPY.git
```

Once slurpy is downloaded (and expanded), change the current directory to the local SLURPY directory and modify the python files as executables. 

```
cd ./SLURPY
chmod +x *.py 

```

## Dependencies 

Slurpy utilizes SLURM and was developed under version 21.08.8-2. The suit of tools in samtools is also required with the minimum version of 1.15.1. 


#### fastp 
For splitting intial input read pairs into subsets for parallele processing we use [fastp](https://github.com/OpenGene/fastp).
The fastp executable needs to be within the "SLURPY" directory. After downloading move this executable to the SLURPY directory:

```
## Change directorys
cd ./SLURPY

## download the latest build, you may need to install wget
conda install wget
wget http://opengene.org/fastp/fastp

## Modify the executable with chmod
chmod a+x ./fastp

## Check that fastp works
./fastp -h

```

#### juicer tools (optional)
Here we currently use Juicer tools (specifically the juicer pre command) to process and format penultimate forms of Hi-C data (.bedpe) into a final .hic file. 
The juicer jar files are hosted on the downloads page of the [juicer github](https://github.com/aidenlab/juicer/wiki/Download). 
Downloading a jar file into the SLURPY directory makes calls to format .hic data easier (shown later). 

```
## Move the jar file into the slurpy directory
mv juicer_tools_1.22.01.jar ./SLURPY
```

### Checking the python environment 
Once within the slrupy directory, run the environment checking script, [modcheck.py](https://github.com/SLUR-m-Py/SLURPY/blob/main/modcheck.py).

```
## Activate the computing environment
conda activate bioenv 
python modcheck.py
```

If the environment was created successfully the script will run to completion and print the following:

```
INFO: Modules loaded correctly.
INFO: We are ready to run slurpy!
```

## Running slurpy, examples
### Setting up a project directory
Currently the ATAC- and ChIP-seq protocols are fully functional. To run slurpy, change the current directory to the target project (in the example below, the project is named 2501_001) and soft link to the path of the slurpy executables.

```
cd /path/to/project/directory/2501_001
ln -s /path/to/SLURPY
```

Within the project directory ensure the paired fastqs are listed or linked within a subdirectory 
“./fastqs”. 

```
ls -l ./fastqs/*.fastq.gz
```

The pairs fastq files must be gzipped (.gz) and the file must end with the extension ".fastq.gz". 
The fastq files with in the fastq directory must have the first and second reads in pair marked with "R1" and "R2" surronded by underlines (see below.)
The first in pair must also be first when listed within the directory.

```
ls -l ./fastqs/*.gz

example_sample_1_R1_0001.fastq.gz
example_sample_2_R2_0001.fastq.gz

```

### Envoking slurpy for Hi-C processing
The help menu (-h) of protocols within slurpy lists all the available arguments and default settings. Be sure to activate the conda environment "bioenv". 
```
## Activate our computing environment
conda activate bioenv 

## Call the help menu for hic.py 
$ ./SLURPY/slurm.py -h
usage: slurm.py [-h] -r ./path/to/reference.fasta [-F 10000000 [10000000 ...]] [-T n] [-P tb gpu fast [tb gpu fast ...]] [-M chrM] [-X chrX, chrY ... [chrX, chrY ... ...]] [-Q 30] [-R step] [-a 666666] [-N n] [-j n] [-f 8] [-t 8] [-b 8] [-B ,-5SMP] [-n name] [-E bp] [-L MboI] [-Z n] [-G ./path/to/list.tsv] [-c ./path/to/control.bam [./path/to/control.bam ...]] [-J ./path/to/juicer.jar] [-x 49152] [-S 25000, 10000, ... [25000, 10000, ... ...]] [-gxg ./path/to/my.gff] [--nodelist NODES [NODES ...]] [--toshort] [--pairs] [--restart] [--force] [--debug] [--skipdedup]
                [--clean] [--count] [--nomerge] [--mcool] [--inter-only] [--atac-seq] [--skipfastp] [--broad] [--skipmacs3] [--dedovetail] [--hicexplorer] [--sam] [--bam]

A SLURM Powered, Pythonic Pipeline, Performing Parallel Processing of Piared-end Sequenced Reads Prepaired from 3D Epigenomic Profiles.

options:
  -h, --help            show this help message and exit
  -r ./path/to/reference.fasta, --refix ./path/to/reference.fasta
                        Path to input reference referecne (in .fasta or .fa format) with an assoicated (.fai) bwa index to use for alignment.
  -F 10000000 [10000000 ...], --fastp-splits 10000000 [10000000 ...]
                        The approximate number of reads per split made by fastp on input fastq files. Default is: 10000000.
  -T n, --threads n     The number of threads used across all applications of the run (fastp, bwa, dask.dataframes). Default is: 8.
  -P tb gpu fast [tb gpu fast ...], --partition tb gpu fast [tb gpu fast ...]
                        The type of partition jobs formatted by slurpy run on. Default is tb.
  -M chrM, --mtDNA chrM
                        Name of the mitochondrial contig. Default is: chrM.
  -X chrX, chrY ... [chrX, chrY ... ...], --exclude chrX, chrY ... [chrX, chrY ... ...]
                        List of chromosomes/contigs to exclude from analysis. Default behavior is to process all within the passed .fasta or .fa file.
  -Q 30, --map-threshold 30
                        Mapping quality threshold to filter alignments. Default is: 30.
  -R step, --rerun-from step
                        Step within the pipeline to re-run from. Options include: fastp, bwa, filter, dedup, concat, gxg, toshort, hic, macs3, sam, count, clean
  -a 666666, --afterok 666666
                        A SLURM job ID, used as a dependency, specifying all jobs in this run to start after succssful termination.
  -N n, --nice n        The SLURM nice parameter, an integer value lowering the job priority of submissions. Default is: 100000000
  -j n, --n-parallel n  Number of bwa and filtering jobs allowed to run in parallel. Default: 33
  -f 8, --fastp-threads 8
                        The number of threads used in fastp to split input fastq files. Default is: 8. Note: must be an even multiple of the number of splits.
  -t 8, --dask-threads 8
                        The number of threads used in calls to functions and calculations with pandas and dask dataframes. Default is: 8.
  -b 8, --bwa-threads 8
                        The number of threads used per bwa alignment on split input fastq files. Default is: 8.
  -B ,-5SMP, --bwa-options ,-5SMP
                        A comma starting and seperated list (no spaces) of options (and their values) for the bwa mem algorithm (for example ,-t,2,-k,10,-y,5,-S). See bwa mem for help and a list of options.
  -n name, --run-name name
                        Run name used to name output files. Default behavior is to take the common name from the input read pairs.
  -E bp, --error-distance bp
                        Minimum fragment size of read pairs scanned for an intersecting restriction fragment site (if passed thru library parameter). Default is 25000. These pairs are also marked for dangling ends and self-circles.
  -L MboI, --library MboI
                        The name of the restriction site enzyme (or library prep) used in Hi-C sample creation. Default is Arima. Options include Arima, MboI, DpnII, Sau3AI, and HindIII. Note: passing none (i.e. Dovetail) is also allowed, but checks for restriction sites and dangling ends will be skipped.
  -Z n, --chunksize n   Number of rows loaded into pandas at a time. Default is: 950000. WARNING: while increasing could speed up pipeline it could also cause memeory issues.
  -G ./path/to/list.tsv, --genomelist ./path/to/list.tsv
                        Path to list of chromosomes (by name) to include in final analysis. Default behavior expects a tab seperated tsv or bed, comma seperated csv, or space seperated txt file with no header.
  -c ./path/to/control.bam [./path/to/control.bam ...], --controls ./path/to/control.bam [./path/to/control.bam ...]
                        Path to control or input bam/bedpe files used in ChIP-seq experiments.
  -J ./path/to/juicer.jar, --jar-path ./path/to/juicer.jar
                        Path to a juicer jar file with the juicer pre command. Required for .hic file creation.
  -x 49152, --Xmemory 49152
                        Amount of Xmx and Xms memory passed to juicer's pre command. Default is: 49152.
  -S 25000, 10000, ... [25000, 10000, ... ...], --bin-sizes 25000, 10000, ... [25000, 10000, ... ...]
                        Space seperated list of chromosome resolutions (i.e. bin sizes) for .hic files. Default: 2500000 2000000 1000000 750000 500000 250000 200000 100000 75000 50000 25000 10000 5000
  -gxg ./path/to/my.gff, --features ./path/to/my.gff
                        The path to a gff or bed file with a feature space (i.e. genes) to count gene x gene interactions. Must have columns named Chrom, Left, and Right specifying genomic coordiantes.
  --nodelist NODES [NODES ...]
                        Space seperated list of nodes to run jobs on.
  --toshort             A boolean flag to convert output bedpe file to short format for hic creation with juicer pre. Defaults to True if the juicer jarpath (-J) is specificed.
  --pairs               Convert final output to pairs format defined by the 4DNucleome consortium.
  --restart             Flag to force the pipeline to reset and run from start.
  --force               Flag to force the overwrite of output files generated from bwa.
  --debug               A flag to run in verbose mode, printing sbatch commands. Default behavior is false.
  --skipdedup           Pass this flag to skip marking and removing duplicates. Default behavior is false (conduct duplicate marking).
  --clean               If included will run clean up script at end of run. The default behavior is false, can be run after pipeline.
  --count               Boolean flag to performe diagnostics on Hi-C and ATAC-seq samples.
  --nomerge             Passing this flag will keep replicates/samples seperate acorss (n) pairs of input fastqs, generating (n) outputs rather than one final output.
  --mcool               Flag to make an mcool file with cooler.
  --inter-only          Flag to return only inter-chromosomal contacts
  --atac-seq            Preset mode to run in ATAC-seq analysis mode.
  --skipfastp           Flag to skip initial quality control and filtering with fastp (i.e. only split reads).
  --broad               Flag to call broad peaks using the --broad-cutoff=0.1 setting in macs3. See macs3 callpeak --help for more details.
  --skipmacs3           A boolean flag to skip peak calling via macs3.
  --dedovetail          Boolean flag to remove dovetailed paired-end reads (paired reads with overlapping mapped coordiantes) from analsyis (Default: is not to remove these).
  --hicexplorer         Flag to run stricter intra-fragment filtering.
  --sam                 Flag to convert output .bedpe file from SLUR(M)-py to .sam format via samtools.
  --bam                 Flag to convert output .bedpe file from SLUR(M)-py to .bam format via samtools.

```
### For ATAC-seq experiments
To call the peaks.py script within the slurpy pipeline to analyze an ATAC-seq experiment and save output to .bam format:

```
./SLURPY/slurm.py -r /path/to/reference/file.fasta --atac-seq --bam
```

### For ChIP-seq experiments 
To run slurpy to analyze a ChIP-seq experiment run:

```
./SLURPY/slurm.py -r /path/to/reference/file.fasta -c /path/to/control/or/input.bam --atac-seq
```