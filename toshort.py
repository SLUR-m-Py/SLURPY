#!/usr/bin/env python
#SBATCH --job-name=hictojuice           ## Name of job
#SBATCH --output=%x.%j.out              ## Name stdout
#SBATCH --error=%x.%j.err               ## Name stderr
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=12              ## Number of tasks to be launched
#SBATCH --partition=tb                  ## Set the partition
## bring in mods
import dask.dataframe as dd 
## Brin in defaults
from gxgcounts import file_end, hicsep
"""
Juicer short format:

        'str1','chr1','pos1','frag1','str2','chr2','pos2','frag2'

str:  the strand of the contact, 0 for positive orientation anything else negtive
chr:  the chormosome of the contact
pos:  the position of the contact
frag: the fragment, if using dummy var must be different for the pair 
"""
## Set use columns
use_cols = ['Seqrev1','Rname1','Pos1','Mapq1','Seqrev2','Rname2','Pos2','Mapq2']
## Generate a pre cursor to a .hic flie
desc = "Converts an input bedpe file (representing Hi-C contacts from SLURPY) to a short formated text file for juicer pre command."
## Set help message
I_help = "Input path to a bedpe file from SLURPY Hi-C pipeline."

## If the script is envoked 
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse
    ## Make the parse
    parser = argparse.ArgumentParser(description=desc)
    ## Add the required arguments
    parser.add_argument("-i", dest="I", type=str, required=True, help=I_help, metavar='./path/to/input.bedpe') 
    ## Set the paresed values as inputs
    inputs = parser.parse_args()

    ## Set input
    input_path  = inputs.I 
    ## Check path
    assert file_end in input_path, "ERROR: We expected an input .bedpe file and didn't find that extension in: %s"%input_path

    ## SEt output path
    output_path = input_path.split(file_end)[0] + '.short' 

    ## Load in data
    df = dd.read_csv(input_path,sep=hicsep,usecols=use_cols)
    ## Remap mapq coulmns as our fragment dummy vars
    df['Mapq1'] = 0
    df['Mapq2'] = 1

    ## Save out data
    df[use_cols].to_csv(output_path,single_file=True,header=False,index=False,sep=hicsep)

    ## Print to log
    print("Finished converting input bedpe file (%s) to short format."%input_path)
## End of file 