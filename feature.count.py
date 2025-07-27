#!/usr/bin/env python
#SBATCH --job-name=bed.count   ## Name of job
#SBATCH --output=%x.%j.out     ## Name stdout
#SBATCH --error=%x.%j.err      ## Name stderr
#SBATCH --nodes=1              ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1    ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=12     ## Number of tasks to be launched
#SBATCH --partition=mpi        ## Set the partition

## Set help 
I_help = "Input path to a bedpe file from SLURPY HI-C pipeline."
B_help = "Path to a bed file with features for cross reference and to count."
O_help = "The path to save out the feature count in bed format (optional)."
## Set description
description = 'Feature counts: counts the number of fragments (from a bedpe file) for each feature within a bed (tsv) file.'

## If the script is envoked 
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse
    ## Make the parse
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument("-i", dest="I", type=str, required=True,  help=I_help, metavar='./path/to/input.bedpe')
    parser.add_argument("-b", dest="B", type=str, required=True,  help=B_help, metavar='./path/to/features.bed')
    parser.add_argument("-o", dest="O", type=str, required=False, help=O_help, metavar='./path/to/save.out.bed', default=None)
    ## Set the paresed values as inputs
    inputs = parser.parse_args()

    ## Set input paths 
    bedpe_path = inputs.I
    inbed_path = inputs.B
    outcu_path = inputs.O 

    ## Reset out path
    outcu_path = outcu_path if outcu_path else bedpe_path.split('.bedpe')[0] + '.counts.bed'

    ## Load in dask
    import dask.dataframe as dd

    ## Bring in the fragments and features 
    fragments = dd.read_csv(bedpe_path,sep=' ',usecols=['Rname1','Pos1','End2'])
    features  = dd.read_csv(inbed_path,sep='\t')

    ## Gather the list of chromosomes
    chrlist = features.Chrom.unique().compute()

    ## Initate total frag counts
    total_fragments = 0
    total_feats     = features.Chrom.count().compute()

    ## Iterate thru the list of chromosomes
    for i,c in enumerate(chrlist):
        ## Gather the features
        chrom_feats = features[(features.Chrom==c)].compute()
        chrom_frags = fragments[(fragments.Rname1==c)].compute()

        ## Iniate feat counts 
        feat_counts = []
        ## Iterate thru the rows
        for k,row in chrom_feats.iterrows():
            rstart,rend = row.Start,row.End
            thecount = chrom_frags[(chrom_frags.Pos1<=rend) & (chrom_frags.End2>=rstart)].Rname1.count()
            ## append the count
            feat_counts.append(thecount)
        
        ## add a column of feat counts
        chrom_feats['Fragments'] = feat_counts
        ## Add to the total counts
        total_fragments += sum(feat_counts)

        ## Save out the counts
        chrom_feats.to_csv(outcu_path,sep='\t',index=False,header=False if i else True,mode='a' if i else 'w')
## Print out a save
print("Finished counting %s fragments across %s features."%(total_fragments,total_feats))
## EOF 