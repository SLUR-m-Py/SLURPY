#!/usr/bin/env python
#SBATCH --job-name=annotator            ## Name of job
#SBATCH --output=%x.%j.out              ## Name stdout
#SBATCH --error=%x.%j.err               ## Name stderr
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=12              ## Number of tasks to be launched
#SBATCH --partition=mpi                 ## Set the partition
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## Set verbosity
verbose = False 
## Load in mods 
import sys, pandas as pd, numpy as np, seaborn as sns 
## Load in current wd 
from os import getcwd
## append path path
sys.path.append(getcwd()) 

## Bring in matplotlib
from matplotlib import pyplot as plt 

## define ftn for turning off spoines
def spinesoff(inax):
    ## SEt the axis
    plt.sca(inax)
    ## turn off the top and right axis 
    [inax.spines[r].set_visible(False) for r in ['top','right']]
    pass 

## Load in gene id 
def geneid(A) -> str:
    return A.split('gene_id "')[-1].split('";')[0]

## Set gff column names
gtfnames = ['Chrom','Source','Feature','Left','Right','Score','Strand','Frame','Attribute']

## Write ftn for loading in a gff file
def loadgtf(inpath:str,slimdown=True,id=True) -> pd.DataFrame:
    """Loads in an input gtf (or gff) file with column names using pandas read csv function."""
    ## Return the gff as a dataframe obj 
    gtf = pd.read_csv(inpath,sep='\t',comment='#',names=gtfnames,header=None)
    ## Drop source
    if slimdown:
        gtf.drop(['Source','Score'],axis=1,inplace=True)
    ## Add gene id
    if id:
        gtf['GeneID'] = gtf.Attribute.apply(geneid) 
    ## REturn gtf obj
    return gtf 

## Load in narrow peak 
from pymacs3 import loadnarrowpeak

## Set the description and help messages 
desc = 'Annotates peaks from input narrowPeak file with gene and genomic regions (taken from a gtf file).'
I_help = 'The path to an input .narrowPeak file (usually from MACS2, MACS3, or ChIPr)'
G_help = 'The path to a .gtf file (representing genomic annotations).'
O_help = 'A given file name to save the output generated here (optional).'
W_help = 'Window size (bp) formed around TSS used to define promoter regions (default: 500)'
P_help = 'Boolean flag to plot results and save a .png file.'

## Set needed vars
not_needed = ['CDS','transcript']
tcols      = ['Chrom','GeneID','Strand','TSS']       
annot_cols = ['Chrom','Feature','Left','Right','Strand','GeneID']

## Set var for plotting 
colors      = ['tab:purple','tab:blue','tab:orange','tab:green','tan','red'][::-1]

## Calculate the distribution of distance from TSS use the folowing distance groups 
"""
0 - 1 kb 
1 - 3
3 - 5
5 - 10
10 - 100
> 100
"""
## Set distribution of bins and labels  
dist_labels = ['0 - 1','1 - 3','3 - 5','5 - 10','10 - 100','> 100'][::-1]
dist_bins = [[0,1000],[1000,3000],[3000,5000],[5000,10000],[10000,100000],[100000,np.inf]]

"""
Visulization and methods sampled from:
https://hbctraining.github.io/Intro-to-ChIPseq/lessons/12_functional_analysis.html
"""

## If the script is envoked 
if __name__ == "__main__":
    ## Bring in argparse and set parser
    import argparse
    ## Make the parse
    parser = argparse.ArgumentParser(description=desc)
    ## Add arguments 
    parser.add_argument("-i", dest="I", type=str,  required=True,  help=I_help, metavar='./path/to/peaks.narrowPeak') 
    parser.add_argument("-g", dest="G", type=str,  required=True,  help=G_help, metavar='./path/to/genomic.gtf') 
    parser.add_argument("-o", dest="O", type=str,  required=False, help=O_help, metavar='./path/to/saveout.csv', default=None)
    parser.add_argument("-w", dest="W", type=int,  required=False, help=W_help, metavar='bp', default=500)
    ## Add boolean vars
    parser.add_argument("--plot", dest="P",  help = P_help, action = 'store_true')

    ## Set the paresed values as inputs
    inputs = parser.parse_args()

    ## Set input vars
    peak_path = inputs.I
    gtfs_path = inputs.G
    outs_path = inputs.O
    prom_oter = inputs.W
    make_plot = inputs.P

    ## Load in gtf
    gtf = loadgtf(gtfs_path)
    ## Drop unneded features
    gtf = gtf[(~gtf.Feature.isin(not_needed))].reset_index(drop=True)
    ## Gather genes only
    genes = gtf[(gtf.Feature=='gene')]

    ## Format gene df into TSS df, gather pos and negative genes 
    pgenes,ngenes = genes[(genes.Strand=='+')].copy(), genes[(genes.Strand=='-')].copy()
    ## Set pos and neg col
    pcols,ncols = tcols[:-1] + ['Left'], tcols[:-1] + ['Right']
    ## Gather positive and negative tss 
    ptss, ntss = pgenes[pcols],ngenes[ncols]
    ## Set new column names 
    ptss.columns = tcols
    ntss.columns = tcols
    ## Concat show head 
    tss = pd.concat([ptss,ntss]).reset_index(drop=True)
    ## Set index 
    tss.index = tss.GeneID

    ## Make a promotor df, we need this seperatly to match the gtf format, set feature to promoter
    pgenes['Feature'] = 'promoter'
    ngenes['Feature'] = 'promoter'
    ## Gather appropriate start of pos and neg gens 
    pstart, nstart = pgenes.Left, ngenes.Right
    ## Add left and right of positive and negative 
    pgenes['Left'],pgenes['Right']  = pstart - prom_oter, pstart + prom_oter
    ## Calc left, start bound 
    ngenes['Left'],ngenes['Right']  = nstart - prom_oter, nstart + prom_oter
    ## Concat to promoters
    promoters = pd.concat([pgenes,ngenes],axis=0)

    ## Concat features from gtf into a single df
    features = pd.concat([gtf,promoters],axis=0).reset_index(drop=True)
    ## Gather unique feitures show score 
    uni_feats = sorted(features.Feature.unique())
    uni_score = [4,5,2,1,3]
    ## SEt feature dict 
    uni_dict = dict(zip(uni_feats,uni_score))
    print(uni_dict) if verbose else None 
    ## Add score column 
    features['Score'] = features.Feature.map(uni_dict)

    ## Load in narrow peak, drop na columns
    tmp = loadnarrowpeak(peak_path)
    tmp.dropna(axis=1,how='all',inplace=True)
    ## Add a position col
    tmp['Pos'] = tmp[['Start','End']].mean(axis=1)
    ## If strand value is not unique, drop it 
    if tmp.Strand.unique()[0] == '.':
        tmp.drop('Strand',axis=1,inplace=True)

    ## Annotate peaks, initilize new peak list 
    new_peaks = []
    ## Group by chrom 
    for c,cdf in tmp.groupby('Chrom'):
        ## gather the features mapping to this chromosome 
        cfeats = features[(features.Chrom==c)]

        ## Iterate thru the peaks for this chormosome 
        for i,k in cdf.iterrows():
            ## Find overlap 
            overlap = cfeats[(cfeats.Left<=k.End) & (cfeats.Right>=k.Start)]
            ## If none was found, continue 
            if not overlap.Chrom.count():
                continue
            
            ## Gather the feature we overlap with by max scores
            oix = overlap.groupby('GeneID').Score.idxmin().tolist()
            annotes = overlap.loc[oix]
            ## Add annotation to list as a df 
            new_peaks.append(pd.DataFrame(k).T.merge(annotes[annot_cols],on='Chrom'))

    ## Concat annotated peaks and replace gene with intron within feature column, calc dist to nearest TSS 
    new_peaks = pd.concat(new_peaks)
    new_peaks['Feature'] = new_peaks.Feature.replace(dict(zip(['gene'],['intron'])))
    new_peaks['Distance'] = new_peaks.Pos - tss.loc[new_peaks.GeneID].TSS.values
    
    ## Gather intergenic peaks from the ones we did not annotate above
    inter_geneic = tmp[(~tmp.Name.isin(new_peaks.Name.tolist()))].copy()
    inter_geneic['Feature'] = 'Intergenic'

    ## Initiate nearest gene search, listing, and index 
    neargene = []
    neardist = []
    neargixs = []

    ## Group by chorm in nearrow peak 
    for c,cdf in inter_geneic.groupby('Chrom'):
        ## Gather the tss by chrom 
        ctss = tss[(tss.Chrom==c)]
        ## Iterate thru the rows of the narrow peak/ by chrom find min
        for i, row in cdf.iterrows():

            ## Calc the raw distance from the tss, then take the absolute L1 norm of the distance 
            raw_dist = ctss.TSS - row.Pos
            abs_dist = abs(raw_dist)
            
            ## Find the min abs distance index, the geneID and the raw distance
            min_abs_dix = abs_dist.idxmin()
            min_gene_id = ctss.loc[min_abs_dix,'GeneID']
            min_raw_dis = raw_dist.loc[min_abs_dix]
            
            ## Append values 
            neardist.append(min_raw_dis)
            neargene.append(min_gene_id)
            neargixs.append(i)

    ## Add to the narrow peak 
    inter_geneic['GeneID']   = pd.Series(neargene,index=neargixs)
    inter_geneic['Distance'] = pd.Series(neardist,index=neargixs)
    
    ## Concatonate results
    annoated_peaks = pd.concat([new_peaks,inter_geneic],axis=0)
    ## Drop the position column
    annoated_peaks.drop('Pos',axis=1,inplace=True)
    ## Save out the peaks, set path
    save_path = outs_path if outs_path else '.'.join(peak_path.split('.')[:-1]) + '.annotated.peaks.csv'
    annoated_peaks.to_csv(save_path,header=True,index=False)

    ## Gather unique peak counts
    l,c = np.unique(annoated_peaks.Name.tolist(),return_counts=True)
    ambiquous_total   = np.sum(c>1)
    peak_total        = len(c)
    ambiquous_percent = np.round(ambiquous_total/peak_total,4)
    ## Print results to screen 
    print('INFO: Annotated %s peaks'%peak_total)
    print('INFO: %s %s of peaks (%s) have more than one annotation.'%(ambiquous_percent,'%',ambiquous_total))

    ## IF we are plotting
    if make_plot:
        ## Calcualte precentages
        percentages = annoated_peaks.groupby('Feature').count()[['Chrom']]/annoated_peaks.shape[0]
        percentages.columns = ['Percentage']
        percentages = percentages.reset_index()
        print(percentages) if verbose else None 

        ## Count the dists that are positive and away from TSS 
        pos_counts = [annoated_peaks[(annoated_peaks.Distance>=a) & (annoated_peaks.Distance<b)].Chrom.count() for (a,b) in dist_bins]
        ## Count the dists that are positive and away from TSS 
        neg_counts = [annoated_peaks[(annoated_peaks.Distance<=-a) & (annoated_peaks.Distance>-b)].Chrom.count() for (a,b) in dist_bins]
        ## Format dist counts
        dis_counts = neg_counts+pos_counts

        ## Calc cum sum 
        pos_sum = np.cumsum(pos_counts)/np.sum(dis_counts)
        neg_sum = np.cumsum(neg_counts)/np.sum(dis_counts)

        ## Call figure, set face color 
        fig,ax = plt.subplots(1,2,figsize=(10,3))
        fig.set_facecolor('w')

        ## Plot a bar plot with percentage of features 
        spinesoff(ax[0])
        sns.barplot(x='Percentage',data=percentages,hue='Feature',alpha=0.9)
        ## SEt xlims
        plt.xlim(0,0.5)

        ## Plot the distirbuiton of distances 
        spinesoff(ax[1])

        for i,c in enumerate(pos_sum[::-1]):
            plt.plot(0,0,'s',color=colors[i],label=dist_labels[i])
            plt.hlines(0,0,c,linewidth=50,color=colors[i])

        for i,c in enumerate(neg_sum[::-1]):
            plt.hlines(0,0,-c,linewidth=50,color=colors[i])

        ## ADd leged, vertical line, and remove yticks 
        plt.legend(title='kb to TSS',frameon=False,bbox_to_anchor=(1.3,0.9))
        plt.vlines(0,-0.01,0.01,color='k')
        plt.yticks([])

        ## REset xticsk 
        new_xticks = [-0.65,-0.45,-0.3,-0.15,0,0.15,0.3,0.45,0.65]
        plt.xticks(new_xticks,np.abs(new_xticks))
        plt.xlabel('Portion of Peaks (%) 5` to 3`')

        ## Save figure
        plt.savefig('.'.join(save_path.split('.')[:-1])+'.png',dpi=300,bbox_inches='tight')
    ## EOF 