#!/usr/bin/env python
#SBATCH --job-name=hictojuice           ## Name of job
#SBATCH --output=%x.%j.out              ## Name stdout
#SBATCH --error=%x.%j.err               ## Name stderr
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=12              ## Number of tasks to be launched
#SBATCH --partition=tb                  ## Set the partition

## bring in mods
import dask.dataframe as dd, pandas as pd, numpy as np    
## Load in hic sep 
from parameters import hicsep
## Bring in diag dir 
from directories import diagdir
## Bring in mat plot lib
from matplotlib import pyplot as plt
## Load in mods
import seaborn as sns 
## Bring in exists 
from os.path import exists

## Ftn for calculating inter score 
def intscore(df,ci,cj,N) -> float:
    ## Interactions between chroms
    N_ij = df[(df.Rname1==ci) & (df.Rname2==cj)].Rname1.count() + 1
    ## All interactions envolving chrom 1
    N_i = df[(df.Rname1==ci) | (df.Rname2==cj)].Rname1.count() + 1
    ## All interaction envolving chrom 2
    N_j = df[(df.Rname1==cj) | (df.Rname2==ci)].Rname1.count() + 1
    ## Return the score
    return round(N_ij/(N *(N_i/N) * (N_j/N)),4)

## Set column names 
colnames = ['Rname1','Rname2','Inter','Distance']

## Set the help messages
i_help = "The path to an input bedpe file from SLURM.py."
w_help = "Window size (bp) for calculating distance decay. Default: 10000 bp."
s_help = "Size of span to calculate rolling median."

## SEt font size
myfs = 12

## If the script is envoked 
if __name__ == "__main__":
    ## Bring in argparse
    import argparse
    ## call a parser
    parser = argparse.ArgumentParser(description='Calcualtes inter-chromosomal scores and intra-chromosomal distance decay.')
    ## Add arguments
    parser.add_argument("-i", dest="I", type=str, required=True,  help=i_help, metavar="./path/to/input.bedpe")
    parser.add_argument("-w", dest="W", type=int, required=False, help=w_help, metavar="n", default=10000)
    parser.add_argument("-s", dest="S", type=int, required=False, help=s_help, metavar="n", default=5)
    ## Set inputs
    inputs = parser.parse_args()
    
    ## Set inputs
    inpath = inputs.I
    window = inputs.W
    span   = inputs.S 
    
    ## Set output path
    basename = '%s.inter.intra'%inpath.split('/')[-1].split('.bedpe')[0]
    ## Set out name 
    outname = (f'./{diagdir}/{basename}') if exists(f'./{diagdir}') else (inpath.split('.bedpe')[0] + '.inter.intra')

    ## Load in a dask dataframe 
    df = dd.read_csv(inpath,sep=hicsep,usecols=colnames)

    ## Set the chromlist
    chrlist = list(df.Rname1.unique().compute())

    ## Format inter and intra contacts 
    inter = df[(df.Inter==1)][['Rname1','Rname2']].compute()
    intra = df[(df.Inter==0)].Distance.compute().values

    ## Count the inter and intra contacts
    ninter = inter.Rname1.count()
    nintra = len(intra)

    ## Calc total and ratio 
    ntotal = ninter + nintra
    rintra = round(nintra/ntotal,4)

    ## Initiate inter score df 
    inter_score_df = pd.DataFrame(index=chrlist,columns=chrlist,dtype=float)

    ## Iterate thru the chromosomes 
    for i,c1 in enumerate(chrlist):
        for j,c2 in enumerate(chrlist):
            ## If the chromosomes are not equal and we are going left to right. 
            if (c1 != c2) and (j>i):
                ## Calc the inter -chrom score 
                ics = intscore(inter,c1,c2,ninter)

                ## Place score in df between c1 and c2 and their mirror pos
                inter_score_df.loc[c1,c2] = ics
                inter_score_df.loc[c2,c1] = ics 

    ## Print outpath
    outpath = outname+'.csv'
    print(outpath)
    ## save the inter score path
    inter_score_df.to_csv(outpath,index=True,header=True)

    ## Set window size 
    binsize   = window 
    distanced = []

    ## While intra has counts 
    while len(intra):
        ## Gatather those less than binsize 
        boolix = (intra<binsize)
        ## Append count 
        distanced.append((binsize,sum(boolix)))

        ## Reset intra list and binsize  binsize
        intra = intra[(~boolix)]
        binsize = binsize + window 

    ## format df and save out
    distanced = pd.DataFrame(distanced,columns=['Distance','Counts'])
    distanced.to_csv(outname.split('.inter.intra')[0]+'.distance.decay.csv',index=False)

    ## Add a pad of one to counts
    distanced['Counts'] = distanced['Counts'] + 1

    ## Convert log to base 10 
    median_dd = np.log10(distanced.Counts.rolling(span).median().values/distanced.Counts.sum())
    log10_dis = np.log10(distanced.Distance.values)

    ## Call figure 
    fig,ax = plt.subplots(1,2,figsize=(11,5))
    fig.set_facecolor('w')

    ## Plot distance decay 
    plt.sca(ax[0])
    plt.plot(log10_dis,median_dd,color='tab:grey')
    plt.xlabel('Genomic Distance (log$_{10}$ bp)',fontsize=myfs)
    plt.ylabel('Contact Frequency (log$_{10}$)',fontsize=myfs)

    ## Add a bar showing percent of contcats
    barax = fig.add_axes([-0.05,0.0975,0.05,0.78])
    plt.sca(barax)

    ## Add fake color
    plt.plot(0,0,'s',color='tan',label='Inter')
    plt.plot(0,0,'s',color='tab:blue',label='Intra')

    plt.vlines(0,0,1,color='tan',linewidth=50)
    plt.vlines(0,0,rintra,color='tab:blue',linewidth=50)
    
    ## Modify x-axis 
    plt.xlim(-0.001,0.001)
    plt.xticks([])

    ## Modify y axis 
    plt.ylim(0,1)
    plt.ylabel('Percent of Contacts',fontsize=myfs)

    ## Add legend
    plt.legend(bbox_to_anchor=(2.5,1),frameon=False)
    
    ## Set a color bar 
    cbar = fig.add_axes([0.92,0.25,0.01,0.5])

    ## Plot interaction scores 
    plt.sca(ax[1])
    sns.heatmap(np.log2(inter_score_df),cbar_ax=cbar,cmap='cividis')
    ## Modify tick marks 
    plt.xticks(fontsize=10);plt.yticks(fontsize=myfs-2)

    ## Set color bar labels 
    plt.sca(cbar)
    plt.ylabel('log$_2$ (Interaction Score)',fontsize=myfs);
    
    ## Saveout the png 
    plt.savefig(outname+'.png',dpi=150,bbox_inches='tight')
## EOF      