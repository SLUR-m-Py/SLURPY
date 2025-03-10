#!/usr/bin/env python
#SBATCH --job-name=inter.count          ## Name of job
#SBATCH --output=%x.%j.out              ## Name stdout
#SBATCH --error=%x.%j.err               ## Name stderr
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=12              ## Number of tasks to be launched
#SBATCH --partition=mpi                 ## Set the partition

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

## Ftn for counting lines
def bedpecount(inpath:str) -> int:
    with open(inpath,'r') as infile:
        for i,l in enumerate(infile):
            pass 
    return i ## NOTE: we are assumign the header is in the line space, so we don't add one

## define ftn for turning off spoines
def spinesoff(inax):
    ## SEt the axis
    plt.sca(inax)
    ## turn off the top and right axis 
    [inax.spines[r].set_visible(False) for r in ['top','right']]
    pass 

## Set the error list
error_list = ['unmapped','dangend','dovetailed','excluded','lowqual','intrafrag','inter']

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

    ## Set paths of duplicates and notused
    not_used_path   = inpath.replace('.valid.','.notused.')
    duplicates_path = inpath.replace('.valid.','.duplicates.') 

    ## Count the number of duplicates
    dup_count = bedpecount(duplicates_path) if exists(duplicates_path) else 0

    ## Load in "not used" list
    if exists(not_used_path):
        ## Read in not used column from df, count by error 
        not_used     = dd.read_csv(not_used_path,sep=hicsep,usecols=['Error'])
        error_counts = [sum(not_used.Error==error) for error in error_list]
        ## Format into dictionary 
        error_dict   = dict(zip(error_list,error_counts))
        ## Appdend duplicate 
        if dup_count:
            error_dict['duplicates'] = dup_count
        #print(error_dict)
    else:
        error_counts = []
    
    ## Set output path
    basename = inpath.split('/')[-1].replace('.bedpe','.inter.intra')
    ## Set out name 
    outname = f'./{diagdir}/{basename}' if exists(f'./{diagdir}') else inpath.replace('.bedpe','.inter.intra')

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

    ## If the dictionary was made 
    if len(error_counts):
        error_dict['Inter'] = ninter
        error_dict['Intra'] = nintra

    ## Calc total and ratio 
    ntotal = ninter + nintra
    rintra = round(nintra/ntotal,4)

    #print(error_dict)

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

    ## Format and rint outpath
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
    distanced.to_csv(outpath.replace('.inter.intra','.distance.decay'),index=False)

    ## Add a pad of one to counts
    distanced['Counts'] = distanced['Counts'] + 1

    ## Convert log to base 10 
    median_dd = np.log10(distanced.Counts.rolling(span).median().values/distanced.Counts.sum())
    log10_dis = np.log10(distanced.Distance.values)

    ## Call figure 
    fig,ax = plt.subplots(1,2,figsize=(11,5))
    fig.set_facecolor('w')

    ## Plot distance decay 
    spinesoff(ax[0])
    plt.plot(log10_dis,median_dd,color='tab:grey')
    plt.xlabel('Genomic Distance (log$_{10}$ bp)',fontsize=myfs)
    plt.ylabel('Contact Frequency (log$_{10}$)',fontsize=myfs)

    ## Add a bar showing percent of contcats
    barax = fig.add_axes([-0.015,0.0975,0.05,0.78])
    spinesoff(barax)

    ## If errors were counted 
    if len(error_counts):
        ## Set keys 
        error_keys  = list(error_dict.keys())
        ## Caluclatte total and preset cumlative percent 
        total_pairs = sum(error_dict.values())
        cumper = 0

        #print(total_pairs)
        for key in error_dict.keys():
            ## GAther value and percent 
            value = error_dict[key]
            vper  = round(value/total_pairs,3)
            ## Plot fake point for label
            m, = plt.plot(-1,0,'s',label=f'{key.capitalize()}: {value} ( {100*vper} %)')
            mcolor = m.get_color()
            ## Plot the percentage 
            plt.vlines(0,cumper,cumper+vper,linewidth=50,color=mcolor)
            ## Add to the percentage 
            cumper = cumper + vper 
            #print(cumper)
        ## Save out the counts
        pd.DataFrame(error_dict.values(),index=error_dict.keys()).T.to_csv(outpath.replace('.inter.intra','.counts'),index=False,header=True)
    
    else:
        ## Add fake color
        plt.plot(0,0,'s',color='tan',label=f'Inter: {ninter}')
        plt.plot(0,0,'s',color='tab:blue',label=f'Intra: {nintra}')

        plt.vlines(0,0,1,color='tan',linewidth=50)
        plt.vlines(0,0,rintra,color='tab:blue',linewidth=50)
    
    ## Modify x-axis 
    plt.xlim(-0.001,0.001)
    plt.xticks([])

    ## Modify y axis 
    plt.ylim(0,1)
    plt.ylabel('Percent of Read Pairs ( %s )'%total_pairs if len(error_counts) else 'Percent of Contacts ( %s )'%ntotal,fontsize=myfs)

    ## Add legend
    plt.legend(bbox_to_anchor=(-1.75,1),frameon=False,fontsize=myfs-2)
    
    ## Set a color bar 
    cbar = fig.add_axes([0.92,0.25,0.01,0.5])

    ## Plot interaction scores 
    spinesoff(ax[1])
    sns.heatmap(np.log2(inter_score_df),cbar_ax=cbar,cmap='cividis')
    ## Modify tick marks 
    plt.xticks(fontsize=10);plt.yticks(fontsize=myfs-2)

    ## Set color bar labels 
    plt.sca(cbar)
    plt.ylabel('log$_2$ (Interaction Score)',fontsize=myfs);

    ## Saveout the png 
    plt.savefig(outname.replace('.inter.intra','contacts')+'.png',dpi=150,bbox_inches='tight')
## EOF      