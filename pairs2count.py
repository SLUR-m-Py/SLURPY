#!/usr/bin/env python
#SBATCH --job-name=pairs2.count         ## Name of job
#SBATCH --output=%x.%j.out              ## Name stdout
#SBATCH --error=%x.%j.err               ## Name stderr
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=12              ## Number of tasks to be launched
#SBATCH --partition=mpi                 ## Set the partition

## Bring in sys, and other mods  and getcwd
import sys, pandas as pd, numpy as np    
## Bring in exists 
from os.path import exists
## Load in current wd 
from os import getcwd
## append path path
sys.path.append(getcwd()) 

## Bring in mat plot lib and other plotting funcions 
from matplotlib import pyplot as plt
## Load in mods
import seaborn as sns 

## Ftn for calculating inter score 
def intscore(df,ci,cj,N,c1='Rname1',c2='Rname2') -> float:
    ## Interactions between chroms
    N_ij = df[(df[c1]==ci) & (df[c2]==cj)].Inter.sum() + 1
    ## All interactions envolving chrom 1
    N_i = df[(df[c1]==ci) | (df[c2]==cj)].Inter.sum()  + 1
    ## All interaction envolving chrom 2
    N_j = df[(df[c1]==cj) | (df[c2]==ci)].Inter.sum()  + 1
    ## Return the score
    return round(N_ij/(N *(N_i/N) * (N_j/N)),4)

## Ftn for counting interchormosomeal scores between chromosomes
def interscoring(df,chrlist) -> pd.DataFrame:
    ## Initiate inter score df 
    inter_score_df = pd.DataFrame(index=chrlist,columns=chrlist,dtype=float)
    ## Iterate thru the chromosomes 
    for i,c1 in enumerate(chrlist):
        for j,c2 in enumerate(chrlist):
            ## If the chromosomes are not equal and we are going left to right. 
            if (c1 != c2) and (j>i):
                ## Calc the inter -chrom score 
                ics = intscore(df,c1,c2,ninter)

                ## Place score in df between c1 and c2 and their mirror pos
                inter_score_df.loc[c1,c2] = ics
                inter_score_df.loc[c2,c1] = ics 
    ## REturn the inter scores across chromosomes 
    return inter_score_df

## Ftn for counting lines
def bedpecount(inpath:str) -> int:
    with open(inpath,'r') as infile:
        for i,l in enumerate(infile):
            pass 
    return i ## NOTE: we are assumign the header is in the line space, so we don't add one

## define ftn for turning off spoines
def spinesoff(inax) -> None:
    ## SEt the axis
    plt.sca(inax)
    ## turn off the top and right axis 
    [inax.spines[r].set_visible(False) for r in ['top','right']]
    pass 

## Set column names 
mycols = ['Rname1','Rname2','Inter','Distance']

## Def fnt for distance count
def countdistance(df:pd.DataFrame,binsize:int) -> pd.DataFrame:
    ## Format into a df, devide by nearest kb 
    k = pd.DataFrame(df[(df.Inter==0)].Distance//binsize)
    ## Modify by binsize, set dummy 
    k['Distance'] = binsize*k 
    k['Contacts']   = 1
    ## Return the count 
    return k.groupby('Distance').count()

## Define ftn for counting inter chromosome
def countinter(df,c1='Rname1',c2='Rname2') -> pd.DataFrame:
    ## Take contacts between chrom, group by chrom's and count, drop the distance column 
    return df[(df.Inter==1)].groupby([c1,c2]).count()

## Ftn for counting errors 
def counterror(df) -> pd.DataFrame:
    ## Return error count 
    return df.groupby(['Error']).count()

## Ftn for counting unused aligments from bedpe file 
def countunused(inpath:str,chunksize:int) -> dict:
    ## Format /iniate dictionary of erros
    with pd.read_csv(inpath,sep=hicsep,chunksize=chunksize,usecols=['Error','Inter']) as chunks:
        for i,chunk in enumerate(chunks):
           ## Recursively define error counts
           error_counts = error_counts.add(counterror(chunk),fill_value=0) if i else counterror(chunk)
    ## Make the last column and inteager 
    error_counts['Inter'] = error_counts['Inter'].apply(int)
    ## REturn error counts 
    return error_counts

## Ftn for formating print
def formatp(p:float) -> str:
    return "{:.2f}".format(100*p)

## Load in ftns from total count 
from totalcount import sumcounts, getjson

## Set the help messages
i_help = "The path to an input bedpe file from SLURM.py."
w_help = "Genomic bin size (bp) for calculating distance decay. Default: 10000 bp."
s_help = "Size of span to calculate rolling median."

## SEt font size
myfs        = 12
split_on    = '.valid.'
inter_save  = '.inter.score'
binsize     = 10000
span        = 5
interc      = mycols[:-1]
debugging   = False

## Load in paraamters from slurpy mods hic sep 
from parameters import hicsep, Z_help, chunksize, atac_help
## Bring in diag dir 
from directories import diagdir

## If the script is envoked 
if __name__ == "__main__":
    ## Bring in argparse
    import argparse
    ## call a parser
    parser = argparse.ArgumentParser(description='Calcualtes inter-chromosomal scores and intra-chromosomal distance decay.')
    ## Add arguments
    parser.add_argument("-i", dest="I", type=str, required=True,  help=i_help, metavar="./path/to/input.bedpe")
    parser.add_argument("-w", dest="W", type=int, required=False, help=w_help, metavar="n", default=binsize)
    parser.add_argument("-s", dest="S", type=int, required=False, help=s_help, metavar="n", default=span)
    parser.add_argument("-c", dest="C", type=int, required=False, help=Z_help, metavar='n', default=chunksize)
    parser.add_argument("--atac-seq",             dest="atac",      help = atac_help,     action = 'store_true')

    ## Set inputs
    inputs = parser.parse_args()
    
    ## Set inputs
    inpath    = inputs.I        ## Set input path
    binsize   = inputs.W        ## Set binsize
    span      = inputs.S        ## Span for rolling median filter 
    chunksize = inputs.C        ## Number of rows per chunk 
    inatac    = inputs.atac     ## ATAC-seq mode 

    ## Set output path
    basename = inpath.split('/')[-1].replace('.bedpe',inter_save)
    ## Set out name 
    outname = f'./{diagdir}/{basename}' if exists(f'./{diagdir}') else inpath.replace('.bedpe',inter_save)
    ## Format the csv output path name 
    outpath = outname+'.csv'

    ## Gather total counts if they exist, set sjon paths and calc totals 
    fastp_paths = getjson(diagdir)
    print(f'WARNING: Unable to find fastp logs in {diagdir}!') if not len(fastp_paths) else None 
    pairs_totals = sumcounts(fastp_paths) if len(fastp_paths) else 0

    ## Set paths of duplicates and notused
    not_used_path   = inpath.replace(split_on,'.notused.')
    duplicates_path = inpath.replace(split_on,'.duplicates.') 

    ## Count errors and duplicates 
    error_counts = countunused(not_used_path,chunksize) if exists(not_used_path) else []
    dup_count    = bedpecount(duplicates_path) if exists(duplicates_path) else 0

    ## If error, format the dictionary
    if len(error_counts):
        error_dict = dict(zip(error_counts.index,error_counts.Inter))
        error_dict['duplicates'] = dup_count

    ## Iniate list of chormosomes 
    chrlist      = []

    ## Iniate chromosome list and append ftn 
    def append_to_chrlist(item):
        if item not in chrlist:
            chrlist.append(item)

    ## If in atac seq mode 
    if inatac:
        with pd.read_csv(inpath,sep=hicsep,chunksize=chunksize,usecols=['Pos1','Pos2','End1','End2']) as chunks:
            for i,chunk in enumerate(chunks):
                ## Clac fragments
                frag_size = pd.DataFrame(chunk.max(axis=1) - chunk.min(axis=1),columns=['Fragsize'])
                frag_size['Counts'] = 1
                ## Count
                frag_dist = frag_dist.add(frag_size.groupby('Fragsize').count(),fill_value=0) if i else frag_size.groupby('Fragsize').count()

        ## Modify counts to an int and reset index
        frag_dist['Counts'] = frag_dist.Counts.apply(int)
        frag_dist.reset_index(inplace=True)

        ## Save out the fragmetn dist
        frag_dist.to_csv(outpath.replace(inter_save,'.frag.dist'),index=False)

        ## Calc totals
        totals = frag_dist.Counts.sum()

        ## Call figure 
        fig,ax = plt.subplots(1,1,figsize=(5,5))
        fig.set_facecolor('w')

        ## If errors were counted 
        if len(error_counts):
            ## Add fragment count
            error_dict['Fragments'] = totals

            ## Add a bar showing percent of contcats
            barax = fig.add_axes([-0.1,0.0975,0.05,0.78])
            spinesoff(barax)

            ## Caluclatte total and preset cumlative percent 
            error_sum = sum(error_dict.values())
            total_pairs = pairs_totals if pairs_totals else error_sum

            ## Calc the those removed by fastp 
            fastp_removed = total_pairs - error_sum
            ## add fratp removed
            if fastp_removed:
                error_dict['fastp'] = fastp_removed

            ## Set keys, iniate cumlative percent 
            error_keys  = list(error_dict.keys())
            cumper = 0

            ## Iterrate over the error keys in reverse orer
            for key in error_keys[::-1]:
                ## GAther value and percent 
                value = error_dict[key]
                vper  = round(value/total_pairs,3)
                ## Plot fake point for label
                m, = plt.plot(-1,0,'s',label=f'{key.capitalize()}: {value} ( {formatp(vper)} %)')
                mcolor = m.get_color()
                ## Plot the percentage 
                plt.vlines(0,cumper,cumper+vper,linewidth=50,color=mcolor)
                ## Add to the percentage 
                cumper = cumper + vper 
                    
            ## Modify x-axis 
            plt.xlim(-0.001,0.001)
            plt.xticks([])

            ## Modify y axis 
            plt.ylim(0,1)
            plt.ylabel('Percent of Read Pairs ( %s )'%total_pairs if len(error_counts) else None,fontsize=myfs)

            ## Add legend
            plt.legend(bbox_to_anchor=(-1.75,1),frameon=False,fontsize=myfs-2)
            
            ## Add the total to error dict
            error_dict['Total'] = total_pairs

            ## Save out the counts
            pd.DataFrame(error_dict.values(),index=error_dict.keys()).T.to_csv(outpath.replace(inter_save,'.readpairs'),index=False,header=True)
        
        ## cald the log 10 distance
        log10_dis = np.log10(frag_dist.Fragsize.values)
        frag_freq = frag_dist.Counts/totals

        ## Plot distance decay 
        spinesoff(ax)
        plt.plot(log10_dis,frag_freq,'.',color='tab:blue')

        ## Modify yticks 
        p,l = plt.yticks()
        yticks = ['%s$^{-3}$'%(int(p*1000)) for p in p[1:-1]]
        plt.yticks(p[1:-1],yticks)

        ## Annotate the xaxis 
        plt.xlabel('Fragment Size (log$_{10}$ bp)',fontsize=myfs)
        plt.ylabel('Size Frequency',fontsize=myfs)
        ## Add ttile 
        plt.title(None if len(error_counts) else 'Fragments: %s'%totals,fontsize=myfs)

        ## Saveout the png 
        plt.savefig(outname.replace(inter_save,'.profile')+'.png',dpi=150,bbox_inches='tight')
    else:
        ## Open with with statment, and iterat thru chunks 
        with pd.read_csv(inpath,sep=hicsep,chunksize=chunksize,usecols=mycols) as chunks:
            for i,chunk in enumerate(chunks):
                ## If this is the first chunk, count otherwise add to existing df 
                inter_counts = inter_counts.add(countinter(chunk[interc]),fill_value=0) if i else countinter(chunk[interc]) 
                distances    = distances.add(countdistance(chunk,binsize),fill_value=0) if i else countdistance(chunk,binsize)
                ## Gather the chromosomes 
                chroms = chunk['Rname1'].drop_duplicates().tolist()
                ## Iterate thru the chroms and append to set, keepign order 
                for c in chroms:
                    append_to_chrlist(c)

        ## Set columns to integers
        inter_counts['Inter'] = inter_counts.Inter.apply(int)
        distances['Contacts'] = distances.Contacts.apply(int)

        ## Reset indect of inter counts
        inter_counts.reset_index(inplace=True)
        distances.reset_index(inplace=True)

        ## Print to sreen
        if debugging:
            print(inter_counts.head())
            print()
            print(distances.head())

        ## Count the inter and intra contacts
        ninter = inter_counts.Inter.sum()
        nintra = distances.Contacts.sum()

        ## Calc total and ratio 
        ntotal = ninter + nintra
        rintra = round(nintra/ntotal,4)

        ## Calcluat the inter chromosomal scores
        inter_score_df = interscoring(inter_counts,chrlist)

        ## save the inter scores, counts, and distance decay to csvs givin outpath 
        inter_score_df.to_csv(outpath,index=True,header=True)
        inter_counts.to_csv(outpath.replace('.score','.counts'),index=False)
        distances.to_csv(outpath.replace(inter_save,'.distance.decay'),index=False)

        ## Add a pad of one to counts
        distances['Contacts'] = distances['Contacts'] + 1

        ## Convert log to base 10 
        median_dd = np.log10(distances.Contacts.rolling(span).median().values/distances.Contacts.sum())
        log10_dis = np.log10(distances.Distance.values+1)

        ## Call figure 
        fig,ax = plt.subplots(1,2,figsize=(11,5))
        fig.set_facecolor('w')

        ## Plot distance decay 
        spinesoff(ax[0])
        plt.plot(log10_dis,median_dd,color='tab:blue')
        plt.xlabel('Genomic Distance (log$_{10}$ bp)',fontsize=myfs)
        plt.ylabel('Contact Frequency (log$_{10}$)',fontsize=myfs)

        ## Add a bar showing percent of contcats
        barax = fig.add_axes([-0.015,0.0975,0.05,0.78])
        spinesoff(barax)

        ## If errors were counted 
        if len(error_counts):
            ## Add inter and intra 
            error_dict['Inter'] = ninter
            error_dict['Intra'] = nintra
 
            ## Caluclatte total and preset cumlative percent 
            error_sum = sum(error_dict.values())
            total_pairs = pairs_totals if pairs_totals else error_sum
            cumper = 0

            ## Set keys 
            error_keys  = list(error_dict.keys())
            ## Iterrate over the error keys in reverse orer
            for key in error_keys[::-1]:
                ## GAther value and percent 
                value = error_dict[key]
                vper  = round(value/total_pairs,3)
                ## Plot fake point for label
                m, = plt.plot(-1,0,'s',label=f'{key.capitalize()}: {value} ( {formatp(vper)} %)')
                mcolor = m.get_color()
                ## Plot the percentage 
                plt.vlines(0,cumper,cumper+vper,linewidth=50,color=mcolor)
                ## Add to the percentage 
                cumper = cumper + vper 

            ## Add the total to error dict
            error_dict['Valid'] = error_dict['Inter'] + error_dict['Intra']
            error_dict['Total'] = total_pairs

            ## Save out the counts
            pd.DataFrame(error_dict.values(),index=error_dict.keys()).T.to_csv(outpath.replace(inter_save,'.contacts'),index=False,header=True)
        
        ## Otherwise, just plot hte inter:intra contact ratio 
        else: ## Add dummby plotting values 
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
        sns.heatmap(np.log2(inter_score_df),cbar_ax=cbar,cmap='cividis',vmin=-2,vmax=2)
        ## Modify tick marks 
        plt.xticks(fontsize=10);plt.yticks(fontsize=myfs-2)

        ## Set color bar labels 
        plt.sca(cbar)
        plt.ylabel('log$_2$ (Interaction Score)',fontsize=myfs)

        ## Saveout the png 
        plt.savefig(outname.replace(inter_save,'.profile')+'.png',dpi=150,bbox_inches='tight')
    ## EOF      