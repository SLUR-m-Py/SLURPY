#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""

## Load in debug dir 
from directories import debugdir, diagdir
## Load in matplot lib
from matplotlib import pyplot as plt 
## Bring in pandas, numpy, seaborn 
import pandas as pd, numpy as np, seaborn as sns 
## Bring in sortglob
from defaults import sortglob, submitter, dictzip

## Set memory sclae
memory_dict = dictzip(['GB','MB','KB'],[1,1/(10**3),1/(10**6)])

## Set job id 
by_jobid = 'Job ID: '

## Ftn for gathering job ids from bwa and filter submitter 
def subjobids(inpath) -> list[int]:
    ## iniate ids
    ids = []
    ## Open the infile 
    with open(inpath,'r') as infile:
        ## Iterate thru lines, if the jobid str in line
        for l in infile:
            if (by_jobid in l):
                ids.append(int(l.split(by_jobid)[-1])) 
        ## Close file
        infile.close()
    ## Return ids 
    return ids 

## Set the memory utlized str 
mem_str = 'Memory Utilized: '

## Ftn for parsing a seff report
def parseseff(inreport) -> str: 
    mem = ''
    ## Iterate thru lines of the reposrt if the mem string is there return it 
    for l in inreport:
        if mem_str in l:
            mem = l.split(mem_str)[-1] 
    ## Return the mem
    return mem 

## Ftn for getting utilized memory from slurm's seff command
def getmemory(jobid:int) -> str:
    ## Gather the used mem 
    return parseseff(submitter(f'seff {jobid}')) if jobid else '0 KB'

## if the script is envoked
if __name__ == "__main__":

    ## Load in command files
    command_file_paths = sortglob(f'./{debugdir}/command.file.*.csv')
    ## Iterate thru the command file
    command_df = pd.concat([pd.read_csv(inpath) for inpath in command_file_paths],axis=1)
    ## Gahter only those commands ran
    command_df = command_df[(command_df.Torun==0)]
    ## Parse the name of the pipeline step 
    command_df['Name'] = [a.split('/')[-1].split('.')[1] for a in command_df.Jobfile]
    ## Gather job memory
    command_df['Memory Utilized'] = [getmemory(int(jid)) for jid in command_df.JobID.replace(np.nan,0).tolist()]
    ## convert memory values to GB
    command_df['GB'] = [ float(k.split(' ')[0])*memory_dict[k.split(' (')[0].split(' ')[-1]] for k in command_df['Memory Utilized'] ]

    ## Gather bwa jobs
    bwa_master_paths = sortglob(f'./{debugdir}/1.bwa.to.bedpe.*.log')
    ## Gather job ides then memory, then convert to GB
    bwa_job_ids = np.concatenate([subjobids(k) for k in bwa_master_paths])
    bwa_job_mem = [getmemory(int(jid)) for jid in bwa_job_ids]
    bwa_job_GB  = pd.DataFrame([('bwa mem',float(k.split(' ')[0])*memory_dict[k.split(' ')[-1]]) 
                                for k in bwa_job_mem],columns=['Name','GB'])
    
    ## Gather filter jobs
    filt_master_paths = sortglob(f'./{debugdir}/2.filter.master.*.log')
    ## Gather job ides
    filt_job_ids = np.concatenate([subjobids(k) for k in  filt_master_paths])
    filt_job_mem = [getmemory(int(jid)) for jid in filt_job_ids]
    filt_job_GB  = pd.DataFrame([('filtering',float(k.split(' ')[0])*memory_dict[k.split(' ')[-1]]) 
                                for k in filt_job_mem],columns=['Name','GB'])
    
    ## concatonate the gb df 
    gb_df = pd.concat([command_df[['Name','GB']],
                       bwa_job_GB[['Name','GB']],
                       filt_job_GB[['Name','GB']]]).reset_index(drop=True)
    
    ## Save out the dataframe
    out_path = f'./{diagdir}/Memory.profile.csv'
    gb_df.to_csv(out_path,index=False,header=True)

    ## Plot the memory profile
    sns.barplot(x='Name',y='GB',data=gb_df)
    plt.ylabel('Gigabytes',fontsize=12)
    plt.xlabel('SLUR(M)-py Pipeline Step',fontsize=12)
    plt.xticks(rotation=90,fontsize=10)
    ## SAve figure 
    plt.savefig(out_path.replace('.csv','.png'),dpi=150,bbox_inches='tight')
    ## Print we have completed this
    print('Finished calculating memory usage.')
## EOF 