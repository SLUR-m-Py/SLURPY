#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
#######################################
##         TIME STAMP (v 0.0.0)      ##
#######################################
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## ----------------------------------------- LOAD IN MODULES ------------------------------------------ ## 
## Load time and sys mods 
import time, sys, pandas as pd, json 
## Load in os path getctime
from os.path import getctime
## Load in date and time
from datetime import datetime
## load in write to file from our pysamtools 
from pysamtools import writetofile
## Load in sort glbo 
from defaults import sortglob
## load in log dir
from directories import debugdir, diagdir

## Ftn for reading in json
def loadjson(jsonpath):
    # Open the JSON file in read mode
    with open(jsonpath, 'r') as file:
        # Load the JSON data from the file
        json_data = json.load(file)
    return json_data

## Get the total data counts 
def gettotal(data) -> int: 
    return int(data['summary']['before_filtering']['total_reads'])

def getfiltered(data) -> int: 
    return int(data['summary']['after_filtering']['total_reads'])

## If the script is envoked 
if __name__ == "__main__":
    ## ------------------------------------------ Input Variables ------------------------------------------ ##
    ## Gather the file name, the start time stamp from input and calculate the end time stamp 
    filename,sstamp,estamp = sys.argv[1], float(sys.argv[2]), time.time()

    ## Format the start time stamp 
    dt0 = datetime.fromtimestamp(sstamp)

    ## Format the end time stamp
    dt2 = datetime.fromtimestamp(estamp)

    ## Gather logs from dir
    fastp_log = sortglob(f'./{debugdir}/*.fastp.*.log')

    ## Format tiem stamp from fastp log if we find one 
    dt1 = datetime.fromtimestamp(getctime(fastp_log[0])) if len(fastp_log) else dt0 

    ## Calculate delta of time stamps and then the approx run time
    runtime = dt2 - dt1 #runtime = datetime.utcfromtimestamp(estamp - sstamp).strftime('%Y-%m-%d %H:%M:%S').split(' ')[-1]

    ## Write the timestamp to file 
    writetofile(filename,['SUBMIT TIME: %s\n'%dt0, 'START TIME: %s\n'%dt1, 'END TIME: %s\n'%dt2, 'RUN TIME: %s\n'%runtime],False)

    ## Gather filtering and duplicate logs
    filter_logs = sortglob(f'./{debugdir}/*.filter.bedpe.*.log')
    duplit_logs = sortglob(f'./{debugdir}/*.dedup.*.log')

    ## Gather counts across logs 
    mappings = pd.concat([pd.read_csv(f,sep='\t',header=None,names=['Mapping','Counts']).dropna() for f in filter_logs]).groupby('Mapping').sum()
    duplicat = pd.concat([pd.read_csv(f,sep='\t',header=None,names=['Mapping','Counts']).dropna() for f in duplit_logs]).groupby('Mapping').sum()

    ## Copy new map
    newmap = pd.concat([mappings,duplicat],axis=0)
    ## Adjust for duplicates
    newmap.loc['INFO: InterHiC','Counts'] = newmap.loc['INFO: InterHiC','Counts'] - duplicat.loc['INFO: InterDuplicates','Counts']
    newmap.loc['INFO: IntraHiC','Counts'] = newmap.loc['INFO: IntraHiC','Counts'] - duplicat.loc['INFO: IntraDuplicates','Counts']
    newmap.loc['INFO: Valid',   'Counts'] = newmap.loc['INFO: Valid',   'Counts'] - (duplicat.loc['INFO: InterDuplicates','Counts'] + duplicat.loc['INFO: IntraDuplicates','Counts'])
    ## Reset the index
    newmap.reset_index(inplace=True)
    ## Reet column name
    newmap['Mapping'] = [v.split('INFO: ')[-1] for v in newmap.Mapping.tolist()]


    ## Check inital (zeroth) fastp logs for total
    initial_fastp_logs = sortglob(f'./{diagdir}/*fastp.*.0.json') 
    ## Calculate total read pairs 
    total_counts = int(sum([gettotal(loadjson(fastlog)) for fastlog in initial_fastp_logs])/2)
    fastp_lost   = int(sum([gettotal(loadjson(fastlog)) - getfiltered(loadjson(fastlog)) for fastlog in initial_fastp_logs])/2)
   
    ## Transpose and add fastp counts
    newmap = newmap.T
    newmap['fastp'] = ['Fastp',fastp_lost]
    newmap['total'] = ['Total',total_counts]

    ## Retranspose newmap
    newmap = newmap.T
    newmap.reset_index(drop=True,inplace=True)
    ## Calc new sum 
    newsum = newmap[(~newmap.Mapping.isin(['Valid','Total']))].Counts.sum()

    ## Check our work 
    if abs(newsum  - total_counts):
        print("WARNING: The total number of counts did not match after accounting for those removed by fastp: %s - %s"%(newsum,total_counts))

    ## Set outfile name
    outfilename = filename.split('.time')[0] + '.slurpy.counts.csv'

    ## Format counts column into integers
    newmap['Counts'] = newmap.Counts.apply(int)
    ## ADd the percent to the new map
    newmap['Percent (Total)'] = [round(r,5) for r in newmap.Counts.values/newmap[(newmap.Mapping=='Total')].Counts.max()]
    newmap['Percent (Valid)'] = [round(r,5) for r in newmap.Counts.values/newmap[(newmap.Mapping=='Valid')].Counts.max()]

    ## Save out new map
    newmap.to_csv(outfilename,index=False)

    ## Print to log
    print("Finished our SLUR(M)-py!")
    ## End of file 