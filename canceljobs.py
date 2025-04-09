#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## Load in pandas 
import pandas as pd
## Bring in submitter
from defaults import submitter
## Load in dir
from directories import debugdir

## Set command file path 
job_file = f'./{debugdir}/command.file.csv'
## Load in job file 
jobs = pd.read_csv(job_file)
## Gather job ids
jobids = jobs[(jobs.Operation!='bwaix')].JobID.tolist()

## Make iterater
i = 0
## Iterate thrut he job ids and canclse them
for job in jobids:
    try:
        submitter(f'scancel {job}')
        i += 1
    except Exception as error:
        print(error)

## Print the nubmer of jobs canned
print("INFO: Cancelled %s jobs assoiated with this directory."%i)
## End of file 