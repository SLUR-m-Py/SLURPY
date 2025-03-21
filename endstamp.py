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
import time, sys
## Load in os path getctime
from os.path import getctime
## Load in date and time
from datetime import datetime
## load in write to file from our pysamtools 
from pysamtools import writetofile
## Load in sort glbo 
from defaults import sortglob
## load in log dir
from directories import debugdir

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

## Print to log
print("Finished our SLUR(M)-py!")
## End of file 