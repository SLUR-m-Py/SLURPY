#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## ----------------------------- v 0.0.0 ------------------------------ ## 
##          CHECKWORK: Checks error logs from runs of SLURPY            ##
## -------------------------------------------------------------------- ##
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## -------------------------------------------------------------------- ##
##      MODULE LOADING 
## Bring in sorted glob 
from defaults import sortglob, getfilesize
from directories import debugdir
## bring in numpy 
import numpy as np 
## -------------------------------------------------------------------- ##

## -------------------------------------------------------------------- ##
##
##      Define Functions 
## Write function for checking for the word error in txt
def checkforerror(inpath:str) -> bool:
    """Checks for the word error within text from a given input path. Returns boolean upon first occurance."""
    with open(inpath,'r') as inhandle:
        for l in inhandle:
            ll = l.lower()
            k = (('error' in ll) or ('no such file' in ll) or ('permission denied' in ll)) and (not 'errors.bedpe' in ll)
            if k:
                break 
            else:
                pass 
    ## Return the check 
    return k

## Write function for checking for the word error in txt
def getwarnings(inpath:str) -> bool:
    """Checks for the word warning within text from a given input path. Returns boolean upon first occurance."""
    with open(inpath,'r') as inhandle:
        for l in inhandle:
            ll = l.lower()
            k = ('warning' in ll)
            if k:
                break 
            else:
                pass 
    ## Return the check 
    return k

## Writing function to check for unfinsihed logs
def unfinished(inpath:str) -> bool:
    with open(inpath,'r') as inhandle:
        for l in inhandle:
            pass 
    ## Return those that have not finsihed
    return not l.lower().startswith('finished')

## Ftn for defining logs name
def logsvslog(n) -> str:
    return 'logs' if n > 1 else 'log'
## -------------------------------------------------------------------- ##

## -------------------------------------------------------------------- ##
##      ANALYSIS of ERROR LOGS
## Bring in the error logs 
all_error_logs = sortglob(f'./{debugdir}/*.log')
## Filter the error logs for those with text only within them
error_logs = np.array([k for k in all_error_logs if getfilesize(k)])
## Gather the sizes 
error_counts = np.array([checkforerror(k) for k in error_logs])
warns_counts = np.array([getwarnings(k) for k in error_logs])
unfin_counts = np.array([unfinished(k) for k in error_logs])
## -------------------------------------------------------------------- ##


## -------------------------------------------------------------------- ##
##      CHECK ERROR STATUS
## Check if we have non-zeros
if np.sum(error_counts):
    ## Gather the logs with erros 
    has_errors = error_logs[error_counts]
    ## calc total erros
    total_erros = len(has_errors)

    ## Set if log vs logs based on error count
    lvls = logsvslog(total_erros)

    ## Print the number of warnings
    print(f'WARNING: Errors were detected in {total_erros} {lvls}!')
    print(f'WARNING: Check the following {lvls}:\n')
    ## Print the logs with erros
    [print("\t%s"%f) for f in has_errors]
    print(f'WARNING: Check the above logs for errors.')

elif np.sum(warns_counts):
    ## Gather lots with warnings 
    has_warnings = error_logs[warns_counts]
    ## Calc total warnigns
    total_warns = len(has_warnings)

    ## Set if log vs logs based on error count
    lvls = logsvslog(total_warns)

    ## Print the number of warnings
    print(f'INFO: Warnings were detected in {total_warns} {lvls}!')
    print(f'INFO: Check the following {lvls}:\n')
    ## Print the logs with erros
    [print("\t%s"%f) for f in has_warnings]
    print(f'INFO: Check the above logs for warnings (may be nothing).')

elif np.sum(unfin_counts):
    ## Agther unfin counts
    is_unfin = error_logs[unfin_counts]
    ## Calc totla number un fin
    total_unfin = len(is_unfin)

    ## Set level of log
    lvls = logsvslog(total_unfin)

    ## Print the number of warnings
    print(f'WARNING: We detected {total_unfin} {lvls} with unfinished marks!')
    print(f'WARNING: Check the following {lvls}:\n')
    ## Print the logs with erros
    [print("\t%s"%f) for f in is_unfin]
    print(f'WARNING: Check the above logs for unfinished processes.')

else: ## Otherwise
    print("INFO: Good news everyone! No errors were detected in this run.\n\t:-)\n\t<3")
## -------------------------------------------------------------------- ##
## End of file 