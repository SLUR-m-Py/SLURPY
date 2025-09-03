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
## Load in rm treat 
from shutil import rmtree, copyfileobj
## Bring in sorted glob 
from defaults import sortglob, fileexists, submitsbatch, readtable, ifprint
## Load in debug dir 
from parameters import debugdir, splitsdir, hicdir, bedtmpdir, checkerdir, comsdir, macs3dir, aligndir, diagdir
## bring in numpy 
import numpy as np, sys, gzip
## laod in remove
from os import remove
## -------------------------------------------------------------------- ##

## -------------------------------------------------------------------- ##
## Write vars
directories = [splitsdir, hicdir, bedtmpdir, checkerdir]
## Set the error message
no_file_error = "WARNING: The input file path -- %s -- could not be found!"
rm_file_error = "ERROR: Unable to gzip and remove input file: %s"

## Define ftn for gzippin 
def gzipfile(infilepath:str,outfilepath:str):
    """
    Given an input file string (INFILEPATH), this functions "gzips" the contents into a clone file (OUTFILEPATH) with the extension .gz."
    """
    ## Check our work
    if fileexists(infilepath): 
        ## Open the input file as binar 
        with open(infilepath, 'rb') as f_in:
            ## open the output gzip file 
            with gzip.open(outfilepath, 'wb') as f_out:
                ## Copy out as zipped file 
                copyfileobj(f_in, f_out)

        ## Clost the files 
        f_in.close()
        f_out.close()

        ## Make sure we have the out file, print an error message if the file exists
        remove(infilepath) if fileexists(outfilepath) else print(rm_file_error%infilepath)
    pass 

##
##      Define Functions 
## Write function for checking for the word error in txt
def checkforerror(inpath:str) -> bool:
    """Checks for the word error within text from a given input path. Returns boolean upon first occurance."""
    with open(inpath,'r') as inhandle:
        for l in inhandle:
            ll = l.lower()
            k = (('error' in ll) or ('no such file' in ll) or ('permission denied' in ll)) and (not 'errors.bedpe' in ll) and (not ll.startswith('(null)'))
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
            k = ('warning' in ll) and (not ('runtimewarning' in ll))
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
## If the script is envoked by name 
if __name__ == "__main__":
    ## Check if the checks and splits match
    ## GAther counts of fastq splits,  bwa checks, bedpe checks
    nsplits       = len(sortglob(f'./{splitsdir}/*_R1_*fastq.gz'))
    nbwa_checks   = len(sortglob(f'./{checkerdir}/*.bwa.log'))
    nbedpe_checks = len(sortglob(f'./{checkerdir}/*.bedpe.log'))

    ## Print error to screen
    assert nbwa_checks >= nsplits,   'ERROR: The number of parallele BWA MEM runs (%s) did not match the nubmer of splits (%s)'%(nbwa_checks,nsplits)
    assert nbedpe_checks >= nsplits, 'ERROR: The number of filtering runs (%s) did not match the nubmer of splits (%s)'%(nbedpe_checks,nsplits)

    ## Check if we have an argument
    if len(sys.argv) > 1:
        ## Gather the input param
        param = sys.argv[1].lower()
        ## decide what we are doing
        ## Cleaning
        if param.startswith('cl'): ## Iterate thru the directoris and remove them if clean was passed
            [rmtree(direct,ignore_errors=True) for direct in directories]
        ## Zipping
        elif param.startswith('g') or ('zip' in param):
            ## gzip the passed files
            [gzipfile(tozip,tozip+'.gz') for tozip in sys.argv[2:]]
        ## Resetting
        elif param.startswith('r') or param.startswith('q') or ('reset' in param) or ('quick' in param):
            [rmtree(d,ignore_errors=True) for d in [debugdir, splitsdir, hicdir, bedtmpdir, checkerdir, comsdir, macs3dir, aligndir, diagdir]]
        ## Cancle the jobs
        elif param.startswith('can') or ('cancel' in param):
            ## Set command file path 
            job_file = f'./{debugdir}/command.file.csv'
            ## Load in job file 
            jobs = readtable(job_file,header=0)
            ## Gather job ids
            jobids = jobs[(jobs.Operation!='bwaix')].JobID.tolist()

            ## Make iterater
            i = 0
            ## Iterate thrut he job ids and canclse them
            for job in jobids:
                try:
                    submitsbatch(f'scancel {job}',returnid=False)
                    i += 1
                except Exception as error:
                    print(error)

            ## Print the nubmer of jobs canned
            print("INFO: Cancelled %s jobs assoiated with this directory."%i)
        ## Otherwise
        else:
            print('WARNING: Arguments were passed but no action was taken.')

    else: ## otherwize check work
        ## -------------------------------------------------------------------- ##
        ##      ANALYSIS of ERROR LOGS
        ## Bring in the error logs 
        all_error_logs = sortglob(f'./{debugdir}/*.log')
        ## Filter the error logs for those with text only within them
        error_logs = np.array([k for k in all_error_logs if fileexists(k)])
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

        if np.sum(warns_counts):
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

        if np.sum(unfin_counts):
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

        else: ## Otherwise, print to screen 
            print("INFO: Good news everyone! No errors were detected in this run.\n\t:-)\n\t<3")
## End of file 