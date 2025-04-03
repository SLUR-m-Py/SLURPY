#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## Bring in shutil and gzip 
import gzip, shutil, sys 

## Remove the old file
## Bring file exists
from os.path import exists
from os import remove

## Set the error message
no_file_error = "WARNING: The input file path -- %s -- could not be found!"
rm_file_error = "ERROR: Unable to gzip and remove input file: %s"

## Define ftn for gzippin 
def gzipfile(infilepath:str,outfilepath:str):
    """
    Given an input file string (INFILEPATH), this functions "gzips" the contents into a clone file (OUTFILEPATH) with the extension .gz."
    """
    ## Check our work
    if exists(infilepath): 
        ## Open the input file as binar 
        with open(infilepath, 'rb') as f_in:
            ## open the output gzip file 
            with gzip.open(outfilepath, 'wb') as f_out:
                ## Copy out as zipped file 
                shutil.copyfileobj(f_in, f_out)

        ## Clost the files 
        f_in.close()
        f_out.close()

        ## Make sure we have the out file, print an error message if the file exists
        remove(infilepath) if exists(outfilepath) else print(rm_file_error%infilepath)
    #else:
        #print(no_file_error%infilepath)
    pass 

## ----------------------------------------------- MAIN EXECUTABLE --------------------------------------------------- ## 
## If the library is called as an executable
if __name__ == "__main__":
    ## Bring in all the files being gzipped, Iterate thru the file pairs, and G-zip the file with our ftn
    [gzipfile(tozip,tozip+'.gz') for tozip in sys.argv[1:]]
## End of file 