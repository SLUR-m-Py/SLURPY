#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## Load in rm treat 
from shutil import rmtree
## Laod in defaults 
from directories import splitsdir, hicdir, bedtmpdir, checkerdir, comsdir

## ----------------------------------------------- MAIN EXECUTABLE --------------------------------------------------- ## 
## If the library is called as an executable
if __name__ == "__main__":
    ## Bring in all the files being gzipped, Iterate thru the file pairs, and G-zip the file with our ftn
    [rmtree(directory,ignore_errors=True) for directory in [splitsdir, hicdir, bedtmpdir, comsdir, checkerdir]]
## End of file 