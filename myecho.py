#!/usr/bin/env python
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
"""
## Import sys 
import sys
## Bring in write to file 
from pysamtools import writetofile

## Set the output message and output file 
outmessage, outfile = ' '.join(sys.argv[1:-1]), sys.argv[-1]

## Write to file
writetofile(outfile,outmessage,debug=False,mode='a')
## End of file 