#!/usr/bin/env python
## Import sys 
import sys
## append path
sys.path.append('./SLURPY/pipeline')
## Bring in write to file 
from defaults.tools.pysamtools import writetofile

## Set the output message and output file 
outmessage, outfile = ' '.join(sys.argv[1:-1]), sys.argv[-1]

## Write to file
writetofile(outfile,outmessage,debug=False,mode='a')
## End of file 