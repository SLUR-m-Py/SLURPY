#!/usr/bin/env python
## List command for canceling all
"""
squeue -u croth | grep 'croth' | awk '{print $1}' | xargs -n 1 scancel
""" 
## Bring in mods 
import sys, os
## Set user nmae
user_name = sys.argv[1]

## Cancel command
can_comand = 'squeue -u %s | grep \'%s\' | awk \'\{print $1\}\' | xargs -n 1 scancel'%user_name
## Print the command 
print(can_comand)
## Submit to os 
os.system(can_comand)
## Eof 