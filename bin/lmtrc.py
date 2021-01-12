#! /usr/bin/env python
#
#   lmtrc.py  rcfile key=val [key2=val2....] > rcfile2
#
#   lmtrc.py  rcfile [rcfile2 ...]  key=val [key2=val2....]
#    

import sys
import math
import numpy as np		
import os
import glob


if len(sys.argv) < 3:
    print("usage: %s rcfile key=val [key=val...]" % sys.argv[0])
    sys.exit(0)


rcfile = sys.argv[1]
keyval = []
for argv in sys.argv[2:]:
    kv = argv.split('=')
    if len(kv) == 2:
        keyval.append([kv[0],kv[1]])


lines = open(rcfile).readlines()
for line in lines:
    line1 = line.strip()
    n=0
    for kv in keyval:
        if line1.find("%s=" % kv[0])==0:
            n = n + 1
            print("%s=%s" % (kv[0],kv[1]))
    if n==0:
        print(line1)
        
        

