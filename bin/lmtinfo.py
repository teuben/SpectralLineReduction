#! /usr/bin/env python
#

import sys
import numpy as np		
import matplotlib.pyplot as pl

import sys, os
import glob
import netCDF4


if len(sys.argv) == 2:
                                                     # mode 1: obsnum or nc_file
    obsnum = sys.argv[1]
    fn = glob.glob('*/ifproc/ifproc*%s*.nc' % obsnum)
    if len(fn) > 0:
        ifproc = fn[0]
    else:
        ifproc = sys.argv[1]
elif len(sys.argv) == 3:
                                                     # mode 2: path and obsnum
    path = sys.argv[1]
    obsnum = sys.argv[2]
    fn = glob.glob('%s/ifproc/ifproc*%s*.nc' % (path,obsnum))
    if len(fn) > 0:
        ifproc = fn[0]
    else:
        print("Warning - no ifproc file found")
        sys.exit(0)        
else:
    print("Usage : %s [path] obsnum" % sys.argv[0])
    sys.exit(0)


nc = netCDF4.Dataset(ifproc)
vlsr = nc.variables['Header.Source.Velocity'][0]
src = b''.join(nc.variables['Header.Source.SourceName'][:]).decode().strip()
skyfreq  = nc.variables['Header.Sequoia.SkyFreq'][0]
restfreq = nc.variables['Header.Sequoia.LineFreq'][0]
nc.close()


print("# <lmtinfo>")
print('vlsr=%g' % vlsr)
print('skyfreq=%g' % skyfreq)
print('restfreq=%g' % restfreq)
print('src=%s' % src)
resolution = 1.15 * 299792458 / skyfreq / 1e9 / 50.0 * 206264.806
print('resolution=%.2f' % resolution)
print('cell=%.2f' % (resolution/2.0))
print("# </lmtinfo>")



