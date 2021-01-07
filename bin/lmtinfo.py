#! /usr/bin/env python
#
#    lmtinfo:    provide some info, in terms of a list of "keyword=value", that could be used
#                by an external pipeline that needs input.
#                This list is both bash and python friendly: we only allow integer/float/string
#

"""
Usage: lmtinfo.py OBSNUM
       lmtinfo.py IFPROCFILE
       lmtinfo.py PATH OBSNUM
       lmtinfo.py PATH

-h --help  This help


"""

import sys
import numpy as np		
import matplotlib.pyplot as pl

import sys, os
import glob
import netCDF4

from docopt import docopt

def summary(ifproc, rc=False):
    """   summary a procnum in a single line
    """
    #   e.g.  M51_data/ifproc/ifproc_2020-02-20_091111_00_0001.nc
    fn  = ifproc.split("/")[-1].replace('.nc','').split('_')
    #   e.g. ['ifproc', '2020-02-20', '091111', '00', '0001']
    
    nc = netCDF4.Dataset(ifproc)
    vlsr = nc.variables['Header.Source.Velocity'][0]
    src = b''.join(nc.variables['Header.Source.SourceName'][:]).decode().strip()
    skyfreq  = nc.variables['Header.Sequoia.SkyFreq'][0]
    restfreq = nc.variables['Header.Sequoia.LineFreq'][0]
    bbtime = nc.variables['Data.IfProc.BasebandTime']
    t0 = float(bbtime[0].data)
    t1 = float(bbtime[-1].data)
    dt = t1-t0
    nc.close()

    if rc:
        print('# <lmtinfo>')
        print('# ifproc="%s"' % ifproc)
        print('# inttime=%g sec' % dt)
        print('vlsr=%g' % vlsr)
        print('skyfreq=%g' % skyfreq)
        print('restfreq=%g' % restfreq)
        print('src="%s"' % src)
        resolution = 1.15 * 299792458 / skyfreq / 1e9 / 50.0 * 206264.806
        print('resolution=%.2f' % resolution)
        print('cell=%.2f' % (resolution/2.0))
        print("# </lmtinfo>")
    else:    
        print("%s %s  %-20s %g %g %g" % (fn[1], fn[2], src, restfreq, vlsr, dt))

arguments = docopt(__doc__,options_first=True, version='0.1')
#print(arguments)

if len(sys.argv) == 2:
                                                     # mode 1: obsnum or nc_file
    obsnum = sys.argv[1]
    fn = glob.glob('*/ifproc/ifproc*%s*.nc' % obsnum)
    if len(fn) > 0:
        ifproc = fn[0]
    else:
        ifproc = sys.argv[1]

    if os.path.isdir(ifproc):
        path = ifproc
        fn = glob.glob('%s/ifproc/ifproc*.nc' % path)
        for f in fn:
            summary(f)
        sys.exit(0)
    elif os.path.exists(ifproc):
        summary(ifproc,rc=True)
               
elif len(sys.argv) == 3:
                                                     # mode 2: path and obsnum
    path = sys.argv[1]
    obsnum = sys.argv[2]
    fn = glob.glob('%s/ifproc/ifproc*%s*.nc' % (path,obsnum))
    if len(fn) > 0:
        ifproc = fn[0]
        summary(ifproc,True)
    else:
        print("Warning - no ifproc file found")
        sys.exit(0)        
else:
                                                     # no other modes
    print("Usage : %s [path] obsnum" % sys.argv[0])
    sys.exit(0)


