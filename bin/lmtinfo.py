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


This routine grabs some useful summary information from the ifproc file, ignoring
the roach files. If one unique OBSNUM is given, it will show this information
in a "rc" style for the pipeline. If more OBSNUM are possible, for example by only
giving a PATH, all possible OBSNUMs are listed with a short summary, one OBSNUM
per line. Example of output:

      #     DATE  OBSNUM   SOURCE     RESTFRQ VLSR INTTIME
      2018-11-16  079447  IRC+10216   115.271  -20       8
      2018-11-16  079448  IRC+10216   115.271  -20     686
      2020-02-18  090910  NGC5194     115.271  463       7
      2020-02-18  090911  NGC5194     115.271  463    3986
      2020-02-20  091111  NGC5194     115.271  463       7
      2020-02-20  091112  NGC5194     115.271  463    6940

"""

import sys
import math
import numpy as np		
import matplotlib.pyplot as pl
import datetime

import sys, os
import glob
import netCDF4

from docopt import docopt

#  ifproc/ifproc_2018-06-29_078085_00_0001.nc
#  spectrometer/roach0/roach0_78085_0_1_CHI-Cyg_2018-06-29_041713.nc
#  RedshiftChassis0/RedshiftChassis0_2015-01-22_033551_00_0001.nc

def slr_summary(ifproc, rc=False):
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
    # Header.Dcs.ObsNum 
    obspgm = b''.join(nc.variables['Header.Dcs.ObsPgm'][:]).decode().strip()
    # the following Map only if obspgm=='Map'
    if obspgm=='Map':
        xlen = nc.variables['Header.Map.XLength'][0] * 206264.806
        ylen = nc.variables['Header.Map.YLength'][0] * 206264.806
        hpbw = nc.variables['Header.Map.HPBW'][0]

    date_obs = nc.variables['Data.TelescopeBackend.TelTime'][0].tolist()
    date_obs = datetime.datetime.fromtimestamp(date_obs).strftime('%Y-%m-%dT%H:%M:%S')
    
    
    t0 = float(bbtime[0].data)
    t1 = float(bbtime[-1].data)
    dt = t1-t0
    nc.close()

    if rc:
        print('# <lmtinfo>')
        print('# ifproc="%s"' % ifproc)
        print('# date-obs="%s"' % date_obs)
        print('# inttime=%g sec' % dt)
        print('# obspgm="%s"' % obspgm)
        print('vlsr=%g' % vlsr)
        print('skyfreq=%g' % skyfreq)
        print('restfreq=%g' % restfreq)
        print('src="%s"' % src)
        resolution = math.ceil(1.0 * 299792458 / skyfreq / 1e9 / 50.0 * 206264.806)
        print('resolution=%g' % resolution)
        print('cell=%g' % (resolution/2.0))
        print('x_extent=%g' % xlen)
        print('y_extent=%g' % ylen)
        
        print("# </lmtinfo>")
    else:    
        print("%s %s  %-5s %-20s %g %g %g" % (date_obs, fn[2], obspgm, src, restfreq, vlsr, dt))


def rsr_summary(rsr_file, rc=False):
    # RedshiftChassis2/RedshiftChassis2_2015-01-22_033551_00_0001.nc
    nc = netCDF4.Dataset(rsr_file)

    # Header.Source.SourceName
    src = b''.join(nc.variables['Header.Source.SourceName'][:]).decode().strip()
    
    # Header.Dcs.ObsNum = 33551 ;
    obsnum = nc.variables['Header.Dcs.ObsNum'][0] 
    
    # Header.Radiometer.UpdateDate = "21/01/2015 23:12:07
    date_obs = b''.join(nc.variables['Header.Radiometer.UpdateDate'][:]).decode().strip()
    
    # Header.Weather.UpdateDate = "22/01/15 0:39:48
    
    nc.close()

    # no rc mode, only one line summary
    print("%s  %d  RSR   %s" %   (date_obs, obsnum, src))


#  although we grab the command line arguments here, they are actually not
#  used in the way most scripts use them. Below there is a more hardcoded
#  parsing of arguments based on how many there are, which are files, and
#  which are directories.
arguments = docopt(__doc__,options_first=True, version='0.1')
#print(arguments)

if len(sys.argv) == 2:
                                                     # mode 1: obsnum or nc_file or path
    obsnum = sys.argv[1]
    fn = glob.glob('*/ifproc/ifproc*%s*.nc' % obsnum)
    if len(fn) > 0:
        ifproc = fn[0]
    else:
        ifproc = sys.argv[1]

    if os.path.isdir(ifproc):
        path = ifproc

        chassis = -1
        if chassis < 0:
            globs = '%s/RedshiftChassis?/RedshiftChassis?_*.nc'  % (path)
        else:
            globs = '%s/RedshiftChassis%d/RedshiftChassis%d_*.nc'  % (path,chassis,chassis)            
        fn = glob.glob(globs)
        for f in fn:
            # print('RSR',f)
            try:
                rsr_summary(f)
            except:
                print("Failing on ",f)
            

        globs = '%s/ifproc/ifproc*.nc' % path
        fn = glob.glob(globs)
        for f in fn:
            try:
                slr_summary(f)
            except:
                try:
                    yyyymmdd = f.split('/')[-1].split('_')[1]
                    obsnum   = f.split('/')[-1].split('_')[2]
                    print("%s          %s   failed on %s" % (yyyymmdd,obsnum,f))
                except:
                    print("1900-00-00       failed on %s" % f)

        sys.exit(0)
    elif os.path.exists(ifproc):
        try:
            slr_summary(ifproc,rc=True)
        except:
            print("%s: failed" % ifproc)
               
elif len(sys.argv) == 3:
                                                     # mode 2: path and obsnum : differentiate between SLR and RSR
    path = sys.argv[1]
    obsnum = sys.argv[2]
    globs = '%s/ifproc/ifproc*%s*.nc' % (path,obsnum)
    fn = glob.glob(globs)
    if len(fn) > 0:
        ifproc = fn[0]
        slr_summary(ifproc,True)
    else:
        globs = '%s/RedshiftChassis?/RedshiftChassis?_*%s*.nc'  % (path,obsnum)
        print("Trying RSR %s" % globs)
        fn = glob.glob(globs)
        if len(fn) > 0:
            for f in fn:
                print(f)
        else:
            print("Warning - no RSR files found")
else:
                                                     # no other modes
    print("Usage : %s [path] obsnum" % sys.argv[0])
    sys.exit(0)


