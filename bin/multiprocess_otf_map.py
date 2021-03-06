#!/usr/bin/env python
'''
Creates a SpecFile from a single OTF mapping observation
'''

# Python Imports	
import sys
import numpy as np		
import matplotlib.pyplot as pl
import netCDF4			 
import multiprocessing

# Line Data Reduction Imports
from lmtslr.spec.spec import *
from lmtslr.ifproc.ifproc import *

from lmtslr.reduction.line_reduction import *
from lmtslr.grid.grid import *

from lmtslr.utils.reader import read_obsnum_otf_multiprocess, count_otf_spectra

#from lmtslr.utils.parser import HandleProcessOptions
from lmtslr.utils.argparser import HandleProcessOptions

import time

def multiprocess_otf_map(Opts, I, ICal):

    # set up the grid geometry
    theGrid = Grid()
    
    # check to see whether output file exists and remove it if it does
    if os.path.isfile(Opts.output_file_name) == True:
        os.remove(Opts.output_file_name) 

    S = read_obsnum_otf_multiprocess(I, ICal, Opts.obsnum,
                          Opts.pix_list,
                          Opts.bank,
                          Opts.use_cal,
                          tsys=Opts.tsys,
                          path=Opts.data_path)
    
    # count the total number of spectra that will be processed and written to file
    total_spectra = count_otf_spectra(S,Opts.pix_list)

    # make a dummy spectrum to count the channels after processing steps
    LD = LineData(I,Opts.bank,S.nchan,S.bandwidth,np.zeros(S.nchan))
    L = LD.vslice(Opts.slice[0],Opts.slice[1])
    nchan_to_save = L.nchan

    # write the netCDF file

    # open Dataset.  If the file exists, we stop here!
    nc = netCDF4.Dataset(Opts.output_file_name, 'w', format='NETCDF4')
    
    # dimension of number of spectra is from total number count
    nc_dimension_nspec = nc.createDimension('nspec',total_spectra)
    
    # dimension of number of channels in spectrum is from trial reduction step
    nc_dimension_nchan = nc.createDimension('nchan',nchan_to_save)
    
    # just doing 20 characters in string
    nc_dimension_nlabel = nc.createDimension('nlabel',20)
    
    # the Observation Header
    nc_obsnum = nc.createVariable('Header.Obs.ObsNum','i4')
    nc.variables['Header.Obs.ObsNum'][0] = S.obsnum
    
    # copy the source name into netCDF header
    nc_source = nc.createVariable('Header.Obs.SourceName','c',('nlabel',))
    if len(S.source) > 19:
        nc_source[0:19] = S.source[0:19]
    else:
        nc_source[0:len(S.source)] = S.source[0:len(S.source)]

    nc_x_position = nc.createVariable('Header.Obs.XPosition','f4')
    nc_y_position = nc.createVariable('Header.Obs.YPosition','f4')
    if S.map_coord == 1:
        nc.variables['Header.Obs.XPosition'][0] = S.ifproc.source_RA/np.pi*180.0
        nc.variables['Header.Obs.YPosition'][0] = S.ifproc.source_Dec/np.pi*180.0
    else:
        nc.variables['Header.Obs.XPosition'][0] = 0.0
        nc.variables['Header.Obs.YPosition'][0] = 0.0

    # using line header information derived from spec bank
    ncl = NetCDFLineHeader(nc)
    ncl.write_line_header_variables(L) # write using the result of trial run 
          
    nc_pix = nc.createVariable('Data.Pixel','i4',('nspec',))
    nc_seq = nc.createVariable('Data.Sequence','i4',('nspec',))
    nc_x = nc.createVariable('Data.XPos','f4',('nspec',))
    nc_x.units = 'arcsec'
    nc_y = nc.createVariable('Data.YPos','f4',('nspec',))
    nc_y.units = 'arcsec'
    nc_rms = nc.createVariable('Data.RMS','f4',('nspec',))
    nc_rms.units = 'K'
    nc_data = nc.createVariable('Data.Spectra','f4',('nspec','nchan'))
    nc_data.units = 'K'

    count = 0
    for ipix in Opts.pix_list:
        i = S.find_pixel_index(ipix)
        n_spectra = len(S.roach[i].xmap[S.roach[i].ons])
        x_spectra = S.roach[i].xmap[S.roach[i].ons] # x coordinate
        y_spectra = S.roach[i].ymap[S.roach[i].ons] # y coordinate
        if I.map_coord == 0:
            gx,gy = theGrid.azel(S.elev/180.*np.pi,I.tracking_beam)
        else:
            parang = np.mean(S.roach[i].pmap[S.roach[i].ons]) # average parang
            gx,gy = theGrid.radec(S.elev/180.*np.pi,parang,I.tracking_beam)

        for j in range(n_spectra):
            # process each spectrum
            L = LineData(I,Opts.bank,S.nchan,S.bandwidth,S.roach[i].reduced_spectra[j])
            LL = L.vslice(Opts.slice[0],Opts.slice[1])
            LL.eliminate(Opts.eliminate_list)
            bbase,nbase = LL.xlist(Opts.b_regions)
            LL.baseline(bbase,nbase,baseline_order=Opts.b_order)
        
            # write the reduced line into the NetCDF file
            nc_data[count,:] = LL.yarray
            nc_rms[count] = LL.rms
            nc_pix[count] = ipix
            nc_seq[count] = j
            nc_x[count] = x_spectra[j]-gx[ipix]
            nc_y[count] = y_spectra[j]-gy[ipix]
            count = count + 1

    nc.close()        
    print('netCDF %s Done'%(Opts.output_file_name))

def main(argv):

    print(time.time(), time.clock())
    
    Opts1 = HandleProcessOptions()
    Opts2 = HandleProcessOptions()
    Opts3 = HandleProcessOptions()
    Opts4 = HandleProcessOptions()

    Opts1.parse_options(argv,'multiprocess_otf_map 1',1,False)
    Opts1.pix_list = [0,1,2,3]
    output_file_name = Opts1.output_file_name + '_1.nc'
    Opts1.output_file_name = output_file_name
    Opts1.print_options()

    Opts2.parse_options(argv,'multiprocess_otf_map 2',1,False)
    Opts2.pix_list = [4,5,6,7]
    output_file_name = Opts2.output_file_name + '_2.nc'
    Opts2.output_file_name = output_file_name
    Opts2.print_options()

    Opts3.parse_options(argv,'multiprocess_otf_map 3',1,False)
    Opts3.pix_list = [8,9,10,11]
    output_file_name = Opts3.output_file_name + '_3.nc'
    Opts3.output_file_name = output_file_name
    Opts3.print_options()

    Opts4.parse_options(argv,'multiprocess_otf_map 4',1,False)
    Opts4.pix_list = [12,13,14,15]
    output_file_name = Opts4.output_file_name + '_4.nc'
    Opts4.output_file_name = output_file_name
    Opts4.print_options()
    
    print('Reading IFProc Files')
    ifproc_file = lookup_ifproc_file(Opts1.obsnum,path=Opts1.data_path+'ifproc/')
    I = IFProcData(ifproc_file)
    calobsnum = I.calobsnum
    ifproc_cal_file = lookup_ifproc_file(calobsnum,path=Opts1.data_path+'ifproc/')
    ICal = IFProcCal(ifproc_cal_file)
    ICal.compute_tsys()

    print('Processing Roach Files in parallel')
    p1 = multiprocessing.Process(target=multiprocess_otf_map, args=(Opts1,I,ICal ))
    p1.start()
    p2 = multiprocessing.Process(target=multiprocess_otf_map, args=(Opts2,I,ICal ))
    p2.start()
    p3 = multiprocessing.Process(target=multiprocess_otf_map, args=(Opts3,I,ICal ))
    p3.start()
    p4 = multiprocessing.Process(target=multiprocess_otf_map, args=(Opts4,I,ICal ))
    p4.start()

    p1.join()
    p2.join()
    p3.join()
    p4.join()
    print('DONE')
    print(time.time(), time.clock())
    
if __name__ == '__main__':
    main(sys.argv[1:])


