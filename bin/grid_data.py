#!/usr/bin/env python
"""
Reads a SpecFile and creates a data cube in FITS format
"""

# Python Imports	
import numpy as np		
import matplotlib.pyplot as pl
import subprocess		 
import netCDF4			 

# Line Data Reduction Imports
from lmtslr.spec.spec import *
from lmtslr.ifproc.ifproc import *

from lmtslr.reduction.line_reduction import *

#from lmtslr.utils.parser import HandleGridOptions
from lmtslr.utils.argparser import HandleGridOptions


def main(argv):

    Opts = HandleGridOptions()
    Opts.parse_options(argv,'grid_data',1,True)
    
    # check to see whether output file exists and remove it if it does
    if os.path.isfile(Opts.output_file_name) == True:
        os.remove(Opts.output_file_name) 

    print(Opts.pix_list)        # default:  [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    print(Opts.program_path)    # default:  'spec_driver_fits' 
    with open('out.txt','w+') as outputfile:
        with open('err.txt','w+') as errorfile:

            #exit_code = subprocess.call(['./test_otf','-i',Opts.input_file_name,'-u',Opts.pix_list],stdout=outputfile,stderr=errorfile)
            exit_code=subprocess.call([Opts.program_path,
                                       '-i',Opts.input_file_name,     # 15 options are passed on
                                       '-o',Opts.output_file_name,
                                       '-l',str(Opts.resolution),     # --resolution
                                       '-c',str(Opts.cell),           # --cell
                                       '-u',str(Opts.pix_list),       # --pix_list 
                                       '-z',str(Opts.rms_cut),        # --rms_cut
                                       '-s',str(Opts.noise_sigma),    # --noise_sigma
                                       '-x',str(Opts.x_extent),       # --x_extent
                                       '-y',str(Opts.y_extent),       # --y_extent
                                       '-f',str(Opts.otf_select),     # --otf_select
                                       '-r',str(Opts.rmax),           # --rmax 
                                       '-n',str(Opts.n_samples),      # --n_samples
                                       '-0',str(Opts.otf_a),          # --otf_a
                                       '-1',str(Opts.otf_b),          # --otf_b
                                       '-2',str(Opts.otf_c),],        # --otf_c
                                      stdout=outputfile,
                                      stderr=errorfile)
            
            # reset stdout file to read from it
            outputfile.seek(0)
            # save output (if any) in variable
            standard_output=outputfile.read()
            print('STDOUT ****************************')
            print(standard_output)

            # reset stderr file to read from it
            errorfile.seek(0) 
            # save errors (if any) in variable
            standard_error = errorfile.read()
            print('STDERR ****************************')
            print(standard_error)

    print('Exit Code: %d'%(exit_code))


if __name__ == '__main__':
    main(sys.argv[1:])

