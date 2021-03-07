#!/usr/bin/env python
#

"""Usage: view_spec_file.py -i INPUT [options]

-i INPUT --input INPUT                 Input SpecFile filename (default: None)
--pix_list PIX_LIST                    list of pixels [Default: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
-p SHOW_PIXEL --show_pixel SHOW_PIXEL  Show one particular pixel (default is all pixels from pix_list)
--rms_cut RMS_CUT                      rms threshold for data [Default: 10]
--plot_range PLOT_RANGE                set plot range for plots [Default: -1,3]
--plot_file PLOT,FILE                  Basename of output plot files, defaults to on screen.
                                       Optional file type extension after a comma, e.g. png or pdf
--skip_tsys                            Skip the tsys plot [Default: False]

-h --help                              show this help


Reads a SpecFile from a single OTF mapping observation and creates
various visualizations

Tsys is only present for SpecFiles created after 6-mar-2021, data produce before this needs --skip_tsys.

Bug:  tsys pixel number can be wrong is not all pixels are in the SpecFile

"""

# Python Imports
import sys
import matplotlib.pyplot as pl

# command line parsing
from docopt import docopt
import lmtslr.utils.convert as acv

#from lmtslr.utils.parser import HandleViewSpecFileOptions
from lmtslr.utils.argparser import HandleViewSpecFileOptions
from lmtslr.viewer.spec_file_viewer import SpecFileViewer
from lmtslr.viewer.plots import Plots

def main(argv):
    av = docopt(__doc__,options_first=True, version='0.1')
    print(av)   # debug

    input_file_name = av['--input']
    pix_list        = acv.listi(av['--pix_list'],  16)
    rms_cut         = acv.listf(av['--rms_cut'],    1)
    plot_range      = acv.listf(av['--plot_range'], 2)
    plot_file       = av['--plot_file']
    skip_tsys       = av['--skip_tsys']
    
    if av['--show_pixel'] != None:
        show_all_pixels = False        
        show_pixel  = acv.listi(av['--show_pixel'], 1)
    else:
        show_all_pixels = True
    
    Plots.init(plot_file)
    
    SV = SpecFileViewer(input_file_name)

    if show_all_pixels:
        SV.sequoia_waterfall_plot(pix_list, rms_cut, plot_range=plot_range)
        SV.sequoia_rms_plot(pix_list, rms_cut, plot_range=[0.,plot_range[1]])
        SV.sequoia_rms_histogram(pix_list, rms_cut)
        SV.sequoia_mean_spectra_plot(pix_list, rms_cut)
        if not skip_tsys: SV.sequoia_tsys_spectra_plot(pix_list)
        SV.xy_position_plot()
        #SV.sx_position_plot()
        #SV.sy_position_plot()        
    else:
        SV.pixel_waterfall_plot(show_pixel, rms_cut, plot_range=plot_range)
        SV.pixel_rms_plot(show_pixel, rms_cut, plot_range=[0.,plot_range[1]])
        SV.pixel_rms_histogram(show_pixel, rms_cut)
        SV.pixel_mean_spectrum_plot(show_pixel, rms_cut)
        if not skip_tsys: SV.pixel_tsys_spectra_plot(show_pixel)
        SV.xy_position_plot(False)
        #SV.sx_position_plot(False)
        #SV.sy_position_plot(False)
    
    Plots.show()
    
if __name__ == '__main__':
    main(sys.argv[1:])
