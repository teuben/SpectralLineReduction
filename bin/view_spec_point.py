#!/usr/bin/env python


"""Usage: view_spec_point.py  -i INPUT [options]

-i INPUT --input INPUT        Input SpecFile (no default)
--pix_list PIX_LIST           Comma separated list of pixels [Default: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
--location X,Y                Location [Default: 0,0]
-r RAD --radius RAD           Radius around location in arcsec [Default: -1] 
-p PIXEL --show_pixel PIXEL   Show a particular pixel flags code (deprecated) [Default: -1]
-z RMS_CUT --rms_cut RMS_CUT  RMS threshold for data, negative allowed for robust MAD method  [Default: 10.0]
--plot_range PLOT_RANGE       Plotting range (deprecated) [Default: -1,4]
--plots METHOD                Plotting style, defaults to interactive.
-h --help                     show this help


VIEW_SPEC_POINT will for selected pixels and a selected location and
radius in the field plot their spectra, a different color for each
pixel. But data are just averages, no weighting as function of
distance.

For convenience it still has an option (using -p) to make it behave
like the old VIEW_SPEC_FILE but only mean_spectra_plot and rms_plot
are plotted.

"""

import sys
from docopt import docopt
import matplotlib.pyplot as pl
from lmtslr.utils.argparser import HandleViewSpecFileOptions
from lmtslr.viewer.spec_file_viewer import SpecFileViewer
from lmtslr.viewer.plots import Plots
import lmtslr.utils.convert as acv

def main(argv):
    av = docopt(__doc__,options_first=True, version='0.1')
    # debug:    show the dictionary
    # print(av)

    nc_file    = av['--input']
    pix_list   = acv.listi(av['--pix_list'],  16)
    location   = acv.listf(av['--location'],   2)
    radius     = acv.listf(av['--radius'],     1)
    show_pixel = acv.listi(av['--show_pixel'], 1)
    rms_cut    = acv.listf(av['--rms_cut'],    1)
    plot_range = acv.listf(av['--plot_range'], 2)
    plots      = av['--plots']

    Plots.init(plots)
    
    SV = SpecFileViewer(nc_file)

    if show_pixel < 0:
        # new routine
        SV.pixel_mean_spectrum_plot2(pix_list, rms_cut, location, radius)
    else:
        # classic one
        SV.pixel_mean_spectrum_plot(show_pixel, rms_cut)
        SV.pixel_rms_plot(show_pixel, rms_cut, plot_range=[0.,plot_range[1]])

    Plots.show()
    
if __name__ == '__main__':
    main(sys.argv[1:])
