#!/usr/bin/env python
'''
Reads a SpecFile from a single OTF mapping observation and creates visualizations
'''

# Python Imports
import sys
import matplotlib.pyplot as pl
#from lmtslr.utils.parser import HandleViewSpecFileOptions
from lmtslr.utils.argparser import HandleViewSpecFileOptions
from lmtslr.viewer.spec_file_viewer import SpecFileViewer

def main(argv):
    
    Opts = HandleViewSpecFileOptions()
    Opts.parse_options(argv, 'view_spec_file', 1, True)

    pl.ion()
    
    SV = SpecFileViewer(Opts.input_file_name)

    if Opts.show_all_pixels:
        SV.sequoia_waterfall_plot(Opts.pix_list, Opts.rms_cut, plot_range=Opts.plot_range)
        SV.sequoia_rms_plot(Opts.pix_list, Opts.rms_cut, plot_range=[0.,Opts.plot_range[1]])
        SV.sequoia_rms_histogram(Opts.pix_list, Opts.rms_cut)
        SV.sequoia_mean_spectra_plot(Opts.pix_list, Opts.rms_cut)
    else:
        SV.pixel_waterfall_plot(Opts.show_pixel, Opts.rms_cut, plot_range=Opts.plot_range)
        SV.pixel_rms_plot(Opts.show_pixel, Opts.rms_cut, plot_range=[0.,Opts.plot_range[1]])
        SV.pixel_rms_histogram(Opts.show_pixel, Opts.rms_cut)
        SV.pixel_mean_spectrum_plot(Opts.show_pixel, Opts.rms_cut)

    SV.xy_position_plot()
    pl.ioff()
    pl.show()
    
if __name__ == '__main__':
    main(sys.argv[1:])
