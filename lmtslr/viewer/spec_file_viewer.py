import numpy as np
import matplotlib.pyplot as pl
import netCDF4
from astropy.stats import mad_std

class SpecFileViewer():
    """
    Class for viewing spectrum files.
    """
    def __init__(self, netcdf_filename):
        """
        Constructor for SpecFileViewer class.
        Args:
            netcdf_filename
        Returns:
            none
        """
        nc = netCDF4.Dataset(netcdf_filename, 'r', format='NETCDF4')
        self.obsnum = nc.variables['Header.Obs.ObsNum'][0]
        self.source_name = netCDF4.chartostring(
            nc.variables['Header.Obs.SourceName'][:])

        self.x_position =  nc.variables['Header.Obs.XPosition'][0]
        self.y_position =  nc.variables['Header.Obs.YPosition'][0]

        self.nchan = nc.variables['Header.Line.NChannels'][0]
        self.chan = nc.variables['Header.Line.ChannelNumber'][:]
        self.cdelt = nc.variables['Header.SpectrumAxis.CDELT'][0]
        self.crpix = nc.variables['Header.SpectrumAxis.CRPIX'][0]
        self.crval = nc.variables['Header.SpectrumAxis.CRVAL'][0]
        self.ctype = netCDF4.chartostring(
            nc.variables['Header.SpectrumAxis.CTYPE'][:])
        self.caxis = nc.variables['Header.SpectrumAxis.CAXIS'][:]

        self.pixel = nc.variables['Data.Pixel'][:]
        self.sequence = nc.variables['Data.Sequence'][:]
        self.xpos = nc.variables['Data.XPos'][:]
        self.ypos = nc.variables['Data.YPos'][:]
        self.rms = nc.variables['Data.RMS'][:]
        self.data = nc.variables['Data.Spectra'][:]
        
        nc.close()

    def sequoia_waterfall_plot(self, pixel_list, rms_cut, plot_range=[-1,1], 
                               figsize=8):
        """
        Makes waterfall plots of spectra from pixels in pixel_list.
        Args:
            pixel_list (list): list of pixel IDs to plot
            rms_cut (float): rms cutoff value for plot
            plot_range (list): range of plot (default is [-1,1])
            figsize (float): size of figure in inches (default is 8)
        Returns:
            none
        """
        fig1, ax1 = pl.subplots(4, 4, sharex='col', sharey='row', 
                                gridspec_kw={'hspace': 0, 'wspace': 0}, 
                                figsize=(figsize, figsize))
        fig1.text(0.02, 0.5, self.ctype, va='center', rotation='vertical')
        fig1.text(0.5, 0.1, 'Sample', ha='center')
        for the_pixel in pixel_list:
            pindex = np.where(self.pixel == the_pixel)[0]
            rindex = np.where(self.rms[pindex] < rms_cut)[0]
            ax1[np.mod(the_pixel, 4), the_pixel // 4].imshow(
                self.data[pindex[rindex]].transpose(), origin='lower', 
                extent=[0, float(len(rindex)), self.caxis[0], self.caxis[-1]],
                clim=plot_range, aspect='auto')
            ax1[np.mod(the_pixel, 4), the_pixel // 4].text(0.05 * len(rindex),
                self.caxis[0] + 0.85 * (self.caxis[-1] - self.caxis[0]), 
                '%d'%(the_pixel))

    def pixel_waterfall_plot(self, the_pixel, rms_cut, plot_range=[-1,1]):
        """
        Makes waterfall plot of spectrum from pixel the_pixel.
        Args:
            the_pixel (int): pixel ID to plot
            rms_cut (float): rms cutoff value for plot
            plot_range (list): range of plot (default is [-1,1])
            figsize (float): size of figure in inches (default is 8)
        Returns:
            none
        """
        pl.figure()
        pindex = np.where(self.pixel == the_pixel)[0]
        rindex = np.where(self.rms[pindex] < rms_cut)[0]
        pl.imshow(self.data[pindex[rindex]].transpose(), origin='lower', 
                  extent=[0, float(len(rindex)), self.caxis[0], 
                  self.caxis[-1]], clim=plot_range, aspect='auto')
        pl.title('PIXEL: %d'%(the_pixel))
        pl.ylabel(self.ctype)
        pl.xlabel('Sample')
        pl.colorbar()

    def sequoia_rms_plot(self, pixel_list, rms_cut, plot_range=[0,10], 
                         figsize=8):
        """
        Makes rms plots of spectra from pixels in pixel_list.
        Args:
            pixel_list (list): list of pixel IDs to plot
            rms_cut (float): rms cutoff value for plot
            plot_range (list): range of plot (default is [0,10])
            figsize (float): size of figure in inches (default is 8)
        Returns:
            none
        """
        fig2, ax2 = pl.subplots(4, 4, sharex='col', sharey='row', 
                                gridspec_kw={'hspace': 0, 'wspace': 0}, 
                                figsize=(figsize,figsize))
        fig2.text(0.02, 0.5, 'RMS', va='center', rotation='vertical')
        fig2.text(0.5, -0.1, 'Sample', ha='center')

        for the_pixel in pixel_list:
            pindex = np.where(self.pixel == the_pixel)[0]
            rindex = np.where(self.rms[pindex] < rms_cut)[0]
            ax2[np.mod(the_pixel,4), the_pixel // 4].plot(
                self.rms[pindex[rindex]], 'k.')
            #ax2[np.mod(the_pixel,4), the_pixel//4].text(0.05*len(rindex),plot_range[0] + 0.9*(plot_range[-1]-plot_range[0]), '%d'%(the_pixel))
            #ax2[np.mod(the_pixel,4), the_pixel//4].ylim(plot_range)

    def pixel_rms_plot(self, the_pixel, rms_cut, plot_range=[0,10]):
        """
        Makes rms plot of spectrum from pixel the_pixel.
        Args:
            the_pixel (int): pixel ID to plot
            rms_cut (float): rms cutoff value for plot
            plot_range (list): range of plot (default is [0,10])
            figsize (float): size of figure in inches (default is 8)
        Returns:
            none
        """
        pl.figure()
        pindex = np.where(self.pixel == the_pixel)[0]
        rindex = np.where(self.rms[pindex] < rms_cut)[0]
        pl.plot(self.rms[pindex[rindex]], 'k.')
        pl.ylim(plot_range)
        pl.ylabel('RMS')
        pl.xlabel('Sample')
        pl.title('PIXEL: %d'%(the_pixel))
        
    def xy_position_plot(self, all=True):
        """
        Makes x-y position plot.
        Args:
            none
        Returns:
            none
        """
        s0 = 0
        s1 = max(self.sequence)+1
        print("Max sequence=%d" % s1)
        pl.figure()
        if all:
            pl.plot(self.xpos,self.ypos, 'k.')
        else:
            pl.plot(self.xpos[s0:s1],self.ypos[s0:s1], 'k.')
        pl.xlabel('X')
        pl.ylabel('Y')

    def sx_position_plot(self, all=True):
        """
        Makes sequence-x plot.
        Args:
            none
        Returns:
            none
        """
        s0 = 0
        s1 = max(self.sequence)+1
        pl.figure()
        if all:
            pl.plot(self.sequence,self.xpos, 'k.')
        else:
            pl.plot(self.sequence[s0:s1],self.xpos[s0:s1], 'k.')
        pl.xlabel('Sequence')
        pl.ylabel('X')

    def sy_position_plot(self, all=True):
        """
        Makes sequence-y position plot.
        Args:
            none
        Returns:
            none
        """
        s0 = 0
        s1 = max(self.sequence)+1
        pl.figure()
        if all:
            pl.plot(self.sequence,self.ypos, 'k.')
        else:
            pl.plot(self.sequence[s0:s1],self.ypos[s0:s1], 'k.')
        pl.xlabel('Sequence')
        pl.ylabel('Y')
        
    def sequoia_rms_histogram(self, pixel_list, rms_cut, figsize=8):
        """
        Makes rms histogram of spectra from pixels in pixel_list.
        Args:
            pixel_list (list): list of pixel IDs to plot
            rms_cut (float): rms cutoff value for plot
            figsize (float): size of figure in inches (default is 8)
        Returns:
            none
        """
        fig3, ax3 = pl.subplots(4, 4, sharex='col', sharey='row', 
            gridspec_kw={'hspace': 0, 'wspace': 0}, figsize=(figsize,figsize))
        fig3.text(0.5, -0.1, 'RMS', ha='center')
        for the_pixel in pixel_list:
            pindex = np.where(self.pixel == the_pixel)[0]
            rindex = np.where(self.rms[pindex] < rms_cut)[0]
            ax3[np.mod(the_pixel, 4), the_pixel // 4].hist(
                self.rms[pindex[rindex]], bins=np.arange(0,3.02,.02))

    def pixel_rms_histogram(self, the_pixel, rms_cut):
        """
        Makes rms histogram of spectra from pixel the_pixel.
        Args:
            the_pixel (int): pixel ID to plot
            rms_cut (float): rms cutoff value for plot
            figsize (float): size of figure in inches (default is 8)
        Returns:
            none
        """
        pl.figure()
        pindex = np.where(self.pixel == the_pixel)[0]
        rindex = np.where(self.rms[pindex] < rms_cut)[0]
        pl.hist(self.rms[pindex[rindex]],bins = np.arange(0,3.02,.02))
        pl.xlabel('RMS')
        pl.ylabel('N')
        pl.title('PIXEL: %d'%(the_pixel))
        
    def sequoia_mean_spectra_plot(self, pixel_list, rms_cut, figsize=8):
        """
        Makes mean spectra plot of spectra from pixels in pixel_list.
        Args:
            pixel_list (list): list of pixel IDs to plot
            rms_cut (float): rms cutoff value for plot
            figsize (float): size of figure in inches (default is 8)
        Returns:
            none
        """
        fig4, ax4 = pl.subplots(4, 4, sharex='col', sharey='row', 
            gridspec_kw={'hspace': 0, 'wspace': 0}, figsize=(figsize,figsize))
        fig4.text(0.5, -0.1, self.ctype, ha='center')
        for the_pixel in pixel_list:
            pindex = np.where(self.pixel == the_pixel)[0]
            rindex = np.where(self.rms[pindex] < rms_cut)[0]
            ax4[np.mod(the_pixel, 4), the_pixel // 4].plot(self.caxis, 
                np.mean(self.data[pindex[rindex]], axis=0))
 

    def pixel_mean_spectrum_plot(self, the_pixel, rms_cut):
        """
        Makes mean spectra plot of spectra from pixel the_pixel.
        Args:
            the_pixel (int): pixel ID to plot
            rms_cut (float): rms cutoff value for plot
            figsize (float): size of figure in inches (default is 8)
        Returns:
            none
        """
        pl.figure()
        pindex = np.where(self.pixel == the_pixel)[0]
        rindex = np.where(self.rms[pindex] < rms_cut)[0]
        pl.plot(self.caxis, np.mean(self.data[pindex[rindex]], axis=0))
        pl.xlabel(self.ctype)
        pl.ylabel('TA*')
        pl.title('PIXEL: %d'%(the_pixel))
                
    def sequoia_mean_spectra_plot2(self, pixel_list, rms_cut, figsize=8):
        """
        Makes mean spectra plot of spectra from pixels in pixel_list.
        Args:
            pixel_list (list): list of pixel IDs to plot
            rms_cut (float): rms cutoff value for plot
            figsize (float): size of figure in inches (default is 8)
        Returns:
            none
        """
        fig4, ax4 = pl.subplots(4, 4, sharex='col', sharey='row', 
            gridspec_kw={'hspace': 0, 'wspace': 0}, figsize=(figsize,figsize))
        fig4.text(0.5, -0.1, self.ctype, ha='center')
        for the_pixel in pixel_list:
            pindex = np.where(self.pixel == the_pixel)[0]
            rindex = np.where(self.rms[pindex] < rms_cut)[0]
            ax4[np.mod(the_pixel, 4), the_pixel // 4].plot(self.caxis, 
                np.mean(self.data[pindex[rindex]], axis=0))
 
    def pixel_mean_spectrum_plot2(self, pixel_list, rms_cut, location, radius):
        """
        Overplots mean spectra plot of spectra from a pixel_list
        Args:
            pixel_list (list of int): pixel IDs to plot
            rms_cut (float): rms cutoff value for plot (negative allowed)
            location (list of 2 floats):   location in grid
            radius (float); radius of circle within which points selected
        Returns:
            none
        """

        if radius > 0:
            dx = self.xpos - location[0]
            dy = self.ypos - location[1]
            r2 = dx*dx+dy*dy
            rad2 = radius*radius

        pl.figure()

        for the_pixel in pixel_list:
            pindex = np.where(self.pixel == the_pixel)[0]

            print("Warning: MAD test",pindex)
            med1 = np.median(self.rms[pindex])
            std1 = mad_std(self.rms[pindex])
            if rms_cut  < 0:
                cut1 = med1 - rms_cut*std1
            else:
                cut1 = rms_cut
            print("MAD:",med1,std1,cut1)
        
            rindex = np.where(self.rms[pindex] < cut1)[0]
            if radius > 0:
                cindex = np.where(r2[pindex[rindex]] < rad2)[0]
                print('Pixel %d %d' % (the_pixel, len(cindex)))
                pl.plot(self.caxis, np.mean(self.data[pindex[rindex[cindex]]], axis=0), label="%d" % the_pixel)
            else:
                print('Pixel %d %d' % (the_pixel, len(rindex)))
                pl.plot(self.caxis, np.mean(self.data[pindex[rindex]], axis=0), label="%d" % the_pixel)
        pl.xlabel(self.ctype)
        pl.ylabel('TA*')
        pl.title('PIXEL: %s   rms_cut: %g'%(str(pixel_list),cut1))
        pl.legend()
                
