"""
Module for maintaining plotting status. This is an easier way to
flip the code between interactive and batch mode.   in LMTSLR
all routines that follow this should initialize via the --plots
argument option


Command usage:
   view_xxx.py  -i M51_91112.nc --plots=M51_91112
and advanced (NYI)
   view_xxx.py  -i M51_91112.nc --plots=M51_91112,pdf,10


Usage:
   import lmtslr.viewer.Plots as Plots

       old                new
    --------------        -----------------------------------------
    pl.ion()              Plots.init() -or- Plots.init("M51_91112") 
    pl.figure(1)          pl.figure()
                          Plots.savefig()
    pl.figure(2)          pl.figure()
                          Plots.savefig()
    pl.ioff()
    pl.show()             Plots.show()

"""

import matplotlib.pyplot as pl

_plnum = 0
_plint = 1
_plsiz = (8,8)
_plnam = "Figure"
_plfmt = ".%d"
_plext = "png"

class Plots(object):
    """ Static class to maining plot characteristics

    """
    @staticmethod
    def init(method=None):
        """
        This should parse the optional commandline argument --plots=XXX
        If not given the plots are interactive, and figures are labeled
        1,2,3...
        If a string is given, this will control the batch mode :
           basename 
           extension type [optional]
           first plot number [optional]

        Plots.init("M51_91112")     batch mode with M51_91112.1.png as first plot
        Plots.init("test,pdf,12")   batch mode with test.12.pdf as first plot
        
        support a comma separated string ?

        """
        global _plnum, _plint, _plnam
        if method == None:
            print("SLR Plots in interactive mode")
        else:
            print("SLR Plots initialized with %s" % str(method))
        _plnum = 0
        _plint = 1
        if method != None:     # batch mode, it needs to be a string   @todo
            _plint = 0
            _plnam = method
            
    @staticmethod
    def debug():
        print("current plot number = %d" % _plnum)
        print("current plot name   = %s" % _plnam)
        print("current interactive = %d" % _plint)
        print("current plot size   = %s" % str(_plsiz))


    @staticmethod
    def figure():
        global _plnum
        _plnum = _plnum + 1
        print("New Figure %d" % _plnum)
        if _plint:
            pl.ion()
        pl.figure(_plnum)
    
    @staticmethod
    def show():
        global _plint
        if _plint:
            pl.ioff()
            pl.show()

    @staticmethod
    def savefig(name=None):
        global _plint, _plnum, _plnam, _plext
        if _plint == 1:
            return
        if name == None:
            fmt = _plfmt % _plnum
            plname = "%s%s.%s" % (_plnam,fmt,_plext)
        print("Plots saving in %s" % plname)
        pl.savefig(plname)



if __name__ == '__main__':

    x = [0, 0.5, 1]
    y = [0, 0.7, 1]

    interactive = False

    if interactive:
        pl.ion()
        
        pl.figure(1)
        pl.plot(x,y)
        pl.title("__main__ Test1 old style")

        pl.figure(2)
        pl.plot(y,x)
        pl.title("__main__ Test2 old style")        

        
        pl.ioff()
        pl.show()
    else:
        print("Old style interactive mode, creating Test1.png and Test2.pdf")
        pl.figure(1)
        pl.plot(x,y)
        pl.title("__main__ Test1 old style")
        pl.savefig("Test1.png")

        pl.figure(2)
        pl.plot(y,x)
        pl.title("__main__ Test2 old style")
        pl.savefig("Test2.pdf")        

    # new proposed style - a common --plots= command line option would set the arguments to init()

    Plots.init("LMTPlots")
    #Plots.init()
    
    Plots.figure()
    pl.plot(x,y)
    pl.title("__main__ Test1 new style")
    Plots.savefig()

    Plots.figure()    
    pl.plot(y,x)
    pl.title("__main__ Test2 new style")
    Plots.savefig()    
    
    Plots.show()

        
