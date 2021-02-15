#    allow a site file, pick from available makedefs.XXX example or write your own
-include makedefs

#   
ifdef CFITSIO_PATH
INCLUDES += -I $(CFITSIO_PATH)/include
LIBS += -L $(CFITSIO_PATH)/lib -Xlinker -rpath -Xlinker $(CFITSIO_PATH)
endif
LIBS += -lcfitsio

ifdef NETCDF_PATH
INCLUDES += -I $(NETCDF_PATH)/include
LIBS += -L $(NETCDF_PATH)/lib -Xlinker -rpath -Xlinker $(NETCDF_PATH)
endif
LIBS += -lnetcdf

LIBS += -lm

CC = gcc
#CFLAGS = -O3 -std=c99 -Wall
CFLAGS = -O3

# an extra library (e.g. -lefence)
ELIBS = 

# just in case
CFLAGS += `pkg-config --cflags cfitsio netcdf`
LIBS   += `pkg-config --libs   cfitsio netcdf`

all:	spec_driver_fits
	@./spec_driver_fits

clean:
	rm -f *.o *~ spec_driver_fits

install: spec_driver_fits
	cp spec_driver_fits ../bin

Cube.o : Cube.h Cube.c 
	$(CC) $(CFLAGS) -c Cube.c $(INCLUDES)

Plane.o : Plane.h Plane.c
	$(CC) $(CFLAGS) -c Plane.c $(INCLUDES)

SpecFile.o : SpecFile.h SpecFile.c Version.h
	$(CC) $(CFLAGS) -c SpecFile.c $(INCLUDES)

ConvolveFunction.o : ConvolveFunction.h ConvolveFunction.c
	$(CC) $(CFLAGS) -c ConvolveFunction.c

OTFParameters.o : OTFParameters.h OTFParameters.c
	$(CC) $(CFLAGS) -c OTFParameters.c

Stats.o : Stats.h Stats.c
	$(CC) $(CFLAGS) -c Stats.c

DEPS = spec_driver_fits.c Cube.o Plane.o SpecFile.o ConvolveFunction.o OTFParameters.o Stats.o


spec_driver_fits : $(DEPS)
	$(CC) $(CFLAGS) -o spec_driver_fits $(DEPS)  $(INCLUDES) $(LIBS)  $(ELIBS)

