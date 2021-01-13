#ifndef _SPECFILE_H_
#define _SPECFILE_H_

#include <netcdf.h>

/** error handling for netcdf reading and writing
 */

#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

typedef struct
{
  int nspec,nchan;
  // data from obs header
  int obsnum;
  char source[20];
  double x_position, y_position;
  double restfreq;
  float vlsr;
  float zsource;
  char date_obs[20];
  // axis parameters from line header
  double CRPIX, CRVAL, CDELT;
  double *CAXIS;
  char CTYPE[20];
  // the data from the file
  float *theData;
  int *Pixel;
  float *RMS_cut;  
  int *Sequence;
  int *use;
  float *XPos;
  float *YPos;
  float *RMS;
  double *Date;        // not used yet
  char history[512];   // cmdline how the specfile was made ( > 10-jan-2021)
} SpecFile;
  

int read_spec_file(SpecFile *S, char *filename);
void make_spec_beam(SpecFile *S);
void free_spec_file(SpecFile *S);
float *get_spectrum(SpecFile *S, int i);

#endif
