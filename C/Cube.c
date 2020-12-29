/** cube.c - methods for handling data cube
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <netcdf.h>
#include "Cube.h"
#include "SpecFile.h"
#include "Version.h"
#include "fitsio.h"

/*************************  CUBE METHODS ***********************************/

/** initialize_cube - allocated data for the cube
 */
void initialize_cube(Cube* C, int *n)
{
  int i;
  for(i=0;i<3;i++)
    C->n[i] = n[i];
  C->ncube  = n[0]*n[1]*n[2];
  C->nplane = n[0]*n[1];
#if defined(MDMAXDIM)
  C->cube = allocate_mdarray3(n[2],n[0],n[1]);  // Fortran style contiguous
#else  
  C->cube = (float*)malloc(C->ncube*sizeof(float));
  if(C->cube == NULL)
    fprintf(stderr,"Cube: Failed to Allocate Cube\n");
  else
    for(i=0;i<C->ncube;i++)
      C->cube[i] = 0.0;
#endif  
}

/** initialize_axis - initializes the axis data 
 */
void initialize_cube_axis(Cube *C, int axis, float crval, float crpix, float cdelt, char *ctype, char *units)
{
  int i, nn;

  C->crval[axis] = crval;
  C->crpix[axis] = crpix;
  C->cdelt[axis] = cdelt;
  for(i=0;i<16;i++)
    {
      C->ctype[axis][i] = '\0';
      C->units[axis][i] = '\0';
    }
  strcpy(C->ctype[axis],ctype);
  strcpy(C->units[axis],units);

  C->caxis[axis] = (float*)malloc(C->n[axis]*sizeof(float));
  if(C->caxis[axis] == NULL)
    fprintf(stderr,"Cube: Failed to allocate axis %d\n",axis);

  for(i=0;i<C->n[axis];i++)
      C->caxis[axis][i] = (i-C->crpix[axis])*C->cdelt[axis]+C->crval[axis];
}

/** axis_index - finds the index in an axis array given a value
 */
int cube_axis_index(Cube *C, int axis, float value)
{
  float result_f;
  int result_i;
  
  result_f = (value-(C->crval[axis]-C->cdelt[axis]/2.))/C->cdelt[axis] + C->crpix[axis];
  result_i = (int)floor(result_f);
  if((result_i<0) || (result_i>=C->n[axis]))
    result_i = -1;
  return(result_i);
}

/** cube_z_index - finds the index of the first element of z array in cube
    given x and y positions
*/
int cube_z_index(Cube *C, float x, float y)
{
  int ix, iy;
  int result;

  ix = cube_axis_index(C, X_AXIS, x);
  iy = cube_axis_index(C, Y_AXIS, y);

  result = iy*C->n[X_AXIS]*C->n[Z_AXIS] + ix*C->n[Z_AXIS];
  return(result);
}

/** write netcdf data cube (actually not used, we write fits cubes)
 */
#if !defined(MDMAXDIM)
void write_netcdf_cube(Cube *C, char *filename)
{
  int ncid;
  int retval;
  int dimid_x, dimid_y, dimid_z, dimid_label, dimid_T;
  int T_dimids[3], l_dimids[3];
  int varid_obsnum, varid_source;
  int varid_NAXIS1, varid_NAXIS2, varid_NAXIS3;
  int varid_CRVAL1, varid_CRPIX1, varid_CDELT1, varid_CTYPE1, varid_CAXIS1; 
  int varid_CRVAL2, varid_CRPIX2, varid_CDELT2, varid_CTYPE2, varid_CAXIS2; 
  int varid_CRVAL3, varid_CRPIX3, varid_CDELT3, varid_CTYPE3, varid_CAXIS3; 

  if ((retval = nc_create(filename, NC_NOCLOBBER, &ncid)) != NC_NOERR)
    ERR(retval);

  if((retval = nc_def_dim(ncid, "x", (size_t)C->n[X_AXIS], &dimid_x)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_dim(ncid, "y", (size_t)C->n[Y_AXIS], &dimid_y)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_dim(ncid, "z", (size_t)C->n[Z_AXIS], &dimid_z)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_dim(ncid, "label", (size_t)16, &dimid_label)) != NC_NOERR)
    ERR(retval);
  l_dimids[0] = dimid_label;
  l_dimids[1] = dimid_label;
  l_dimids[2] = dimid_label;


  if((retval = nc_def_var(ncid, "Header.Obs.ObsNum", NC_INT, 0, &dimid_y, &varid_obsnum)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "Header.Obs.SourceName", NC_CHAR, 1, &dimid_label, &varid_source)) != NC_NOERR)
    ERR(retval);

  if((retval = nc_def_var(ncid, "NAXIS1", NC_INT, 0, &dimid_z, &varid_NAXIS1)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CRVAL1", NC_FLOAT, 0, &dimid_z, &varid_CRVAL1)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CRPIX1", NC_FLOAT, 0, &dimid_z, &varid_CRPIX1)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CDELT1", NC_FLOAT, 0, &dimid_z, &varid_CDELT1)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CTYPE1", NC_CHAR, 1, &dimid_label, &varid_CTYPE1)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CAXIS1", NC_FLOAT, 1, &dimid_z, &varid_CAXIS1)) != NC_NOERR)
    ERR(retval);

  if((retval = nc_def_var(ncid, "NAXIS2", NC_INT, 0, &dimid_x, &varid_NAXIS2)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CRVAL2", NC_FLOAT, 0, &dimid_x, &varid_CRVAL2)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CRPIX2", NC_FLOAT, 0, &dimid_x, &varid_CRPIX2)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CDELT2", NC_FLOAT, 0, &dimid_x, &varid_CDELT2)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CTYPE2", NC_CHAR, 1, &dimid_label, &varid_CTYPE2)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CAXIS2", NC_FLOAT, 1, &dimid_x, &varid_CAXIS2)) != NC_NOERR)
    ERR(retval);

  if((retval = nc_def_var(ncid, "NAXIS3", NC_INT, 0, &dimid_y, &varid_NAXIS3)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CRVAL3", NC_FLOAT, 0, &dimid_y, &varid_CRVAL3)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CRPIX3", NC_FLOAT, 0, &dimid_y, &varid_CRPIX3)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CDELT3", NC_FLOAT, 0, &dimid_y, &varid_CDELT3)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CTYPE3", NC_CHAR, 1, &dimid_label, &varid_CTYPE3)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CAXIS3", NC_FLOAT, 1, &dimid_y, &varid_CAXIS3)) != NC_NOERR)
    ERR(retval);


  T_dimids[0] = dimid_y;
  T_dimids[1] = dimid_x;
  T_dimids[2] = dimid_z;
  if((retval = nc_def_var(ncid, "T", NC_FLOAT, 3, T_dimids, &dimid_T)) != NC_NOERR)
    ERR(retval);

  if((retval = nc_enddef(ncid)) != NC_NOERR)
    ERR(retval);

  if((retval = nc_put_var_int(ncid, varid_obsnum, &C->obsnum)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var(ncid, varid_source, C->source)) != NC_NOERR)
    ERR(retval);



  if((retval = nc_put_var_int(ncid, varid_NAXIS1, &C->n[Z_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CRVAL1, &C->crval[Z_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CRPIX1, &C->crpix[Z_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CDELT1, &C->cdelt[Z_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var(ncid, varid_CTYPE1, C->ctype[Z_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CAXIS1, C->caxis[Z_AXIS])) != NC_NOERR)
    ERR(retval);

  if((retval = nc_put_var_int(ncid, varid_NAXIS2, &C->n[X_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CRVAL2, &C->crval[X_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CRPIX2, &C->crpix[X_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CDELT2, &C->cdelt[X_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var(ncid, varid_CTYPE2, C->ctype[X_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CAXIS2, C->caxis[X_AXIS])) != NC_NOERR)
    ERR(retval);

  if((retval = nc_put_var_int(ncid, varid_NAXIS3, &C->n[Y_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CRVAL3, &C->crval[Y_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CRPIX3, &C->crpix[Y_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CDELT3, &C->cdelt[Y_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var(ncid, varid_CTYPE3, C->ctype[Y_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CAXIS3, C->caxis[Y_AXIS])) != NC_NOERR)
    ERR(retval);

  if((retval = nc_put_var_float(ncid, dimid_T, C->cube)) != NC_NOERR)
    ERR(retval);
	   
  nc_close(ncid);  
}
#endif

void write_fits_cube(Cube *C, char *filename)
{
  int i,j,k,ii,ic;
  int retval, status;
  int naxis;
  long naxes[3], obsnum;
  float equinox;
  char radesys[20];
  float *buffer;
  fitsfile *fptr;

  char ctype[20], cunit[20], comment[60];
  float crval, cdelt, crpix, bmaj, bpa;
  
  naxis = 3;
  naxes[0] = C->n[X_AXIS];
  naxes[1] = C->n[Y_AXIS];
  naxes[2] = C->n[Z_AXIS];

  // create the buffer to reorder the cube to FITS standard
  // buffer(ix,iy,iz) in fortran, or buffer[iz][iy][ix]
#if defined(MDMAXDIM)
  // no reordering needed if we picked the Fortran ordering
  buffer = &C->cube[0][0][0];
#else  
  buffer = (float*)malloc(C->ncube*sizeof(float));

  ic = 0;
  for(i=0;i<C->n[Z_AXIS];i++)
      for(j=0;j<C->n[Y_AXIS];j++)
	  for(k=0;k<C->n[X_AXIS];k++)
	    {
	      ii = i + j*C->n[Y_AXIS]*C->n[Z_AXIS] + (C->n[X_AXIS]-k-1)*C->n[Z_AXIS];
	      buffer[ic] = C->cube[ii];
	      ic++;
	    }
#endif
  // you MUST initialize status
  status = 0;
  if((retval=fits_create_file(&fptr, filename, &status)) != 0)
    print_fits_error(status);

  if((retval=fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status)) != 0)
    print_fits_error(status);

  // Write the header variables
  if((retval=fits_update_key(fptr, TSTRING, "TELESCOP", "LMT", " ", &status)) != 0)
    {
      printf("TELESCOP\n");
      print_fits_error(status);
    }
  // INSTRUME is not in https://heasarc.gsfc.nasa.gov/docs/fcg/common_dict.html
  strcpy(comment,"or WARES ?");  
  if((retval=fits_update_key(fptr, TSTRING, "INSTRUME", "SEQUOIA", comment, &status)) != 0)
    {
      printf("INSTRUME\n");
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TSTRING, "OBJECT  ", C->source, " ", &status)) != 0)
    {
      printf("OBJECT\n");
      print_fits_error(status);
    }
  strcpy(comment,"Data.TelescopeBackend.TelTime");
  if((retval=fits_update_key(fptr, TSTRING, "DATE-OBS", C->date_obs, comment, &status)) != 0)
    {
      printf("DATE-OBS\n");
      print_fits_error(status);
    }
  strcpy(comment,"Software version");
  if((retval=fits_update_key(fptr, TSTRING, "ORIGIN  ", LMTSLR_VERSION, comment, &status)) != 0)
    {
      printf("ORIGIN\n");
      print_fits_error(status);
    }
  strcpy(comment,"LMT observing number");  
  obsnum = (long)C->obsnum;
  if((retval=fits_update_key(fptr, TLONG,   "OBSNUM  ", &obsnum, comment, &status)) != 0)
    {
      printf("OBSNUM\n");
      print_fits_error(status);
    }

  // BUNIT
  strcpy(cunit,"K       ");
  if((retval=fits_update_key(fptr, TSTRING, "BUNIT   ", cunit, " ", &status)) != 0)
    {
      printf("BUNIT %s\n",cunit);
      print_fits_error(status);
    }

  // scale axes to standards
  strcpy(ctype,"RA---SFL");          // nominal projection Sanson-Flamsteed
  crval = C->x_position;             // degrees
  cdelt = -C->cdelt[X_AXIS] / 3600.;  // degrees - we flipped the RA axis
  crpix = C->crpix[X_AXIS];
  strcpy(cunit,"deg     ");
  if((retval=fits_update_key(fptr, TSTRING, "CTYPE1  ", ctype, " ", &status)) != 0)
    {
      printf("CTYPE1 %s\n",ctype);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TFLOAT,  "CRVAL1  ", &crval, cunit, &status)) != 0)
    {
      printf("CRVAL1 %f\n",crval);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TFLOAT,  "CDELT1  ", &cdelt, cunit, &status)) != 0)
    {
      printf("CDELT1 %f\n",cdelt);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TFLOAT,  "CRPIX1  ", &crpix, " ", &status)) != 0)
    {
      printf("CRPIX1 %f\n",crpix);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TSTRING, "CUNIT1  ", cunit, " ", &status)) != 0)
    {
      printf("CUNIT1 %s\n",cunit);
      print_fits_error(status);
    }

  strcpy(ctype,"DEC--SFL");          // nominal projection Sanson-Flamsteed
  crval = C->y_position;             // degrees
  cdelt = C->cdelt[Y_AXIS] / 3600.;  // degrees 
  crpix = C->crpix[Y_AXIS];
  strcpy(cunit,"deg     ");
  if((retval=fits_update_key(fptr, TSTRING, "CTYPE2  ", ctype, " ", &status)) != 0)
    {
      printf("CTYPE2 %s\n",ctype);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TFLOAT,  "CRVAL2  ", &crval, cunit, &status)) != 0)
    {
      printf("CRVAL2 %f\n",crval);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TFLOAT,  "CDELT2  ", &cdelt, cunit, &status)) != 0)
    {
      printf("CDELT2 %f\n",cdelt);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TFLOAT,  "CRPIX2  ", &crpix, " ", &status)) != 0)
    {
      printf("CRPIX2 %f\n",crpix);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TSTRING, "CUNIT2  ", cunit, " ", &status)) != 0)
    {
      printf("CUNIT2 %s\n",cunit);
      print_fits_error(status);
    }

  strcpy(ctype,"VELO-LSR");          // nominal projection  
  crval = C->crval[Z_AXIS]*1000.;    // m/s
  cdelt = C->cdelt[Z_AXIS]*1000.;    // m/s
  crpix = C->crpix[Z_AXIS];
  strcpy(cunit,"m/s     ");
  
  if((retval=fits_update_key(fptr, TSTRING, "CTYPE3  ", ctype, " ", &status)) != 0)
    {
      printf("CTYPE3 %s\n",ctype);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TFLOAT,  "CRVAL3  ", &crval, cunit, &status)) != 0)
    {
      printf("CRVAL3 %f\n",crval);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TFLOAT,  "CDELT3  ", &cdelt, cunit, &status)) != 0)
    {
      printf("CDELT3 %f\n",cdelt);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TFLOAT,  "CRPIX3  ", &crpix, " ", &status)) != 0)
    {
      printf("CRPIX3 %f\n",crpix);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TSTRING, "CUNIT3  ", cunit, " ", &status)) != 0)
    {
      printf("CUNIT3 %s\n",cunit);
      print_fits_error(status);
    }

  strcpy(cunit,"deg     ");
  bmaj = C->resolution_size / 3600.;
  if((retval=fits_update_key(fptr, TFLOAT,  "BMAJ    ", &bmaj, cunit, &status)) != 0)
    {
      printf("BMAJ %f\n",bmaj);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TFLOAT,  "BMIN    ", &bmaj, cunit, &status)) != 0)
    {
      printf("BMIN %f\n",bmaj);
      print_fits_error(status);
    }
  bpa = 0.0;
  if((retval=fits_update_key(fptr, TFLOAT,  "BPA     ", &bpa, cunit, &status)) != 0)
    {
      printf("BPA %f\n",bmaj);
      print_fits_error(status);
    }
      
  strcpy(radesys,"FK5     ");
  equinox = 2000.;
  if((retval=fits_update_key(fptr, TFLOAT,  "EQUINOX ", &equinox, " ", &status)) != 0)
    {
      printf("EQUINOX %f\n", equinox);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TSTRING, "RADESYS ", radesys, " ", &status)) != 0)
    {
      printf("RADESYS %s\n", radesys);
      print_fits_error(status);
    }

  strcpy(comment, "Header.LineData.VSource");
  if((retval=fits_update_key(fptr, TFLOAT,  "VLSR", &C->vlsr, comment, &status)) != 0)
    {
      printf("VLSR %f\n", C->vlsr);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TFLOAT,  "VSOURCE", &C->vlsr, comment, &status)) != 0)
    {
      printf("VSOURCE %f\n", C->vlsr);
      print_fits_error(status);
    }
  // alma: ZSOURCE  (even astropy claims this is the preferred one)
  // alma: LINTRN

  strcpy(comment, "Header.LineData.LineRestFrequency");
  if((retval=fits_update_key(fptr, TDOUBLE,  "RESTFRQ", &C->restfreq, comment, &status)) != 0)
    {
      printf("RESTFRQ %f\n", C->restfreq);
      print_fits_error(status);
    }

  // @todo this is assumed here, but should be fixed if upstream not VELO-LSR was choosen
  strcpy(cunit,"LSRK");
  strcpy(comment,"could be wrong");
  if((retval=fits_update_key(fptr, TSTRING,  "SPECSYS", cunit, comment, &status)) != 0)
    {
      printf("SPECSYS %s\n", cunit);
      print_fits_error(status);
    }

  // write the data cube
  if((retval=fits_write_img(fptr, TFLOAT, 1, C->ncube, buffer, &status)) != 0)
    print_fits_error(status);

  // close the file
  if((retval=fits_close_file(fptr, &status)) != 0)
    print_fits_error(status);

  //printf("PJT all done\n");
} 

void print_fits_error(int status)
{
  char e[81];
  if(status)
    {
      fits_get_errstatus(status, e);
      fprintf(stderr, "%s\n",e);
      exit(status);
    }
  return;
}
