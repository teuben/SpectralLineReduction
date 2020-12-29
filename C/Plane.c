/** plane.c - methods for handling data plane
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Plane.h"
#include "fitsio.h"

/*************************  PLANE METHODS ***********************************/

/** initialize_plane - allocated data for the plane
 */
void initialize_plane(Plane* P, int *n)
{
  int i;
  for(i=0;i<2;i++)
    P->n[i] = n[i];
  P->nplane = n[0]*n[1];
#if defined(MDMAXDIM)
  //P->plane = allocate_mdarray2(n[1],n[0]);  // Fortran style contiguous
  P->plane = allocate_mdarray2(n[0],n[1]);  // Fortran style contiguous
#else  
  P->plane = (float*)malloc(P->nplane*sizeof(float));
  if(P->plane == NULL)
    fprintf(stderr,"Plane: Failed to Allocate Plane\n");
  else
    for(i=0;i<P->nplane;i++)
      P->plane[i] = 0.0;
#endif  
}

/** initialize_axis - initializes the axis data 
 */
void initialize_plane_axis(Plane *P, int axis, float crval, float crpix, float cdelt, char *ctype, char *cunit)
{
  int i;

  P->crval[axis] = crval;
  P->crpix[axis] = crpix;
  P->cdelt[axis] = cdelt;
  strcpy(P->ctype[axis],ctype);
  strcpy(P->cunit[axis],cunit);

  P->caxis[axis] = (float*)malloc(P->n[axis]*sizeof(float));
  if(P->caxis[axis] == NULL)
    fprintf(stderr,"Plane: Failed to allocate axis %d\n",axis);

  for(i=0;i<P->n[axis];i++)
      P->caxis[axis][i] = (i-P->crpix[axis])*P->cdelt[axis]+P->crval[axis];
}

/** axis_index - finds the index in an axis array given a value
 */
int plane_axis_index(Plane *P, int axis, float value)
{
  float result_f;
  int result_i;
  
  result_f = (value-P->crval[axis])/P->cdelt[axis] + P->crpix[axis];
  result_i = (int)floor(result_f);
  if((result_i<0) || (result_i>=P->n[axis]))
    result_i = -1;
  return(result_i);
}

/** plane_index - finds the index of an element in the plane according to 
    given x and y positions
*/
int plane_index(Plane *P, float x, float y)
{
  int ix, iy;
  int result;

  ix = plane_axis_index(P, PLANE_X_AXIS, x);
  iy = plane_axis_index(P, PLANE_Y_AXIS, y);

  result = iy*P->n[PLANE_X_AXIS] + ix;
  return(result);
}


void write_fits_plane(Plane *P, char *filename)
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
  
  naxis = 2;
  naxes[0] = P->n[PLANE_X_AXIS];
  naxes[1] = P->n[PLANE_Y_AXIS];

  // create the buffer to reorder the cube to FITS standard
#if defined(MDMAXDIM)
  buffer = &P->plane[0][0];
#else  
  buffer = (float*)malloc(P->nplane*sizeof(float));

  ic = 0;
  i=0;
  for(j=0;j<P->n[PLANE_Y_AXIS];j++)
    for(k=0;k<P->n[PLANE_X_AXIS];k++)
      {
	ii = i + j*P->n[PLANE_X_AXIS] + (P->n[PLANE_X_AXIS]-k-1);
	buffer[ic] = P->plane[ii];
	ic++;
      }
#endif
  
  // you MUST initialize status
  status = 0;
  if((retval=fits_create_file(&fptr, filename, &status)) != 0)
    print_fits_error(status);

  if((retval=fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status)) != 0)
    print_fits_error(status);

  // scale axes to standards
  strcpy(ctype,"RA---SFL");          // nominal projection Sanson-Flamsteed
  crval = 0.0;
  cdelt = -P->cdelt[PLANE_X_AXIS] / 3600.;  // degrees - we flipped the RA axis
  crpix = P->crpix[PLANE_X_AXIS];
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
  crval = 0.0;
  cdelt = P->cdelt[PLANE_Y_AXIS] / 3600.;  // degrees 
  crpix = P->crpix[PLANE_Y_AXIS];
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


  // write the data cube
  if((retval=fits_write_img(fptr, TFLOAT, 1, P->nplane, buffer, &status)) != 0)
    print_fits_error(status);

  // close the file
  if((retval=fits_close_file(fptr, &status)) != 0)
    print_fits_error(status);

  //printf("PJT all done\n");
} 
