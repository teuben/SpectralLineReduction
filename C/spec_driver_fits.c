#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "Cube.h"
#include "Plane.h"
#include "ConvolveFunction.h"
#include "OTFParameters.h"
#include "SpecFile.h"
#include "Version.h"

int main(int argc, char *argv[])
{
  Cube C;
  Plane Weight, Mask;
  SpecFile S;
  ConvolveFunction CF;
  OTFParameters OTF;

  float *spectrum;
  int ifile;
  int i,j,k;
  int ii,jj;
  int ix,iy,iz,ixp,iyp,izp;
  int ngood=0;
  float x,y,X,Y,distance,weight,rmsweight;
  int n[3];

  printf("%s %s\n", argv[0], LMTSLR_VERSION);

  //printf("1\n");
  // initialize
  initialize_otf_parameters(&OTF, argc, argv);
  //printf("2\n");

  // read the first SpecFile 
  read_spec_file(&S, OTF.i_filename[0]);
  //printf("3\n");
  // copy over obs header variables
  C.obsnum = S.obsnum;
  printf("%d\n",C.obsnum);
  //printf("%s\n",S.source);
  strncpy(C.source,S.source,18);  // 18, seriously?   - there's 20, 32 and now 18?
  printf("%s\n",C.source);
  strncpy(C.date_obs,S.date_obs,20); 
  printf("DATE-OBS %s\n",C.date_obs);  
  C.x_position = S.x_position;
  C.y_position = S.y_position;
  // 
  C.resolution_size = OTF.resolution_size;
  C.restfreq = S.restfreq;
  C.vlsr = S.vlsr;
  // set up convolution array for the gridding process.
  initialize_convolve_function(&CF, OTF.resolution_size, OTF.cell_size, OTF.rmax, OTF.nsamples);
  //printf("4\n");
  if(OTF.otf_select == 1)
    initialize_jinc_filter(&CF, OTF.otf_jinc_a, OTF.otf_jinc_b, OTF.otf_jinc_c);
  else if(OTF.otf_select == 2)
    initialize_gauss_filter(&CF, OTF.otf_jinc_b);
  else if(OTF.otf_select == 3)
    initialize_triangle_filter(&CF, OTF.resolution_size);
  else
    initialize_box_filter(&CF, OTF.cell_size/2.);

  // prints the convolution function ; n_cells denotes how much we will use?
  // @todo the scaling of delta is wrong, but irrelevant
  printf("CF.n_cells= %d cell=%g  oft_select=%d\n",CF.n_cells,CF.delta,OTF.otf_select);
#if 0
  printf("r(arcsec)  conv.array\n");
  for(i=0;i< CF.npts;i++)
    printf("%5.2f %8.4f\n",i*CF.delta, CF.array[i]);
#endif

  // initialize cube and axes
  n[0] = 2 * (int)(floor((OTF.x_extent+OTF.cell_size/2.)/OTF.cell_size)) + 1;
  n[1] = 2 * (int)(floor((OTF.y_extent+OTF.cell_size/2.)/OTF.cell_size)) + 1;
  n[2] = S.nchan;

  //printf("5\n");
  initialize_cube(&C, n);
  // note that we add one to crpix's per fits convention
  initialize_cube_axis(&C, Z_AXIS, S.CRVAL, S.CRPIX+1.,      S.CDELT,       S.CTYPE, "km/s");
  initialize_cube_axis(&C, X_AXIS, 0.0,     (n[0]-1.)/2.+1., OTF.cell_size, "X",     "arcsec");
  initialize_cube_axis(&C, Y_AXIS, 0.0,     (n[1]-1.)/2.+1., OTF.cell_size, "Y",     "arcsec");

  initialize_plane(&Weight, n);
  initialize_plane_axis(&Weight, PLANE_X_AXIS, 0.0, (n[0]-1.)/2.+1., OTF.cell_size, "X", "arcsec");
  initialize_plane_axis(&Weight, PLANE_Y_AXIS, 0.0, (n[1]-1.)/2.+1., OTF.cell_size, "Y", "arcsec");
  //printf("6\n");

  initialize_plane(&Mask, n);
  initialize_plane_axis(&Mask, PLANE_X_AXIS, 0.0, (n[0]-1.)/2.+1., OTF.cell_size, "X", "arcsec");
  initialize_plane_axis(&Mask, PLANE_Y_AXIS, 0.0, (n[1]-1.)/2.+1., OTF.cell_size, "Y", "arcsec");
  

  free_spec_file(&S);
  printf("axes initialized\n");

  for(i=0;i<16;i++)
    printf("%d %d\n",i,OTF.use_pixels[i]);

  for(ifile=0;ifile<OTF.nfiles;ifile++)
    {
      // read the new specfile for gridding
      printf("file %d %s\n",ifile,OTF.i_filename[ifile]);
      read_spec_file(&S, OTF.i_filename[ifile]);
      int nout = 0;

      // now we do the gridding
      for(i=0;i<S.nspec;i++)
	{
	  //if (i != 10000) continue;  // PJT test
	  if(OTF.use_pixels[S.Pixel[i]] == 1)
	    {
	      if(S.RMS[i] < OTF.rms_cutoff)
		{
		  ngood++;
		  spectrum = get_spectrum(&S,i);
		  ix = cube_axis_index(&C, X_AXIS, S.XPos[i]);
		  iy = cube_axis_index(&C, Y_AXIS, S.YPos[i]);
		  if( (ix>=0) && (iy>=0) )
		    {
		      for(ii=-CF.n_cells; ii<=CF.n_cells; ii++)
			for(jj=-CF.n_cells; jj<=CF.n_cells; jj++)
			  {
			    if (ix+ii < 0 || iy+jj<0 || ix+ii >= C.n[X_AXIS] || iy+jj >= C.n[Y_AXIS]) {
			      nout++;
			      continue;
			    }
			    x = S.XPos[i]-C.caxis[X_AXIS][ix+ii];
			    y = S.YPos[i]-C.caxis[Y_AXIS][iy+jj];
			    distance = sqrt(x*x+y*y);
			    // @todo what if S.RMS[i] == 0.0
			    if ((S.RMS[i] != 0.0) && (OTF.noise_sigma > 0.0)) {
			      rmsweight = 1.0 /(S.RMS[i] * S.RMS[i]);
			    } else {
			      rmsweight = 1.0;
			    }
			    weight = get_weight(&CF, distance) * rmsweight;
#if defined(MDMAXDIM)
			    for(k=0;k<C.n[Z_AXIS];k++)
			      C.cube[iz+k][iy+jj][ix+ii] += weight * spectrum[k];
			    Weight.plane[iy+jj][ix+ii] += weight;
			    if (ii==0 && jj==0)
			      Mask.plane[iy+jj][ix+ii] = 1;
#else			    
			    iz = cube_z_index(&C, C.caxis[X_AXIS][ix+ii], C.caxis[Y_AXIS][iy+jj]);
			    for(k=0;k<C.n[Z_AXIS];k++)
			      C.cube[iz+k] = C.cube[iz+k] + weight * spectrum[k];
			    izp = plane_index(&Weight, C.caxis[X_AXIS][ix+ii], C.caxis[Y_AXIS][iy+jj]);
			    Weight.plane[izp] = Weight.plane[izp] + weight;
			    if (ii==0 && jj==0)
			      Mask.plane[izp] = 1;
#endif			    
			  }		    
		    }
		}
	    }
	}
      free_spec_file(&S);
      printf("Found %d points outside convolving array size +/-%d\n",nout,CF.n_cells);
    }

  printf("Cube Completed, %d spectra accepted\n",ngood);

  // compute averages for each map point; if no data assign NAN
  for(i=0;i<C.n[X_AXIS];i++)
    {
      x = C.caxis[X_AXIS][i];
      for(j=0;j<C.n[Y_AXIS];j++)
	{
	  y = C.caxis[Y_AXIS][j];
#if defined(MDMAXDIM)
	  if(Mask.plane[j][i] > 0.0)
	    for(k=0;k<C.n[Z_AXIS];k++)
	      C.cube[k][j][i] /= Weight.plane[j][i];
	  else
	      for(k=0;k<C.n[Z_AXIS];k++)
		C.cube[k][j][i] = NAN;
#else	  
	  izp = plane_index(&Weight, x, y);
	  iz = cube_z_index(&C, x, y);
	  //if(Weight.plane[izp] > 0.0)     //  classic method with fuzzy edges
	  if(Mask.plane[izp] > 0.0)         //  only expose cells if it had a pixel
	    {
	      for(k=0;k<C.n[Z_AXIS];k++)
		C.cube[iz+k] = C.cube[iz+k] / Weight.plane[izp];
	    }
	  else
	    {
	      for(k=0;k<C.n[Z_AXIS];k++)
		C.cube[iz+k] = NAN;
	    }
#endif	  
	}
    }

  printf("Weighting Completed\n");

  // dumping the spectrum at 0,0 for fun...
#if defined(MDMAXDIM)
  printf("TBD code for mdarray\n");
#else  
  izp = plane_index(&Weight, 0.0, 0.0);
  printf("Weight of %f %f is %f\n",0.0,0.0,Weight.plane[izp]);
  if (Weight.plane[izp] == 0.0)
    printf("*** Warning: zero weights\n");
  iz = cube_z_index(&C, 0.0, 0.0);
#if 1
  // debug
  FILE *fp = fopen("spec-tmp.tab","w!");
  for(i=0;i<S.nchan;i++)
    fprintf(fp,"%d %8.3f %6.2f\n ",i, C.caxis[Z_AXIS][i],C.cube[iz+i]);
  fclose(fp);
#else
  for(i=0;i<S.nchan;i++)
    printf("%d %8.3f %6.2f\n ",i, C.caxis[Z_AXIS][i],C.cube[iz+i]);
#endif
#endif  
  
  printf("write to %s\n",OTF.o_filename);
  // finally write the data cube as FITS file
  write_fits_cube(&C, OTF.o_filename);
  // and the weight plane @todo need a flag for this, 7 times
  if (strlen(OTF.w_filename) > 0) {
    unlink(OTF.w_filename);
#if 1
    printf("write weights to %s\n",OTF.w_filename);
    write_fits_plane(&Weight, OTF.w_filename);
#else
    printf("write 0/1 mask to %s\n",OTF.w_filename);    
    write_fits_plane(&Mask, OTF.w_filename);
#endif    
  }
}
