#ifndef _CUBE_H_
#define _CUBE_H_

#define X_AXIS 1
#define Y_AXIS 0
#define Z_AXIS 2

typedef struct 
{
  int obsnum;
  char source[32];
  float x_position, y_position;
  float *cube;
  float *caxis[3];
  float crval[3], crpix[3], cdelt[3];
  char ctype[3][16], units[3][16];
  int n[3],ncube,nplane;
} Cube;

void initialize_cube(Cube* C, int *n);

void initialize_cube_axis(Cube* C, int axis, float crval, float crpix, float cdelt, char *ctype, char *units);

int cube_axis_index(Cube* C, int axis, float value);

int cube_z_index(Cube* C, float x, float y);

void write_netcdf_cube(Cube* C, char *filename); 

void write_fits_cube(Cube *C, char *filename);
void print_fits_error(int);

#endif
