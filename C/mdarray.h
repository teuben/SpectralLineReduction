
/*
 * simple multi-dimensional arrays  - sorry, no template programming
 * (real = double or float, depending on the setting in NEMO)
 *
 * MDMAXDIM: maximum number of dimensions we support here;
 * it's hardcoded.
 *
 * See: NEMO/src/kernel/misc/mdarray.c for the MDMAXDIM=8 version
 *
 */

#define _MDARRAY_
#define MDMAXDIM   3

typedef float real, *realptr;     // a NEMO-ism
#define allocate malloc           // no error checking

typedef real        *mdarray1;    /* v1[n1]                              */
typedef mdarray1    *mdarray2;    /* v2[n2][n1]                          */
typedef mdarray2    *mdarray3;    /* v3[n3][n2][n1]                      */

mdarray1 allocate_mdarray1(int n1);
mdarray2 allocate_mdarray2(int n2, int n1);
mdarray3 allocate_mdarray3(int n3, int n2, int n1);

void free_mdarray1(mdarray1 x, int n1);
void free_mdarray2(mdarray2 x, int n2, int n1);
void free_mdarray3(mdarray3 x, int n3, int n2, int n1);

