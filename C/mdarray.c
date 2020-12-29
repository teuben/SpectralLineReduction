/*
 * mdarray: a simple multi-dimensional array (allocator)
 *          cannot be used in a multi-threaded environment
 *
 * See also: Iliffe Vector
 *           Shortridge ADASS 2019 proceedings - https://github.com/KnaveAndVarlet/ADASS2019
 *
 * Note:  MDMAXDIM defines the highest dimension we've implemented (currently 7)
 *        see mdarray.h
 *
 *  5-may-2003   Created                                 Peter Teuben
 * 18-feb-2006   use sequential memory for the data      PJT
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "mdarray.h"

static int top = 0;    /* highest mdarray will have this at 0, each call increases top */
static real *p = 0;    /* keep pointer to the big alloc incrementing by n1 */


mdarray1 allocate_mdarray1(int n1)
{
  real *p1;
  if (top) {
    //dprintf(7,"  top @ 0x%x\n",p);
    p1 = p;
    p += n1;
    return p1;
  } else
    return (mdarray1 )allocate(n1*sizeof(real));
}

mdarray2 allocate_mdarray2(int n2,int n1)
{
  mdarray2 x = (mdarray2)allocate(sizeof(mdarray2)*n2);
  int i;
  if (top == 0) 
    p = (real *) allocate(sizeof(real)*n1*n2);
  top++;
  for (i=0;i<n2;i++)
    x[i] = allocate_mdarray1(n1);
  top--;
  return x;
}

mdarray3 allocate_mdarray3(int n3,int n2,int n1)
{
  mdarray3 x = (mdarray3)allocate(sizeof(mdarray3)*n3);
  int i;
  if (top == 0)
    p = (real *) allocate(sizeof(real)*n1*n2*n3);
  top++;
  for (i=0;i<n3;i++)
    x[i] = allocate_mdarray2(n2,n1);
  top--;
  return x;
}

void free_mdarray1(mdarray1 x, int n1)
{
  if (top==0) free(x);
}

void free_mdarray2(mdarray2 x, int n2, int n1)
{
  int i;
  if (top==0) free(x[0]);
  top++;
  for (i=0;i<n2;i++)
    free_mdarray1(x[i],n1);
  top--;
  free(x);
}

void free_mdarray3(mdarray3 x, int n3, int n2, int n1)
{
  int i;
  if (top==0) free(x[0][0]);
  top++;
  for (i=0;i<n3;i++)
    free_mdarray2(x[i],n2,n1); 
  top--;
  free(x);
}



