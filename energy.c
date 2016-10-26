/* energy.c: Compute energetic properties of quantum systems

   Copyright 2013 Hendrik Weimer

   This file is part of libquantum

   libquantum is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation; either version 3 of the License,
   or (at your option) any later version.

   libquantum is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with libquantum; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
   MA 02110-1301, USA

*/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "energy.h"
#include "qureg.h"
#include "qtime.h"
#include "complex.h"

extern void dstevd_(char *jobz, int *n, double *d, double *e, double *z, 
		    int *ldz, double *work, int *lwork, int *iwork, int *liwork,
		    int *info);


/* Modified Lanczos algorithm that iterates over a series of 2x2
   matrix diagonalzations [E. Dagotto & A. Moreo, Phys. Rev. D 31, 865
   (1985)] */

double 
quantum_lanczos_modified(quantum_reg H(MAX_UNSIGNED, double), double epsilon, 
			 quantum_reg *reg)
{
  double E0=DBL_MAX, Eold=DBL_MAX, E1, E2, t;
  quantum_reg tmp, tmp2;
  int i;
  COMPLEX_FLOAT h01;
  double h00, h11;

  for(i=0; i<reg->size; i++)
    {
      quantum_normalize(reg);

      tmp = quantum_matrix_qureg(H, 0, reg, QUANTUM_RK4_NODELETE);

      h00 = quantum_real(quantum_dot_product(&tmp, reg));

      E0 = h00;

      if(fabs(E0-Eold)<epsilon)
	return E0;

      Eold = E0;

      quantum_copy_qureg(reg, &tmp2);
      quantum_scalar_qureg(-h00, &tmp2);

      quantum_vectoradd_inplace(&tmp, &tmp2);

      quantum_normalize(&tmp);

      quantum_delete_qureg(&tmp2);

      tmp2 = quantum_matrix_qureg(H, 0, &tmp, QUANTUM_RK4_NODELETE);

      h11 = quantum_real(quantum_dot_product(&tmp2, &tmp));
      h01 = quantum_dot_product(&tmp2, reg);

      t = sqrt(h11*h11-2*h00*h11+4*h01*quantum_conj(h01)+h00*h00);

      E1 = -(t-h11-h00)/2.;
      E2 = (t+h11+h00)/2.;

      if(E1<E2)
	{
	  quantum_scalar_qureg(-(t-h11+h00)/2./h01, &tmp);
	  quantum_vectoradd_inplace(reg, &tmp);
	}
      else
	{
	  quantum_scalar_qureg((t+h11-h00)/2./h01, &tmp);
	  quantum_vectoradd_inplace(reg, &tmp);
	}
      
      quantum_delete_qureg(&tmp);
      quantum_delete_qureg(&tmp2);
	
    
    }

  quantum_error(QUANTUM_ENOCONVERGE);  
  return nan("0");
  
}

/* Standard Lanczos algorithm without reorthogonalization (see, e.g.,
   [E. Dagotto, Rev. Mod. Phys. 66, 763 (1994)]. */

double 
quantum_lanczos(quantum_reg H(MAX_UNSIGNED, double), double epsilon, 
		quantum_reg *reg)
{
#ifdef HAVE_LIBLAPACK
  double E0=DBL_MAX, Eold=DBL_MAX, *a, *b, *d, *e, norm, *eig, *work;
  quantum_reg *phi, tmp;
  int n, i, j;
  char jobz = 'V';
  int lwork, *iwork, liwork, info;

  phi = calloc(2, sizeof(quantum_reg));
  a = calloc(2, sizeof(double));
  b = calloc(2, sizeof(double));

  work = malloc(sizeof(double));
  iwork = malloc(sizeof(int));

  eig = malloc(sizeof(double));
  d = malloc(sizeof(double));
  e = malloc(sizeof(double));

  if(!(phi && a && b && work && iwork && eig && d && e))
    quantum_error(QUANTUM_ENOMEM);

  quantum_memman(2*sizeof(quantum_reg)+4*sizeof(double));

  quantum_copy_qureg(reg, &phi[0]);
  quantum_normalize(&phi[0]);

  tmp = quantum_matrix_qureg(H, 0, &phi[0], QUANTUM_RK4_NODELETE);

  a[0] = quantum_dot_product(&tmp, &phi[0]);
  
  quantum_copy_qureg(&phi[0], &phi[1]);
  quantum_scalar_qureg(-a[0], &phi[1]);
  quantum_vectoradd_inplace(&phi[1], &tmp);

  quantum_delete_qureg(&tmp);

  tmp = quantum_matrix_qureg(H, 0, &phi[1], QUANTUM_RK4_NODELETE);
  
  norm = quantum_dot_product(&phi[1], &phi[1]);

  a[1] = quantum_dot_product(&tmp, &phi[1]) / norm;
  b[0] = norm / quantum_dot_product(&phi[0], &phi[0]);

  for(n=2; n<reg->size; n++)
    {
      lwork = n*n+4*n+1;
      work = realloc(work, lwork*sizeof(double));

      liwork = 5*n+3;
      iwork = realloc(iwork, lwork*sizeof(int));

      eig = realloc(eig, n*n*sizeof(double));
      d = realloc(d, n*sizeof(double));
      e = realloc(e, n*sizeof(double));

      if(!(work && iwork && eig && d && e))
	quantum_error(QUANTUM_ENOMEM);

      memcpy(d, a, n*sizeof(double));

      for(i=0; i<n; i++)
	e[i] = sqrt(b[i]);

      dstevd_(&jobz, &n, d, e, eig, &n, work, &lwork, iwork, &liwork, &info);

      if(info < 0)
	quantum_error(QUANTUM_ELAPACKARG);

      else if(info > 0)
	quantum_error(QUANTUM_ELAPACKCONV);

      E0 = d[0];

      if(fabs(E0-Eold) < epsilon)
	break;

      Eold = E0;

      phi = realloc(phi, (n+1)*sizeof(quantum_reg));
      a = realloc(a, (n+1)*sizeof(double));
      b = realloc(b, (n+1)*sizeof(double));

      if(!(phi && a && b))
	quantum_error(QUANTUM_ENOMEM);

      quantum_memman(sizeof(quantum_reg)+2*sizeof(double));

      quantum_copy_qureg(&phi[n-1], &phi[n]);
      quantum_scalar_qureg(-a[n-1], &phi[n]);
      quantum_vectoradd_inplace(&phi[n], &tmp);

      quantum_delete_qureg(&tmp);

      quantum_copy_qureg(&phi[n-2], &tmp);
      quantum_scalar_qureg(-b[n-2], &tmp);
      quantum_vectoradd_inplace(&phi[n], &tmp);

      /*      printf("%i %f\n", n, quantum_prob(quantum_dot_product(&phi[n], 
	      &phi[0]))); */

      quantum_delete_qureg(&tmp);

      tmp = quantum_matrix_qureg(H, 0, &phi[n], QUANTUM_RK4_NODELETE);  
      
      norm = quantum_dot_product(&phi[n], &phi[n]);

      a[n] = quantum_dot_product(&tmp, &phi[n]) / norm;
      b[n-1] = norm / quantum_dot_product(&phi[n-1], &phi[n-1]);

    }

  if(n == reg->size)
    {
      quantum_error(QUANTUM_ENOCONVERGE);  
      return nan("0");
    }

  for(i=0; i<n; i++)
    quantum_normalize(&phi[i]);

  for(i=0; i<reg->size; i++)
    {
      reg->amplitude[i] = 0;
      for(j=0; j<n; j++)
	reg->amplitude[i] += eig[j]*phi[j].amplitude[i];
    }

  quantum_delete_qureg(&tmp);

  for(i=0; i<n; i++)
    quantum_delete_qureg(&phi[i]);

  free(phi);
  free(a);
  free(b);
  free(d);
  free(e);
  free(eig);
  free(work);
  free(iwork);

  return E0;


#else
  quantum_error(QUANTUM_ENOLAPACK);

#endif /* HAVE_LIBLAPACK */

}

/* Imaginary time evolution algorithm */

double 
quantum_imaginary_time(quantum_reg H(MAX_UNSIGNED, double),  double epsilon, 
		       double dt, quantum_reg *reg)
{
  double E0=DBL_MAX, Eold=DBL_MAX;
  quantum_reg reg2;
  int i;

  for(i=0; i<reg->size; i++)
    {
      quantum_rk4(reg, 0, dt, H, QUANTUM_RK4_IMAGINARY | QUANTUM_RK4_NODELETE);
      reg2 = quantum_matrix_qureg(H, 0, reg, QUANTUM_RK4_NODELETE);

      E0 =  quantum_real(quantum_dot_product(&reg2, reg));

      quantum_delete_qureg(&reg2);

      if(fabs(Eold-E0)<epsilon)
      	break;
      Eold = E0;
    }

  if(i == reg->size)
    {
      quantum_error(QUANTUM_ENOCONVERGE);
      return nan("0");
    }
  
  else
    return E0;
}

/* Wrapper around the various solver functions */

double 
quantum_groundstate(quantum_reg *reg, double epsilon, 
		    quantum_reg H(MAX_UNSIGNED, double), int solver,
		    double stepsize)
{
  switch(solver)
    {
    case QUANTUM_SOLVER_LANCZOS:
      return quantum_lanczos(H, epsilon, reg);
    case QUANTUM_SOLVER_LANCZOS_MODIFIED:
      return quantum_lanczos_modified(H, epsilon, reg);
    case QUANTUM_SOLVER_IMAGINARY_TIME:
      return quantum_imaginary_time(H, epsilon, stepsize, reg);
    default:
      quantum_error(QUANTUM_ENOSOLVER);
      return nan("0");
    }
}
