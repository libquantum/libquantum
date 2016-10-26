/* lapack.c: LAPACK interface

   Copyright 2008-2013 Hendrik Weimer

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

#include <stdlib.h>
#include <math.h>

#include "lapack.h"
#include "matrix.h"
#include "complex.h"
#include "qureg.h"
#include "error.h"
#include "config.h"

extern void cheev_(char *jobz, char *uplo, int *n, float _Complex *A, int *lda,
		   float *w, float _Complex *work, int *lwork, float *rwork,
		   int *info);

extern void zheev_(char *jobz, char *uplo, int *n, double _Complex *A, int *lda,
		   double *w, double _Complex *work, int *lwork, double *rwork,
		   int *info);

void 
quantum_diag_time(double t, quantum_reg *reg0, quantum_reg *regt, 
		  quantum_reg *tmp1, quantum_reg *tmp2, quantum_matrix H, 
		  REAL_FLOAT **w)
{
#ifdef HAVE_LIBLAPACK
  char jobz = 'V';
  char uplo = 'U';
  int dim = H.cols;
  COMPLEX_FLOAT *work;
  int lwork = -1;
  REAL_FLOAT rwork[3*dim-2];
  int info;
  int i, j;
  void *p;
  
  if(tmp2->size != reg0->size)
    {
      /* perform diagonalization */

      for(i=0; i<dim; i++)
	{
	  for(j=0; j<dim; j++)
	    {
	      if(sqrt(quantum_prob(M(H, i, j) - quantum_conj(M(H, j, i)))) 
		 > 1e-6)
		quantum_error(QUANTUM_EHERMITIAN);
	    }
	}

      p = regt->amplitude;
      *regt = *reg0;
      regt->amplitude = realloc(p, regt->size*sizeof(COMPLEX_FLOAT));
      
      p = tmp1->amplitude;
      *tmp1 = *reg0;
      tmp1->amplitude = realloc(p, regt->size*sizeof(COMPLEX_FLOAT));

      p = tmp2->amplitude;
      *tmp2 = *reg0;
      tmp2->amplitude = realloc(p, regt->size*sizeof(COMPLEX_FLOAT));

      if(!(regt->amplitude && tmp1->amplitude && tmp2->amplitude))
	quantum_error(QUANTUM_ENOMEM);

      *w = malloc(dim*sizeof(float));

      if(!*w)
	quantum_error(QUANTUM_ENOMEM);

      work = malloc(sizeof(COMPLEX_FLOAT));

      if(!work)
	quantum_error(QUANTUM_ENOMEM);

      QUANTUM_LAPACK_SOLVER(&jobz, &uplo, &dim, H.t, &dim, *w, work, &lwork, 
			    rwork, &info);

      if(info < 0)
	quantum_error(QUANTUM_ELAPACKARG);

      else if(info > 0)
	quantum_error(QUANTUM_ELAPACKCONV);
      
      lwork = (int) work[0];
      work = realloc(work, lwork*sizeof(COMPLEX_FLOAT));

      if(!work)
	quantum_error(QUANTUM_ENOMEM);

      QUANTUM_LAPACK_SOLVER(&jobz, &uplo, &dim, H.t, &dim, *w, work, &lwork, 
			    rwork, &info);

      if(info < 0)
	quantum_error(QUANTUM_ELAPACKARG);

      else if(info > 0)
	quantum_error(QUANTUM_ELAPACKCONV);
      
      free(work);

      quantum_mvmult(tmp1, H, reg0);

      quantum_adjoint(&H);
    }

  if(tmp1->size != reg0->size)
    {
      p = regt->amplitude;
      *regt = *reg0;
      regt->amplitude = realloc(p, regt->size*sizeof(COMPLEX_FLOAT));

      p = tmp1->amplitude;
      *tmp1 = *reg0;
      tmp1->amplitude = realloc(p, regt->size*sizeof(COMPLEX_FLOAT));

      quantum_adjoint(&H);
      
      quantum_mvmult(tmp1, H, reg0);

      quantum_adjoint(&H);
    }

  for(i=0; i<dim; i++)
    tmp2->amplitude[i] = quantum_cexp(-(*w)[i]*t)*tmp1->amplitude[i];

  quantum_mvmult(regt, H, tmp2);

#else
  quantum_error(QUANTUM_ENOLAPACK);

#endif /* HAVE_LIBLAPACK */
  
}


