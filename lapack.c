/* lapack.c: LAPACK interface

   Copyright 2008 Bjoern Butscher, Hendrik Weimer

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

#include "lapack.h"
#include "matrix.h"
#include "complex.h"
#include "qureg.h"
#include "error.h"
#include "config.h"

extern void cheev_(char *jobz, char *uplo, int *n, float _Complex *A, int *lda,
		   float *w, float _Complex *work, int *lwork, float *rwork,
		   int *info);

void 
quantum_diag_time(float t, quantum_reg *reg0, quantum_reg *regt, 
		  quantum_reg *tmp1, quantum_reg *tmp2, quantum_matrix H, 
		  float **w)
{
#ifdef HAVE_LIBLAPACK
  char jobz = 'V';
  char uplo = 'U';
  int dim = H.cols;
  COMPLEX_FLOAT *work;
  int lwork = -1;
  float rwork[3*dim-2];
  int info;
  int i;
  void *p;
  
  if(tmp2->size != reg0->size)
    {
      /* perform diagonalization */

      p = regt->node;
      *regt = *reg0;
      regt->node = realloc(p, regt->size*sizeof(quantum_reg_node));
      for(i=0; i<reg0->size; i++)
	regt->node[i].state = i;
      
      p = tmp1->node;
      *tmp1 = *reg0;
      tmp1->node = realloc(p, regt->size*sizeof(quantum_reg_node));
      for(i=0; i<reg0->size; i++)
	tmp1->node[i].state = i;

      p = tmp2->node;
      *tmp2 = *reg0;
      tmp2->node = realloc(p, regt->size*sizeof(quantum_reg_node));
      for(i=0; i<reg0->size; i++)
	tmp2->node[i].state = i;

      *w = malloc(dim*sizeof(float));

      if(!*w)
	quantum_error(QUANTUM_ENOMEM);

      work = malloc(sizeof(COMPLEX_FLOAT));

      if(!work)
	quantum_error(QUANTUM_ENOMEM);

      cheev_(&jobz, &uplo, &dim, H.t, &dim, *w, work, &lwork, rwork, &info);

      if(info < 0)
	quantum_error(QUANTUM_ELAPACKARG);

      else if(info > 0)
	quantum_error(QUANTUM_ELAPACKCHEEV);
      
      lwork = (int) work[0];
      work = realloc(work, lwork*sizeof(COMPLEX_FLOAT));

      if(!work)
	quantum_error(QUANTUM_ENOMEM);

      cheev_(&jobz, &uplo, &dim, H.t, &dim, *w, work, &lwork, rwork, &info);

      if(info < 0)
	quantum_error(QUANTUM_ELAPACKARG);

      else if(info > 0)
	quantum_error(QUANTUM_ELAPACKCHEEV);
      
      free(work);

      quantum_mvmult(tmp1, H, reg0);

      quantum_adjoint(&H);
    }

  if(tmp1->size != reg0->size)
    {
      p = regt->node;
      *regt = *reg0;
      regt->node = realloc(p, regt->size*sizeof(quantum_reg_node));
      for(i=0; i<reg0->size; i++)
	regt->node[i].state = i;

      p = tmp1->node;
      *tmp1 = *reg0;
      tmp1->node = realloc(p, regt->size*sizeof(quantum_reg_node));
      for(i=0; i<reg0->size; i++)
	tmp1->node[i].state = i;

      quantum_adjoint(&H);
      
      quantum_mvmult(tmp1, H, reg0);

      quantum_adjoint(&H);
    }

  for(i=0; i<dim; i++)
    tmp2->node[i].amplitude = quantum_cexp(-(*w)[i]*t)*tmp1->node[i].amplitude;

  quantum_mvmult(regt, H, tmp2);

#else
  quantum_error(QUANTUM_ENOLAPACK);

#endif /* HAVE_LIBLAPACK */
  
}


