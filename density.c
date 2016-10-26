/* density.c: Density operator formalism

   Copyright 2004 Bjoern Butscher, Hendrik Weimer

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
#include <stdio.h>

#include "density.h"
#include "qureg.h"
#include "config.h"
#include "matrix.h"
#include "complex.h"
#include "error.h"

/* Build a new density operator from multiple state vectors */

quantum_density_op
quantum_new_density_op(int num, float *prob, quantum_reg *reg)
{
  int i;
  quantum_density_op rho;
  int *phash;
  int hashw;

  rho.num = num;
  
  rho.prob = calloc(num, sizeof(float));

  if(!rho.prob)
    quantum_error(QUANTUM_ENOMEM);
  
  rho.reg = calloc(num, sizeof(quantum_reg));

  if(!rho.reg)
    quantum_error(QUANTUM_ENOMEM);

  quantum_memman(num * (sizeof(float) + sizeof(quantum_reg)));

  /* Take the hash table from the first quantum register */

  rho.prob[0] = prob[0];
  phash = reg[0].hash;
  hashw = reg[0].hashw;
  rho.reg[0] = reg[0];

  /* Destroy the quantum register */

  reg[0].size = 0;
  reg[0].width = 0;
  reg[0].node = 0;
  reg[0].hash = 0;

  for(i=1; i<num; i++)
    {
      rho.prob[i] = prob[i];
      rho.reg[i] = reg[i];
      rho.reg[i].hash = phash;
      rho.reg[i].hashw = hashw;

      reg[i].size = 0;
      reg[i].width = 0;
      reg[i].node = 0;
      reg[i].hash = 0;
    }

  return rho;

}
    
/* Convert a state vector to a density operator */

quantum_density_op
quantum_qureg2density_op(quantum_reg *reg)
{
  float f = 1;

  return quantum_new_density_op(1, &f, reg);
  
}

/* Compute the reduced density operator of a system. Bit POS will be
   traced out. */

void
quantum_reduced_density_op(int pos, quantum_density_op *rho)
{
  int i, j;
  double p0=0, ptmp;
  MAX_UNSIGNED pos2;
  quantum_reg rtmp;

  rho->prob = realloc(rho->prob, 2*rho->num*sizeof(float));

  if(!rho->prob)
    quantum_error(QUANTUM_ENOMEM);

  rho->reg = realloc(rho->reg, 2*rho->num*sizeof(quantum_reg));

  if(!rho->reg)
    quantum_error(QUANTUM_ENOMEM);

  quantum_memman(rho->num * (sizeof(float) + sizeof(quantum_reg)));

  pos2 = (MAX_UNSIGNED) 1 << pos;

  for(i=0; i<rho->num; i++)
    {
      ptmp = rho->prob[i];
      rtmp = rho->reg[i];
      p0 = 0;

      /* Sum up the probability for 0 being the result for this state
	 vector */
  
      for(j=0; j<rho->reg[i].size; j++)
	{
	  if(!(rho->reg[i].node[j].state & pos2))
	    p0 += quantum_prob_inline(rho->reg[i].node[j].amplitude);
	}

      rho->prob[i] = ptmp * p0;
      rho->prob[rho->num + i] = ptmp * (1-p0);

      rho->reg[i] = quantum_state_collapse(pos, 0, rtmp);
      rho->reg[rho->num + i] = quantum_state_collapse(pos, 1, rtmp);

      quantum_delete_qureg_hashpreserve(&rtmp); 
    }

  rho->num *= 2;
  
}

/* Convert the density operator to a full density matrix */

quantum_matrix
quantum_density_matrix(quantum_density_op *rho)
{
  int i, j, k, l1, l2, dim;
  quantum_matrix m;

  dim = 1 << rho->reg[0].width;

  if(dim < 0)
    quantum_error(QUANTUM_EMLARGE);

  m = quantum_new_matrix(dim, dim);

  for(k=0; k<rho->num; k++)
    {
      quantum_reconstruct_hash(&rho->reg[k]);

      for(i=0; i<dim; i++)
	{
	  for(j=0; j<dim; j++)
	    {
	      l1 = quantum_get_state(i, rho->reg[k]);
	      l2 = quantum_get_state(j, rho->reg[k]);
	      if((l1 > -1) && (l2 > -1))
		M(m, i, j) += rho->prob[k] * rho->reg[k].node[l2].amplitude
		  * quantum_conj(rho->reg[k].node[l1].amplitude);
	    }
	}
    }

  return m;
}

/* Print the whole density matrix. */

void
quantum_print_density_matrix(quantum_density_op *rho)
{
  quantum_matrix m;

  m = quantum_density_matrix(rho);
  quantum_print_matrix(m);
  quantum_delete_matrix(&m);

}

/* Delete a density operator */

void
quantum_delete_density_op(quantum_density_op *rho)
{
  int i;

  /* Destroy hash table only once */

  quantum_destroy_hash(&rho->reg[0]);

  for(i=0; i<rho->num; i++)
    quantum_delete_qureg_hashpreserve(&rho->reg[i]);

  free(rho->prob);
  free(rho->reg);
  quantum_memman(-rho->num * (sizeof(float) + sizeof(quantum_reg)));

  rho->prob = 0;
  rho->reg = 0;

}

/* Compute the purity of a density operator */

float
quantum_purity(quantum_density_op *rho)
{
  int i, j , k, l;
  float f = 0;
  COMPLEX_FLOAT g, dp;
  
  /* Diagonal elements */

  for(i=0; i<rho->num; i++)
    f += rho->prob[i]*rho->prob[i];

  for(i=0; i<rho->num; i++)
    {
      for(j=0; j<i; j++)
	{
	  dp = quantum_dot_product(&rho->reg[i], &rho->reg[j]);

	  for(k=0; k<rho->reg[i].size; k++)
	    {
	      /* quantum_dot_product makes sure that rho->reg[j] has a
		 correct hash table */

	      l = quantum_get_state(rho->reg[i].node[k].state, rho->reg[j]);

	      /* Compute p_i p_j <k|\psi_iX\psi_i|\psi_jX\psi_j|k> */
	      
	      if(l > -1)
		g = rho->prob[i] * rho->prob[j] * dp 
		  * rho->reg[i].node[k].amplitude
		  * quantum_conj(rho->reg[j].node[l].amplitude);
	      else
		g = 0;

	      f += 2 * quantum_real(g);

	    }
	}
    }

  return f;

}
