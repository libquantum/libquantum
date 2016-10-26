/* density.c: Density operator formalism

   Copyright 2004 Bjoern Butscher, Hendrik Weimer

   This file is part of libquantum

   libquantum is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation; either version 2 of the License,
   or (at your option) any later version.

   libquantum is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with libquantum; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA

*/

#include <stdlib.h>
#include <stdio.h>

#include "density.h"
#include "qureg.h"
#include "config.h"
#include "matrix.h"
#include "complex.h"

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
    {
      printf("Error allocating probability array!\n");
      exit(1);
    }
  
  rho.reg = calloc(num, sizeof(quantum_reg));
  if(!rho.reg)
    {
      printf("Error allocating state vector array!\n");
      exit(1);
    }

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
    {
      printf("Error re-allocating probability array!\n");
      exit(1);
    }

  rho->reg = realloc(rho->reg, 2*rho->num*sizeof(quantum_reg));
  if(!rho->reg)
    {
      printf("Error re-allocating state vector array!\n");
      exit(1);
    }

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

/* Print the whole density matrix. */

void
quantum_print_density_matrix(quantum_density_op *rho)
{
  int i, j, k, dim;
  quantum_matrix m;
  COMPLEX_FLOAT f;

  dim = 1 << rho->reg[0].width;

  if(dim < 0)
    {
      printf("Density matrix is too big!\n");
      exit(1);
    }

  m = quantum_new_matrix(dim, dim);

  /* \rho_ij = \sum_k p_k <i|\psi_kX\psi_k|j> */

  for(i=0; i<rho->num; i++)
    {
      for(j=0; j<rho->reg[i].size; j++)
	{
	  M(m, rho->reg[i].node[j].state, rho->reg[i].node[j].state) 
	    += rho->prob[i] 
	    * quantum_prob_inline(rho->reg[i].node[j].amplitude);
	  
	  for(k=0; k<j; k++)
	    {
	      f = rho->prob[i] * quantum_conj(rho->reg[i].node[j].amplitude)
		* rho->reg[i].node[k].amplitude;
	      M(m, rho->reg[i].node[j].state , rho->reg[i].node[k].state) += f;
	      M(m, rho->reg[i].node[k].state , rho->reg[i].node[j].state) 
		+= quantum_conj(f);
	    }
	}
    }

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
