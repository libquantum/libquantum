/* decoherence.c: Simulation of decoherence effects

   Copyright 2003 Bjoern Butscher, Hendrik Weimer

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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "measure.h"
#include "qureg.h"
#include "gates.h"
#include "complex.h"
#include "error.h"

/* Status of the decoherence simulation. Non-zero means enabled and
   decoherence effects will be simulated. */

int quantum_status = 0;

/* Decoherence parameter. The higher the value, the greater the
   decoherence impact. */

float quantum_lambda = 0;

float
quantum_get_decoherence()
{
  return quantum_lambda;
}

/* Initialize the decoherence simulation and set the decoherence
   parameter. */

void 
quantum_set_decoherence(float l)
{
  if(l)
    {
      quantum_status = 1;
      quantum_lambda = l;
    }
  else
    quantum_status = 0;
}

/* Perform the actual decoherence of a quantum register for a single
   step of time. This is done by applying a phase shift by a normal
   distributed angle with the variance LAMBDA. */

void
quantum_decohere(quantum_reg *reg)
{
  float u, v, s, x;
  float *nrands;
  float angle;
  int i, j;

  /* Increase the gate counter */

  quantum_gate_counter(1);

  if(quantum_status)
    {
      
      nrands = calloc(reg->width, sizeof(float));

      if(!nrands)
	quantum_error(QUANTUM_ENOMEM);

      quantum_memman(reg->width * sizeof(float));

      for(i=0; i<reg->width; i++)
	{
	  /* Generate normal distributed random numbers */
	  
     	  do {
	    u = 2 * quantum_frand() - 1;
	    v = 2 * quantum_frand() - 1;
	    s = u * u + v * v;
	  } while (s >= 1);

	  x = u * sqrt(-2 * log(s) / s);

	  x *= sqrt(2 * quantum_lambda);

	  nrands[i] = x/2;
	}

  
      /* Apply the phase shifts for decoherence simulation */

      for(i=0; i<reg->size; i++)
	{
	  angle = 0;

	  for(j=0; j<reg->width; j++)
	    {
	      if(reg->state[i] & ((MAX_UNSIGNED) 1 << j))
		angle += nrands[j];
	      else
		angle -= nrands[j];
	    }

	  reg->amplitude[i] *= quantum_cexp(angle);
	  
	}
      free(nrands);
      quantum_memman(-reg->width * sizeof(float));  
  
    }
}
