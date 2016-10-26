/* decoherence.c: Simulation of decoherence effects

   Copyright 2003 Bjoern Butscher, Hendrik Weimer

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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "measure.h"
#include "qureg.h"
#include "gates.h"
#include "complex.h"

/* Status of the decoherence simulation. Non-zero means enabled and
   decoherence effects will be simulated. */

int status = 0;

/* Decoherence parameter. The higher the value, the greater the
   decoherence impact. */

float lambda = 0;

/* Initialize the decoherence simulation and set the decoherence
   parameter. */

void 
quantum_set_decoherence(float l)
{
  if(l)
    {
      status = 1;
      lambda = l;
    }
  else
    status = 0;
}

/* Perform the actual decoherence of a quantum register for a single
   step of time. This is done by applying a phase shift by a normal
   distributed angle with the variance LAMBDA. */

void
quantum_decohere(quantum_reg *reg)
{
  float u, v, s, x;
  COMPLEX_FLOAT c0, c1;
  int i, j;

  /* Increase the gate counter */

  quantum_gate_counter(1);

  if(status)
    {
      for(j=0; j<reg->width; j++)
	{
	  /* Generate a normal distributed random number */
	  
     	  do {
	    u = 2 * quantum_frand() - 1;
	    v = 2 * quantum_frand() - 1;
	    s = u * u + v * v;
	  } while (s >= 1);

	  x = u * sqrt(-2 * log(s) / s);

	  x *= sqrt(2 * lambda);
	  
	  /* Apply the phase shift gate for decoherence simulation */

	  c0 = quantum_cexp(-x / 2);
	  c1 = quantum_cexp(x / 2);

	  for(i=0; i<reg->size; i++)
	    {
	      if(reg->node[i].state & ((MAX_UNSIGNED) 1 << j))
		reg->node[i].amplitude *= c1;

	      else
		reg->node[i].amplitude *= c0;
	    }
	}
    }
}