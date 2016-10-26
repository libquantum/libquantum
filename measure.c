/* measure.c: Quantum register measurement

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

#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <stdio.h>

#include "qureg.h"
#include "complex.h"
#include "config.h"

/* Generate a uniformly distributed random number between 0 and 1 */

double 
quantum_frand()
{
  return (double) rand() / RAND_MAX;
}

/* Measure the contents of a quantum register */

MAX_UNSIGNED
quantum_measure(quantum_reg reg)
{
  double r;
  int i;

  /* Get a random number between 0 and 1 */
  
  r = quantum_frand();

  for (i=0; i<reg.size; i++)
    {
      /* If the random number is less than the probability of the
	 given base state - r, return the base state as the
	 result. Otherwise, continue with the next base state. */

      r -= quantum_prob_inline(reg.node[i].amplitude);
      if(quantum_prob_inline(reg.node[i].amplitude) >= r)
	return reg.node[i].state;
    }

  /* The sum of all probabilities is less than 1. Usually, the cause
     for this is the application of a non-normalized matrix, but there
     is a slim chance that rounding errors may lead to this as
     well. */

  return -1;
}

/* Measure a single bit of a quantum register. The bit measured is
   indicated by its position POS, starting with 0 as the least
   significant bit. The new state of the quantum register depends on
   the result of the measurement. */

int
quantum_bmeasure(int pos, quantum_reg *reg)
{
  int i, j, k;
  int size=0, result=0;
  double d=0, pa=0, r;
  MAX_UNSIGNED lpat=0, rpat=0, pos2;
  quantum_reg out;

  pos2 = (MAX_UNSIGNED) 1 << pos;

  /* Sum up the probability for 0 being the result */

  for(i=0; i<reg->size; i++)
    {
      if(!(reg->node[i].state & pos2))
	pa += quantum_prob_inline(reg->node[i].amplitude);
    }

  /* Compare the probability for 0 with a random number and determine
     the result of the measurement */

  r = quantum_frand();
  
  if (r > pa)
    result = 1;

  /* Eradicate all amplitudes of base states which have been ruled out
     by the measurement and get the absolute of the new register */

  for(i=0;i<reg->size;i++)
    {
      if(reg->node[i].state & pos2)
	{
	  if(!result)
	    reg->node[i].amplitude = 0;
	  else
	    {
	      d += quantum_prob_inline(reg->node[i].amplitude);
	      size++;
	    }
	}
      else
	{
	  if(result)
	    reg->node[i].amplitude = 0;
	  else
	    {
	      d += quantum_prob_inline(reg->node[i].amplitude);
	      size++;
	    }
	}
    }

  /* Build the new quantum register */

  out.width = reg->width-1;
  out.size = size;
  out.node = calloc(size, sizeof(quantum_reg_node));
  if(!out.node)
    {
      printf("Not enough memory for %i-sized qubit!\n", size);
      exit(1);
    }
  quantum_memman(size * sizeof(quantum_reg_node));
  out.hashw = reg->hashw;
  out.hash = reg->hash;

  /* Determine the numbers of the new base states and norm the quantum
     register */
  
  for(i=0, j=0; i<reg->size; i++)
    {
      if(reg->node[i].amplitude)
	{
	  for(k=0, rpat=0; k<pos; k++)
	  rpat += (MAX_UNSIGNED) 1 << k;

	  rpat &= reg->node[i].state;

	  for(k=sizeof(MAX_UNSIGNED)*8-1, lpat=0; k>pos; k--)
	    lpat += (MAX_UNSIGNED) 1 << k;

	  lpat &= reg->node[i].state;
  
	  out.node[j].state = (lpat >> 1) | rpat;
	  out.node[j].amplitude = reg->node[i].amplitude * 1 / (float) sqrt(d);
	
	  j++;
	}
    }

  quantum_delete_qureg_hashpreserve(reg);
  *reg = out;
  return result;
}

/* Measure a single bit, but do not remove it from the quantum
   register */

int
quantum_bmeasure_bitpreserve(int pos, quantum_reg *reg)
{
  int i, j;
  int size=0, result=0;
  double d=0, pa=0, r;
  MAX_UNSIGNED pos2;
  quantum_reg out;

  pos2 = (MAX_UNSIGNED) 1 << pos;

  /* Sum up the probability for 0 being the result */

  for(i=0; i<reg->size; i++)
    {
      if(!(reg->node[i].state & pos2))
	pa += quantum_prob_inline(reg->node[i].amplitude);
    }

  /* Compare the probability for 0 with a random number and determine
     the result of the measurement */

  r = quantum_frand();
  
  if (r > pa)
    result = 1;

  /* Eradicate all amplitudes of base states which have been ruled out
     by the measurement and get the absolute of the new register */

  for(i=0;i<reg->size;i++)
    {
      if(reg->node[i].state & pos2)
	{
	  if(!result)
	    reg->node[i].amplitude = 0;
	  else
	    {
	      d += quantum_prob_inline(reg->node[i].amplitude);
	      size++;
	    }
	}
      else
	{
	  if(result)
	    reg->node[i].amplitude = 0;
	  else
	    {
	      d += quantum_prob_inline(reg->node[i].amplitude);
	      size++;
	    }
	}
    }

  /* Build the new quantum register */

  out.size = size;
  out.node = calloc(size, sizeof(quantum_reg_node));
  if(!out.node)
    {
      printf("Not enough memory for %i-sized qubit!\n", size);
      exit(1);
    }
  quantum_memman(size * sizeof(quantum_reg_node));
  out.hashw = reg->hashw;
  out.hash = reg->hash;
  out.width = reg->width;

  /* Determine the numbers of the new base states and norm the quantum
     register */
  
  for(i=0, j=0; i<reg->size; i++)
    {
      if(reg->node[i].amplitude)
	{
	  out.node[j].state = reg->node[i].state;
	  out.node[j].amplitude = reg->node[i].amplitude * 1 / (float) sqrt(d);
	
	  j++;
	}
    }

  quantum_delete_qureg_hashpreserve(reg);
  *reg = out;
  return result;
}
