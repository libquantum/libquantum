 /* gates.c: Basic gates for quantum register manipulation

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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "matrix.h"
#include "defs.h"
#include "complex.h"
#include "qureg.h"
#include "decoherence.h"
#include "qec.h"

/* Apply a controlled-not gate */

void
quantum_cnot(int control, int target, quantum_reg *reg)
{
  int i;
  int qec;

  quantum_qec_get_status(&qec, NULL);

  if(qec)
    quantum_cnot_ft(control, target, reg);
  else
    {
      for(i=0; i<reg->size; i++)
	{
	  /* Flip the target bit of a base state if the control bit is set */
      
	  if((reg->node[i].state & ((MAX_UNSIGNED) 1 << control)))
	    reg->node[i].state ^= ((MAX_UNSIGNED) 1 << target);
	}
      quantum_decohere(reg);
    }
}

/* Apply a toffoli (or controlled-controlled-not) gate */

void
quantum_toffoli(int control1, int control2, int target, quantum_reg *reg)
{
  int i;
  int qec;

  quantum_qec_get_status(&qec, NULL);

  if(qec)
    quantum_toffoli_ft(control1, control2, target, reg);
  else
    {
      for(i=0; i<reg->size; i++)
	{
	  /* Flip the target bit of a base state if both control bits are
	     set */

	  if(reg->node[i].state & ((MAX_UNSIGNED) 1 << control1))
	    {
	      if(reg->node[i].state & ((MAX_UNSIGNED) 1 << control2))
		{
		  reg->node[i].state ^= ((MAX_UNSIGNED) 1 << target);
		}
	    }
	}
      quantum_decohere(reg);
    }
}

/* Apply a sigma_x (or not) gate */

void
quantum_sigma_x(int target, quantum_reg *reg)
{
  int i;
  int qec;

  quantum_qec_get_status(&qec, NULL);

  if(qec)
    quantum_sigma_x_ft(target, reg);
  else
    {
      for(i=0; i<reg->size; i++)
	{
	  /* Flip the target bit of each base state */

	  reg->node[i].state ^= ((MAX_UNSIGNED) 1 << target);
	} 
      quantum_decohere(reg);
    }
}

/* Apply a sigma_y gate */

void
quantum_sigma_y(int target, quantum_reg *reg)
{
  int i;
  
  for(i=0; i<reg->size;i++)
    {
      /* Flip the target bit of each base state and multiply with 
	 +/- i */

      reg->node[i].state ^= ((MAX_UNSIGNED) 1 << target);
      
      if(reg->node[i].state & ((MAX_UNSIGNED) 1 << target))
	reg->node[i].amplitude *= IMAGINARY;
      else
	reg->node[i].amplitude *= -IMAGINARY;
    }

  quantum_decohere(reg);
}

/* Apply a sigma_y gate */

void
quantum_sigma_z(int target, quantum_reg *reg)
{
  int i;

  for(i=0; i<reg->size; i++)
    {
      /* Multiply with -1 if the target bit is set */

      if(reg->node[i].state & ((MAX_UNSIGNED) 1 << target))
	reg->node[i].amplitude *= -1;
    }
  quantum_decohere(reg);
}

/* Swap the first WIDTH bits of the quantum register. This is done
   classically by renaming the bits, unless QEC is enabled. */

void
quantum_swaptheleads(int width, quantum_reg *reg)
{
  int i, j;
  int pat1, pat2;
  int qec;
  MAX_UNSIGNED l;

  quantum_qec_get_status(&qec, NULL);

  if(qec)
    {
      for(i=0; i<width; i++)
	{
	  quantum_cnot(i, width+i, reg);
	  quantum_cnot(width+i, i, reg);
	  quantum_cnot(i, width+i, reg);
	}
    }
  else
    {
      for(i=0; i<reg->size; i++)
	{
	  /* calculate left bit pattern */
	  
	  pat1 = reg->node[i].state % ((MAX_UNSIGNED) 1 << width);
	  
	  /*calculate right but pattern */
	  
	  pat2 = 0;

	  for(j=0; j<width; j++)
	    pat2 += reg->node[i].state & ((MAX_UNSIGNED) 1 << (width + j));
	  
	  /* construct the new basis state */
	  
	  l = reg->node[i].state - (pat1 + pat2);
	  l += (pat1 << width);
	  l += (pat2 >> width);
	  reg->node[i].state = l;
	}
    }
}

/* Swap WIDTH bits starting at WIDTH and 2*WIDTH+2 controlled by
   CONTROL */

void
quantum_swaptheleads_omuln_controlled(int control, int width, quantum_reg *reg)
{
   int i;

  for(i=0; i<width; i++)
    {
      quantum_toffoli(control, width+i, 2*width+i+2, reg);
      quantum_toffoli(control, 2*width+i+2, width+i, reg);
      quantum_toffoli(control, width+i, 2*width+i+2, reg);
    }
}

/* Apply the 2x2 matrix M to the target bit. M should be unitary and
   having a determinant of 1. */

void 
quantum_gate1(int target, quantum_matrix m, quantum_reg *reg)
{
  int i, j, k;
  int addsize=0, decsize=0;
  COMPLEX_FLOAT t, tnot=0;
  char *done;
  quantum_reg out;
  quantum_reg_hash *p;

  if((m.cols != 2) || (m.rows != 2))
    {
      printf("Matrix is not a 2x2 matrix!\n");
      exit(1);
    }

  /* Build hash table */

  for(i=0; i<(1 << reg->hashw); i++)
    {
      while(reg->hash[i])
	{
	  p = reg->hash[i]->next;
	  free(reg->hash[i]);
	  quantum_memman(-sizeof(quantum_reg_hash));
	  reg->hash[i] = p;
	}
    }
      
  for(i=0; i<reg->size; i++)
    quantum_add_hash(reg->node[i].state, i, reg);

  /* calculate the number of base states to be added */

  for(i=0; i<reg->size; i++)
    {
      j = quantum_get_state(reg->node[i].state ^ ((MAX_UNSIGNED) 1 << target),
			    *reg);
      if(j == -1)
	{
	  if((m.t[1] != 0) && (reg->node[i].state 
			       & ((MAX_UNSIGNED) 1 << target)))
	    addsize++;
	  if((m.t[2] != 0) && !(reg->node[i].state 
				& ((MAX_UNSIGNED) 1 << target)))
	    addsize++;
	}
    }

  /* allocate memory for the new base states */
  
  reg->node = realloc(reg->node, 
		      (reg->size + addsize) * sizeof(quantum_reg_node));
  if(!reg->node) 
    {
      printf("Not enough memory for %i-sized qubit!\n", reg->size + addsize);
      exit(1);
    }
  quantum_memman(addsize*sizeof(quantum_reg_node));

  for(i=0; i<addsize; i++)
    {
      reg->node[i+reg->size].state = 0;
      reg->node[i+reg->size].amplitude = 0;
    }

  done = calloc(reg->size + addsize, sizeof(char));
  if(!done)
    {
      printf("Not enough memory for %i bytes array!\n", 
	     (reg->size + addsize) * sizeof(char));
      exit(1);
    }
  quantum_memman(reg->size + addsize * sizeof(char));

  k = reg->size;

  /* perform the actual matrix multiplication */

  for(i=0; i<reg->size; i++)
    {
      if(!done[i])
	{
	  tnot = 0;
	  j = quantum_get_state(reg->node[i].state 
				^ ((MAX_UNSIGNED) 1<<target), *reg);
	  t = reg->node[i].amplitude;

	  if(j >= 0)
	    tnot = reg->node[j].amplitude;

	  if(reg->node[i].state & ((MAX_UNSIGNED) 1 << target))
	    reg->node[i].amplitude = m.t[2] * tnot + m.t[3] * t;

	  else
	    reg->node[i].amplitude = m.t[0] * t + m.t[1] * tnot;

	  if(j >= 0)
	    {
	      if(reg->node[i].state & ((MAX_UNSIGNED) 1 << target))
		reg->node[j].amplitude = m.t[0] * tnot + m.t[1] * t;

	      else
		reg->node[j].amplitude = m.t[2] * t + m.t[3] * tnot;
	    }

	  
	  else /* new base state will be created */
	    {
	      
	      if((m.t[1] == 0) 
		 && (reg->node[i].state & ((MAX_UNSIGNED) 1 << target)))
		break;
	      if((m.t[2] == 0)
		 && !(reg->node[i].state & ((MAX_UNSIGNED) 1 << target)))
		 break; 

	      reg->node[k].state = reg->node[i].state 
		^ ((MAX_UNSIGNED) 1 << target);

	      if(reg->node[i].state & ((MAX_UNSIGNED) 1 << target))
		reg->node[k].amplitude = m.t[1] * t;

	      else
		reg->node[k].amplitude = m.t[2] * t;

	      k++;
	    }

	  if(j >= 0)
	    done[j] = 1;

	  if(!reg->node[i].amplitude)
	    decsize++;

	  if(j >= 0)
	    {
	      if(!reg->node[j].amplitude)
		decsize++;
	    }
	}
    }

  reg->size += addsize;

  free(done);
  quantum_memman(-reg->size * sizeof(char));

  /* remove base states with zero amplitude */

  if(decsize)
    {
      out.width = reg->width;
      out.size = reg->size - decsize;
      out.node = calloc(out.size, sizeof(quantum_reg_node));
      if(!out.node)
	{
	  printf("Not enough memory for %i-sized qubit!\n", out.size);
	  exit(1);
	}
      quantum_memman(out.size * sizeof(quantum_reg_node));
      out.hashw = reg->hashw;
      out.hash = reg->hash;

      for(i=0, j=0; i<reg->size; i++)
	{
	  if(reg->node[i].amplitude)
	    {
	      out.node[j].state = reg->node[i].state;
	      out.node[j].amplitude = reg->node[i].amplitude;
	      j++;
	    }
	}

      quantum_delete_qureg_hashpreserve(reg);
      *reg = out;
    }

  quantum_decohere(reg);
}

/* Apply a hadamard gate */

void
quantum_hadamard(int target, quantum_reg *reg)
{
  quantum_matrix m;

  m = quantum_new_matrix(2, 2);

  m.t[0] = sqrt(1.0/2);  m.t[1] = sqrt(1.0/2);
  m.t[2] = sqrt(1.0/2);  m.t[3] = -sqrt(1.0/2);

  quantum_gate1(target, m, reg);
  
  quantum_delete_matrix(&m);

}

/* Apply a walsh-hadamard gate */

void
quantum_walsh(int target, quantum_reg *reg)
{
  quantum_matrix m;

  m = quantum_new_matrix(2, 2);

  m.t[0] = IMAGINARY*sqrt(1.0/2);  m.t[1] = IMAGINARY*sqrt(1.0/2);
  m.t[2] = IMAGINARY*sqrt(1.0/2);  m.t[3] = -IMAGINARY*sqrt(1.0/2);

  quantum_gate1(target, m, reg);
  
  quantum_delete_matrix(&m);
  
}

/* Apply a rotation about the x-axis by the angle GAMMA */

void
quantum_r_x(int target, float gamma, quantum_reg *reg)
{
  quantum_matrix m;

  m = quantum_new_matrix(2, 2);

  m.t[0] = cos(gamma / 2);              m.t[1] = -IMAGINARY * sin(gamma / 2);
  m.t[2] = -IMAGINARY * sin(gamma / 2); m.t[3] = cos(gamma / 2);

  quantum_gate1(target, m, reg);

  quantum_delete_matrix(&m);

}

/* Apply a rotation about the y-axis by the angle GAMMA */

void
quantum_r_y(int target, float gamma, quantum_reg *reg)
{
  quantum_matrix m;

  m = quantum_new_matrix(2, 2);

  m.t[0] = cos(gamma / 2);  m.t[1] = -sin(gamma / 2);
  m.t[2] = sin(gamma / 2);  m.t[3] = cos(gamma / 2);

  quantum_gate1(target, m, reg);

  quantum_delete_matrix(&m);

}

/* Apply a rotation about the z-axis by the angle GAMMA */

void
quantum_r_z(int target, float gamma, quantum_reg *reg)
{
  int i;
  
  for(i=0; i<reg->size; i++)
    {
      if(reg->node[i].state & ((MAX_UNSIGNED) 1 << target))
	reg->node[i].amplitude *= quantum_cexp(gamma/2);
      else
	reg->node[i].amplitude *= quantum_cexp(-gamma/2);
    }

  quantum_decohere(reg);
}

/* Scale the phase of qubit */

void
quantum_phase_scale(int target, float gamma, quantum_reg *reg)
{
  int i;
  
  for(i=0; i<reg->size; i++)
    {
      reg->node[i].amplitude *= quantum_cexp(gamma);
    }

  quantum_decohere(reg);
}


/* Apply a phase kick by the angle GAMMA */

void
quantum_phase_kick(int target, float gamma, quantum_reg *reg)
{
  int i;
  
  for(i=0; i<reg->size; i++)
    {
      if(reg->node[i].state & ((MAX_UNSIGNED) 1 << target))
	reg->node[i].amplitude *= quantum_cexp(gamma);
    }

  quantum_decohere(reg);
}

/* Apply a conditional phase shift by PI / 2^(CONTROL - TARGET) */

void
quantum_cond_phase(int control, int target, quantum_reg *reg)
{
  int i;

  for(i=0; i<reg->size; i++)
    {
      if(reg->node[i].state & ((MAX_UNSIGNED) 1 << control))
	{
	  if(reg->node[i].state & ((MAX_UNSIGNED) 1 << target))
	    reg->node[i].amplitude 
	      *= quantum_cexp(pi / ((MAX_UNSIGNED) 1 << (control - target)));
	}
    }

  quantum_decohere(reg);
}


void
quantum_cond_phase_inv(int control, int target, quantum_reg *reg)
{
  int i;

  for(i=0; i<reg->size; i++)
    {
      if(reg->node[i].state & ((MAX_UNSIGNED) 1 << control))
	{
	  if(reg->node[i].state & ((MAX_UNSIGNED) 1 << target))
	    reg->node[i].amplitude 
	      *= quantum_cexp(-pi / ((MAX_UNSIGNED) 1 << (control - target)));
	}
    }

  quantum_decohere(reg);
}



void
quantum_cond_phase_kick(int control, int target, float gamma, quantum_reg *reg)
{
  int i;

  for(i=0; i<reg->size; i++)
    {
      if(reg->node[i].state & ((MAX_UNSIGNED) 1 << control))
	{
	  if(reg->node[i].state & ((MAX_UNSIGNED) 1 << target))
	    reg->node[i].amplitude *= quantum_cexp(gamma);
	}
     }
  quantum_decohere(reg);
}



/* Increase the gate counter by INC steps or reset it if INC < 0. The
   current value of the counter is returned. */

int
quantum_gate_counter(int inc)
{
  static int counter = 0;

  if(inc > 0)
    counter += inc;
  else if(inc < 0)
    counter = 0;

  return counter;
}
