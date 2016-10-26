/* qureg.c: Quantum register management

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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "matrix.h"
#include "qureg.h"
#include "config.h"
#include "complex.h"

/* Convert a vector to a quantum register */

quantum_reg
quantum_matrix2qureg(quantum_matrix *m, int width)
{
  quantum_reg reg;
  int i, j, size=0;

  if(m->cols != 1)
    {
      printf("Error! Cannot convert a multi-column-matrix (%i)!\n", m->cols);
      exit(1);
    }

  reg.width = width;

  /* Determine the size of the quantum register */

  for(i=0; i<m->rows; i++)
    {
      if(m->t[i])
	size++;
    }

  /* Allocate the required memory */

  reg.size = size;
  reg.hashw = width + 2;

  reg.node = calloc(size, sizeof(quantum_reg_node));
  if(!reg.node)
    {
      printf("Not enough memory for %i-sized qubit!\n", size);
      exit(1);
    }
  quantum_memman(size * sizeof(quantum_reg_node));

  /* Allocate the hash table */

  reg.hash = calloc(1 << reg.hashw, sizeof(quantum_reg_hash *));
  if(!reg.hash)
    {
      printf("Not enough memory for %i-sized hash!\n", 1 << reg.hashw);
      exit(1);
    }
  quantum_memman((1 << reg.hashw) * sizeof(quantum_reg_hash *));

  /* Copy the nonzero amplitudes of the vector into the quantum
     register */

  for(i=0, j=0; i<m->rows; i++)
    {
      if(m->t[i])
	{
	  reg.node[j].state = i;
	  reg.node[j].amplitude = m->t[i];
	  j++;
	}
    }

  /* Initialize the PRNG */

  srandom(time(0));

  return reg;
}

/* Create a new quantum register from scratch */

quantum_reg
quantum_new_qureg(MAX_UNSIGNED initval, int width)
{
  quantum_reg reg;

  reg.width = width;
  reg.size = 1;
  reg.hashw = width + 2;

  /* Allocate memory for 1 base state */

  reg.node = calloc(1, sizeof(quantum_reg_node));
  if(!reg.node)
    {
      printf("Not enough memory for %i-sized qubit!\n", 1);
      exit(1);
    }
  quantum_memman(sizeof(quantum_reg_node));

  /* Allocate the hash table */

  reg.hash = calloc(1 << reg.hashw, sizeof(quantum_reg_hash *));
  if(!reg.hash)
    {
      printf("Not enough memory for %i-sized hash!\n", 1 << reg.hashw);
      exit(1);
    }
  quantum_memman((1 << reg.hashw) * sizeof(quantum_reg_hash *));

  /* Initialize the quantum register */
  
  reg.node[0].state = initval;
  reg.node[0].amplitude = 1;

  /* Initialize the PRNG */

  srandom(time(0));

  return reg;
}

/* Convert a quantum register to a vector */

quantum_matrix
quantum_qureg2matrix(quantum_reg reg)
{
  quantum_matrix m;
  int i;

  m = quantum_new_matrix(1, 1 << reg.width);
  
  for(i=0; i<reg.size; i++)
    m.t[reg.node[i].state] = reg.node[i].amplitude;

  return m;
}

/* Delete a quantum register. Note that the hash table is still
   alive. */

void
quantum_delete_qureg(quantum_reg *reg)
{
  free(reg->node);
  quantum_memman(-reg->size * sizeof(quantum_reg_node));
  reg->node = 0;
}
    
/* Print the contents of a quantum register to stdout */

void
quantum_print_qureg(quantum_reg reg)
{
  int i;
  
  for(i=0; i<reg.size; i++)
    {
      printf("%f %+fi|%lli>\n", quantum_real(reg.node[i].amplitude),
	     quantum_imag(reg.node[i].amplitude), reg.node[i].state);
    }

  printf("\n");
}

/* Print the output of the modular exponentation algorithm */

void
quantum_print_expn(quantum_reg reg)
{
  int i;
  
  for(i=0; i<reg.size; i++)
    {
      printf("%i: %lli\n", i, reg.node[i].state - i * (1 << (reg.width / 2)));
    }
}

/* Add additional space to a qureg. It is initialized to zero and can
   be used as scratch space. Note that the space gets added at the LSB */

void
quantum_addscratch(int bits, quantum_reg *reg)
{
  int i, oldwidth;
  MAX_UNSIGNED l;
  
  oldwidth = reg->width;

  reg->width += bits;

  for(i=0; i<reg->size; i++)
    {
      l = reg->node[i].state << bits;
      reg->node[i].state = l;
    }
}

/* Add an element to the hash table */

void
quantum_add_hash(MAX_UNSIGNED a, int pos, quantum_reg *reg)
{
  int i;
  quantum_reg_hash *p;

  i = quantum_hash64(a, reg->hashw);

  p = calloc(1, sizeof(quantum_reg_hash));
  if(!p)
    {
      printf("Not enough memory for hash element!\n");
      exit(1);
    }

  quantum_memman(sizeof(quantum_reg_hash));

  p->i = pos;
  p->next = reg->hash[i];

  reg->hash[i] = p;
  
}

/* Remove an element from the hash table */

void
delete_hash(MAX_UNSIGNED a, int pos, quantum_reg *reg)
{
  int i;
  quantum_reg_hash *p, *q=0;

  i = quantum_hash64(a, reg->hashw);

  p = reg->hash[i];

  while(p->i != pos)
    {
      q = p;
      p = p->next;
    }

  if(q)
    {
      q->next = p->next;
      free(p);
    }
  else
    {
      if(p)
	{
	  reg->hash[i] = p->next;
	  free(p);
	}
      else
	{
	  reg->hash[i] = 0;
	  free(reg->hash[i]);
	}
    }
}

/* Print the hash table to stdout and test if the hash table is
   corrupted */

void
print_hash(quantum_reg reg)
{
  int i;
  quantum_reg_hash *p;
  char *done;

  done = calloc(reg.size, sizeof(char));
  if(!done)
    {
      printf("Not enough memory for %i bytes array!\n", 
	     (reg.size)*sizeof(char));
      exit(1);
    }

  for(i=0; i < (1 << reg.hashw); i++)
    {
      printf("%i: ", i);

      if(reg.hash[i])
	{
	  p = reg.hash[i];

	  while(p->next)
	    {
	      printf("%llu ", reg.node[p->i].state);
	      if(quantum_hash64(reg.node[p->i].state, reg.hashw) != i)
		printf(" Corrupted hash table!\n... ");
	      done[p->i] = 1;
	      p = p->next;
	    }

	  printf("%llu", reg.node[p->i].state);
	  if(quantum_hash64(reg.node[p->i].state, reg.hashw) != i)
	    printf(" Corrupted hash table!");
	  done[p->i] = 1;
	}

      printf("\n");
    }

  /* Test if there are elements in the quantum register which are not
     in the hash table */

  for(i=0; i<reg.size; i++)
    {
      if(!done[i])
	printf("Corrupted hash table: %llu is detached!\n", reg.node[i].state);
    }

  free(done);
}
