 /* grover.c: Implementation of Grover's search algorithm

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

#include <quantum.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#ifdef M_PI
#define pi M_PI
#else
#define pi 3.141592654
#endif

void oracle(int state, quantum_reg *reg)
{
  int i;

  for(i=0;i<reg->width;i++)
    {
      if(!(state & (1 << i)))
	{
	  quantum_sigma_x(i, reg);
	}
    }

  quantum_toffoli(0, 1, reg->width+1, reg);

  for(i=1;i<reg->width;i++)
    {
      quantum_toffoli(i, reg->width+i, reg->width+i+1, reg);
    }

  quantum_cnot(reg->width+i, reg->width, reg);
  
  for(i=reg->width-1;i>0;i--)
    {
      quantum_toffoli(i, reg->width+i, reg->width+i+1, reg);
    }
  
  quantum_toffoli(0, 1, reg->width+1, reg);

  for(i=0;i<reg->width;i++)
    {
      if(!(state & (1 << i)))
	quantum_sigma_x(i, reg);
    }

}

void inversion(quantum_reg *reg)
{
  int i;
  
  for(i=0;i<reg->width;i++)
    quantum_sigma_x(i, reg);

  quantum_hadamard(reg->width-1, reg);

  if(reg->width==3)
    quantum_toffoli(0, 1, 2, reg);

  else
    {
      quantum_toffoli(0, 1, reg->width+1, reg);

      for(i=1;i<reg->width-1;i++)
	{
	  quantum_toffoli(i, reg->width+i, reg->width+i+1, reg);
	}

      quantum_cnot(reg->width+i, reg->width-1, reg);
  
      for(i=reg->width-2;i>0;i--)
	{
	  quantum_toffoli(i, reg->width+i, reg->width+i+1, reg);
	}
  
      quantum_toffoli(0, 1, reg->width+1, reg);
    }
  
  quantum_hadamard(reg->width-1, reg);

  for(i=0;i<reg->width;i++)
    quantum_sigma_x(i, reg);
}
  

void grover(int target, quantum_reg *reg)
{
  int i;

  oracle(target, reg);

  for(i=0;i<reg->width;i++)
    quantum_hadamard(i, reg);

  inversion(reg);

  for(i=0;i<reg->width;i++)
    quantum_hadamard(i, reg);

}

int main(int argc, char **argv)
{
  quantum_reg reg;
  int i, N, width=0;

  srand(time(0));

  if(argc==1)
    {
      printf("Usage: grover [number] [[qubits]]\n\n");
      return 3;
    }

  N=atoi(argv[1]);

  if(argc > 2)
    width = atoi(argv[2]);

  if(width < quantum_getwidth(N+1))
    width = quantum_getwidth(N+1);

  reg = quantum_new_qureg(0, width);
  //  reg.width--;

  quantum_sigma_x(reg.width, &reg);

  for(i=0;i<reg.width;i++)
    quantum_hadamard(i, &reg);

  quantum_hadamard(reg.width, &reg);

  printf("Iterating %i times\n", (int) (pi/4*sqrt(1<<reg.width)));

  for(i=1; i<=pi/4*sqrt(1 << reg.width); i++)
    {
      printf("Iteration #%i\n", i);
      grover(N, &reg);
    }

  quantum_hadamard(reg.width, &reg);

  reg.width++;

  quantum_bmeasure(reg.width-1, &reg);
  
  for(i=0; i<reg.size; i++)
    {
      if(reg.state[i] == N)
	printf("\nFound %i with a probability of %f\n\n", N, 
	       quantum_prob(reg.amplitude[i]));
    }

  quantum_delete_qureg(&reg);

  return 0;
}
