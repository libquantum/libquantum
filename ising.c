/* ising.c: Calculate the ground state of the transverse field Ising model

   Copyright 2013 Bjoern Butscher, Hendrik Weimer

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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <quantum.h>

quantum_reg *hreg;
int N;
double g;
int *V;

quantum_reg H(MAX_UNSIGNED i, double t)
{
  quantum_reg reg;
  int j;

  reg = quantum_new_qureg_sparse(N+1, N);

  /* Transverse field part */

  for(j=0; j<N; j++)
    {
      reg.state[j] = i^(1 << j);
      reg.amplitude[j] = g;
    }

  reg.state[N] = i;

  /* Interaction part */

  reg.amplitude[N] = V[i];

  return reg;
}

quantum_reg H2(MAX_UNSIGNED i, double t)
{
  return hreg[i];
}

int main()
{
  quantum_reg reg;
  int i, j, k;
  double E0, m, m2;

  printf("# Ground state properties of the transverse Ising chain\n");
  printf("# g: Transverse field in units of the Ising interaction\n");
  printf("# N: Number of spins\n");
  printf("# E_0: Ground state energy\n");
  printf("# m: Spontaneous magnetization\n");
  printf("# x: Spin susceptibility\n");
  printf("# g\t\tN\tE_0\t\tm\t\tx\n");

  for(N=8; N<=18; N+=2)
    {

      /* Precompute interaction energies */

      V = calloc(1<<N, sizeof(int));

      assert(V);

      for(i=0; i<(1<<N); i++)
	{
	  k = 0;
	  for(j=0; j<N-1; j++)
	    {
	      if(i & (1<<j))
		{
		  if(i & (1<<(j+1)))
		    k--;
		  else
		    k++;
		}
	      else
		{
		  if(i & (1<<(j+1)))
		    k++;
		  else
		    k--;
		}
	    }

	  /* Periodic boundary conditions */

	  if(i & (1<<(N-1)))
	    {
	      if(i & 1)
		k--;
		  else
		    k++;
	    }
	    else
	      {
		if(i & 1)
		  k++;
		else
		  k--;
	      }
	    
	  V[i] = k;
	}

      for(g=0.9; g<1.1; g+=0.01)
	{

	  reg = quantum_new_qureg_size(1<<N, N);

	  for(i=0; i<(1<<N); i++)
	    reg.amplitude[i] = rand();

	  hreg = calloc(1<<N, sizeof(quantum_reg));

	  assert(hreg);

	  for(i=0; i<(1<<N); i++)
	      hreg[i] = H(i, 0);

	  E0 = quantum_groundstate(&reg, 1e-12, H2, QUANTUM_SOLVER_LANCZOS, 0);

	  m = 0;
	  m2 = 0;

	  for(i=0; i<(1<<N); i++)
	    {
	      k = 0;
	      for(j=0; j<N; j++)
		{
		  if(i & (1<<j))
		    k--;
		  else
		    k++;
		}
	      m += quantum_prob(reg.amplitude[i])*abs(k);
	      m2 += quantum_prob(reg.amplitude[i])*k*k;
	    }

	  m /= N;
	  m2 /= N;

	  printf("%f\t%i\t%f\t%f\t%f\n", g, N, E0, m, m2-m*m);

	  quantum_delete_qureg(&reg);

	  for(i=0; i<(1<<N); i++)
	    quantum_delete_qureg(&hreg[i]);

	  free(hreg);
	}

      free(V);
    
    }

  return 0;

}
