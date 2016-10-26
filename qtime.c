/* qtime.c: Time evolution of a quantum system

   Copyright 2006,2007 Bjoern Butscher, Hendrik Weimer

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
#include <string.h>
#include <stdio.h>

#include "qtime.h"
#include "qureg.h"
#include "complex.h"
#include "config.h"

/* Forth-order Runge-Kutta */

void
quantum_rk4(quantum_reg *reg, double t, double dt, 
	    quantum_reg H(MAX_UNSIGNED, double))
{
  quantum_reg k, out, tmp;
  double r = 0;
  int i;
  void *hash;
  int hashw;

  hash = reg->hash;
  reg->hash = 0;

  hashw = reg->hashw;
  reg->hashw = 0;

  /* k1 */
  k = quantum_matrix_qureg(H, t, reg);
  quantum_scalar_qureg(-IMAGINARY*dt/2.0, &k);
  tmp = quantum_vectoradd(reg, &k);
  quantum_scalar_qureg(1.0/3.0, &k);
  out = quantum_vectoradd(reg, &k);
  quantum_delete_qureg(&k);

  /* k2 */
  k = quantum_matrix_qureg(H, t+dt/2.0, &tmp);
  quantum_delete_qureg(&tmp);
  quantum_scalar_qureg(-IMAGINARY*dt/2.0, &k);
  tmp = quantum_vectoradd(reg, &k);
  quantum_scalar_qureg(2.0/3.0, &k);
  quantum_vectoradd_inplace(&out, &k);
  quantum_delete_qureg(&k);

  /* k3 */
  k = quantum_matrix_qureg(H, t+dt/2.0, &tmp);
  quantum_delete_qureg(&tmp);
  quantum_scalar_qureg(-IMAGINARY*dt, &k);
  tmp = quantum_vectoradd(reg, &k);
  quantum_scalar_qureg(1.0/3.0, &k);
  quantum_vectoradd_inplace(&out, &k);
  quantum_delete_qureg(&k);

  /* k4 */
  k = quantum_matrix_qureg(H, t+dt, &tmp);
  quantum_delete_qureg(&tmp);
  quantum_scalar_qureg(-IMAGINARY*dt/6.0, &k);
  quantum_vectoradd_inplace(&out, &k);
  
  quantum_delete_qureg(&k);
  quantum_delete_qureg(reg);

  /* Normalize quantum register */

  for(i=0; i<out.size; i++)
    r += quantum_prob(out.node[i].amplitude);

  //  quantum_scalar_qureg(sqrt(1.0/r), &out);

  out.hash = hash;
  out.hashw = hashw;

  *reg = out;
  
}

/* Adaptive Runge-Kutta. Stores the new stepsize in dt and returns the
   stepsize actually used. */

double
quantum_rk4a(quantum_reg *reg, double t, double *dt, double epsilon, 
	     quantum_reg H(MAX_UNSIGNED, double))
{
  quantum_reg reg2, old, tmp;
  double delta, r, dtused;
  int i;
  void *hash;
  int hashw;

  hash = reg->hash;
  reg->hash = 0;

  hashw = reg->hashw;
  reg->hashw = 0;

  quantum_copy_qureg(reg, &old);
  quantum_copy_qureg(reg, &reg2);

  do
    {
      quantum_rk4(reg, t, *dt, H);
      quantum_rk4(&reg2, t, *dt/2.0, H);
      quantum_rk4(&reg2, t, *dt/2.0, H);

      delta = 0;

      for(i=0;i<reg->size;i++)
	{
	  if(quantum_real(reg->node[i].amplitude - reg2.node[i].amplitude)
	     > quantum_imag(reg->node[i].amplitude - reg2.node[i].amplitude))
	    r = 2*quantum_real(reg->node[i].amplitude - reg2.node[i].amplitude)
	      / quantum_real(reg->node[i].amplitude + reg2.node[i].amplitude);
	  else
	    r = 2*quantum_imag(reg->node[i].amplitude - reg2.node[i].amplitude)
	      / quantum_imag(reg->node[i].amplitude + reg2.node[i].amplitude);
	  
	  if(r > delta)
	    delta = r;
	}
      
      dtused = *dt;
      *dt *= pow(epsilon/delta, 0.2);

      if(delta > epsilon)
	{
	  tmp = *reg;
	  *reg = old;
	  old = tmp;
	  memcpy(reg2.node, reg->node, reg->size*sizeof(quantum_reg_node));
	  memcpy(old.node, reg->node, reg->size*sizeof(quantum_reg_node));
	}
      
    } while(delta > epsilon);

  reg->hash = hash;
  reg->hashw = hashw;
  
  return dtused;
}
