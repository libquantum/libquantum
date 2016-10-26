/* qtime.c: Time evolution of a quantum system

   Copyright 2006-2013 Hendrik Weimer

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

/* Forth-order Runge-Kutta

Flags: QUANTUM_RK4_NODELETE:  Do not delete quantum_reg returned by H
       QUANTUM_RK4_IMAGINARY: Imaginary time evolution */

void
quantum_rk4(quantum_reg *reg, double t, double dt, 
	    quantum_reg H(MAX_UNSIGNED, double), int flags)
{
  quantum_reg k, out, tmp;
  double r = 0;
  int i;
  void *hash;
  int hashw;
  COMPLEX_FLOAT step = dt;

  hash = reg->hash;
  reg->hash = 0;

  hashw = reg->hashw;
  reg->hashw = 0;

  if(!(flags & QUANTUM_RK4_IMAGINARY))
    step *= IMAGINARY;

  /* k1 */
  k = quantum_matrix_qureg(H, t, reg, flags & QUANTUM_RK4_NODELETE);
  quantum_scalar_qureg(-step/2.0, &k);
  tmp = quantum_vectoradd(reg, &k);
  quantum_scalar_qureg(1.0/3.0, &k);
  out = quantum_vectoradd(reg, &k);
  quantum_delete_qureg(&k);

  /* k2 */
  k = quantum_matrix_qureg(H, t+dt/2.0, &tmp, flags & QUANTUM_RK4_NODELETE);
  quantum_delete_qureg(&tmp);
  quantum_scalar_qureg(-step/2.0, &k);
  tmp = quantum_vectoradd(reg, &k);
  quantum_scalar_qureg(2.0/3.0, &k);
  quantum_vectoradd_inplace(&out, &k);
  quantum_delete_qureg(&k);

  /* k3 */
  k = quantum_matrix_qureg(H, t+dt/2.0, &tmp, flags & QUANTUM_RK4_NODELETE);
  quantum_delete_qureg(&tmp);
  quantum_scalar_qureg(-step, &k);
  tmp = quantum_vectoradd(reg, &k);
  quantum_scalar_qureg(1.0/3.0, &k);
  quantum_vectoradd_inplace(&out, &k);
  quantum_delete_qureg(&k);

  /* k4 */
  k = quantum_matrix_qureg(H, t+dt, &tmp, flags & QUANTUM_RK4_NODELETE);
  quantum_delete_qureg(&tmp);
  quantum_scalar_qureg(-step/6.0, &k);
  quantum_vectoradd_inplace(&out, &k);
  
  quantum_delete_qureg(&k);
  quantum_delete_qureg(reg);

  /* Normalize quantum register */

  if(flags & QUANTUM_RK4_IMAGINARY)
    {

      for(i=0; i<out.size; i++)
	r += quantum_prob(out.amplitude[i]);

      quantum_scalar_qureg(sqrt(1.0/r), &out);
    }

  out.hash = hash;
  out.hashw = hashw;

  *reg = out;
  
}

/* Adaptive Runge-Kutta. Stores the new stepsize in dt and returns the
   stepsize actually used. For further details, see Press et al.,
   Numerical Recipes in C (Second Edition, CUP, 1992), Sec. 16.3 */

double
quantum_rk4a(quantum_reg *reg, double t, double *dt, double epsilon, 
	     quantum_reg H(MAX_UNSIGNED, double), int flags)
{
  quantum_reg reg2, old;
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
      quantum_rk4(reg, t, *dt, H, flags);
      quantum_rk4(&reg2, t, *dt/2.0, H, flags);
      quantum_rk4(&reg2, t+*dt/2.0, *dt/2.0, H, flags);

      delta = 0;

      for(i=0;i<reg->size;i++)
	{
	  r = 2*sqrt(quantum_prob(reg->amplitude[i] - reg2.amplitude[i])/
		     quantum_prob(reg->amplitude[i] + reg2.amplitude[i]));
	  
	  if(r > delta)
	    delta = r;
		  
	}
      
      dtused = *dt;

      if(delta < epsilon)
	*dt *= 0.9*pow(epsilon/delta, 0.2);
      else
	*dt *= 0.9*pow(epsilon/delta, 0.25);

      if(*dt > 4*dtused)
	*dt = 4*dtused;
      
      else if(*dt < 0.25*dtused)
	*dt = 0.25*dtused;

      if(delta > epsilon)
	{
	  memcpy(reg->amplitude, old.amplitude, 
		 reg->size*sizeof(COMPLEX_FLOAT));
	  memcpy(reg2.amplitude, old.amplitude, 
		 reg->size*sizeof(COMPLEX_FLOAT));
	  if(reg->state && old.state)
	    memcpy(reg->state, old.state, reg->size*sizeof(MAX_UNSIGNED));
	  if(reg2.state && old.state)
	    memcpy(reg2.state, old.state, reg->size*sizeof(MAX_UNSIGNED));
	}
      
    } while(delta > epsilon);

  reg->hash = hash;
  reg->hashw = hashw;

  quantum_delete_qureg(&old);
  quantum_delete_qureg(&reg2);
  
  return dtused;
}
