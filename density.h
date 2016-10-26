/* density.h: Declarations for density.c

   Copyright 2004, 2005 Bjoern Butscher, Hendrik Weimer

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

#ifndef __DENSITY_H

#define __DENSITY_H

#include "qureg.h"

#define quantum_density_operation(function, rho, ...) \
do{ \
  int quantum_int; \
  for(quantum_int=0; quantum_int < rho.num; quantum_int++) \
    function(__VA_ARGS__, &rho.reg[quantum_int]); \
} while(0)


struct quantum_density_op_struct
{
  int num;          /* total number of state vectors */
  float *prob;      /* probabilities of the state vectors */
  quantum_reg *reg; /* state vectors */
};

typedef struct quantum_density_op_struct quantum_density_op;

extern quantum_density_op quantum_new_density_op(int num, float *prob,
						 quantum_reg *reg);
extern quantum_density_op quantum_qureg2density_op(quantum_reg *reg);
extern void quantum_reduced_density_op(int pos, quantum_density_op *rho);

extern void quantum_print_density_matrix(quantum_density_op *rho);
extern void quantum_delete_density_op(quantum_density_op *rho);

extern float quantum_purity(quantum_density_op *rho);

#endif
