/* energy.h: Declarations for energy.c

   Copyright 2013 Hendrik Weimer

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

#ifndef __ENERGY_H

#define __ENERGY_H

#include "config.h"
#include "qureg.h"

enum {
  QUANTUM_SOLVER_LANCZOS,
  QUANTUM_SOLVER_LANCZOS_MODIFIED,
  QUANTUM_SOLVER_IMAGINARY_TIME
};

extern double quantum_groundstate(quantum_reg *reg, double epsilon, 
				  quantum_reg H(MAX_UNSIGNED, double), 
				  int solver, double stepsize);

#endif
