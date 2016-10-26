/* qtime.h: Declarations for qtime.c

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

#ifndef __QTIME_H

#define __QTIME_H

#include "qureg.h"
#include "config.h"

extern void quantum_rk4(quantum_reg *reg, double t, double dt, 
			quantum_reg H(MAX_UNSIGNED, double), int flags);
extern double quantum_rk4a(quantum_reg *reg, double t, double *dt, 
			   double epsilon, 
			   quantum_reg H(MAX_UNSIGNED, double), int flags);

#endif
