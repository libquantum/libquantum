/* time.c: Declarations for lapack.h

   Copyright 2006 Bjoern Butscher, Hendrik Weimer

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

#ifndef __LAPACK_H

#define __LAPACK_H

#include "config.h"
#include "matrix.h"
#include "qureg.h"

extern void quantum_diag_time(float t, quantum_reg *reg0, quantum_reg *regt, 
			      quantum_reg *tmp1, quantum_reg *tmp2, 
			      quantum_matrix H, float **w);

#endif
