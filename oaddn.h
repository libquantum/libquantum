#ifndef __OADDN_H

#define __OADDN_H

#include "qureg.h"

extern void test_sum(int, int, quantum_reg *);

extern void muxfa(int, int, int, int, int, int, int, quantum_reg *);

extern void muxfa_inv(int, int, int, int, int, int, int, quantum_reg *);

extern void muxha(int, int, int, int, int, int, quantum_reg *);

extern void muxha_inv(int, int, int, int, int, int, quantum_reg *);

extern void madd(int, int, int, quantum_reg *);

extern void madd_inv(int, int, int,quantum_reg *);

extern void addn(int,int,int, quantum_reg *);

extern void addn_inv(int, int, int, quantum_reg *);

extern void add_mod_n(int, int, int, quantum_reg *);
#endif
