#ifndef __OMULN_H

#define __OMULN_H

#include "qureg.h"

void emul(int, int, int, quantum_reg *);

void muln(int, int, int, int, quantum_reg *);

void muln_inv(int, int, int, int, quantum_reg *);

void mul_mod_n(int, int, int, int, quantum_reg *);

#endif
