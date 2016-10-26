#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "defs.h"
#include "gates.h"
#include "omuln.h"
#include "qureg.h"

void 
quantum_exp_mod_n(int N, int x, int width_input, int width, quantum_reg *reg)
{
	
	int i, j, f;
	
	

	quantum_sigma_x(2*width+2, reg);
	for (i=1; i<=width_input;i++){
		f=x%N;			//compute
		for (j=1;j<i;j++)
		  { 
		    f*=f;	//x^2^(i-1)
		    f= f%N;
		  }
		mul_mod_n(N,f,3*width+1+i, width, reg);
		}
	}
