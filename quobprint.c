/* quobprint.c: Examine quantum object code file

   Copyright 2003 Bjoern Butscher, Hendrik Weimer

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
#include <string.h>
#include <stdlib.h>

#include "objcode.h"
#include "types.h"

int main(int argc, char **argv)
{
  int i, j, k, l;
  FILE *fhd;
  unsigned char operation;
  unsigned char buf[OBJBUF_SIZE];
  char output[OBJBUF_SIZE];
  char opname[256][25];
  MAX_UNSIGNED mu;
  double d;

  strncpy(opname[INIT], "init", 24);
  strncpy(opname[CNOT], "cnot", 24);
  strncpy(opname[TOFFOLI], "toffoli", 24);
  strncpy(opname[SIGMA_X], "sigma_x", 24);
  strncpy(opname[SIGMA_Y], "sigma_y", 24);
  strncpy(opname[SIGMA_Z], "sigma_z", 24);
  strncpy(opname[HADAMARD], "hadamard", 24);
  strncpy(opname[ROT_X], "rotate_x", 24);
  strncpy(opname[ROT_Y], "rotate_y", 24);
  strncpy(opname[ROT_Z], "rotate_z", 24);
  strncpy(opname[PHASE_KICK], "phase_kick", 24);
  strncpy(opname[PHASE_SCALE], "phase_scale", 24);
  strncpy(opname[COND_PHASE], "cond_phase", 24);
  strncpy(opname[CPHASE_KICK], "cond_phase_kick", 24);
  strncpy(opname[MEASURE], "measure", 24);
  strncpy(opname[BMEASURE], "bmeasure", 24);
  strncpy(opname[BMEASURE_P], "bmeasure_preserve", 24);
  strncpy(opname[SWAPLEADS], "swaptheleads", 24);
  strncpy(opname[NOP], "nop", 24);
  if(argc != 2)
    {
      printf("Usage: quobprint [file]\n\n");
      return 1;
    }

  fhd = fopen(argv[1], "r");

  if(!fhd)
    {
      fprintf(stderr, "Could not open %s: ", argv[1]);
      perror(0);
      return 1;
    }

  for(i=0; !feof(fhd); i++)
    {
      for(j=0; j<OBJBUF_SIZE; j++)
	{
	  buf[j] = 0;
	  output[j] = 0;
	}
      operation = fgetc(fhd);

      if(feof(fhd))
	break;

      switch(operation)
	{
	case INIT:
	  fread(buf, sizeof(MAX_UNSIGNED), 1, fhd);
	  mu = quantum_char2mu(buf);
	  printf("%5i: %s %llu\n", i, opname[INIT], mu);
	  break;
	case CNOT:
	case COND_PHASE:
	  fread(buf, sizeof(int), 1, fhd);
	  j = quantum_char2int(buf);
	  fread(buf, sizeof(int), 1, fhd);
	  k = quantum_char2int(buf);
	  printf("%5i: %s %i, %i\n", i, opname[operation], j, k);
	  break;
	case TOFFOLI:
	  fread(buf, sizeof(int), 1, fhd);
	  j = quantum_char2int(buf);
	  fread(buf, sizeof(int), 1, fhd);
	  k = quantum_char2int(buf);
	  fread(buf, sizeof(int), 1, fhd);
	  l = quantum_char2int(buf);
	  printf("%5i: %s %i, %i, %i\n", i, opname[TOFFOLI], j, k, l);
	  break;
	case SIGMA_X:
	case SIGMA_Y:
	case SIGMA_Z:
	case HADAMARD:
	case BMEASURE:
	case BMEASURE_P:
	case SWAPLEADS:
	  fread(buf, sizeof(int), 1, fhd);
	  j = quantum_char2int(buf);
	  printf("%5i: %s %i\n", i, opname[operation], j);
	  break;
	case ROT_X:
	case ROT_Y:
	case ROT_Z:
	case PHASE_KICK:
	case PHASE_SCALE:
	  fread(buf, sizeof(int), 1, fhd);
	  j = quantum_char2int(buf);
	  fread(buf, sizeof(double), 1, fhd);
	  d = quantum_char2double(buf);
	  printf("%5i: %s %i, %f\n", i, opname[operation], j, d);
	  break;
	case CPHASE_KICK:
	  fread(buf, sizeof(int), 1, fhd);
	  j = quantum_char2int(buf);
	  fread(buf, sizeof(int), 1, fhd);
	  k = quantum_char2int(buf);
	  fread(buf, sizeof(double), 1, fhd);
	  d = quantum_char2double(buf);
	  printf("%5i: %s %i, %i, %f\n", i, opname[operation], j, k, d);
	  break;
	case MEASURE:
	case NOP:
	  printf("%5i: %s\n", i, opname[operation]);
	  break;
	default:
	  printf("%i: Unknown opcode 0x(%X)!\n", i, operation);
	  exit(EXIT_FAILURE);
	}

    }

  fclose(fhd);

  return 0;
}
