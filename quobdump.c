/* quopdump.c: Generate quantum opcode from a program

   Copyright 2003 Bjoern Butscher, Hendrik Weimer

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

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include "objcode.h"

int main(int argc, char **argv)
{
  char *envstr;
  extern char **environ;

  if(argc < 3)
    {
      printf("Usage: quopdump [file] [program] [[args]]\n\n");
      return 1;
    }

  envstr = malloc(strlen(argv[1]) + 20);

  snprintf(envstr, strlen(argv[1]) + 19, "QUOBFILE=%s", argv[1]);

  putenv(envstr);

  execve(argv[2], &argv[2], environ);

  fprintf(stderr, "Unable to execute %s: ", argv[2]);
  perror(NULL);

  return 1;
}
