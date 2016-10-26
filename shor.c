#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <quantum.h>

int main(int argc, char **argv) {

  quantum_reg qr;
  int i;
  int width, swidth;
  int x;
  int N;
  int c,q,a,b, factor;

  srandom(time(0));

  if(argc==1)
    {
      printf("Usage: sim [number]\n\n");
      return 3;
    }

  N=atoi(argv[1]);

  if(N<15)
    {
      printf("Invalid number\n\n");
      return 3;
    }

  width=quantum_getwidth(N*N);
  swidth=quantum_getwidth(N);

  printf("N = %i, %i qubits required\n", N, width+3*swidth+2);

  do
    {
      x = random() % N;
    } while((quantum_gcd(N, x) > 1) || (x < 2));

  printf("Random seed: %i\n", x);

  qr=quantum_new_qureg(0, width);

  for(i=0;i<width;i++)
    quantum_hadamard(i, &qr);

  quantum_addscratch(3*swidth+2, &qr);

  quantum_exp_mod_n(N, x, width, swidth, &qr);

  for(i=0;i<3*swidth+2;i++)
    {
      quantum_bmeasure(0, &qr);
    }

  quantum_qft(width, &qr); 

  //  quantum_print_qureg(qr);

  c=quantum_measure(qr);

  if(c==-1)
    {
      printf("Impossible Measurement!\n");
      exit(1);
    }

  if(c==0)
    {
      printf("Measured zero, try again.\n");
      exit(2);
    }

  q = 1<<(width);

  printf("Measured %i (%f), ", c, (float)c/q);

  quantum_frac_approx(&c, &q, width);

  printf("fractional approximation is %i/%i.\n", c, q);

  if((q % 2 == 1) && (2*q<(1<<width)))
    {
      printf("Odd denominator, trying to expand by 2.\n");
      q *= 2;
    }
    
  if(q % 2 == 1)
    {
      printf("Odd period, try again.\n");
      exit(2);
    }

  printf("Possible period is %i.\n", q);
  
  a = quantum_ipow(x, q/2) + 1 % N;
  b = quantum_ipow(x, q/2) - 1 % N;
  
  a = quantum_gcd(N, a);
  b = quantum_gcd(N, b);
  
  if(a>b)
    factor=a;
  else
    factor=b;

  if((factor < N) && (factor > 1))
    {
      printf("%i = %i * %i\n", N, factor, N/factor);
    }
  else
    {
      printf("Unable to determine factors, try again.\n");
      exit(2);
    }
    
  quantum_delete_qureg(&qr);

  //  printf("Memory usage: %i bytes\n", (int )memman(0));

  return 0;
}
