#include "parallel.h"

void ShowResult(mpz_t testNum, bool b)
{
  if (b)
  {
    gmp_printf("Test Number:%Zd may be a Prime Number.\n", testNum);
  }
  else
  {
    gmp_printf("Test Number:%Zd is not a Prime Number.\n", testNum);
  }
}

int main(int argc, char *argv[])
{
  double t;
  mpz_t testNum;
  bool isPrime;

  mpz_init(testNum);

  printf("Test Number:");
  gmp_scanf("%Zd", testNum);

  debug_gmp_printf("Input Number : %Zd\n", testNum);

  t = omp_get_wtime();
  isPrime = ParallelMillerTest(testNum);
  t = omp_get_wtime() - t;
  ShowResult(testNum, isPrime);

  printf("Process Time : %f\n", t);

  mpz_clear(testNum);

  return 0;
}