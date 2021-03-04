#include "single.h"
#include <omp.h>

#define TEST_NUM 1024

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
  int myID, rank;
  double t;
  mpz_t testNum;
  bool isPrime;

  mpz_init(testNum);

  printf("Test Number:");
  gmp_scanf("%Zd", testNum);

  debug_gmp_printf("Input Number : %Zd\n", testNum);

  t = omp_get_wtime();
  // isPrime = MillerRabinTest(testNum);
  // isPrime = SingleMillerTest(testNum);
  isPrime = NaivePrimeTest(testNum);
  t = omp_get_wtime() - t;
  ShowResult(testNum, isPrime);

  printf("Processing Time : %f\n", t);

  mpz_clear(testNum);

  return 0;
}