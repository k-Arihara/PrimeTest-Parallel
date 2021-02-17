#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define DEBUG

/* printf for debuggin */
#ifdef DEBUG
#define debug_printf printf
#else
#define debug_printf 1 ? (void) 0 : printf
#endif

void ShowResult(long testNum, bool b)
{
  if (b)
  {
    printf("Test Number:%d is Prime.\n", testNum);
  }
  else
  {
    printf("Test Number:%d is not Prime.\n", testNum);
  }
}
/******************************************************
 * Miller-Rabin 
 * input: Number 'num' to test a prime number
 * output: True when the number may be a prime number
 *        ,False when it is not a prime number.
 *
 * Reference: https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test
 * 
 *****************************************************/
bool MillerRabin(long testNum)
{
  if (testNum == 2)
  {
    return true;
  }

  if (2 < testNum && (testNum & 1) == 0)
  {
    return false;
  }

  long s = 0;
  long t = testNum - 1;
  while ((t & 1) == 0)
  {
    t = t >> 1;
    s++;
  }
  debug_printf("s:%ld, t:%ld\n", s, t);

  return true;
}

int main(int argc, char *argv[])
{
  int myID, rank;
  long testNum = 19;
  bool result;

  result = MillerRabin(testNum);
  ShowResult(testNum, result);

  testNum = 93;
  result = MillerRabin(testNum);
  ShowResult(testNum, result);

  return 0;
}