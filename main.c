#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#ifdef _OPENMP
#include <omp.h>
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

  return true;
}

int main(int argc, char *argv[])
{
  int myID, rank;
  long testNum = 3;
  bool result;

  result = MillerRabin(testNum);
  ShowResult(testNum, result);

  testNum = 6;
  result = MillerRabin(testNum);
  ShowResult(testNum, result);

  return 0;
}