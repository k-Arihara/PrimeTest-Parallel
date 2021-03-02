#include "prime.h"

/* Trial division */
bool NaivePrimeTest(mpz_t testNum)
{
  debug_printf("Test : Naive Test\n");
  bool isPrime;
  if (mpz_cmp_ui(testNum, 2) < 0)
  {
    gmp_printf("%Zd is incorrect input.\n", testNum);
    return false;
  }

  if (mpz_cmp_ui(testNum, 2) == 0)
    return true;

  if (mpz_even_p(testNum))
  {
    return false;
  }

  mpz_t op_index, op_limit, result;
  mpz_init(op_index);
  mpz_init(op_limit);
  mpz_init(result);
  mpz_sqrt(op_limit, testNum);
  mpz_set_ui(op_index, 3);

  /* op_index <= op_limit */
  while (!(mpz_cmp(op_index, op_limit) > 0))
  {
    mpz_mod(result, testNum, op_index);
    if (mpz_cmp_ui(result, 0) == 0)
    {
      isPrime = false;
      goto end_of_func;
    }
    mpz_add_ui(op_index, op_index, 2);
  }
  isPrime = true;

end_of_func:
  mpz_clear(result);
  mpz_clear(op_index);
  mpz_clear(op_limit);
  mpz_clear(result);
  return isPrime;
}