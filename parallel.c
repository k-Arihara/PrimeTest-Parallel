#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <gmp.h>
#include <stdbool.h>
#include <omp.h>

// #define DEBUG

/* printf for debuggin */
#ifdef DEBUG
#define debug_printf printf
#define debug_gmp_printf gmp_printf
#else
#define debug_printf 1 ? (void)0 : printf
#define debug_gmp_printf 1 ? (void)0 : gmp_printf
#endif

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

/******************************************************
 * Miller-Rabin 
 * input: Number 'testNum' to test a prime number
 * output: True when the number may be a prime number
 *        ,False when it is not a prime number.
 *
 * Reference: https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test
 * 
 *****************************************************/
bool MillerTest(mpz_t testNum)
{
  bool isPrime;
  mpz_t result;
  mpz_init(result);

  /* testNum < 2 that return false. */
  if (mpz_cmp_ui(testNum, 2) < 0)
  {
    gmp_printf("%Zd is incorrect input.\n", testNum);
    isPrime = false;
    goto end_of_func;
  }

  /* testNum == 2 */
  if (mpz_cmp_ui(testNum, 2) == 0)
  {
    isPrime = true;
    goto end_of_func;
  }

  /* result = testNum % 2 */
  mpz_mod_ui(result, testNum, 2);
  /*  2 < testNum && result == 0 */
  if (mpz_cmp_ui(testNum, 2) > 0 && mpz_cmp_ui(result, 0) == 0)
  {
    isPrime = false;
    goto end_of_func;
  }

  mpz_t op_s, op_t, op_u;
  mpz_init(op_s);
  mpz_init(op_t);
  mpz_init(op_u);

  /* op_s = 0 */
  mpz_set_ui(op_s, 0);
  /* op_t = testNum -1 */
  mpz_sub_ui(op_t, testNum, 1);
  mpz_sub_ui(op_u, testNum, 1);

  debug_gmp_printf("s:%Zd, t:%Zd\n", op_s, op_t);

  /* result = op_t % 2 */
  mpz_mod_ui(result, op_t, 2);
  while (mpz_cmp_ui(result, 0) == 0)
  {
    /* t = t >> 1 */
    mpz_tdiv_q_ui(op_t, op_t, 2);
    /* s++ */
    mpz_add_ui(op_s, op_s, 1);
    /* result = op_t % 2 */
    mpz_mod_ui(result, op_t, 2);
  }
  debug_gmp_printf("s:%Zd, t:%Zd\n", op_s, op_t);

  mpz_t op_a;
  mpz_init(op_a);
  mpz_set_ui(op_a, 1);

  bool isLoop = true;
  isPrime = true;
#pragma omp parallel shared(isLoop, isPrime, op_a) firstprivate(testNum)
  {
    mpz_t op_private_a, op_private_s, op_private_t, op_private_u, private_result;
    mpz_init(op_private_a);
    mpz_init(op_private_s);
    mpz_init(op_private_t);
    mpz_init(op_private_u);
    mpz_init(private_result);
    mpz_set(op_private_s, op_s);
    mpz_set(op_private_t, op_t);
    mpz_set(op_private_u, op_u);

#pragma omp critical(init_op_a)
    {
      mpz_add_ui(op_a, op_a, 1);
      mpz_set(op_private_a, op_a);
      debug_printf("thread num : %d\n", omp_get_thread_num());
      debug_gmp_printf("%Zd\n", op_private_a);
      debug_gmp_printf("%Zd\n", op_private_s);
      debug_gmp_printf("%Zd\n", op_private_t);
      debug_gmp_printf("%Zd\n", op_private_u);
      debug_gmp_printf("\n");
    }

    if (mpz_cmp(op_private_a, op_private_u) > 0)
      goto end_of_loop;

    while (isLoop)
    {
      /* Felmat Test */
      /* result = a^t % testNum */
      mpz_powm(private_result, op_private_a, op_private_t, testNum);
      /* result == 1 */
      if (mpz_cmp_ui(private_result, 1) == 0)
        goto end_of_loop;

      /* while(s > i) */
      for (unsigned int i = 0; mpz_cmp_ui(op_private_s, (unsigned long)i) > 0; i++)
      {
        /* result = 2^i */
        mpz_ui_pow_ui(private_result, 2, (unsigned long)i);
        /* result = result * t */
        mpz_mul(private_result, private_result, op_t);
        /* result = a^result % testNum */
        mpz_powm(private_result, op_private_a, private_result, testNum);
        /* a^(2^i * t) % testNum == testNum - 1 */
        if (mpz_cmp(private_result, op_private_u) == 0)
          goto end_of_loop;
      }
      isLoop = false;
      isPrime = false;
      break;

    end_of_loop:
#pragma omp critical(increment_op_a)
      mpz_add_ui(op_a, op_a, 1);
      mpz_set(op_private_a, op_a);
      if (mpz_cmp(op_private_a, op_private_u) > 0)
        isLoop = false;
      else
        mpz_set(op_private_a, op_a);
    }
    mpz_clear(op_private_a);
    mpz_clear(op_private_s);
    mpz_clear(op_private_t);
    mpz_clear(op_private_u);
    mpz_clear(private_result);
#pragma omp barrier
  }

  mpz_clear(op_a);
  mpz_clear(op_s);
  mpz_clear(op_t);
  mpz_clear(op_u);
end_of_func:
  mpz_clear(result);
  return isPrime;
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
  isPrime = MillerTest(testNum);
  t = omp_get_wtime() - t;
  ShowResult(testNum, isPrime);

  printf("Process Time : %f\n", t);

  mpz_clear(testNum);

  return 0;
}