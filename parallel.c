#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <stdbool.h>
#include <omp.h>

#define DEBUG

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

bool gmp2mpfr(mpz_t op_in, mpfr_t op_out)
{
  char str[1024];

  mpz_get_str(str, 10, op_in);
  printf("str:%s\n", str);
  mpfr_set_str(op_out, str, 10, MPFR_RNDN);
  mpfr_ceil(op_out, op_out);
  mpfr_out_str(stdout, 10, 0, op_out, MPFR_RNDD);
  putchar('\n');

  return true;
}

void GetLimitNum(mpz_t op_limit, mpz_t testNum, mpz_t testNum_1)
{
  mpfr_t op_result;
  mpfr_init(op_result);

  if (!gmp2mpfr(testNum, op_result))
  {
    mpz_set(op_limit, testNum);
    return;
  }

  /* op_result = loge(op_result) */
  mpfr_log(op_result, op_result, MPFR_RNDN);
  /* op_result = op_result * op_result */
  mpfr_sqr(op_result, op_result, MPFR_RNDN);
  /* op_result = 2 * op_result */
  mpfr_mul_ui(op_result, op_result, 2, MPFR_RNDN);
  mpfr_floor(op_result, op_result);
  mpfr_out_str(stdout, 10, 0, op_result, MPFR_RNDN);
  putchar('\n');

  mpfr_get_z(op_limit, op_result, MPFR_RNDN);

  debug_gmp_printf("floor( 2(loge(n)^2) ) : %Zd,  testNum - 1 : %Zd\n", op_limit, testNum_1);
  /* op_limit = min(op_limit, testNum) */
  if (mpz_cmp(op_limit, testNum_1) > 0)
  {
    mpz_set(op_limit, testNum_1);
  }
  debug_gmp_printf("limit number : %Zd\n", op_limit);
  return;
}

/******************************************************
 * Miller Test 
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

  mpz_t op_s, op_t, op_u, op_limit;
  mpz_init(op_s);
  mpz_init(op_t);
  mpz_init(op_u);
  mpz_init(op_limit);

  /* op_s = 0 */
  mpz_set_ui(op_s, 0);
  /* op_t = testNum -1 */
  mpz_sub_ui(op_t, testNum, 1);
  mpz_sub_ui(op_u, testNum, 1);

  GetLimitNum(op_limit, testNum, op_u);

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

  bool isLoop = true;
  isPrime = true;
#pragma omp parallel shared(isLoop, isPrime) firstprivate(testNum)
  {
    mpz_t op_private_a, op_private_s, op_private_t, op_private_u, op_private_limit, private_result, private_iter;
    mpz_init(op_private_a);
    mpz_init(op_private_s);
    mpz_init(op_private_t);
    mpz_init(op_private_u);
    mpz_init(op_private_limit);
    mpz_init(private_result);
    mpz_init(private_iter);
    mpz_set(op_private_s, op_s);
    mpz_set(op_private_t, op_t);
    mpz_set(op_private_u, op_u);
    mpz_set(op_private_limit, op_limit);

    unsigned int numThread, myID;
    numThread = omp_get_num_threads();
    myID = omp_get_thread_num();

    mpz_cdiv_q_ui(private_iter, op_private_limit, numThread);
    mpz_mul_ui(op_private_a, private_iter, myID);

    if (mpz_cmp_ui(op_private_a, 0) == 0)
      mpz_set_ui(op_private_a, 1);

#pragma omp critical(init_op_a)
    {
      debug_printf("thread num : %d\n", omp_get_thread_num());
      debug_gmp_printf("%Zd\n", op_private_a);
      debug_gmp_printf("%Zd\n", op_private_s);
      debug_gmp_printf("%Zd\n", op_private_t);
      debug_gmp_printf("%Zd\n", op_private_limit);
      debug_gmp_printf("\n");
    }

    if (mpz_cmp(op_private_a, op_private_limit) > 0)
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
#pragma omp critical
      debug_gmp_printf("Finished : %Zd\n", op_private_a);
      isLoop = false;
      isPrime = false;
      break;

    end_of_loop:
      mpz_add_ui(op_private_a, op_private_a, 1);
      mpz_sub_ui(private_iter, private_iter, 1);
      if (mpz_cmp_ui(private_iter, 0) == 0)
        isLoop = false;
    }
    mpz_clear(op_private_a);
    mpz_clear(op_private_s);
    mpz_clear(op_private_t);
    mpz_clear(op_private_u);
    mpz_clear(op_private_limit);
    mpz_clear(private_result);
    mpz_clear(private_iter);
#pragma omp barrier
  }

  mpz_clear(op_s);
  mpz_clear(op_t);
  mpz_clear(op_limit);
end_of_func:
  mpz_clear(result);
  return isPrime;
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
  isPrime = MillerTest(testNum);
  t = omp_get_wtime() - t;
  ShowResult(testNum, isPrime);

  printf("Process Time : %f\n", t);

  mpz_clear(testNum);

  return 0;
}