#include "single.h"

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

bool SingleMillerTest(mpz_t testNum)
{
  debug_printf("Test : Miller Test\n");

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
  /* op_u = testNum -1 */
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

  mpz_t op_a, op_limit;
  mpz_init(op_a);
  mpz_init(op_limit);
  mpz_set_ui(op_a, 2);

  GetLimitNum(op_limit, testNum, op_u);

  isPrime = true;
  while (true)
  {
    /* Felmat Test */
    /* result = a^t % testNum */
    mpz_powm(result, op_a, op_t, testNum);
    /* result == 1 */
    if (mpz_cmp_ui(result, 1) == 0)
      goto end_of_loop;

    /* while(i < s) */
    for (unsigned int i = 0; mpz_cmp_ui(op_s, (unsigned long)i) > 0; i++)
    {
      /* result = 2^i */
      mpz_ui_pow_ui(result, 2, (unsigned long)i);
      /* result = result * t */
      mpz_mul(result, result, op_t);
      /* result = a^result % testNum */
      mpz_powm(result, op_a, result, testNum);
      /* a^(2^i * t) % testNum == testNum - 1 */
      if (mpz_cmp(result, op_u) == 0)
        goto end_of_loop;
    }
    isPrime = false;
    break;
  end_of_loop:
    mpz_add_ui(op_a, op_a, 1);
    if (mpz_cmp(op_a, op_limit) > 0)
      break;
  }

  mpz_clear(op_a);
  mpz_clear(op_s);
  mpz_clear(op_t);
  mpz_clear(op_u);
end_of_func:
  mpz_clear(result);
  return isPrime;
}