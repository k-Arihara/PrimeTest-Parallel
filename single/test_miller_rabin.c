#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>
#include "single.h"

void Test001(void)
{
  mpz_t op_test;
  mpz_init(op_test);
  mpz_set_ui(op_test, 2);
  CU_ASSERT_TRUE(MillerRabinTest(op_test));
  mpz_set_ui(op_test, 3);
  CU_ASSERT_TRUE(MillerRabinTest(op_test));
  mpz_set_ui(op_test, 5);
  CU_ASSERT_TRUE(MillerRabinTest(op_test));
  mpz_clear(op_test);
}

void Test002(void)
{
  mpz_t op_test;
  mpz_init(op_test);
  mpz_set_ui(op_test, 6700417);
  CU_ASSERT_TRUE(MillerRabinTest(op_test));
  mpz_set_ui(op_test, 2147483647);
  CU_ASSERT_TRUE(MillerRabinTest(op_test));
  mpz_clear(op_test);
}

void Test003(void)
{
  mpz_t op_test;
  mpz_init(op_test);
  mpz_set_ui(op_test, 4);
  CU_ASSERT_FALSE(MillerRabinTest(op_test));
  mpz_set_ui(op_test, 8187);
  CU_ASSERT_FALSE(MillerRabinTest(op_test));
  mpz_clear(op_test);
}

int main(void)
{
  CU_pSuite suite;

  CU_initialize_registry();
  suite = CU_add_suite("Miller-Rabin Test", NULL, NULL);
  CU_ADD_TEST(suite, Test001);
  CU_ADD_TEST(suite, Test002);
  CU_ADD_TEST(suite, Test003);
  CU_basic_run_tests();
  CU_cleanup_registry();

  return 0;
}