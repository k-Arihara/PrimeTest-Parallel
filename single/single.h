#include <stdio.h>
#include <stdbool.h>
#include <gmp.h>
#include <mpfr.h>

#ifndef _SINGLE_H_
#define _SINGLE_H_

bool SingleMillerTest(mpz_t testNum);
bool MillerRabinTest(mpz_t testNum);
bool NaivePrimeTest(mpz_t testNum);

// #define DEBUG

/* printf for debuggin */
#ifdef DEBUG
#define debug_printf printf
#define debug_gmp_printf gmp_printf
#else
#define debug_printf 1 ? (void)0 : printf
#define debug_gmp_printf 1 ? (void)0 : gmp_printf
#endif

#endif