#include <stdio.h>
#include <stdbool.h>
#include <gmp.h>
#include <mpfr.h>
#include <omp.h>

#ifndef _PARALLEL_H_
#define _PARALLEL_H_

bool ParallelMillerTest(mpz_t testNum);

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