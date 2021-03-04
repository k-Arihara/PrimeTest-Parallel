single:
	icc -Wall -O2 -qopenmp -lmpfr -lgmp -o single.out single.c miller_rabin.c miller.c naive.c
	
test_naive:
	icc -Wall -O2 -qopenmp -lmpfr -lgmp -lcunit -o test_naive.out test_naive.c miller_rabin.c miller.c naive.c

test_miller:
	icc -Wall -O2 -qopenmp -lmpfr -lgmp -lcunit -o test_miller.out test_miller.c miller_rabin.c miller.c naive.c

test_miller_rabin:
	icc -Wall -O2 -qopenmp -lmpfr -lgmp -lcunit -o test_miller_rabin.out test_miller_rabin.c miller_rabin.c miller.c naive.c