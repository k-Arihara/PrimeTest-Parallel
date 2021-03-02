parallel:
	icc -Wall -O2 -qopenmp -lmpfr -lgmp -o parallel.out parallel.c miller.c

test_miller:
	icc -Wall -O2 -qopenmp -lmpfr -lgmp -lcunit -o test_miller.out test_miller.c miller.c