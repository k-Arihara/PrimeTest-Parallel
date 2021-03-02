single:
	icc -Wall -O2 -qopenmp -lmpfr -lgmp -o single.out single.c miller_rabin.c miller.c naive_prime.c

parallel:
	icc -Wall -O2 -qopenmp -lmpfr -lgmp -o parallel.out parallel.c