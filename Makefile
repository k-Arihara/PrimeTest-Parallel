single:
	icc -Wall -O2 -qopenmp -lgmp -o single.out single.c

parallel:
	icc -Wall -O2 -qopenmp -lgmp -o parallel.out parallel.c
