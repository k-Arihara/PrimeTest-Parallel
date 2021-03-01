single:
	icc -Wall -O2 -qopenmp -lmpfr -lgmp -o single.out single.c

parallel:
	icc -Wall -O2 -qopenmp -lmpfr -lgmp -o parallel.out parallel.c