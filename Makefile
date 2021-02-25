single:
	icc -Wall -O2 -qopenmp -lgmp single.c

parallel:
	icc -Wall -O2 -qopenmp -lgmp parallel.c