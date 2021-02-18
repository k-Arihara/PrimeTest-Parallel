main:
	icc -Wall -O2 -qopenmp -lgmp main.c

sample:
	icc -Wall -O2 -qopenmp -lgmp sample.c