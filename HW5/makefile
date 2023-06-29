all: gaussElim gaussElimLT

gaussElim: gaussElim.C gaussElim.h matrix.C matrix.h
	c++ -o gaussElim -DTEST gaussElim.C matrix.C

gaussElimLT: gaussElim.C gaussElim.h matrix.C matrix.h
	c++ -o gaussElimLT -DLARGE_TEST gaussElim.C matrix.C

clean:
	rm -f gaussElim gaussElimLT


