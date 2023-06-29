PROGS = fixedpt gaussElim gaussElimLT newton iterativeNL

all: $(PROGS)

iterativeNL: main.cpp iterativeNL.cpp iterativeNL.h matrix.cpp matrix.h gaussElim.cpp gaussElim.h
	c++ -o iterativeNL main.cpp iterativeNL.cpp matrix.cpp gaussElim.cpp 

fixedpt: fixedpt.cpp matrix.cpp matrix.h
	c++ -o fixedpt fixedpt.cpp matrix.cpp

gaussElim: gaussElim.cpp gaussElim.h matrix.cpp matrix.h
	c++ -o gaussElim -DTEST gaussElim.cpp matrix.cpp

gaussElimLT: gaussElim.cpp gaussElim.h matrix.cpp matrix.h
	c++ -o gaussElimLT -DLARGE_TEST gaussElim.cpp matrix.cpp

newton: newton.cpp
	c++ -o newton newton.cpp

test: iterativeNL
	printf '\nTesting Example 1 Init 0 Broyden 2\n'
	printf '1 1\n1e-8 100\n3\n0\n1\n' | ./iterativeNL

	printf '\nTesting Example 1 Init 1 Broyden 2\n'
	printf '1 1\n1e-8 100\n3\n1\n1\n' | ./iterativeNL

	printf '\nTesting Example 2 Init 0 Broyden 2\n'
	printf '1 1\n1e-8 100\n3\n0\n2\n' | ./iterativeNL

	printf '\nTesting Example 2 Init 1 Broyden 2\n'
	printf '1 1\n1e-8 100\n3\n1\n2\n' | ./iterativeNL

	printf '\nTesting Example 3 Init 0 Broyden 2\n'
	printf '1 1\n1e-8 100\n3\n0\n3\n' | ./iterativeNL

	printf '\nTesting Example 3 Init 1 Broyden 2\n'
	printf '1 1\n1e-8 100\n3\n1\n3\n' | ./iterativeNL


	printf '\nTesting Example 1 Init 0 Broyden 1\n'
	printf '1 1\n1e-8 100\n2\n0\n1\n' | ./iterativeNL

	printf '\nTesting Example 1 Init 1 Broyden 1\n'
	printf '1 1\n1e-8 100\n2\n1\n1\n' | ./iterativeNL

	printf '\nTesting Example 2 Init 0 Broyden 1\n'
	printf '1 1\n1e-8 100\n2\n0\n2\n' | ./iterativeNL

	printf '\nTesting Example 2 Init 1 Broyden 1\n'
	printf '1 1\n1e-8 100\n2\n1\n2\n' | ./iterativeNL

	printf '\nTesting Example 3 Init 0 Broyden 1\n'
	printf '1 1\n1e-8 200000\n2\n0\n3\n' | ./iterativeNL

	printf '\nTesting Example 3 Init 1 Broyden 1\n'
	printf '1 1\n1e-8 100\n2\n1\n3\n' | ./iterativeNL
clean:
	rm -f *.o *~ $(PROGS)

