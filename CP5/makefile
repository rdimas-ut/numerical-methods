PROGS = fixedpt gaussElim gaussElimLT newton iterativeNL leastSquares

all: $(PROGS)

leastSquares: main.cpp matrix.cpp matrix.h gaussElim.cpp gaussElim.h leastSquares.cpp leastSquares.h
	c++ -o leastSquares main.cpp iterativeNL.cpp matrix.cpp gaussElim.cpp leastSquares.cpp

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

test: leastSquares
	printf '\nTesting CO2 Data with Linear \n'
	printf 'y\n0\n2 50\ny\ncotwo.txt\n' | ./leastSquares

	printf '\nTesting CO2 Data with Quadratic \n'
	printf 'y\n1\n2 50\ny\ncotwo.txt\n' | ./leastSquares

	printf '\nTesting CO2 Data with Cubic \n'
	printf 'y\n2\n2 50\ny\ncotwo.txt\n' | ./leastSquares

	printf '\nTesting CO2 Data with Exponential \n'
	printf 'y\n3\n2 50\ny\ncotwo.txt\n' | ./leastSquares

	printf '\nTesting CO2 Data with Power \n'
	printf 'y\n4\n2 50\ny\ncotwo.txt\n' | ./leastSquares

	printf '\nTesting HW421a \n'
	printf 'n\ny\n3 4\ny\n412aA.txt\ny\n412ab.txt\n' | ./leastSquares
	printf '\nTesting HW421b \n'
	printf 'n\ny\n3 4\ny\n412bA.txt\ny\n412bb.txt\n' | ./leastSquares

clean:
	rm -f *.o *~ $(PROGS)

