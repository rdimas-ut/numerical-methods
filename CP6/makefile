CXXFLAGS=-std=c++11

fem1d: fem1d.cpp fem1d.h matrix.cpp matrix.h matlabPlot.cpp matlabPlot.h
	c++ -o fem1d fem1d.cpp matrix.cpp matlabPlot.cpp

test: fem1d
	printf '\nTesting fem1d with 11pts on CP-2a \n'
	printf '1\n11\n' | ./fem1d

	printf '\nTesting fem1d with 21pts on CP-2a \n'
	printf '1\n21\n' | ./fem1d

	printf '\nTesting fem1d with 41pts on CP-2a \n'
	printf '1\n41\n' | ./fem1d

	printf '\nTesting fem1d with 11pts on CP-2b \n'
	printf '2\n11\n' | ./fem1d

	printf '\nTesting fem1d with 21pts on CP-2b \n'
	printf '2\n21\n' | ./fem1d

	printf '\nTesting fem1d with 41pts on CP-2b \n'
	printf '2\n41\n' | ./fem1d

clean:
	rm -f *.o *~ fem1d
