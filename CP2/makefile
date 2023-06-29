CXXFLAGS=-std=c++11

iterativeLA: main.cpp iterativeLA.cpp iterativeLA.h matrix.cpp matrix.h matlabPlot.cpp matlabPlot.h
	c++ -o ${CXXFLAGS} -o iterativeLA main.cpp iterativeLA.cpp matrix.cpp matlabPlot.cpp

test1_6: iterativeLA
	printf '\nTesting 1:Jacobi\n'
	printf '6\ny\nA1_6.txt\ny\nb1_6.txt\n50 1e-6\nJacobi' | ./iterativeLA
	printf '\nTesting 1:Gauss-Seidel\n'
	printf '6\ny\nA1_6.txt\ny\nb1_6.txt\n50 1e-6\nGauss-Seidel' | ./iterativeLA
	printf '\nTesting 1:SOR_1.25\n'
	printf '6\ny\nA1_6.txt\ny\nb1_6.txt\n50 1e-6\nSOR\n1.25' | ./iterativeLA
	printf '\nTesting 1:CCG\n'
	printf '6\ny\nA1_6.txt\ny\nb1_6.txt\n50 1e-6\nCCG\nn' | ./iterativeLA
	printf '\nTesting 1:CCG_Jacobi\n'
	printf '6\ny\nA1_6.txt\ny\nb1_6.txt\n50 1e-6\nCCG\ny\nJacobi' | ./iterativeLA
	printf '\nTesting 1:CCG_Gauss-Seidel\n'
	printf '6\ny\nA1_6.txt\ny\nb1_6.txt\n50 1e-6\nCCG\ny\nGauss-Seidel' | ./iterativeLA
	printf '\nTesting 1:CCG_SOR\n'
	printf '6\ny\nA1_6.txt\ny\nb1_6.txt\n50 1e-6\nCCG\ny\nSOR\n1.25' | ./iterativeLA

test1_12: iterativeLA
	printf '\nTesting 1:Jacobi\n'
	printf '12\ny\nA1_12.txt\ny\nb1_12.txt\n1000 1e-6\nJacobi' | ./iterativeLA
	printf '\nTesting 1:Gauss-Seidel\n'
	printf '12\ny\nA1_12.txt\ny\nb1_12.txt\n1000 1e-6\nGauss-Seidel' | ./iterativeLA
	printf '\nTesting 1:SOR_1.25\n'
	printf '12\ny\nA1_12.txt\ny\nb1_12.txt\n1000 1e-6\nSOR\n1.25' | ./iterativeLA
	printf '\nTesting 1:CCG\n'
	printf '12\ny\nA1_12.txt\ny\nb1_12.txt\n1000 1e-6\nCCG\nn' | ./iterativeLA
	printf '\nTesting 1:CCG_Jacobi\n'
	printf '12\ny\nA1_12.txt\ny\nb1_12.txt\n1000 1e-6\nCCG\ny\nJacobi' | ./iterativeLA
	printf '\nTesting 1:CCG_Gauss-Seidel\n'
	printf '12\ny\nA1_12.txt\ny\nb1_12.txt\n1000 1e-6\nCCG\ny\nGauss-Seidel' | ./iterativeLA
	printf '\nTesting 1:CCG_SOR\n'
	printf '12\ny\nA1_12.txt\ny\nb1_12.txt\n1000 1e-6\nCCG\ny\nSOR\n1.25' | ./iterativeLA

test2: iterativeLA
	printf '\nTesting 2:Jacobi\n'
	printf '100\ny\nA2.txt\ny\nb2.txt\n2000 1e-10\nJacobi' | ./iterativeLA
	printf '\nTesting 2:Gauss-Seidel\n'
	printf '100\ny\nA2.txt\ny\nb2.txt\n2000 1e-10\nGauss-Seidel' | ./iterativeLA
	printf '\nTesting 2:SOR_1.1\n'
	printf '100\ny\nA2.txt\ny\nb2.txt\n2000 1e-10\nSOR\n1.1' | ./iterativeLA
	printf '\nTesting 2:SOR_1.25\n'
	printf '100\ny\nA2.txt\ny\nb2.txt\n2000 1e-10\nSOR\n1.25' | ./iterativeLA
	printf '\nTesting 2:SOR_1.67\n'
	printf '100\ny\nA2.txt\ny\nb2.txt\n2000 1e-10\nSOR\n1.67' | ./iterativeLA
	printf '\nTesting 2:CCG\n'
	printf '100\ny\nA2.txt\ny\nb2.txt\n2000 1e-10\nCCG\nn' | ./iterativeLA
	printf '\nTesting 2:CCG_Jacobi\n'
	printf '100\ny\nA2.txt\ny\nb2.txt\n2000 1e-10\nCCG\ny\nJacobi' | ./iterativeLA
	printf '\nTesting 2:CCG_Gauss-Seidel\n'
	printf '100\ny\nA2.txt\ny\nb2.txt\n2000 1e-10\nCCG\ny\nGauss-Seidel' | ./iterativeLA
	printf '\nTesting 2:CCG_SOR_1.1\n'
	printf '100\ny\nA2.txt\ny\nb2.txt\n2000 1e-10-\nCCG\ny\nSOR\n1.1' | ./iterativeLA
	printf '\nTesting 2:CCG_SOR_1.25\n'
	printf '100\ny\nA2.txt\ny\nb2.txt\n2000 1e-10-\nCCG\ny\nSOR\n1.25' | ./iterativeLA
	printf '\nTesting 2:CCG_SOR_1.67\n'
	printf '100\ny\nA2.txt\ny\nb2.txt\n2000 1e-10-\nCCG\ny\nSOR\n1.67' | ./iterativeLA

test3: iterativeLA
	printf '\nTesting 2:CCG\n'
	printf '100\ny\nA3.txt\ny\nb3.txt\n5000 1e-16\nCCG\nn' | ./iterativeLA
	printf '\nTesting 2:CCG_Jacobi\n'
	printf '100\ny\nA3.txt\ny\nb3.txt\n5000 1e-16\nCCG\ny\nJacobi' | ./iterativeLA
	printf '\nTesting 2:CCG_Gauss-Seidel\n'
	printf '100\ny\nA3.txt\ny\nb3.txt\n5000 1e-16\nCCG\ny\nGauss-Seidel' | ./iterativeLA

clean:
	rm -f *.o *~ iterativeLA
