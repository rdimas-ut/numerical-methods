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

clean:
	rm -f *.o *~ $(PROGS)

