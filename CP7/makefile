
CXXFLAGS=-std=c++11

heat: main_heat.cpp heat.cpp heat.h matrix.cpp matrix.h 
	c++ -o heat main_heat.cpp heat.cpp matrix.cpp 
 
#heat: main_heat.cpp heat.cpp heat.h matrix.cpp matrix.h \
      matlabPlot.cpp matlabPlot.h 
#	c++ -o heat main_heat.cpp heat.cpp matrix.cpp matlabPlot.cpp

clean:
	rm -f *.o *~ heat
