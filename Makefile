EIGENPATH = /usr/local/include/eigen3
GSLPATH = /usr/local/Cellar/gsl/2.2.1/include
INCPATH = -I/$(EIGENPATH) -I$(GSLPATH)
LIBPATH = -L/usr/local/lib -lgsl
CXXFLAGS = -c -O3
CXX=g++

tightCollimated : tightCollimated.o RNG.o
	$(CXX) -Wall $(LIBPATH) tightCollimated.o RNG.o -o tightCollimated

tightCollimated.o : tightCollimated.cpp 
	$(CXX) $(CXXFLAGS) $(INCPATH) tightCollimated.cpp

RNG.o : RNG.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) RNG.cpp

clean :
	rm *.o tightCollimated