CFLAGS=--ansi --pedantic -O2
#CFLAGS=--ansi --pedantic -g 
INCLUDE= -I../src/
LDFLAGS=  ../src/*.o

all: src gsl  main

src: 
	$(MAKE) -C ../src/

gsl:
	$(MAKE) -C ../gsl_subset/


main: main.cpp
	g++ -std=c++0x   $(CFLAGS) main.cpp $(INCLUDE) $(LDFLAGS) -o run -larmadillo

clean:
	rm -f main
