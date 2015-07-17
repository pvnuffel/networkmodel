CFLAGS=--ansi --pedantic -O2 -g  

 # compiler flags:
  #  -g    adds debugging information to the executable file
  #  -Wall turns on most, but not all, compiler warnings

INCLUDE= -I../src/
LDFLAGS=  ../src/*.o

# include ${PETSC_DIR}/lib/petsc/conf/variables
# include ${PETSC_DIR}/lib/petsc/conf/rules

all:  src gsl  main    

src: 
	$(MAKE) -C ../src/

gsl:
	$(MAKE) -C ../gsl_subset/


main: main.cpp 
	g++ -std=c++0x   $(CFLAGS) main.cpp $(INCLUDE) $(LDFLAGS) -o run -larmadillo 

# myprog: myprog.o chkopts
# 	-${CLINKER} -o myprog myprog.o ${PETSC_LIB}
# 	 ${RM} myprog.o

cleanmain:	rm -f main

# runinput:
# 	-@${MPIEXEC} -n 1 run <input.txt