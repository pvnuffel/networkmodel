# PRACE 2IP PETSc video courses
# by Vaclav Hapla <vaclav.hapla@vsb.cz>
#
# Needs PETSC_DIR & PETSC_ARCH environment variables to be set.
#
# Usage:
# make ex2; ./ex2
# make ex3; ./ex3
# make ex4; ./ex4
# make ex5; ./ex5
# make ex-clean

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
CLINKER=g++

ex1: petsc_1_ex.o  chkopts
	${CLINKER} -w  -o ex2 petsc_1_ex.o ${PETSC_LIB}
	${RM} -f ex1.o
	./ex1

ex2: petsc_2_ex.o  chkopts
	-${CLINKER} -o ex2 petsc_2_ex.o ${PETSC_SYS_LIB}
	${RM} -f ex2.o

ex3: petsc_3_ex.o  chkopts
	-${CLINKER} -o ex3 petsc_3_ex.o ${PETSC_VEC_LIB}
	${RM} -f ex3.o

ex4: petsc_4_ex.o  chkopts
	-${CLINKER} -o ex4 petsc_4_ex.o ${PETSC_MAT_LIB}
	${RM} -f ex4.o

ex5: petsc_5_ex.o  chkopts
	-${CLINKER} -o ex5 petsc_5_ex.o ${PETSC_KSP_LIB}
	${RM} -f ex5.o

ex-clean:

	${RM} -f ex2 ex3 ex4 ex5 *.o #PRACE 2IP PETSc video courses

runex2:
	-@${MPIEXEC} -n 1 ex2