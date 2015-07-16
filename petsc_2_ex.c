#include <petsc.h>

static char help[] = "PETSc Hello world program.\n\n";

int main(int argc, char **argv)
{
  PetscMPIInt rank;
  PetscErrorCode ierr;
  PetscInt myint = -1;
  
  PetscInitialize(&argc, &argv, (char*) 0, help);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Synchronized Hello from %d\n", rank);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nHello from ZERO\n\n");CHKERRQ(ierr);
  
  ierr = PetscOptionsGetInt(PETSC_NULL,"-myint",&myint,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "myint = %d\n", myint);CHKERRQ(ierr);
  
  PetscFinalize();
  return 0;
}
