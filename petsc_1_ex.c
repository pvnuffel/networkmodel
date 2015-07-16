// EX1.CPP
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>

Vec Arr2Vec(double *arr2, int SIZE);

// MAIN FUNCTION
int main(int argc,char **argv)
{
    // Initialize PetSc
    PetscInitialize(&argc,&argv,(char*)0,"Testing a program!");

    // Initialize parameters
    int SIZE = 3;
    PetscErrorCode ierr;

    // **** Create a regular arary and set it with random numbers
    double  * arr2;
    arr2 = new double [SIZE];

    arr2[0] = 0.1;
    arr2[1] = 0.4;
    arr2[2] = 0.2;

    // Convert regular arary to PETSc vector [Note that this must do the same effect as the two-step process of conversion from regular array to PETSc arary and that to PETSc vector as listed above]
    Vec x = Arr2Vec(arr2, SIZE);

    printf("Reciprocal Vector : \n"); VecReciprocal(x);
    VecView(x,PETSC_VIEWER_STDOUT_WORLD);

    //Cleanup
    ierr = VecDestroy(&x);
    CHKERRQ(ierr);
    PetscFinalize();

    return 0;
}

Vec Arr2Vec(double *arr2, int SIZE)
{
  PetscScalar *array1;
  PetscMalloc(SIZE*sizeof(PetscScalar),&array1);

  for(int i=0;i<SIZE;i++)
    array1[i]=arr2[i];

    // Setup vector
  Vec x;
  VecCreate(PETSC_COMM_WORLD,&x);
  VecSetSizes(x,PETSC_DECIDE,SIZE);
  VecSetFromOptions(x);

  // Place PetSc array as Vector
  VecPlaceArray(x,array1);

  return x;

}
