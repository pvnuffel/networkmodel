#include "MersenneTwister.h"
#include <iostream>


using namespace std;
MTRand mtrand;

double method (MTRand rand)
{

  cout << rand.rand() << endl;
}

double method2 ()
{

  cout << mtrand.rand() << endl;
}


int main (int argc, char** argv) {

  MTRand mtrand;
  mtrand.seed(15);
  cout <<  mtrand.rand() << endl;
  cout <<  mtrand.rand() << endl;
  mtrand.seed(15);
  cout <<  mtrand.rand() << endl << endl;

  for (int i=0; i< 3; i++)
    {
       cout << mtrand.rand() << endl;
    }


  for (int i=0; i< 3; i++)
    {
      method(mtrand);
    }

  for (int i=0; i< 3; i++)
    {
      method2();
    }

 mtrand.seed(15);
  for (int i=0; i< 3; i++)
    {
      method2();
    }

 mtrand.seed(15);

  for (int i=0; i< 3; i++)
    {
      method2();
    }


  // MTRand* mtrand;
  // mtrand->seed(15);
  // cout <<  mtrand->rand() << endl; 
  // cout <<  mtrand->rand() << endl;  
  //gives segmenation fault
  
    return 0;
}


