#include "MersenneTwister.h"
#include <iostream>
#include <fstream>

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
  // mtrand.seed(15);
  //  cout <<  mtrand.rand() << endl;
  //cout <<  mtrand.rand() << endl;

  MTRand::uint32 randState[ MTRand::SAVE ];
  mtrand.save( randState ); 

  cout <<" Print two random numbers using the generator after the state is saved" << endl;
  cout<<mtrand.rand()<<endl;
  cout<<mtrand.rand() <<endl;
	// We might want to save the state of the generator at
	// an arbitrary point after seeding so a sequence
	// could be replicated	
	// The array must be of type uint32 and length SAVE.
  mtrand.load( randState );
 
  cout << " Print the number from the generator after it is loaded and verify it is restarted at the number already printed " << endl;
  cout << mtrand.rand() << endl;
  cout<<mtrand.rand() <<endl;




  ofstream stateOut( "state.data" );
  if( stateOut )
	{
	  stateOut << mtrand;
	  stateOut.close();
	}
	
  cout <<" Print two random numbers using the generator after the state is saved" << endl;
  cout << mtrand.rand() << endl;
  cout<<mtrand.rand() <<endl;
  

  ifstream stateIn( "state.data" );
  if( stateIn )
    {
      stateIn >> mtrand;
      stateIn.close();
    }
  cout << " Print the number from the generator after it is loaded and verify it is restarted at the number already printed " << endl;
  cout << mtrand.rand() << endl;
  cout<<mtrand.rand() <<endl;

  
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


