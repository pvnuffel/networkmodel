#ifndef OPINION_FORMATION_H
#define OPINION_FORMATION_H

#include "Simulator.h"


class Opinion_formation: public Simulator
{
  // private:   FILE *file;
 protected:   double  beta=10;    // let op: non-static data member initializers only available with -std=c++11 or -std=gnu++11 [enabled by default]
  
 public:    ofstream outfile; 
  
        Opinion_formation():Simulator() {};
  // Opinion_formation(Network* net,string &filename ):Simulator(net, filename) {};
 Opinion_formation(Network* net):Simulator(net) {};
  //Opinion_formation(vector<Network*> net_list,string &filename):Simulator(net_list, filename) {};
        ~Opinion_formation() { };

       
       void time_step() {
            time++;
	    //    double ensemble_average=0;
	    // double M= net_list.size();
	    //	    for(int j= 0; j < M; j++){
      	    //vector<Node*> nodes = net_list[j]->get_nodes();
	   vector<Node*> nodes = net->get_nodes();
	      for (unsigned int i = 0; i < nodes.size(); i++) { //agents first simultaneously inspect their neighbourhoods, and then select their products.
		nodes[i]->update_utility_function();
	      } 
	      for (unsigned int i = 0; i < nodes.size(); i++) {  
		Node* inode = nodes[i];
		double delta_utility = inode->get_utility();
		//	if (i==10) cout << "utility: " << delta_utility << "\n\n";
		double p1 =  exp(beta*delta_utility)/(exp(beta*delta_utility)+exp(-beta*delta_utility));	      //  p0=1-p1;
		//	cerr << "p1: " << p1 << "\n\n";
		if (mtrand->rand() < p1) inode->set_state(1);   
		else inode->set_state(0);
	      }	 //   initial_average += net_list[j]->get_coarse_state();
	 // };
	 //	 fprintf(file, "%d %f\n", time,net->get_coarse_state()  );
	      coarse_states.push_back(net->get_coarse_state()); //saves coarse_state on each timestep as last entry in a vector
       }

       void run_simulation(int max_time)
       {
	 coarse_states.push_back(net->get_coarse_state()); 
	 while (time <max_time) {
	   time_step() ;
	    //  cerr <<"random number test output " << time << " :  " << mtrand->rand() << " \n\n";
	 }	 
       }
       
 
	 



  /* void newFile() */
  /* { */
  /* 	  char temp[4096];  */
  /* 	  if (file != NULL) fclose(file); */
  /* 	  sprintf(temp, "%s%s_%s", "data/", fname.c_str(), "dat");  */
  /* 	  file = fopen(temp, "w"); */
  /* 	  // filecount2++;  */
  /* 	  // 	      ostringstream name; */
  /* 	  // 	      name << "data/data_"   << setfill('0') << setw(4)   << filecount   << ".dat"; */
  /* } */


	

    	/* void ensemble_average(int number_of_simulations) */
	/* {     */
	/*   // open a file in write mode. */
	/*   // outfile.open("densities.dat"); */
    
	/*   // starttime = time(NULL); // now */
	/*   for (int i = 0; i < number_of_simulations; i++){ */
	/*     run_simulation(); */
	/*   } */
	/*   //  outfile.close(); */
	/*   //  SimulationTime(stdout); */
	/* } */


        /* int count_infected() { */
        /*     return state1.size(); */
        /* } */

        /* int epidemic_size() { */
        /*     return state0.size(); */
        /* } */

        void reset() {
            reset_time();
        }

        void summary() {
	  // cerr << "Network size: " << net->size();
	    // cerr << "\tCoarse state " << net->get_coarse_state() << "\n\n";
        }
};




/* class Ensemble */
/* { */


/*  protected:  */
/*   time_t starttime; */
/*   Network* net; */

/*  public:  */
/*   ofstream outfile; */
/*   Ensemble(Network* net) { this->net = net; };  */
  
/*   void SimulationTime(FILE* f)  //  Write out the duration of the whole simulation */
/*   { */
/*     int64_t seconds = (int64_t) difftime(time(NULL),starttime); */
/*     int64_t days = seconds / 86400; seconds -= days * 86400; */
/*     int64_t hours = seconds / 3600; seconds -= hours * 3600; */
/*     int64_t minutes = seconds / 60; seconds -= minutes * 60; */
/*     fprintf(f,"Simulation time: %d days %d hours %d minutes %d seconds\n",(int)days,(int)hours,(int)minutes,(int)seconds); */
/*     } */
    
/*   void ensemble_average(int number_of_simulations)            //beter niet in klasse opinion_formation omdat we steeds met de initiele isntantie van netwerk willen beginnen */
/*   {     */
/*     // open a file in write mode. */
/*     outfile.open("densities.dat"); */
    
/*     starttime = time(NULL); // now */
/*     for (int i = 0; i < number_of_simulations; i++){ */
/*       // Choose and run simulation */
/*       Opinion_formation sim(net); */
/*       sim.run_simulation(); */
   
/*       sim.summary(); */
/*       sim.reset(); */
/*     } */
/*     outfile.close(); */
/*     SimulationTime(stdout); */
/*   } */
  
/* }; */

#endif
