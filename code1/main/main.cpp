#include <Opinion_formation.h>

#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "sys/times.h"
#include "sys/vtimes.h"


static clock_t lastCPU, lastSysCPU, lastUserCPU;
static int numProcessors;


void init(){
  FILE* file;
  struct tms timeSample;
  char line[128];
  
  
  lastCPU = times(&timeSample);
  lastSysCPU = timeSample.tms_stime;
  lastUserCPU = timeSample.tms_utime;
  

  file = fopen("/proc/cpuinfo", "r");
  numProcessors = 0;
  while(fgets(line, 128, file) != NULL){
    if (strncmp(line, "processor", 9) == 0) numProcessors++;
  }
  fclose(file);
}


double getCurrentValue(){
  struct tms timeSample;
  clock_t now;
  double percent;
  
  
  now = times(&timeSample);
  if (now <= lastCPU || timeSample.tms_stime < lastSysCPU ||
      timeSample.tms_utime < lastUserCPU){
    //Overflow detection. Just skip this value.
    percent = -1.0;
  }
  else{
    percent = (timeSample.tms_stime - lastSysCPU) +
      (timeSample.tms_utime - lastUserCPU);
    percent /= (now - lastCPU);
    percent /= numProcessors;
    percent *= 100;
  }
  lastCPU = now;
  lastSysCPU = timeSample.tms_stime;
  lastUserCPU = timeSample.tms_utime;
  
  
  return percent;
}


int main(int argc, char* argv[])  {

  vector<Network*> net_list;  
  int network_size, M, time_steps, lambda, rg_seed;
  double mean_preference, var_preference, mean_coupling, var_coupling, average_state, p_rewiring;
  string output_name, network_type;
  FILE *file_s0, *file_s1, *file_m;
  time_t starttime; 

  //MTRand* mtrand;
  // cout<<  "argc: "<< argc <<"\n" ;
    // We print argv[0] assuming it is the program name.
      
      cout<<  "usage: "<< argv[0] ;          
      cout << "Interactive mode. Enter 0 for default value\n";
      cout << "Enter the number of time steps "; cin >> time_steps;     
      cout << "Enter number of Monte Carlo steps: "; cin >> M;
      cout << "Enter number of nodes: "; cin >> network_size;      
      cout << "Enter netwerk type: 'fully_connected', 'random', 'small_world', 'square_lattice' or 'ring'"; cin >> network_type;
      cout << "Enter average number of edges: "; cin >> lambda;                 
      cout << "Enter initial average state (in [0,1] "; cin >> average_state;          
      cout << "Enter mean preference in [-1,1]"; cin >> mean_preference;          
      cout << "Enter preference variance in [0,1]"; cin >> var_preference;         
      cout << "Enter mean interaction strength in [0,1]"; cin >> mean_coupling;         
      cout << "Enter interaction strength variance in [0,1]"; cin >> var_coupling;           
      cout << "Enter output-file name (without .dat)"; cin >> output_name;      
      cout << "Enter seed if you want repeatable data"; cin >> rg_seed;
      if (network_type.compare("small_world")==0)
	{ 
	  cout << "Enter rewiring probability for creating small-world network with the Watts-Strogatz algorithm"; cin >> p_rewiring;
	}
      cout << endl;
      //other options (add later): network type, directed?, output_filename, etc.
 
 //   // Correct values in case of invalid input interactively or through the command line ( 0 or - ).
 // //cout << "time" << time << endl;
 if (time_steps == 0) time_steps = 30;
 if (M== 0) M = 100;
 if (network_size== 0) network_size = 400;
 if ( lambda== 0) lambda =(int) network_size/40;  
 if (average_state == 0) average_state = 0.5;
 if (var_preference == 0) var_preference = 0.236;    
 if (mean_coupling == 0) mean_coupling = 0.5; 
 if (var_coupling ==0)var_coupling=0.167;   
 if (rg_seed !=0) cout <<" fdv " << endl; // mtrand->seed(rg_seed);     //else: random seed
 if (strlen(output_name.c_str())<2) {
   char temp[1024];
   sprintf(temp,"Sim_N_%d_M_%d",(int) network_size, (int) M);
   output_name = string(temp);
 }
 if(p_rewiring ==0) p_rewiring = 0.1;

 //mtrand->seed(42); geeft segmentation fault
      //mtrand->seed(42);          //set seed if you want repeatable data



 cout << "Time = " << time_steps << " steps" << endl;
 cout << "M = " << M << endl;
 cout << "Network size = "<< network_size<< " nodes" << endl;

/*!
 * constructing and initializing the networks..
 */


 try
   {
     vec lockin_states_0 = zeros<vec>(time_steps+1);  // (don't forget to make space for the initial state t=0)
     vec lockin_states_1 = zeros<vec>(time_steps+1); 
     vec mixed_states = zeros<vec>(time_steps+1); 
     int state_1_counter =0;
     int state_0_counter =0;
     int mixed_counter =0;
          
     starttime = time(NULL); //now */
     init();
                                                        
	  //	  Network net("consumers", Network::Undirected);
       Network* net = new Network("consumers", Network::Undirected); // Construct Networks

       if (network_type.compare("random")==0)                           //evaluatie input_string binnen  lus over MC-steps niet-ideaal
	 {
	   //	   net->erdos_renyi(lambda);      
	   net->populate(network_size);
	   net->fast_random_graph(lambda); //uses erdos-renyi-algorithm except for sparse networks (edges/nodes < 1 %).
	   net->get_deg_series();
	   //for (int x = 0; x != example.size(); ++x)
	   //     cout << example[x] << "- subscripting" << endl;
	   cout << "mean deg = " << net->mean_deg();
	 }
       else if (network_type.compare("fully_connected")==0)
	 {
	   net->populate(network_size);
	   net->all_to_all_coupling();
	 }
       else if (network_type.compare("small_world")==0)
	 {
	   net->small_world(network_size, lambda, p_rewiring);      //p_rewiring=  , p_rewiring=1 -> erdos random network, p_rewiring=0 -> lattice-model
	   net->get_deg_series();
	   //for (int x = 0; x != example.size(); ++x)
	   //     cout << example[x] << "- subscripting" << endl;
	   cout << "mean deg = " << net->mean_deg();
	 }
       else if(network_type.compare("square_lattice")==0)
	 {
	   int number_of_rows= (int) round(sqrt(network_size));
	   cout << number_of_rows << endl;
	   int number_of_columns= number_of_rows;
	   int new_network_size =  number_of_rows*number_of_columns;
	   fprintf(stdout,"used %d nodes instead of %d nodes to generate the square lattice \n", new_network_size, network_size); 
	   //	   net->square_lattice(number_of_rows, number_of_rows, 1); //no diagonal interactions chosen: lambda =4
	   net->square_lattice(number_of_rows, number_of_columns, 1);
	 }
       else if(network_type.compare("ring")==0)
	 {
	   net->ring_lattice(network_size, lambda);
	   cout << "mean deg = " << net->mean_deg();
	 }
       else if(network_type.compare("poisson")==0)
	 {
	   net->populate(network_size);
	   net->rand_connect_poisson(lambda);
	 }
       else { cerr << "No valid network type!" << endl;}

     for(int j= 0; j < M; j++)
       {       
	 net->initialize(mean_coupling,var_coupling, mean_preference, var_preference,average_state); 
	 Opinion_formation sim(net); 
	 sim.run_simulation(time_steps);
       //vector<double> coarse_states = sim.get_coarse_states();
       vec coarse_states = conv_to< vec >::from(sim.get_coarse_states());  // Conversion from std::vector to Armadillo vector
       if ( coarse_states[time_steps] > 0.7 )  { lockin_states_1 += coarse_states; state_1_counter++;   }       //0.7 is nog willekeurig gekozen, kies grenswaarde later ifv aantal nodes bijv.?
       else if (coarse_states[time_steps] < 0.3 ) { lockin_states_0 += coarse_states; state_0_counter++; }
       else { mixed_states += coarse_states; mixed_counter++;} 
       }
     cout << "mixed states = " << mixed_counter << endl;
     cout << "lockin_states product 1 = " << state_1_counter << endl;
     cout << "lockin_states product 0 = " << state_0_counter << endl;

     if (state_0_counter>0) 
       {
	 char temp[4096];
	 sprintf(temp, "%s%s%s%s", "data/", output_name.c_str(),"_state0", ".dat");
	 file_s0 = fopen(temp, "w");
	 if (file_s0 == 0) { cerr << "Error: can't open file \n" << output_name << endl; }
	 for (int t=0; t < time_steps+1; t++) {
	   fprintf(file_s0, "%d %f \n", t, lockin_states_0[t]/(state_0_counter)  );
	 }
	  }
     if (state_1_counter>0) 
       { 
	 char temp[4096];
	 sprintf(temp, "%s%s%s%s", "data/", output_name.c_str(),"_state1", ".dat");
	 file_s1 = fopen(temp, "w");
	 if (file_s1 == 0) { cerr << "Error: can't open file \n" << output_name << endl; }
	 for (int t=0; t < time_steps+1; t++) {
	   fprintf(file_s1, "%d %f  \n", t, lockin_states_1[t]/(state_1_counter)  );
	 }
       }
     if (mixed_counter>0) 
       { 
	 char temp[4096];
	 sprintf(temp, "%s%s%s%s", "data/", output_name.c_str(),"_mixed", ".dat");
	 file_m = fopen(temp, "w");
	 if (file_m == 0) { cerr << "Error: can't open file \n" << output_name << endl; }
	 for (int t=0; t < time_steps+1; t++) {
	   fprintf(file_m, "%d %f  \n", t, mixed_states[t]/(mixed_counter)  );
	 }
       }
     
 
	int64_t seconds = (int64_t) difftime(time(NULL),starttime);
	int64_t days = seconds / 86400; seconds -= days * 86400;
	int64_t hours = seconds / 3600; seconds -= hours * 3600;
	int64_t minutes = seconds / 60; seconds -= minutes * 60;
	fprintf(stdout,"Simulation time: %d days %d hours %d minutes %d seconds\n",(int)days,(int)hours,(int)minutes,(int)seconds);
	fprintf(stdout,"CPU currently used by current process: %f % \n ", getCurrentValue());
	
	return 0;
   }
 catch (const std::bad_alloc& ex)
   {   //Observe that the local variable `data` no longer exists here.
     cerr << "Oops. Looks like you need to use a 64 bit system or get a bigger hard disk for that calculation!"<< endl;
     return -1;
   }
 



       	//   Opinion_formation sim(net_list[5]);   
	//   net.write_edgelist("output_edgelist.csv"); 
	// SimulationTime(stdout);
	//	Ensemble ensemble(&net);
	//	ensemble.ensemble_average(7);

}





