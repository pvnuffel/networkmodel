#include <Opinion_formation.h>

#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "sys/times.h"
#include "sys/vtimes.h"
//#include <math.h>    

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




Network* construct_network(string network_type,int network_size, int lambda,  double p_rewiring)
{

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
	   net->small_world(network_size, lambda, p_rewiring);      // p_rewiring=1 -> erdos random network, p_rewiring=0 -> lattice-model
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
       
 
       return net;

}


double analytical_coarse_step(double U_guess, double mean_coupling, double mean_preference, double var_preference)
{
  return 0.5*erfc( (mean_coupling*(1-2*U_guess)/(1-mean_coupling) - mean_preference ) /(sqrt(2)*var_preference));               //Avitabile14 eq 3.6
}

double anal_derivative_coarse_step(double U_guess, double mean_coupling, double mean_preference, double var_preference)
{
  const double pi = 3.14159265358979323846;
  return   exp( - pow( (mean_coupling*(1-2*U_guess)/(1-mean_coupling) - mean_preference ) /(sqrt(2)*var_preference) , 2) ) *  mean_coupling/((1-mean_coupling)*var_preference*sqrt(2*pi));         //derivative of erfc(U)*0.5
}

double newton_raphson_analytical(double U_guess, double mean_coupling, double mean_preference, double var_preference)
{ 
  double tol=1e-6;
  double err=tol+1;
  double it = 0;
  double maxit=50;
  double U_coarse;
  while(  err > tol && it < maxit )
    {
      U_coarse=  U_guess-  ( analytical_coarse_step(U_guess, mean_coupling, mean_preference, var_preference)  - U_guess)/ (anal_derivative_coarse_step(U_guess, mean_coupling, mean_preference, var_preference) -1);  //newton-step  U_n+1 = U_n - f(U_n)/f'(u_n)
      err = fabs( U_coarse - U_guess );
      U_guess= U_coarse;
      it++;
      cout << "it = " << it <<  "   u_coarse= " <<  U_coarse << endl;
    }
  if( err <= tol )
    cout << "The Newton-Raphson root U_coarse = " << U_coarse <<  " after " << it <<  " iterations" << endl;
  else
    cerr << "Error, no convergence in Newton Raphson \n";
}

double coarse_step(double U_guess, double mean_coupling, double var_coupling, double mean_preference, double var_preference, Network* network, int M, int time_horizon)
{
  double U_coarse =0;
       for(int j= 0; j < M; j++)
	 {   
	   network->initialize(mean_coupling,var_coupling, mean_preference, var_preference, U_guess);        //=lifting step
	   Opinion_formation sim(network); 
	   sim.run_simulation(time_horizon);
	   U_coarse += network->get_coarse_state();   //what with realizations cancelling each other out?
	 }
    return   U_coarse /= M;    
}

double coarse_step_weighted(double U_guess, double mean_coupling, double var_coupling, double mean_preference, double var_preference, Network* network, int M, int time_horizon)
{
  double U_coarse =0;
       for(int j= 0; j < M; j++)
	 {   
	   network->initialize(mean_coupling,var_coupling, mean_preference, var_preference, U_guess);        //=lifting step
	   Opinion_formation sim(network); 
	   sim.run_simulation(time_horizon);
	   U_coarse += network->get_coarse_state();   //what with realizations cancelling each other out?
	 }
    return   U_coarse /= M;    
}


double derivative_coarse_step(double U_guess, double mean_coupling, double var_coupling, double mean_preference, double var_preference, Network* network, int M, int time_horizon) //nog implementeren
{
  double epsilon = 1e-5;
  return (coarse_step(U_guess + epsilon, mean_couling,var_coupling, mean_reference, var_preference, network, M, time_horizon)  - coarse_step(U_guess, mean_couling,var_coupling, mean_reference, var_preference, network, M, time_horizon ) )/epsilon;
}




 double newton_raphson(double U_guess, double mean_coupling, double var_coupling,  double mean_preference, double var_preference, Network* network, int M, int time_horizon)
{ 
  double tol=1e-6;
  double err=tol+1;
  double it = 0;
  double maxit=50;
  double U_coarse;
  while(  err > tol && it < maxit )
    {
      U_coarse=  U_guess-  ( coarse_step(U_guess, mean_coupling, var_coupling, mean_preference, var_preference, network, M, time_horizon)  - U_guess)/ (derivative_coarse_step(U_guess, mean_coupling,var_coupling,  mean_preference, var_preference, network, M, time_horizon) -1);  //newton-step  U_n+1 = U_n - f(U_n)/f'(u_n)
      cout << "it = " << it <<  "   u_coarse= " <<  U_coarse << endl;
      err = fabs( U_coarse - U_guess );
      U_guess= U_coarse;
      it++;
    }
  if( err <= tol )
    cout << "The Newton-Raphson root U_coarse = " << U_coarse <<  " after " << it <<  " iterations" << endl;
  else
    cerr << "Error, no convergence in Newton Raphson \n";
}




// double coarse_step_vector(double U_guess, double mean_coupling, double mean_preference, double var_preference, Network* network, int M, int time_horizon)
// {
//        for(int j= 0; j < M; j++)
// 	 {   
// 	   network->lift(mean_coupling,var_coupling, mean_preference, var_preference, U_guess);
// 	   Opinion_formation sim(network); 
// 	   sim.run_simulation(time_horizon);
// 	   network->get_node_states();
// 	   vec node_states_j = conv_to< vec >::from(net_j->get_node_states());  // Conversion from std::vector to Armadillo vector  (dimension: time_steps)
// 	   node_states +=node_states_j;
// 	       //net_list.push_back(net_j);
// 	 }
//        node_states /= M;    
//        cout << "State N=10 after " << node_states[10] << endl;
//        cout << "State N=11 after " << node_states[11] << endl;

// }

double fixpoint_iteration_analytical(double U_guess, double mean_coupling, double mean_preference, double var_preference)
{ 
  double tol=1e-6;
  double err=tol+1;
  double it = 0;
  double maxit=50;
  double U_coarse;
  
  cout << "it = " << it <<  "   u_coarse= " <<  U_coarse << endl;
  while( err> tol  && it < maxit) 
    {
      U_coarse = analytical_coarse_step(U_guess, mean_coupling, mean_preference, var_preference); 
      err = fabs( U_coarse - U_guess );
      U_guess= U_coarse;
      it++;
      cout << "it = " << it <<  "   u_coarse= " <<  U_coarse << endl;
    }
  if( err <= tol )
    cout << "The fixed point U_coarse = " << U_coarse <<  " after " << it <<  " iterations" << endl;
  else
    cerr << "Error, no convergence in fixed point iteration  \n";
}   

double fixpoint_iteration(double U_guess, double mean_coupling, double var_coupling, double mean_preference, double var_preference,Network* network, int M, int time_horizon)
{ 
  double tol=1e-6;
  double err=tol+1;
  double it = 0;
  double maxit=50;
  double U_coarse;
  cout << "it = " << it <<  "   u_coarse= " <<  U_coarse << endl;
  while( err> tol  && it < maxit) 
    {
      U_coarse = coarse_step(U_guess, mean_coupling, var_coupling, mean_preference, var_preference, network, M, time_horizon); 
      err = fabs( U_coarse - U_guess );
      U_guess= U_coarse;
      it++;
      cout << "it = " << it <<  "   u_coarse= " <<  U_coarse << endl;
    }
  if( err <= tol )
    cout << "The fixed point U_coarse = " << U_coarse <<  " after " << it <<  " iterations" << endl;
  else
    cerr << "Error, no convergence in fixed point iteration  \n";
}   





int main(int argc, char* argv[])  {

  vector<Network*> net_list;  
  int network_size, M, time_steps, degree, rg_seed;
  double mean_preference, var_preference, mean_coupling, var_coupling, average_state, p_rewiring;
  string output_name, network_type;
  FILE *file_s0, *file_s1, *file_m, *file_mean, *file_var, *file;
  time_t starttime; 
  bool separate_coarse_states = false;

  // cout << "armadillotest"  << endl;
  // // vec y = zeros<vec>(5);  
  // //mat A = randu<mat>(10,10);
  // vec z;
  // z << 2 << 4 << 9 <<1; 
  // vec q = square(z);
  // cout << q[0] << endl;

  //MTRand* mtrand;
  // cout<<  "argc: "<< argc <<"\n" ;
    // We print argv[0] assuming it is the program name.
      
      cout<<  "usage: "<< argv[0] ;          
      cout << "Interactive mode. Enter 0 for default value\n";
      cout << "Enter the number of time steps "; cin >> time_steps;     
      cout << "Enter number of Monte Carlo steps: "; cin >> M;
      cout << "Enter number of nodes: "; cin >> network_size;      
      cout << "Enter netwerk type: 'fully_connected', 'random', 'small_world', 'square_lattice' or 'ring'"; cin >> network_type;
      cout << "Enter average number of edges: "; cin >> degree;                 
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
 if ( degree== 0) degree =(int) network_size/40;  
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
     vec  averages= zeros<vec>(time_steps+1); 
     vec averages2 = zeros<vec>(time_steps+1); 
     vec variances  = zeros<vec>(time_steps+1); 
 
     vec node_states  = zeros<vec>(network_size);
     bool use_weights = 0;
 
     starttime = time(NULL); //now */
     init();
                           
     Network* net= construct_network(network_type, network_size, degree, p_rewiring);                                                     
                    
       for(int j= 0; j < M; j++)
       {       
	 net->initialize(mean_coupling,var_coupling, mean_preference, var_preference,average_state); 
	 Opinion_formation sim(net); 

	 sim.run_simulation(time_steps);
       //vector<double> coarse_states = sim.get_coarse_states();

	     vec coarse_states = conv_to< vec >::from(sim.get_coarse_states());  // Conversion from std::vector to Armadillo vector  (dimension: time_steps)

	     if (separate_coarse_states) //set this boolean true if you want to look at time evolution (without realizations cancelling each other out)
	       {
		 if ( coarse_states[time_steps] > 0.7 )  { lockin_states_1 += coarse_states; state_1_counter++;   }       //0.7 is nog willekeurig gekozen, kies grenswaarde later ifv aantal nodes bijv.?
		 else if (coarse_states[time_steps] < 0.3 ) { lockin_states_0 += coarse_states; state_0_counter++; }
		 else { mixed_states += coarse_states; mixed_counter++;}
	       }
	     else // no separation of coarse states
	       {
		 averages += coarse_states;
		 averages2 += square(coarse_states);
	       } 	     

       }
       averages /= M;
       averages2 /=M;
       variances= averages2 - square(averages);

       // solve for equilibrium state
       double U_guess=average_state;
       fixpoint_iteration_analytical(U_guess, mean_coupling, mean_preference, var_preference);
       newton_raphson_analytical(U_guess, mean_coupling, mean_preference, var_preference);


       // LIFTING/restriction
       int time_horizon = 1; // 20;
       Network*  network =  construct_network(network_type, network_size, degree, p_rewiring);
       //   newton_raphson(U_guess, mean_coupling, var_coupling, mean_preference, var_preference, network, M, time_horizon);  
       fixpoint_iteration(U_guess, mean_coupling, var_coupling, mean_preference, var_preference, network, M, time_horizon);  
       // vector<double> U;
       // for (int i=0; i < network_size; i++)
       // 	 { U.push_back(0.4);}  //Suppose U is a N-dim testvector     zoals  U = net->get_node_states(); hoewel in het algemeen is elke toestand reëel getal is ipv (0,1)
 

       
       



       char temp[4096];
       sprintf(temp, "%s%s%s%s", "data/", output_name.c_str(),"_o", ".dat");
       file= fopen(temp, "w");
       if (file == 0) { cerr << "Error: can't open file \n" << output_name << endl; }
       for (int t=0; t < time_steps+1; t++) {
       	 fprintf(file, "%d %f %f \n", t, averages[t], sqrt(variances[t])  );
            }

     int64_t seconds = (int64_t) difftime(time(NULL),starttime);
     int64_t days = seconds / 86400; seconds -= days * 86400;
     int64_t hours = seconds / 3600; seconds -= hours * 3600;
     int64_t minutes = seconds / 60; seconds -= minutes * 60;
     cout << "Simulation " << output_name  << " done " <<  endl;
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


