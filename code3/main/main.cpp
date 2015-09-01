#include <Opinion_formation.h>

#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "sys/times.h"
#include "sys/vtimes.h"
#include <iostream>
#include <fstream>


// #include <petsc.h>
// #include <petsctao.h>

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

vec calculate_weights (rowvec u_realizations, double U)
{
  int M = u_realizations.n_elem;
  // u_realizations.print();
  // cout << "M= "  << M << endl;
  mat I =  eye<mat>(M,M);
  mat C =   join_cols(u_realizations, ones<rowvec>(M));                          //u_coarse_states.insert_rows(1,ones(M) ) ;
  mat D;
  //D << 1 << 3 << 5 << endr     << 2 << 4 << 6 << endr;
  // cout << "check (1,2)=3 : " <<  D[1,2] << endl;
  //cout << "sampled check (1,2): " <<  C[1,2] << endl;
  mat A_l = join_cols(I, C);
  mat A_r =  join_cols( C.t(), zeros<mat>(2,2));
  mat A = join_rows(A_l,A_r);
  //A.print();
  // mat R2 = chol(A);  failed to converge
  // R2.print();
  vec b;
  b << U*M << M;
  vec g = ones<vec>(M);
  vec B= join_cols(g,b);
  vec X = solve(A,B);  //should be replaced by cholevsky
  vec W =  zeros<vec>(M);
  for (int m=0; m < M; m++)
    {
      if ( X[m] <0 ) cerr << "error: negative weight!" << endl;
      else W[m]= X[m];
    }
  //  X.print();
  double delta = 1e-9; 
  if ( (sum(W) -M > delta)  || (  ( dot(u_realizations,W )/M -U) > delta ) )
      {
  	cerr << "Constraints not fulfilled in calculated weights!" << endl;
      }
  FILE  *file_weights;
  char temp[4096];
  sprintf(temp, "%s%s%s", "data/", "weights", ".dat");
  file_weights= fopen(temp, "w");
  if (file_weights == 0) { cerr << "Error: can't open weights file \n" << endl; }
  for (int m=0; m < M; m++) {
    fprintf(file_weights, "%d %f \n", m, W(m) );
  }  
  return W;
}

vec calculate_weights2 (mat u_realizations, vec U)
{

 
  int M_red = u_realizations.n_cols;
  int N= u_realizations.n_rows -1;        //omdat je begint te labelen bij 0 (en omdat later de laatste rij verwijderd wordt)
  // u_realizations.print();
  // cout << "M= "  << M << endl;
  mat I =  eye<mat>(M_red,M_red);
  
  rowvec rowg =  u_realizations.row(N);                        // get the cardinalities (from the last row of u_realizations-matrix)
  cout << "here0" << endl;
  vec g = rowg.t();
  g.print();
  u_realizations.shed_row( N );                     // remove that row

  mat C =   join_cols(u_realizations, ones<rowvec>(M_red));                          //u_coarse_states.insert_rows(1,ones(M) ) ;
  mat A_l = join_cols(I, C);
  cout << "here1" << endl;
  mat A_r =  join_cols( C.t(), zeros<mat>(N+1,N+1));
  mat A = join_rows(A_l,A_r);
  cout << "print A" << endl; A.print();
  // mat R2 = chol(A);  failed to converge
  // R2.print();


  vec norm; norm << M_red;
  vec b= join_cols(U*M_red, norm);

  //b.print();
  vec B= join_cols(g,b);

  cout << "print B" << endl;  B.print();  cout << "and now solve the equation!" <<  endl;
  vec X = solve(A,B);  //should be replaced by cholevsky
  vec W =  zeros<vec>(M_red);
  for (int m=0; m < M_red; m++)
    {
      if ( X[m] <0 ) cerr << "error: negative weight!" << endl;
      else W[m]= X[m];
    }
   X.print();
  double delta = 1e-9; 
  if ( (sum(W) -M_red) > delta )
      {
  	cerr << "Constraint not fulfilled in calculated weights: bad norm" << endl;
      }
  u_realizations.print();
  vec check =  (u_realizations*W )/M_red -U;
for (int i=0; i < N; i++)
  { 
    if( check(i)> delta ) 
       {
  	cerr << "Constraints not fulfilled in calculated weights: bad average" << endl;
      }
    
  }
  FILE  *file_weights;
  char temp[4096];
  sprintf(temp, "%s%s%s", "data/", "weights", ".dat");
  file_weights= fopen(temp, "w");
  if (file_weights == 0) { cerr << "Error: can't open weights file \n" << endl; }
  for (int m=0; m < M_red; m++) {
    fprintf(file_weights, "%d %f \n", m, W(m) );
  }  
  return W;
}

double coarse_step(double U_guess, double mean_coupling, double var_coupling, double mean_preference, double var_preference, Network* network, int M, int time_horizon)
{
  double U_coarse =0;
       for(int j= 0; j < M; j++)
	 {   
	   network->lift(mean_coupling,var_coupling, mean_preference, var_preference, U_guess);        //=lifting step
	   double state_before = network->get_node(7)->get_state();
	   //	   cout <<  "state_before= " << state_before << endl;
	   Opinion_formation sim(network); 
	   sim.run_simulation(time_horizon);	  
	   // U_coarse += network->get_coarse_state();   //what with realizations cancelling each other out?
	   double state = network->get_node(7)->get_state();
	   //	   cout << "state_after= " << state << endl;
	   U_coarse +=state;   //what with realizations cancelling each other out?
	 }
       
    return   U_coarse /= M;    
}

vec coarse_stepper(double U_guess, double mean_coupling, double var_coupling, double mean_preference, double var_preference, Network* network, int M, int time_horizon)
{
  //double U_coarse =0;
  vec U_coarse = zeros<vec>(network->size());
       for(int j= 0; j < M; j++)
	 {   
	   //vec states_before = zeros<vec>(N); 
	   network->lift(mean_coupling,var_coupling, mean_preference, var_preference, U_guess);   //voorlopig laat ik U_guess nog scalair     //=lifting step
	   //  vec states_before = network->get_states();
	   vec states_before = conv_to< vec >::from(network->get_states());  // Conversion from std::vector to Armadillo vector  (dimension: N
	   //	   cout << "sampled states in simple coarse stepper " << endl; 	   states_before.print();
	   Opinion_formation sim(network); 
	   sim.run_simulation(time_horizon);	  
	   vec states = conv_to< vec >::from(network->get_states()); 
	   U_coarse +=states;   
	 }
       //states_before.print();
 
    return   U_coarse /= M;    
}



double lift_restrict(double U_guess, double mean_coupling, double var_coupling, double mean_preference, double var_preference, Network* network, int M, int time_horizon)
{
  // network->set_seed(1);
  double U_coarse =0;
       for(int j= 0; j < M; j++)
	 {   
	   network->lift(mean_coupling,var_coupling, mean_preference, var_preference, U_guess);        //=lifting step
	   // U_coarse += network->get_coarse_state();   //what with realizations cancelling each other out?
	   double state = network->get_node(7)->get_state();
	   //	   cout << state << endl;
	   U_coarse +=state;   //what with realizations cancelling each other out?
	 }
       
    return   U_coarse /= M;    
}

double lift_restrict_weighted(double U_guess, double mean_coupling, double var_coupling, double mean_preference, double var_preference, Network* network, int M, int time_horizon)
{
  double U_coarse =0;
  rowvec sampled_states = zeros<rowvec>(M); 
  //mat sampled_states = zeros<mat>(N,M);

       for(int j= 0; j < M; j++)
	 {   
	   network->lift(mean_coupling,var_coupling, mean_preference, var_preference, U_guess);        //=lifting step
	   sampled_states(j)= network->get_node(7)->get_state(); //omdat we eendimensionaal werken    
	   // Opinion_formation sim(network); 
	   //	   sim.run_simulation(time_horizon);
	 }
       //    sampled_states.print();
       vec w= calculate_weights(sampled_states, U_guess);
       return dot(sampled_states,w)/M;
    
}


double coarse_step_weighted(double U_guess, double mean_coupling, double var_coupling, double mean_preference, double var_preference, Network* network, int M, int time_horizon)
{
  double U_coarse =0;
  rowvec sampled_states = zeros<rowvec>(M); 
 rowvec new_states = zeros<rowvec>(M); 
  //mat sampled_states = zeros<mat>(N,M);

       for(int j= 0; j < M; j++)
	 {   
	   network->lift(mean_coupling,var_coupling, mean_preference, var_preference, U_guess);        //=lifting step 
	   sampled_states(j)= network->get_node(7)->get_state();  //omdat we eendimensionaal werken  
	   Opinion_formation sim(network); 
	   sim.run_simulation(time_horizon);
	   new_states(j)= network->get_node(7)->get_state();
	 }
       vec w= calculate_weights(sampled_states, U_guess);
       // cout << "before simulation" << endl;
       // sampled_states.print();       
       // cout << "after simulation" << endl;
       // new_states.print();
       return dot(new_states,w)/M;
    
}

vec coarse_stepper_weighted(vec U_guess, double mean_coupling, double var_coupling, double mean_preference, double var_preference, Network* network, int M, int time_horizon)
{
  int N= network->size();
  vec U_coarse = zeros<vec>(N);
  mat sampled_states = zeros<mat>(N,0);
  //mat sampled_states =0;
  mat new_states = zeros<mat>(N,0);
  
       for(int j= 0; j < M; j++)
	 {   
	   network->lift(mean_coupling,var_coupling, mean_preference, var_preference, U_guess(1));        //=lifting step
	   sampled_states = join_rows ( sampled_states, conv_to< vec >::from(network->get_states()) ) ;  // Conversion from std::vector to Armadillo vector  (dimension: N)
	   Opinion_formation sim(network); 
	   sim.run_simulation(time_horizon);	  
	   new_states = join_rows (new_states, conv_to< vec >::from(network->get_states()) ) ; //add new realization of N states as a new column
	 }       
       sampled_states.print();
       cout << "new states " << endl;
       new_states.print();
       //vec w = ones(M);
       vec w = calculate_weights2(sampled_states, U_guess);
       return  new_states*w/M;    
}


vec coarse_stepper_weighted_cardin(vec U_guess, double mean_coupling, double var_coupling, double mean_preference, double var_preference, Network* network, int M, int time_horizon)
{
  int N= network->size();
  //vec card = zeros<vec>(N);
  mat sampled_states_red = zeros<mat>(N+1,0);
  mat sampled_states_orig = zeros<mat>(N,0);
  //mat sampled_states =0;
  mat new_states = zeros<mat>(N,0);
  // rowvec back_to_orginal_representation = zeros<rowvec>(M);
  
 
       for(int j= 0; j < M; j++)
	 {   
	   network->lift(mean_coupling,var_coupling, mean_preference, var_preference, U_guess(1));        //=lifting step   //voorlopig homogeneous
	   //  vec states_before = network->get_states();
	   sampled_states_orig = join_rows (   sampled_states_orig, conv_to< vec >::from(network->get_states()) ) ;
	   vec card1  = ones<vec>(1); //vector filled with cardinality= 1 (initially)
	   vec new_realization = join_cols(conv_to< vec >::from(network->get_states()), card1); //N+1
	   bool unique_realization= true;
	   // cout << "number of added realizations (= number of columns) =" << sampled_states_red.n_cols  << endl;

	   int col =0;
	   while (  col < sampled_states_red.n_cols && unique_realization) //before adding realization, check whether this realization doesn't already exists (check for identical columns in constraint matrix)
	     {
	       int row= 0;
	       //    sampled_states_red.print();   cout << endl;
	       //   cout << sampled_states(0,0) << endl;	      
	       while (sampled_states_red(row,col) == new_realization(row) && row <N)
		 { row++;}
	    
	       if (row>=N-1)       //identical column -> realization already exists -> cardinality++
		 {   
		   sampled_states_red(N,col) = sampled_states_red(N,col) +1;
		   unique_realization = false;
		 }
	       col++;
	     } 
	   if (unique_realization)  { 	        cout <<"unique realization" << endl;   sampled_states_red = join_rows ( sampled_states_red, new_realization) ; cout <<"added" << endl;	      }  
	   
      
	   Opinion_formation sim(network); 
	   sim.run_simulation(time_horizon);	  
	   new_states = join_rows (new_states, conv_to< vec >::from(network->get_states()) ) ;
	 }       
       cout << "M=  " << M << " realizations" << endl;
       cout << "M_red'=  " << sampled_states_red.n_cols << " unique realizations" << endl;
       cout << "Original sampled states:" << endl;
       sampled_states_orig.print();
       cout << "Reduced sampled states:" << endl;
       sampled_states_red.print();
       cout << "new states " << endl;
       new_states.print();
       //vec w = ones(M);
       vec w_red = calculate_weights2(sampled_states_red, U_guess);
       vec w =  zeros<vec>(M);
       
       int col =0;
       for (int col =0; col  < sampled_states_orig.n_cols ; col++)       //transform weights back to M-vector: get the weight if the column of the original representation equals the column in the reduced representation 
	 {
	   for (int col_red =0; col_red  < sampled_states_red.n_cols ; col_red++)
	     {
	       int row= 0;	      
	       while (sampled_states_orig(row,col) == sampled_states_red(row,col_red) && row < N-1 )
		 { row++;}
	    
	       if (row==N-1) // same realization       
		 {   
		   w(col) = w_red(col_red)/ sampled_states_red(N, col_red);
		 }
	       
	     } 
	 }
       cout << "w_red=  " << endl; w_red.print() ;
       cout << "w =  " << endl; w.print();
       return  new_states*w/M;    
}



rowvec sample_states (double U_guess, int M)   //ter illustratie
{ 
  MTRand mtrand;  
  rowvec sampled_states = zeros<rowvec>(M); 

        for(int m= 0; m < M; m++)
	  {   
	    if (mtrand.rand() < U_guess)     
	      sampled_states(m) = 1;
	    else sampled_states(m)=0; 
	 }
	cout << "sampled s0: " <<  sampled_states[0] << endl;
		cout << "sampled s1: " <<  sampled_states[1] << endl;
		cout << "sampled sM: " <<  sampled_states[M] << endl;
		cout << "sampled sM: " <<  sampled_states[M+1] << endl;
	return sampled_states;

}





double derivative_coarse_step(double U_guess, double mean_coupling, double var_coupling, double mean_preference, double var_preference, Network* network, int M, int time_horizon) 
{
  double epsilon = 1e-20;
  double Ueps ;
  double U;
  FILE  *file_eps;
  char temp[4096];
  sprintf(temp, "%s%s%s", "data/", "epsilon", ".dat");
  file_eps = fopen(temp, "w");
  if (file_eps == 0) { cerr << "Error: can't open epsilon file \n" << endl; }
  while (epsilon < 0.1)
      { 
	//	network->save_generator_state();
 	Ueps = coarse_step(U_guess+epsilon, mean_coupling,var_coupling, mean_preference, var_preference, network, M, time_horizon);
	//	network->load_generator_state(); //reproduce the same random numbers
	U= coarse_step(U_guess, mean_coupling,var_coupling, mean_preference, var_preference, network, M, time_horizon ) ;
	//	fprintf(file_eps, "%e %e \n", epsilon, abs(Ueps-U) );
	fprintf(file_eps, "%e %e \n", epsilon, abs(epsilon - Ueps +U) ); //F(U+e)-F(U)
	cout << "epsilon= " << epsilon << " difference:|Ueps-U| = " << abs(Ueps  - U) << endl ; //  " U= " << U << endl;
	epsilon*=10;
      }    
 return (Ueps-U)/epsilon;
 
}


double derivative_coarse_step_weighted(double U_guess, double mean_coupling, double var_coupling, double mean_preference, double var_preference, Network* network, int M, int time_horizon) 
{
  double epsilon = 1e-20;
  double Ueps ;
  double U;
  FILE  *file_eps_w;
  char temp[4096];
  sprintf(temp, "%s%s%s", "data/", "epsilon_weighted", ".dat");
  file_eps_w = fopen(temp, "w");
  if (file_eps_w == 0) { cerr << "Error: can't open epsilon file \n" << endl; }
  while (epsilon < 0.1)
      { 
	network->save_generator_state();
	Ueps = coarse_step_weighted(U_guess+epsilon, mean_coupling,var_coupling, mean_preference, var_preference, network, M, time_horizon);
	network->load_generator_state(); 
	U= coarse_step_weighted(U_guess, mean_coupling,var_coupling, mean_preference, var_preference, network, M, time_horizon ) ;
	//	fprintf(file_eps_w, "%e %e \n", epsilon, abs(Ueps-U) );
	fprintf(file_eps_w, "%e %e \n", epsilon, abs( epsilon - Ueps +U) ); //F(U+e)-F(U)
	cout << "epsilon= " << epsilon << " difference (weighted): |Ueps-U| = " << abs(Ueps  - U) << endl ; //  " U= " << U << endl;
	epsilon*=10;
      }
 return (Ueps-U)/epsilon;
}

double derivative_coarse_stepper_weighted(vec U_guess, double mean_coupling, double var_coupling, double mean_preference, double var_preference, Network* network, int M, int time_horizon) 
{
  double epsilon=1e-20;
  vec U= zeros<vec>(U_guess.n_elem);
  vec Ueps= zeros<vec>(U_guess.n_elem); 
  vec epsilons= zeros<vec>(U_guess.n_elem);
  FILE  *file_eps_ws;
  char temp[4096];
  sprintf(temp, "%s%s%s", "data/", "epsilon_weighted_cardin", ".dat");
  file_eps_ws = fopen(temp, "w");
  if (file_eps_ws == 0) { cerr << "Error: can't open epsilon file \n" << endl; }
  while (epsilon < 0.1)
      { 
	epsilons.fill(epsilon);
	network->save_generator_state();
	Ueps = coarse_stepper_weighted_cardin(U_guess+epsilons, mean_coupling,var_coupling, mean_preference, var_preference, network, M, time_horizon);
	network->load_generator_state(); 
	U= coarse_stepper_weighted_cardin(U_guess, mean_coupling,var_coupling, mean_preference, var_preference, network, M, time_horizon ) ;
	//	fprintf(file_eps_w, "%e %e \n", epsilon, abs(Ueps-U) );
	fprintf(file_eps_ws, "%e %e \n", epsilon, norm( epsilon - Ueps +U,2) ) ; //F(U+e)-F(U)
	cout << "epsilon= " << epsilon << " 2-norm of difference (weighted): |Ueps-U| = " << abs(Ueps  - U) << endl ; //  " U= " << U << endl;
	epsilon*=10;
      }
 return  norm( epsilon - Ueps +U,2);
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

double newton_raphson_weighted(double U_guess, double mean_coupling, double var_coupling,  double mean_preference, double var_preference, Network* network, int M, int time_horizon)
{ 
  double tol=1e-6;
  double err=tol+1;
  double it = 0;
  double maxit=50;
  double U_coarse;
  while(  err > tol && it < maxit )
    {
 
      U_coarse=  U_guess-  ( coarse_step_weighted(U_guess, mean_coupling, var_coupling, mean_preference, var_preference, network, M, time_horizon)  - U_guess)/ (derivative_coarse_step_weighted(U_guess, mean_coupling,var_coupling,  mean_preference, var_preference, network, M, time_horizon) -1);  //newton-step  U_n+1 = U_n - f(U_n)/f'(u_n)
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
  int network_size, M, time_steps, degree;
  double mean_preference, var_preference, mean_coupling, var_coupling, average_state, p_rewiring ,rg_seed;
  string output_name, network_type;
  FILE *file_s0, *file_s1, *file_m, *file_mean, *file_var, *file;
  time_t starttime; 
  bool separate_coarse_states = false;
  MTRand mtrand;

  // mat D;
  // D << 1 << 3 << 5 << endr     << 2 << 4 << 6 << endr;
  // cout << "check (1,2)=3 : " <<  D[1,2] << endl; //starts counting from 1!
  //  cout << "check (1,2)=6 : " <<  D(1,2) << endl; //starts counting from 0!
//   PetscErrorCode ierr; /* used to check for functions returning nonzeros */
//   Vec x; /* solution vector */
//   Mat H; /* Hessian matrix */
//   Tao tao; /* Tao context */
//   AppCtx user; /* user-defined application context */
//   PetscInitialize(&argc,&argv,(char*)0,0);
//   ierr = PetscPrintf( PETSC_COMM_WORLD, "\nHello World\n" );CHKERRQ(ierr);
//   /* Initialize problem parameters */
//   user.n = 2; user.alpha = 99.0;
// /* Allocate vectors for the solution and gradient */
//   ierr = VecCreateSeq(PETSC_COMM_SELF,user.n,&x); CHKERRQ(ierr);
//   ierr = MatCreateSeqBAIJ(PETSC_COMM_SELF,2,user.n,user.n,1,NULL,&H);
//   /* Create TAO solver with desired solution method */
//   ierr = TaoCreate(PETSC_COMM_SELF,&tao); CHKERRQ(ierr);
//   ierr = TaoSetType(tao,TAOLMVM); CHKERRQ(ierr);
//   /* Set solution vec and an initial guess */
//   ierr = VecSet(x, 0); CHKERRQ(ierr);
//   ierr = TaoSetInitialVector(tao,x); CHKERRQ(ierr);
// /* Set routines for function, gradient, hessian evaluation */
//   ierr = TaoSetObjectiveAndGradientRoutine(tao,FormFunctionGradient,&user);
//   ierr = TaoSetHessianRoutine(tao,H,H,FormHessian,&user); CHKERRQ(ierr);
//   /* Check for TAO command line options */
//   ierr = TaoSetFromOptions(tao); CHKERRQ(ierr);
// /* Solve the application */
//   ierr = TaoSolve(tao); CHKERRQ(ierr);
//   /* Free data structures */
//   ierr = TaoDestroy(&tao); CHKERRQ(ierr);
//   ierr = VecDestroy(&x); CHKERRQ(ierr);
//   ierr = MatDestroy(&H); CHKERRQ(ierr);
//   PetscFinalize();

 

  // cout << "armadillotest"  << endl;
  // // vec y = zeros<vec>(5);  
  // //mat A = randu<mat>(10,10);
  // vec z;
  // z << 2 << 4 << 9 <<1; 
  // vec q = square(z);
  // cout << q[0] << endl;


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
 if (network_type.compare("0")==0) network_type="random";
 if ( degree== 0) degree =(int) network_size/40;  
 if (average_state == 0) average_state = 0.5;
 if (var_preference == 0) var_preference = 0.236;    
 if (mean_coupling == 0) mean_coupling = 0.5; 
 if (var_coupling ==0)var_coupling=0.167;   
 if (rg_seed ==0)  rg_seed=time(0);        // mtrand->seed(rg_seed);  
 if (strlen(output_name.c_str())<2) {
   char temp[1024];
   sprintf(temp,"Sim_N_%d_M_%d",(int) network_size, (int) M);
   output_name = string(temp);
 }
 if(p_rewiring ==0) p_rewiring = 0.1;


 //mtrand.seed(rg_seed);   //set seed if you want repeatable data
      //mtrand->seed(42);        //segmenation fault



 cout << "Time = " << time_steps << " steps" << endl;
 cout << "M = " << M << endl;
 cout << "Network size = "<< network_size<< " nodes" << endl;
 //cout << "RG seed  = "<< rg_seed <<  endl;

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
       // double U_guess=average_state;
       // fixpoint_iteration_analytical(U_guess, mean_coupling, mean_preference, var_preference);
       //  newton_raphson_analytical(U_guess, mean_coupling, mean_preference, var_preference);
  

       // LIFTING/restriction
        int time_horizon =20; // 20;
	Network*  network =  construct_network(network_type, network_size, degree, p_rewiring);
	network->set_seed(rg_seed);
	vec U_guess= zeros<vec>(network_size);
	U_guess.fill(average_state);
       //   mtrand = network->get_rng();
       //network->setSeed(42);
       // newton_raphson(U_guess, mean_coupling, var_coupling, mean_preference, var_preference, network, M, time_horizon, mtrand);  
       // fixpoint_iteration(U_guess, mean_coupling, var_coupling, mean_preference, var_preference, network, M, time_horizon);  
       // vector<double> U;
       // for (int i=0; i < network_size; i++)
       // 	 { U.push_back(0.4);}  //Suppose U is a N-dim testvector     zoals  U = net->get_node_states(); hoewel in het algemeen is elke toestand reëel getal is ipv (0,1)
 
       //double U= 0.7; //average_state;
	//       double epsilon = 1e-3;
       // rowvec u_realizations = sample_states(U,M);
       // calculate_weights(u_realizations,U);
       //       calculate_weights(u_realizations, U+epsilon);
  
       //  newton_raphson_weighted(U_guess, mean_coupling, var_coupling, mean_preference, var_preference, network, M, time_horizon, mtrand);  
        // network->set_seed(rg_seed);
        // cout << "check identity relation: U: " << lift_restrict(U_guess(1), mean_coupling, var_coupling, mean_preference, var_preference, network, M, time_horizon) << endl;
        // network->set_seed(rg_seed);
        // cout << "check identity relation using weights" << endl;
        // cout <<   "U: " <<   lift_restrict_weighted(U_guess(1), mean_coupling, var_coupling, mean_preference, var_preference, network, M, time_horizon) << endl;  
        // network->set_seed(rg_seed);
        // cout << "course step" << endl;
        // network->set_seed(rg_seed);
	// cout << "U =" <<  coarse_stepper(average_state, mean_coupling, var_coupling, mean_preference, var_preference, network, M, time_horizon) << endl;
	//  cout << "weighted course step" << endl;
        //network->set_seed(rg_seed);
	//    cout <<   "U: " <<   coarse_stepper_weighted(U_guess, mean_coupling, var_coupling, mean_preference, var_preference, network, M, time_horizon) << endl;  
	// cout << "checksame random values " << endl;
	//	network->set_seed(rg_seed);
	//	cout << "only save unique realization " << endl;
	//        cout <<   "U: " <<   coarse_stepper_weighted_cardin(U_guess, mean_coupling, var_coupling, mean_preference, var_preference, network, M, time_horizon) << endl;  
	cout <<  "derivative" <<   derivative_coarse_step(U_guess(1), mean_coupling, var_coupling, mean_preference, var_preference, network, M, time_horizon) << endl;  
	cout <<  "derivative weighted " <<   derivative_coarse_step_weighted(U_guess(1), mean_coupling, var_coupling, mean_preference, var_preference, network, M, time_horizon) << endl;  
	cout <<  "derivative cardin weighted " <<   derivative_coarse_stepper_weighted(U_guess, mean_coupling, var_coupling, mean_preference, var_preference, network, M, time_horizon) << endl; 

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


