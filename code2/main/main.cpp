#include <Opinion_formation.h>

int main(int argc, char* argv[])  {

  vector<Network*> net_list;  
  int network_size, M, time_steps, lambda, rg_seed;
  double mean_preference, var_preference, mean_coupling, var_coupling, average_state;
  string output_name, network_type;
  time_t starttime; 
  //MTRand* mtrand;

      cout<<  "usage: "<< argv[0] ;          
      cout << "Interactive mode. Enter 0 for default value\n";
      cout << "Enter the number of time steps "; cin >> time_steps;     
      cout << "Enter number of Monte Carlo steps: "; cin >> M;
      cout << "Enter number of nodes: "; cin >> network_size;      
      cout << "Enter netwerk type: erdos or mean-field "; cin >> network_type;
      cout << "Enter average number of edges: "; cin >> lambda;                 
      cout << "Enter initial average state (in [0,1] "; cin >> average_state;          
      cout << "Enter mean preference in [-1,1]"; cin >> mean_preference;          
      cout << "Enter preference variance in [0,1]"; cin >> var_preference;         
      cout << "Enter mean interaction strength in [0,1]"; cin >> mean_coupling;         
      cout << "Enter interaction strength variance in [0,1]"; cin >> var_coupling;           
      cout << "Enter output-file name (without .dat)"; cin >> output_name;      
      cout << "Enter seed if you want repeatable data"; cin >> rg_seed;
      cout << endl;
      //other options (add later): network type, directed?, output_filename, etc.
 

 //   // Correct values in case of invalid input interactively or through the command line ( 0 or - ).
 // //cout << "time" << time << endl;
 if (time_steps == 0) time_steps = 20;
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
 //mtrand->seed(42); geeft segmentation fault
      //mtrand->seed(42);          //set seed if you want repeatable data



 cout << "Time = " << time_steps << " steps" << endl;
 cout << "M = " << M << endl;
 cout << "Network size = "<< network_size<< " nodes" << endl;



  	for(int j= 0; j < M; j++){                                                               
	  //	  Network net("consumers", Network::Undirected);
	  Network* net = new Network("consumers", Network::Undirected); // Construct Networks
	  net->populate(network_size);
	  
	    if (network_type.compare("erdos")==0)
	      {
		net->erdos_renyi(lambda);      //evaluatie input_string binnen  lus over MC-steps niet-ideaal
	      }
	    else if (network_type.compare("mean_field")==0)
	      {
		net->all_to_all_coupling();
	      }
	    else if (network_type.compare("small_world")==0)
	      {
		net->small_world(network_size, lambda, 0.5);      //beta=0.3   , beta=1 -> erdos random network, beta=0 -> latice-model
	      }
	    else { cerr << "No valid network type!" << endl;}

	  net->initialize(mean_coupling,var_coupling, mean_preference, var_preference,average_state); 
	  net_list.push_back(net); 
	}
       
	//   Opinion_formation sim(net_list[5]);   
	Opinion_formation sim(net_list, output_name);   
	//cout << net_list[5]-> << rss << endl;
 try
   {
       starttime = time(NULL); //now */
        sim.run_simulation(time_steps);

	int64_t seconds = (int64_t) difftime(time(NULL),starttime);
	int64_t days = seconds / 86400; seconds -= days * 86400;
	int64_t hours = seconds / 3600; seconds -= hours * 3600;
	int64_t minutes = seconds / 60; seconds -= minutes * 60;
	fprintf(stdout,"Simulation time: %d days %d hours %d minutes %d seconds\n",(int)days,(int)hours,(int)minutes,(int)seconds);

        return 0;
    }


    catch (const bad_alloc& ex)
    {   //Observe that the local variable `data` no longer exists here.
      cerr << "Oops. Looks like you need to use a 64 bit system or get a bigger hard disk for that calculation!"<< endl;
        return -1;
    }






	//   net.write_edgelist("output_edgelist.csv"); 
	// SimulationTime(stdout);
	//	Ensemble ensemble(&net);
	//	ensemble.ensemble_average(7);

	// double vm, rss;
	// process_mem_usage(vm, rss);
	// cout << "VM: " << vm << "; RSS: " << rss << endl;

}


// void process_mem_usage(double& vm_usage, double& resident_set)  // http://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-run-time-in-c
// {
//    using std::ios_base;
//    using std::ifstream;
//    using std::string;

//    vm_usage     = 0.0;
//    resident_set = 0.0;

//    // 'file' stat seems to give the most reliable results
//    //
//    ifstream stat_stream("/proc/self/stat",ios_base::in);

//    // dummy vars for leading entries in stat that we don't care about
//    //
//    string pid, comm, state, ppid, pgrp, session, tty_nr;
//    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
//    string utime, stime, cutime, cstime, priority, nice;
//    string O, itrealvalue, starttime;

//    // the two fields we want
//    //
//    unsigned long vsize;
//    long rss;

//    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
//                >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
//                >> utime >> stime >> cutime >> cstime >> priority >> nice
//                >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

//    stat_stream.close();

//    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
//    vm_usage     = vsize / 1024.0;
//    resident_set = rss * page_size_kb;
// }
