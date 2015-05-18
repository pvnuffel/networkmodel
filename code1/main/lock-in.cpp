#include <Opinion_formation.h>
#include <armadillo> //linear algebra library
using namespace arma;

int main(int argc, char* argv[])  {
cout << "test"  << endl;
  vec y = zeros<vec>(5);  
  //mat A = randu<mat>(10,10);
  vec z = 1 << 2 << 4 << 9 <<1; 
  vex q = z+y
    cout << q[3] << endl;
  vector<Network*> net_list;  
  int network_size, M, time, lambda, rg_seed;
  double mean_preference, var_preference, mean_coupling, var_coupling, average_state;
  string output_name, network_type;
 
  //MTRand* mtrand;
  // cout<<  "argc: "<< argc <<"\n" ;
  //if ( argc !=2  ) //  at least 8 parameters for correct execution 
    // We print argv[0] assuming it is the program name.
      // {
      cout<<  "usage: "<< argv[0] ;          
      cout << "Interactive mode. Enter 0 for default value\n";
      cout << "Enter the number of time steps "; cin >> time;     
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
      //	}
 //  else {
 //   // We assume argv[1] is a filename to open
 //    ifstream the_file ( argv[1] );
 //    char input[15]; 
 //    // Always check to see if file opening succeeded
 //    if ( !the_file.is_open() )
 //      cout<<"Could not open file\n";
 //    else {
 //      char x; char input[15]; int i=0;
 //      // the_file.get ( x ) returns false if the end of the file
 //      //  is reached or an error occurs
 //      while ( the_file.get ( x ) ) {
 //         cout<< x;
 //      i++;
 //      input[i] = x;
 //      }
 //    }
 //    time = (int) input[1];                    
 //    //   M = atoi(input[2]);                  //Convert string to int 
 //    //   network_size = atoi(argv[3]);
 //    // lambda =  atoi(argv[4]);
 //    // average_state = atoi(argv[5]);                  
 //    // mean_preference = atoi(argv[6]);
 //    // var_preference = atoi(argv[7]);                  
 //    // mean_coupling = atoi(argv[8]);
 //    // var_coupling = atoi(argv[9]);                  
 //  
 //  }
 //   // Correct values in case of invalid input interactively or through the command line ( 0 or - ).
 // //cout << "time" << time << endl;
 if (time == 0) time = 20;
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



 cout << "Time = " << time << " steps" << endl;
 cout << "M = " << M << endl;
 cout << "Network size = "<< network_size<< " nodes" << endl;

/*!
 * constructing and initializing the networks..
 */
// vector<vector<double> > lockin_states_0;
// vector<vector<double> > lockin_states_1;
// vector<vector<double> > mixed_states;

//double*  lockin_states_0 = ;
//colvec lockin_states_0(time, 0.0);
    // file = fopen(string(output_name).append(".dat").c_str(),"w");
    // char temp[4096];
    // sprintf(temp, "%s%s%s", "data/", fname.c_str(), ".dat");
    // file = fopen(temp, "w");
    // if (file == 0) { cerr << "Error: can't open file \n" << fname << endl; }

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
            else { cerr << "No valid network type!" << endl;}

	  net->initialize(mean_coupling,var_coupling, mean_preference, var_preference,average_state); 
	  Opinion_formation sim(net); 
	  sim.run_simulation(time);
	  vector<double> coarse_states = sim.get_coarse_states();
	  //	  colvec coarse_states = conv_to< colvec >::from(sim.get_coarse_states());  // Conversion from std::vector to Armadillo vector
          if ( coarse_states[time] > 0.7 )  { cout << "lockin "  << endl;   }// lockin_states_0 = lockin_states_0 + coarse_states;  }
	  else { cout << "no lockin "  << endl;    }  
	  
	  // net_list.push_back(net); 
	}
       	//   Opinion_formation sim(net_list[5]);   
	//   net.write_edgelist("output_edgelist.csv"); 
	// SimulationTime(stdout);
	//	Ensemble ensemble(&net);
	//	ensemble.ensemble_average(7);

	// double vm, rss;
	// process_mem_usage(vm, rss);
	// cout << "VM: " << vm << "; RSS: " << rss << endl;
	return 0;
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




