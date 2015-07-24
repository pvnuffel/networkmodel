#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <string.h>
#include "Network.h"
#include <assert.h>
#include <queue>
#include <ctime>
#include <time.h>
#include <armadillo> //linear algebra library
using namespace arma;

class Simulator
{

    public:

        int time;
        Network* net;
        MTRand* mtrand;
	vector<double> coarse_states;
        Simulator() { time = 0; net=NULL; mtrand=NULL; };
	//   Simulator(Network* net) { this->net = net; this->time = 0; this->mtrand = net->get_rng();};	//   Simulator(Network* net) { this->net = net; this->time = 0; this->mtrand = net->get_rng()};
        Simulator(Network* net)  {this->net = net; this->mtrand = net->get_rng(); this->time = 0; this->coarse_states = get_coarse_states(); };
	//    	Simulator(Network* net, MTRand mtrand)  {this->net = net; this->time = 0; this.mtrand= mtrand, this->coarse_states = get_coarse_states(); };

        void set_network( Network* net ) { this->net = net; this->mtrand = net->get_rng(); };
        Network* network() { return(net); };

        int get_time() { return(time); };

        vector<double> get_coarse_states() { return(coarse_states); };

        void reset_time() { time = 0; };

        void set_all_nodes_to_state ( stateType s ) {
            vector<Node*> nodes = net->get_nodes();
            set_these_nodes_to_state(nodes, s);
        };

        void set_these_nodes_to_state (vector<Node*> nodes, stateType s) {
            for (unsigned int i = 0; i < nodes.size(); i++) nodes[i]->set_state(s);
        }

        // choose n nodes without replacement
        vector<Node*> rand_choose_nodes (int n) {
            assert(n > -1 and n <= net->size());
            vector<Node*> nodes = net->get_nodes();
            vector<Node*> sample(n);
            vector<int> sample_ids(n);
            rand_nchoosek(net->size(), sample_ids, mtrand);
            Node* node;
            for (unsigned int i = 0; i < sample_ids.size(); i++) {
                node = nodes[ sample_ids[i] ];
                sample[i] = node;
            };
            return sample;
        }

        // change n random nodes to state s (e.g. vaccinate them or infect them randomly)
        vector<Node*> rand_set_nodes_to_state (int n, stateType state) {
            vector<Node*> sample = rand_choose_nodes(n);
            for (unsigned int i = 0; i < sample.size(); i++) {
                sample[i]->set_state(state);
            };
            return sample;
        }

        /* vector<Node*> set_node_to_random_states (stateType state) {               // verplaatst naar network.cpp */
	/*     vector<Node*> nodes = net->get_nodes(); */
        /*     for (unsigned int i = 0; i < net->size(); i++) { */
	/*       if (mtrand->rand() < 0.5)  */
        /*         nodes[i]->set_state(state); */
        /*     }; */
        /*     return sample; */
        /* } */




        //these functions must be derived in child class
        virtual void time_step() {};
        virtual void run_simulation() {};
	// virtual vector<Node*> rand_infect(int) {vector<Node*>x; return x;};
                                 // cumulative infected
	// virtual int epidemic_size() = 0;
                                 // current infected
        virtual int count_infected() {
            return 0;
        };
        virtual void reset() {};

        protected:
	//	time_t starttime;
            //~Simulator() { };

};

/*
class Derived_Example_Simulator: public Simulator {

    public:
        Derived_Example_Simulator();
        Derived_Example_Simulator(Network* net);
        ~Derived_Example_Simulator();

        void step_simulation() { time++; do_stuff_to_network;  };
        void run_simulation() { while(conditional) { step_simulation() } };

};
*/
#endif
