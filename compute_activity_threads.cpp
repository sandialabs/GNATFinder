/* 
 * GNATFinder
 *
 * Copyright 2023 National Technology & Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, 
 * this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote 
 * products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

/*
 * First order activity graph computation
 * Brad Theilman 2023
 */

#include <set>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <stdlib.h>

#define TICKS_PER_MS 1000000
#define GNATS 1
#define CDH   2

typedef unsigned long tstamp_t; // spike timestamp 
typedef unsigned long idx_t;
typedef double real_t;

double gamma(tstamp_t t1, tstamp_t t2, double weight, double delay, double tau) {
    double res;
    int theta;
    theta = (t2 - t1) >= delay ? 1 : 0;
    res = -log(weight * theta * exp(-(t2 - t1 - delay) / tau));
    return res;
}


/**********************************************************/
class SpikeRaster {
    public:
        SpikeRaster(idx_t n_neurons);
        ~SpikeRaster() { };

        idx_t n_neurons;
        int read_event_file(std::string fname);
        void get_spikes_in_range(std::list<tstamp_t>& res, idx_t neuron_idx, tstamp_t low, tstamp_t high);
        std::vector< std::set<tstamp_t> > evtlist; // Vector of sets of spikes, one set for each neuron
};

// Constructor: Initialize the event list by pushing back an empty list for each neuron
SpikeRaster::SpikeRaster(idx_t N) {

    n_neurons = 0;
    for (idx_t i = 0; i < N; ++i) {
        std::set<tstamp_t> s;
        evtlist.push_back(s);
        n_neurons++;
    }
}

// This function populates the list res with spike times from neuron <neuron_idx> that fall within the time range
// [low, high]
void SpikeRaster::get_spikes_in_range(std::list<tstamp_t>& res, idx_t neuron_idx, tstamp_t low, tstamp_t high) {
    
    std::set<tstamp_t>::iterator it, itlow, ithigh;

    itlow = evtlist[neuron_idx].lower_bound(low);
    ithigh = evtlist[neuron_idx].upper_bound(high);

    for (it = itlow; it != ithigh; ++it) {
        res.push_back(*it);
    }

}

// Reads spikes from a text file
// Each line in the file corresponds to a spike
// Each line has the format <event_type> <timestamp> <neuron_index>
// event_type = 0 for spikes
// timestamp is specified in a hexadecimal string
int SpikeRaster::read_event_file(std::string fname) {

    int evttype;
    tstamp_t evtstamp;
    idx_t evtidx;
    
    // Check that file exists
    // try opening file
    std::ifstream infile;
    std::string line;

    infile.open(fname.c_str(), std::ifstream::in);
    if (!infile.is_open()) {
        std::cout << "Error opening event file\n";
        exit(EXIT_FAILURE);
        return -1;
    } else { 
        std::cout << "Opened file: " << fname << "\n";
        while(std::getline(infile, line)) {
            std::istringstream iss(line);
            //std::cout << "Processing line: " << line << "\n";
            iss >> evttype >> std::hex >> evtstamp >> std::dec >> evtidx;
            if (evtidx > n_neurons) {
                std::cout << "Neuron index of event greater than number of neurons; ignoring..\n";
                break;
            }
            if (evttype == 0) {
                evtlist[evtidx].insert(evtstamp);
            }
        }
    }
    return 0;
}
/**********************************************************/
/**********************************************************/

// Each edge in the network has a source index, weight, and delay
struct edge {
    idx_t  idx;    // SOURCE index
    real_t weight;
    real_t delay;
};

class Network {
    public:
        Network(idx_t _n_neurons) {this->n_neurons = _n_neurons;};
        ~Network() {};

        int read_connectivity_csr(std::string fname);
        int read_connectivity(std::string fname);
        int compute_activity_threads(SpikeRaster& raster, std::string fname, double gamma_thresh, double temporal_radius, double tau, int func);

    private:
        // For each neuron, we have a list of presynaptic edges
        idx_t n_targets;
        idx_t n_neurons;
        std::vector<std::vector<struct edge> > presynaptic_edges;

        void emit_causal_neighbors(SpikeRaster& sr, idx_t neuron_idx, double gamma_thresh, double temporal_radius, double tau, int func, std::ofstream& outfile);
};


// Reads connectivity information from a file
// Each line specifies the presynaptic connectivity of a target neuron (target index is line number)
// Each line has the format: N_edges <edge 0 idx> <edge 0 weight> <edge 0 delay> <edge 1 idx> ...
int Network::read_connectivity_csr(std::string fname) {
    // Check that file exists
    std::ifstream infile;
    std::string line;
    idx_t line_idx = 0;
    idx_t n_edges;

    // try opening file
    infile.open(fname.c_str(), std::ifstream::in);
    if (!infile.is_open()) {
        std::cout << "Error opening connectivity file\n";
        exit(EXIT_FAILURE);
        return -1;
    } else { 
        std::cout << "Opened connectivity file: " << fname << "\n";
        while(std::getline(infile, line)) {
            // line format is : N_edges <edge 0 idx> <edge 0 weight> <edge 0 delay> <edge 1 idx> ...
            std::istringstream iss(line);
            std::vector<struct edge> edge_list;
            // read number of edges
            iss >> n_edges;
            for (idx_t edg_idx = 0; edg_idx < n_edges; ++edg_idx) {
                struct edge edg;
                iss >> edg.idx >> edg.weight >> edg.delay;
                edge_list.push_back(edg);
            }
            presynaptic_edges.push_back(edge_list);
            line_idx++;
        }
        n_targets = line_idx;
    }
    infile.close();
    return 0;
}

// Reads connectivity information from a file
// Each line specifies a synapse
// Each line has the format: src_idx tgt_idx rel_w delay
int Network::read_connectivity(std::string fname) {
    // Check that file exists
    std::ifstream infile;
    std::string line;
    idx_t line_idx = 0;
    idx_t n_edges;

    // try opening file
    infile.open(fname.c_str(), std::ifstream::in);
    if (!infile.is_open()) {
        std::cout << "Error opening connectivity file\n";
        exit(EXIT_FAILURE);
        return -1;
    } else { 
        std::cout << "Opened connectivity file: " << fname << "\n";

        // Initialize edge lists
        for (idx_t cell_idx = 0; cell_idx < n_neurons; ++cell_idx) {
            std::vector<struct edge> _edg;
            presynaptic_edges.push_back(_edg);
        }

        while(std::getline(infile, line)) {
            // line format is : N_edges <edge 0 idx> <edge 0 weight> <edge 0 delay> <edge 1 idx> ...
            std::istringstream iss(line);

            struct edge edg;
            idx_t tgt_idx;
            iss >> edg.idx >> tgt_idx >> edg.weight >> edg.delay;
            presynaptic_edges[tgt_idx].push_back(edg);

            line_idx++;
        }
    }
    infile.close();
    return 0;
}

// For each neuron in the network, compute the causal neighbors and write these to the file
// specified by filename
int Network::compute_activity_threads(SpikeRaster& raster, std::string fname, double gamma_thresh, double temporal_radius, double tau, int func) {

    if (n_neurons < raster.n_neurons) {
        std::cout << "Number of neurons in connectivity file is less than the number of neurons in the raster\n";
        exit(EXIT_FAILURE);
       return -1;
    } 
    std::ofstream outfile;
    outfile.open(fname.c_str(), std::ofstream::out | std::ofstream::trunc);
    if (!outfile.is_open()) {
        std::cout << "Error opening activity thread output file\n";
        exit(EXIT_FAILURE);
        return -1;
    }

    // if event list is empty, do nothing
    idx_t neuron_idx;
    if (raster.n_neurons > 0) {
        for (neuron_idx = 0; neuron_idx < raster.n_neurons; ++neuron_idx) {
            emit_causal_neighbors(raster, neuron_idx, gamma_thresh, temporal_radius, tau, func, outfile);
        }
    }
    outfile.close();
    return 0;
}

// computes all the causal neighbors of each spike emitted by neuron_idx
// outputs directed edges as lines in the ofstream outfile
// each line consists of <presynaptic_neuron_idx> <presynaptic_spike_time> <postsynaptic_neuron_idx> <postsynaptic_spike_time>
void Network::emit_causal_neighbors(SpikeRaster& sr, idx_t neuron_idx, double gamma_thresh, double temporal_radius, double tau, int func, std::ofstream& outfile) {

    // get spike train from this neuron
    std::set<tstamp_t> postsyn_neuron_spikes = sr.evtlist[neuron_idx];

    // for each spike, find all spikes from presynaptic neurons within temporal radius
    //
    std::set<tstamp_t>::iterator curr_spike;
    for (curr_spike = postsyn_neuron_spikes.begin(); curr_spike != postsyn_neuron_spikes.end(); ++curr_spike) {

        std::list<tstamp_t> pre_spikes;

        // past_limit is the earliest spike time we consider
        // clamp to zero so that we don't go past start of recording
        tstamp_t past_limit;
        past_limit = (*curr_spike > temporal_radius) ? (*curr_spike - temporal_radius) : 0;

        // loop over presynaptic neurons
        for (idx_t presyn_idx = 0; presyn_idx < presynaptic_edges[neuron_idx].size(); ++presyn_idx) { 

            real_t weight, delay;
            idx_t presyn_neuron_idx;

            presyn_neuron_idx = presynaptic_edges[neuron_idx][presyn_idx].idx;
            weight = presynaptic_edges[neuron_idx][presyn_idx].weight;
            delay = presynaptic_edges[neuron_idx][presyn_idx].delay; 

            // get all spikes within temporal radius from this presynaptic neuron
            pre_spikes.clear();
            sr.get_spikes_in_range(pre_spikes, presyn_neuron_idx, past_limit, *curr_spike);

            std::list<tstamp_t>::iterator pre_spike;
            if (!pre_spikes.empty()) {

                // loop over all presynaptic spikes from this presynaptic neuron
                for (pre_spike = pre_spikes.begin(); pre_spike != pre_spikes.end(); ++pre_spike) {

                    double g;

                    g = gamma(*pre_spike, *curr_spike, weight, delay, tau);
                    // check for causality
                    if (g <= gamma_thresh && func == GNATS) {
                        // emit edge
                        outfile << presyn_neuron_idx << " " << (tstamp_t)*pre_spike << " " << neuron_idx << " " << (tstamp_t)*curr_spike << "\n";
                    } else if (func == CDH) {
                        outfile << g << "\n";
                    }
                }
            }
        }
    }
}



/**********************************************************/

int main(int argc, char *argv[]) {

    double gamma_thresh = 4;
    double temporal_radius = 100;
    double tau = 5;
    if (argc != 9) {
        std::cout << "usage: " << argv[0] << " <n_neurons> <connection_file> <spike_file> <func> <out_file> <tau> <thresh> <causal_radius>\n";
        std::cout << "func = 1 | Compute GNATS\nfunc = 2 | Compute causal distances\n";
    } else {

        tau = std::stod(argv[6]);
        gamma_thresh = std::stod(argv[7]);
        temporal_radius = std::stod(argv[8]); 

        std::cout << "Reading event file...\n";
        SpikeRaster raster = SpikeRaster(std::stoi(argv[1]));
        raster.read_event_file(argv[3]);

        std::cout << "Reading connectivity file...\n";
        Network net = Network(std::stoi(argv[1]));
        net.read_connectivity(argv[2]);

        std::cout << "Computing activity threads...\n";
        net.compute_activity_threads(raster, argv[5], gamma_thresh, temporal_radius, tau, std::stoi(argv[4]));
        std::cout << "Done\n";
    }

    return 0;
}
