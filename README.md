gnatfinder generates the second-order causal activity graph from a list of spike times and network connectivity.

Weakly connected components of the second-order activity graph correspond to recurrances of causal patterns in the spiking activity.

## Usage:
`gnatfinder <N cells> <activity file> <network file>`

The activity file is a text file containing spikes sorted in time.

Each line is a spike and has the format:
`<type> <timestamp> <neuron id>`

`type` is currently ignored, but could be used to tag different classes of events in the future.

`timestamp` is the time of the spike in hexadecimal.  Usually expressed in units of milliseconds. Integers only!

`neuron id` is the integer identifying the neuron producing the spike.

The network file is a text file listing connectivity.

Each line is a synapse.  The line format is 
`<src_id> <tgt_id> <rel_w> <delay>`

`src_id` is the id of the presynaptic neuron

`tgt_id` is the id of the postsynaptic neuron

`rel_w` is the (float) relative weight of the synaptic connection

`delay` is the conduction delay along this connection, in the same units as timestamps in the activity file.

The output file is by default "./gnat_output.txt".  It is a text file. Each line is an edge int he second order graph.  

The line format is:
`<presyn_id> <time_1> <time_2> <postsyn_id> <time_1> <time_2>`

`presyn_id` is the id of the presynaptic neuron for this spike pair.

`postsyn_id` is the id of the postsynaptic neuron for this spike pair.

`time_1` is the timestamp of the earlier presynaptic spike of the repeated interaction

`time_2` is the timestamp of the later pre/post synaptic spike of the recurring interaction


## Compilation
To compile gnatfinder, use the command:
`gcc -o gnatfinder gnatfinder.c quadtree.c raster.c network.c gnats.c -Wall -Wextra -g -lm`

To compile the first order gnatfinder:
`g++ -std=c++11 -o gnat1 compute_activity_threads.cpp`

No other libraries besides the math library are needed for now. 

First order GNAT invocation:
`<progname> <n_neurons> <connection file> <activity file> <function> <output file>`

function = 1 to compute GNATs

function = 2 to compute causal distances (for histogram)
