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


#ifndef NETWORK_H
#define NETWORK_H

struct Synapse {

    unsigned long src_id; /* presynaptic id */
    unsigned long tgt_id; /* postsynaptic id */

    float rel_w; /* relative weight */
    float neg_log_rel_w; /* negative log of the relative weight */
    float delay; /* axonal conduction delay */

    struct Synapse *next;
};

struct PhysNetwork {

    unsigned long n_cells;
    struct Synapse **presyns; /* list of lists of presynaptic partners */

};

/* Network api */

int PhysNetworkInit(struct PhysNetwork *pn, unsigned long _n_cells);
void PhysNetworkAddSynapse(struct PhysNetwork *pn, struct Synapse *edg);
void PhysNetworkReadFile(struct PhysNetwork *pn, char *fname);
void PhysNetworkPrint(struct PhysNetwork *pn);

struct Synapse *SynapseCreate(unsigned long _src, unsigned long _tgt, float _rel_w, float delay);

void SynapsePrint(struct Synapse *syn);



#endif
