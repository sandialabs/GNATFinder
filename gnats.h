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


#ifndef GNATS_H
#define GNATS_H
#define N_EDGBUF 8192 /* size of the GNAT edge buffer */

struct GNATEdge {

    struct SpikePair *spp_pre;
    struct SpikePair *spp_post;
    float cd_ratio; /* causal distance ratio */
};

void finalize_edge_buffer();
int initialize_edge_buffer(char* fname);
void QTreeMapGNATEdge(struct QuadTree *qt, struct BoundingBox *r, struct SpikePair *spp_post, struct Synapse *syn, float tau, float theta);
int GNAT_test_for_edge(struct SpikePair *spp_pre, struct SpikePair *spp_post, struct Synapse *edg, float tau, float thresh);
float compute_gamma(struct Spike *sp_pre, struct Spike *sp_post, struct Synapse *edg, float tau);
float compute_omega(struct Spike *sp_pre, struct Spike *sp_post, struct Synapse *edg, float tau);
void flush_edge_buffer();
void fprint_GNAT_edge(FILE *fp, struct GNATEdge *edg);

#endif
