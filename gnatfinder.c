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
 * GNATFinder main 
 * Brad Theilman 2023
 */

#include <stdlib.h>
#include <stdio.h>

#include "quadtree.h"
#include "raster.h"
#include "network.h"
#include "gnats.h"

/* global spike raster */
struct SpikeRaster g_raster;

/* global network */
struct PhysNetwork g_network;

/* global neuron quadtree array */
struct QuadTree **g_qtarray;

void compute_gnat_edges(float tau, float thresh, float c_radius) {

    struct BoundingBox query_bbox;
    struct QuadTree *presyn_qtree;
    struct Synapse *presyn;
    struct Spike *sp_a, *sp_b;
    struct SpikePair *spp_post;

    unsigned long tgt_id;

    unsigned int post_idx;
    for (post_idx = 0; post_idx < g_network.n_cells; ++post_idx) {

        /* print status */
        if ((post_idx % 10) == 0) {
            printf("Cell %d of %d\n", post_idx, g_network.n_cells);
        }

        /* iterate over spike pairs in post qtree */
        sp_a = g_raster.sp_lists[post_idx];
        while (sp_a) {
            sp_b = sp_a->next;
            while(sp_b) {
                if (!spike_equals(sp_a, sp_b)) {
                    spp_post = create_spike_pair(sp_a, sp_b);
                    //print_spike_pair(spp_post);
                    tgt_id = post_idx;
                    /* list of presynaptic partners */
                    presyn = g_network.presyns[tgt_id];

                    while (presyn) {
                        /* quadtree associated to presynaptic neuron */
                        presyn_qtree = g_qtarray[presyn->src_id];

                        /* set query bounding box */
                        query_bbox.c_x = spp_post->sp1->ts;
                        query_bbox.c_y = spp_post->sp2->ts;
                        query_bbox.w2  = c_radius;

                        /* apply edge test to queried range */
                        QTreeMapGNATEdge(presyn_qtree, &query_bbox, spp_post, presyn, tau, thresh);
                        presyn = presyn->next;

                    } 
                }
                sp_b = sp_b->next;
            }
            sp_a = sp_a->next;
        }
    }
}


void insert_spike_pairs (struct QuadTree *qt, struct Spike *list_head) {

    struct Spike *sp_a, *sp_b;
    struct SpikePair *spp;

    sp_a = list_head;

    while(sp_a) {

        sp_b = sp_a->next;
        while(sp_b) {
            if (!spike_equals(sp_a, sp_b)) {
                spp = create_spike_pair(sp_a, sp_b);
                QTreeInsert(qt, spp);
            }
            sp_b = sp_b->next;
        }
        sp_a = sp_a->next;
    }
}

int main(int argc, char **argv) {

    float _cx, _cy, _hw;
    float tau, thresh, c_radius;
    unsigned long _n_cells;
    struct BoundingBox *bbox_top_level;

    /* check num args */
    if (argc < 7) {
        printf("Usage: %s <N cells> <spike file> <network file> <tau> <thresh> <causal_radius>\n", argv[0]);
        exit(-1);
    }

    _n_cells = strtol(argv[1], NULL, 0);
    tau = strtof(argv[4], NULL);
    thresh = strtof(argv[5], NULL);
    c_radius = strtof(argv[6], NULL);

    if (RasterInit(&g_raster, _n_cells)) {
        printf("Problem initializing raster\n");
    }

    if (PhysNetworkInit(&g_network, _n_cells)) {
        printf("Problem initializing network\n");
    }

    /* Attempt to allocate space for each quadtree */
    g_qtarray = malloc(_n_cells * sizeof(struct QuadTree *));
    if (!g_qtarray) {
        printf("FATAL: Unable to allocate space for neuron quadtrees\n");
        exit(-1);
    }

    /* Read spikes from file into global raster */
    RasterReadFile(&g_raster, argv[2]);

    /* Attempt to read network connectivity file */
    PhysNetworkReadFile(&g_network, argv[3]);
    //PhysNetworkPrint(&g_network);

    /* build top-level bouding box */
    _cx = (float)(g_raster.t_max + g_raster.t_min)/2;
    _cy = _cx;
    _hw = (float)(g_raster.t_max - g_raster.t_min)/2;
    bbox_top_level = BBoxCreate(_cx, _cy, _hw);

    /* build quadtrees for each cell */
    unsigned int idx;
    for (idx = 0; idx < _n_cells; ++idx) {
        g_qtarray[idx] = QTreeCreate(bbox_top_level);
        insert_spike_pairs(g_qtarray[idx], g_raster.sp_lists[idx]);
#ifdef SPDEBUG
        printf("-------- QuadTree --------\n");
        QTreePrint(g_qtarray[idx]);
        printf("-------- End QuadTree --------\n");
#endif
    }

    /* initialize output file */
    initialize_edge_buffer("gnat2_out.txt");

    /* compute gnats here */
    compute_gnat_edges(tau, thresh, c_radius);

    /* clean up */
    finalize_edge_buffer();

    return 0;
}
