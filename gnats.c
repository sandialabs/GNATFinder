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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "quadtree.h"
#include "network.h"
#include "gnats.h"


#define LARGE_GAMMA 999999

/*
 * Second-order activity graph computation 
 * Brad Theilman 2022-11-02
 */


/* 
 * GNAT Edge buffer
 * Buffers activity graph edges before writing to disk
 */

static struct GNATEdge g_edgbuf[N_EDGBUF];
static unsigned long edgbuf_sz = 0; /* number of edges currently in buffer */

static FILE *fp_edgbuf;

void fprint_GNAT_edge(FILE *fp, struct GNATEdge *edg) {


    unsigned long n_id_1, n_id_2;
    long t_11, t_12, t_21, t_22;

    n_id_1 = edg->spp_pre->sp1->n_id;
    n_id_2 = edg->spp_post->sp1->n_id;

    t_11 = edg->spp_pre->sp1->ts;
    t_12 = edg->spp_pre->sp2->ts;

    t_21 = edg->spp_post->sp1->ts;
    t_22 = edg->spp_post->sp2->ts;

    /* format is <pre neuron id> <spike time 1> <spike time 2> <post neuron id> <spike time 1> <spike time 2>*/
    fprintf(fp, "%ld %ld %ld %ld %ld %ld\n", n_id_1, t_11, t_12, n_id_2, t_21, t_22);
}

/*
 * Opens file fname for writing
 * Zeros edge buffer
 */
int initialize_edge_buffer(char* fname) {

    /* attempt to open file */
    fp_edgbuf = fopen(fname, "w");
    if (!fp_edgbuf) {
        printf("FATAL: unable to open output file %s\n", fname);
        exit(-1);
    }

    /* clear edge buffer */
    edgbuf_sz = 0;
    memset(g_edgbuf, 0, sizeof(g_edgbuf));
    return 0;

}

void finalize_edge_buffer() {

    flush_edge_buffer();
    fclose(fp_edgbuf);
}

void flush_edge_buffer() {

    unsigned int idx;

    if (!fp_edgbuf) {
        printf("FATAL: output file not initialized\n");
        exit(-1);
    }

    if (edgbuf_sz == 0) return;

    for (idx = 0; idx < edgbuf_sz; ++idx) {
        fprint_GNAT_edge(fp_edgbuf, &g_edgbuf[idx]);
    }
    edgbuf_sz = 0;
}


static void GNAT_add_edge(struct SpikePair *_spp_pre, struct SpikePair *_spp_post, float _cd_ratio) {


    /* Check if buffer is full */
    if (edgbuf_sz >= N_EDGBUF) {
        flush_edge_buffer();
        edgbuf_sz = 0;
    }

    g_edgbuf[edgbuf_sz].spp_pre = _spp_pre;
    g_edgbuf[edgbuf_sz].spp_post = _spp_post;
    g_edgbuf[edgbuf_sz].cd_ratio = _cd_ratio;
    edgbuf_sz++;

}

float compute_omega(struct Spike *sp_pre, struct Spike *sp_post, struct Synapse *edg, float tau) {

    float omega, theta;
    float delta_t;

    delta_t = (float)(sp_post->ts - sp_pre->ts);

    /* heaviside function */
    theta = (delta_t >= edg->delay) ? 1 : 0;

    omega = theta * (edg->rel_w) * exp(-(delta_t - edg->delay) / tau);
    return omega;
}

/*
 * Computes gamma for a pair of spikes from neurons connected 
 * via synapse edg.
 *
 * The synapse struct stores the precomputed negative log of the relative weight.
 * That way this function requires only a couple adds and one division.
 * No hefty logs or exps needed.
 */
float compute_gamma(struct Spike *sp_pre, struct Spike *sp_post, struct Synapse *edg, float tau) {

    float gamma, theta;
    float delta_t;

    delta_t = (float)(sp_post->ts - sp_pre->ts);

    /* heaviside function */
    theta = (delta_t >= edg->delay) ? 0 : 1;

    gamma = (theta*LARGE_GAMMA) + (edg->neg_log_rel_w + (delta_t - edg->delay)/tau);
    return gamma;
}

int GNAT_test_for_edge(struct SpikePair *spp_pre, struct SpikePair *spp_post, struct Synapse *edg, float tau, float thresh) {

    float gamma_1, gamma_2;

    gamma_1 = compute_gamma(spp_pre->sp1, spp_post->sp1, edg, tau);
    gamma_2 = compute_gamma(spp_pre->sp2, spp_post->sp2, edg, tau);

    return (gamma_1 <= thresh) && (gamma_2 <= thresh);
}


void QTreeMapGNATEdge(struct QuadTree *qt, struct BoundingBox *r, struct SpikePair *spp_post, struct Synapse *syn, float tau, float theta) {

    /* 
     * Maps the function func to all elements in the QuadTree qt
     * that fall within the bounding box r
     */

    /* If the region does not intersect our BBox, return */
    if (!BBoxIntersects(qt->bdry, r)) return;

    struct SpikePair *spp_pre = qt->pairs;

    while (spp_pre) {

        /* apply func to the spike pair */
        if (GNAT_test_for_edge(spp_pre, spp_post, syn, tau, theta)) {
            /* add edge */
            GNAT_add_edge(spp_pre, spp_post, 1);
        }
        spp_pre = spp_pre->next;
    }

    if (!qt->NW) return;

    QTreeMapGNATEdge(qt->NW, r, spp_post, syn, tau, theta);
    QTreeMapGNATEdge(qt->SW, r, spp_post, syn, tau, theta);
    QTreeMapGNATEdge(qt->NE, r, spp_post, syn, tau, theta);
    QTreeMapGNATEdge(qt->SE, r, spp_post, syn, tau, theta);

}
