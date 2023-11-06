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
#include <math.h>
#include <string.h>

#include "network.h"

#define SYN_LINE_MAX 512

/*
 * Network implementation for GNATFinder
 * Brad Theilman 2022-11-02
 *
 */

/*
 * Initialize a PhysNetwork structure by allocating memory for the 
 * synapse lists.
 */
int PhysNetworkInit(struct PhysNetwork *pn, unsigned long _n_cells) {

    pn->n_cells = _n_cells;
    pn->presyns = malloc(_n_cells * sizeof(struct Synapse *));
    if (!pn->presyns) {
        printf("FATAL: Unable to allocate space for synapse lists\n");
        exit(-1);
    }
    return 0;

}

void PhysNetworkAddSynapse(struct PhysNetwork *pn, struct Synapse *edg) {

    if (edg->tgt_id >= pn->n_cells) {
        printf("FATAL: Trying to add synapse onto a cell outside of the network population.\n");
        exit(-1);
    }

    if (!pn->presyns[edg->tgt_id]) {
        pn->presyns[edg->tgt_id] = edg;
    } else {
        edg->next = pn->presyns[edg->tgt_id];
        pn->presyns[edg->tgt_id] = edg;
    }
}

struct Synapse *SynapseCreate(unsigned long _src, unsigned long _tgt, float _rel_w, float delay) {

    float _nl_rel_w;

    struct Synapse *res = malloc(sizeof(struct Synapse));
    if (!res) {
        printf("FATAL: Cannot allocate synapse\n");
        exit(-1);
    }

    res->src_id = _src;
    res->tgt_id = _tgt;
    res->rel_w = _rel_w;
    res->delay = delay;
    res->next = (struct Synapse *)NULL;

    _nl_rel_w = -1*log(_rel_w);
    res->neg_log_rel_w = _nl_rel_w;

    return res;
}

void PhysNetworkReadFile(struct PhysNetwork *pn, char *fname) {

    /* 
     * File format:
     * <src_id> <tgt_id> <rel_w> <delay>
     */

    /* Attempt to open file */
    FILE *syn_file;
    syn_file = fopen(fname, "r");
    if (!syn_file) {
        printf("FATAL: Unable to open synapse file %s\n", fname);
        exit(-1);
    }

    /* for each line create a synapse */
    char line[SYN_LINE_MAX];
    char *field, *end;

    unsigned long src_id, tgt_id;
    float rel_w, delay;

    struct Synapse *syn;

    while (fgets(line, SYN_LINE_MAX, syn_file)) {
        field = strtok(line, " ");
        src_id = strtol(field, &end, 0);
        if (end == field) {
            printf("FATAL: Unable to parse source neuron\n");
            exit(-1);
        }

        field = strtok(NULL, " ");
        tgt_id = strtol(field, &end, 0);
        if (end == field) {
            printf("FATAL: Unable to parse target neuron\n");
            exit(-1);
        }

        field = strtok(NULL, " ");
        rel_w = strtof(field, &end);
        if (end == field) {
            printf("FATAL: Unable to parse relative weight\n");
            exit(-1);
        }

        field = strtok(NULL, " ");
        delay = strtof(field, &end);
        if (end == field) {
            printf("FATAL: Unable to parse delay\n");
            exit(-1);
        }

        syn = SynapseCreate(src_id, tgt_id, rel_w, delay);
        PhysNetworkAddSynapse(pn, syn);
    }

    fclose(syn_file);
}

void SynapsePrint(struct Synapse *syn) {

    if (!syn) return;

    printf("%ld --> %ld [%.2f, %.2f]\n", syn->src_id, syn->tgt_id, syn->rel_w, syn->delay);

}

void PhysNetworkPrint(struct PhysNetwork *pn) {

    unsigned long idx;
    struct Synapse *syn;

    if (!pn) return;

    for (idx=0; idx < pn->n_cells; ++idx) {

        syn = pn->presyns[idx];
        while (syn) {
            SynapsePrint(syn);
            syn = syn->next;
        }
    }

}

