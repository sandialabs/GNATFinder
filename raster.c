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


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "quadtree.h"
#include "raster.h"

#define LINE_MAX 256

/*
 * Raster Routines
 * Brad Theilman 2022-10-31
 */

/*
 * Initializes storage of spike lists for a raster containing _n_cells
 */
int RasterInit(struct SpikeRaster *sr, const unsigned int _n_cells) {

    sr->n_cells = _n_cells;
    sr->sp_lists = (struct Spike **)malloc(_n_cells * sizeof(struct Spike *));
    if (!sr->sp_lists) {
        printf("FATAL:  Unable to Allocate Spike Lists for Raster\n");
        exit(-1);
    }
    sr->t_min = 0;
    sr->t_max = 0;
    sr->n_spikes = 0;
    return 0;
}

/*
 * Adds a spike to the head of the cell's spike list in the raster
 */
void RasterHeadAppend(struct SpikeRaster *sr, struct Spike *sp) {

    if (sp->n_id >= sr->n_cells) {
        printf("FATAL: Attempting to add spike from neuron outside of raster population\n");
        exit(-1);
    }

    sp->next = sr->sp_lists[sp->n_id];
    sr->sp_lists[sp->n_id] = sp;
    
    /* update t_max and t_min */
    long min_diff, max_diff;
    if (!sr->n_spikes) {
        sr->t_min = sp->ts;
        sr->t_max = sp->ts;
    } else {
        min_diff = sp->ts - sr->t_min;
        max_diff = sp->ts - sr->t_max;
        sr->t_min += (min_diff > 0) ? 0 : min_diff;
        sr->t_max += (max_diff < 0) ? 0 : max_diff;
    }
    sr->n_spikes++;


}

/*
 * Reverses a list of spikes; returns the new list head
 */
struct Spike *SpikeListReverse(struct Spike *list_head) {

    struct Spike *prev, *cur, *next;

    if (!list_head) return list_head;

    prev = (struct Spike *)NULL;
    cur = list_head;
    next = list_head->next;

    while (next) {

        cur->next = prev;
        prev = cur;
        cur = next;
        next = next->next;
    }
    cur->next = prev;

    return cur; /* new list head */
        
}

/*
 * Reverses all the spike lists in a raster
 */
void RasterReverse(struct SpikeRaster *sr) {

    unsigned int idx;
    struct Spike *new_head;
    for (idx = 0; idx < sr->n_cells; ++idx) {
        new_head = SpikeListReverse(sr->sp_lists[idx]);
        sr->sp_lists[idx] = new_head;
    }
}

/*
 * Reads spikes from a file into raster sr 
 * Assumes raster has been initialized
 * Also assumes that the file contains spikes in time sorted order
 */
void RasterReadFile(struct SpikeRaster *sr, const char *fname) {

    /* open file */
    FILE *sp_file;

    sp_file = fopen(fname, "r");
    if (!sp_file) {
        printf("FATAL: Could not open spike file %s\n", fname);
        exit(-1);
    }

    /* for each line create a spike */
    char line[LINE_MAX];
    char *field, *end;
    long _ts, _n_id, _sp_type;
    struct Spike *sp;
    while (fgets(line, LINE_MAX, sp_file)) {
        /*
         * Line format is <type> <timestamp> <neuron_id>
         */
        field = strtok(line, " ");
        _sp_type = strtol(field, &end, 0);
        if (end == field) {
            printf("FATAL: Unable to parse spike type\n");
            exit(-1);
        }

        field = strtok(NULL, " ");
        _ts = strtol(field, &end, 16);
        if (end == field) {
            printf("FATAL: Unable to parse timestamp\n");
            exit(-1);
        }

        field = strtok(NULL, " ");
        _n_id = strtol(field, &end, 0);
        if (end == field) {
            printf("FATAL: Unable to parse neuron id\n");
            exit(-1);
        }

        sp = create_spike(_n_id, _ts);
        /* add spike to raster */
        RasterHeadAppend(sr, sp);
    }

    /* done with the file */
    fclose(sp_file);

    /* after all spikes added, reverse the raster */
    RasterReverse(sr);
}

void RasterPrint(struct SpikeRaster *sr) {

    /*
     * For each cell, print the cell's spike train stored in 
     * the raster
     */

    unsigned int idx;
    struct Spike *sp;

    if (!sr) return;

    printf("------ Spike Raster ------\n");
    for (idx = 0; idx < sr->n_cells; ++idx) {
        sp = sr->sp_lists[idx];
        printf("Cell %d\n", idx);
        while (sp) {
            print_spike(sp);
            printf("\n");
            sp = sp->next;
        }
    }
    printf("------ End Spike Raster ------\n");

}
