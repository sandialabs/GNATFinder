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


#ifndef SPIKE_QT_H
#define SPIKE_QT_H

#include <stdint.h>

#define QT_MAX_CAP 4
#define FALSE 0
#define TRUE  1 
/* Data structure definitions */

struct Spike {

    uint32_t n_id; /* neuron id */
    long     ts;   /* timestamp */

    struct Spike *next;

};

struct SpikePair {
    struct Spike *sp1;
    struct Spike *sp2;

    struct SpikePair *prev;
    struct SpikePair *next;
};


struct BoundingBox {

    float   c_x;   /* center x */
    float   c_y;   /* center y */
    float   w2;    /* half-width */

};

struct QuadTree {

    int capacity;
    struct BoundingBox *bdry; /* boundary box */

    struct SpikePair *pairs;

    struct QuadTree *NW;
    struct QuadTree *SW;
    struct QuadTree *NE;
    struct QuadTree *SE;

};

/* prototypes */

struct Spike *create_spike(uint32_t neuron_id, long timestamp);
void          destroy_spike(struct Spike *sp);
int           spike_equals(struct Spike *sp1, struct Spike *sp2);
void          print_spike(struct Spike *sp);


struct SpikePair *create_spike_pair(struct Spike *sp1, struct Spike *sp2);
void              destroy_spike_pair(struct SpikePair *spp);
void              print_spike_pair(struct SpikePair *spp);

struct BoundingBox *BBoxCreate(float center_x, float center_y, float half_width);
void                BBoxDestroy(struct BoundingBox *bb);
int                 BBoxContainsPoint(struct BoundingBox *bb, struct SpikePair *spp);
int                 BBoxIntersects(struct BoundingBox *bb1, struct BoundingBox *bb2);


struct QuadTree *QTreeCreate(struct BoundingBox *bb);
int QTreeInsert(struct QuadTree *qt, struct SpikePair *spp);
void QTreeMapQueryRange(struct QuadTree *qt, struct BoundingBox *r, void (*func)(struct SpikePair *));
void QTreePrint(struct QuadTree *qt);

#endif
