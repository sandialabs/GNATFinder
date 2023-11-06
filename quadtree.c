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
 * Quadtree routines
 * Brad Theilman 2023
 */

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "quadtree.h"

struct Spike *create_spike(uint32_t neuron_id, long timestamp) {

    struct Spike *res = malloc(sizeof(struct Spike));
    if (res == NULL) {
        printf("FATAL: Unable to allocate Spike\n");
        exit(1);
    }

    res->n_id = neuron_id;
    res->ts   = timestamp;
    res->next = (struct Spike *) NULL;
    return res;

}

void destroy_spike(struct Spike *sp) {

    free(sp);

}

int spike_equals(struct Spike *sp1, struct Spike *sp2) {

    int p1, p2;
    p1 = sp1->n_id == sp2->n_id;
    p2 = sp1->ts == sp2->ts;
    return p1 && p2;
}

void print_spike(struct Spike *sp) {

    if (!sp) return;

    printf("Spike[%d, %ld]", sp->n_id, sp->ts);
}

/* SpikePair routines */

struct SpikePair *create_spike_pair(struct Spike *_sp1, struct Spike *_sp2) {

    struct SpikePair *res = malloc(sizeof(struct SpikePair));
    if (res == NULL) {
        printf("FATAL: Unable to allocate SpikePair\n");
        exit(1);
    }

    /* check that the spikes belong to the same cell */
    if (_sp1->n_id != _sp2->n_id) {
        printf("WARNING: Creating spike pair from spikes from two different cells!\n");
    }

    /* check if the spikes are equal */
    /* at this point we know the cell ids are the same */
    /* so just check the timestamps */
    if (_sp1->ts == _sp2->ts) {
        printf("WARNING: Creating spike pair from identical spikes!\n");
    }

    /* go ahead and fill in the information */
    res->sp1 = _sp1;
    res->sp2 = _sp2;
    res->next = (struct SpikePair *) NULL;
    res->prev = (struct SpikePair *) NULL;
    return res;

}

void destroy_spike_pair(struct SpikePair *spp) {

    /* also destroys the spikes contained */
    if (spp->sp1) {
        destroy_spike(spp->sp1);
    }

    if (spp->sp2) {
        destroy_spike(spp->sp2);
    }

    free(spp);
}

void print_spike_pair(struct SpikePair *spp) {

    if (spp->sp1) {
        print_spike(spp->sp1);
    } else {
        printf("[NULL SPIKE]");
    }

    printf (" <---> ");

    if (spp->sp2) {
        print_spike(spp->sp2);
    } else {
        printf("[NULL SPIKE]");
    }

    printf("\n");

}


/* BoundingBox routines */

struct BoundingBox *BBoxCreate(float center_x, float center_y, float half_width) {

    struct BoundingBox *res = malloc(sizeof(struct BoundingBox));
    if (res == NULL) {
        printf("FATAL: Unable to allocate Bounding Box\n");
        exit(1);
    }

    res->c_x = center_x;
    res->c_y = center_y;
    res->w2  = half_width;

    return res;

}

void BBoxDestroy(struct BoundingBox *bb) {

    free(bb);
}

int BBoxContainsPoint(struct BoundingBox *bb, struct SpikePair *spp) {

    int s1, s2;

    s1 = (abs(spp->sp1->ts - bb->c_x) < bb->w2);
    s2 = (abs(spp->sp2->ts - bb->c_y) < bb->w2);
    return s1 && s2;

}

int BBoxIntersects(struct BoundingBox *bb1, struct BoundingBox *bb2) {

    float d;
    int b1, b2;

    d = bb1->w2 + bb2->w2;
    b1 = (abs(bb2->c_x - bb1->c_x) <= d);
    b2 = (abs(bb2->c_y - bb1->c_y) <= d);
    return b1 && b2;

}



/* QuadTree routines */
static void QTreeSubdivide(struct QuadTree *qt);

struct QuadTree *QTreeCreate(struct BoundingBox *bbox) {

    struct QuadTree *res = malloc(sizeof(struct QuadTree));

    if (res == NULL) {
        printf("FATAL: Unable to allocate QuadTree\n");
        exit(1);
    }

    res->capacity = 0;
    res->bdry = bbox;
    res->pairs = (struct SpikePair *) NULL;

    res->NW = (struct QuadTree *) NULL;
    res->SW = (struct QuadTree *) NULL;
    res->NE = (struct QuadTree *) NULL;
    res->SE = (struct QuadTree *) NULL;

    return res;
}

void QTreeSubdivide(struct QuadTree *qt) {

    /* build four new bounding boxes */
    float d2;
    struct BoundingBox *bbNW;
    struct BoundingBox *bbSW;
    struct BoundingBox *bbNE;
    struct BoundingBox *bbSE;

    d2 = qt->bdry->w2 / 2;

    bbNW = BBoxCreate(qt->bdry->c_x - d2, qt->bdry->c_y + d2, d2);
    bbSW = BBoxCreate(qt->bdry->c_x - d2, qt->bdry->c_y - d2, d2);
    bbNE = BBoxCreate(qt->bdry->c_x + d2, qt->bdry->c_y + d2, d2);
    bbSE = BBoxCreate(qt->bdry->c_x + d2, qt->bdry->c_y - d2, d2);

    qt->NW = QTreeCreate(bbNW);
    qt->SW = QTreeCreate(bbSW);
    qt->NE = QTreeCreate(bbNE);
    qt->SE = QTreeCreate(bbSE);

    while (qt->capacity > 0) {

        /* pop off the top of the spike pair list */
        struct SpikePair *spp = qt->pairs;
        qt->pairs = spp->next;
        if (qt->pairs) {
            qt->pairs->prev = (struct SpikePair *) NULL;
        }
        qt->capacity--;

        spp->next = (struct SpikePair *) NULL;
        spp->prev = (struct SpikePair *) NULL;

        if (QTreeInsert(qt->NW, spp)) continue;
        if (QTreeInsert(qt->SW, spp)) continue;
        if (QTreeInsert(qt->NE, spp)) continue;
        if (QTreeInsert(qt->SE, spp)) continue;

    }

}


int QTreeInsert(struct QuadTree *qt, struct SpikePair *spp) {

    /* if spike pair is not in our bounding box, return FALSE */
    if (!BBoxContainsPoint(qt->bdry, spp)) {
        return FALSE;
    }

    /* if we are not full, add spike pair to list */
    if ((qt->capacity < QT_MAX_CAP) && (qt->NW == NULL)) {
        if (!(qt->capacity)) {
            /* list empty */
            spp->next = (struct SpikePair *) NULL;
            spp->prev = (struct SpikePair *) NULL;
            qt->pairs = spp;
        }

        struct SpikePair *ll = qt->pairs;
        while(ll->next != NULL) {
            ll = ll->next;
        }
        ll->next = spp;
        spp->prev = ll;
        spp->next = (struct SpikePair *) NULL;
        qt->capacity++;
        return TRUE;
    }

    /* If we haven't been subdivided yet, do subdivision */
    if (!qt->NW) {
        QTreeSubdivide(qt);
    }

    /* Insert new spike pair into quad tree */
    if (QTreeInsert(qt->NW, spp)) {
        return TRUE;
    }

    if (QTreeInsert(qt->SW, spp)) {
        return TRUE;
    }

    if (QTreeInsert(qt->NE, spp)) {
        return TRUE;
    }

    if (QTreeInsert(qt->SE, spp)) {
        return TRUE;
    }

    /* shouldn't reach here */
    return FALSE;

}

void QTreeMapQueryRange(struct QuadTree *qt, struct BoundingBox *r, void (*func)(struct SpikePair *)) {

    /* 
     * Maps the function func to all elements in the QuadTree qt
     * that fall within the bounding box r
     */

    /* If the region does not intersect our BBox, return */
    if (!BBoxIntersects(qt->bdry, r)) return;

    struct SpikePair *spp = qt->pairs;

    while (spp) {

        /* apply func to the spike pair */
        func(spp);
        spp = spp->next;
    }

    if (!qt->NW) return;

    QTreeMapQueryRange(qt->NW, r, func);
    QTreeMapQueryRange(qt->SW, r, func);
    QTreeMapQueryRange(qt->NE, r, func);
    QTreeMapQueryRange(qt->SE, r, func);

}


void QTreePrint(struct QuadTree *qt) {

    if (!qt) return;

    if (qt->capacity > 0) {
        struct SpikePair *spp;
        spp = qt->pairs;
        while (spp) {
            print_spike_pair(spp);
            spp = spp->next;
        }
    }

    if (qt->NW) {

        QTreePrint(qt->NW);
        QTreePrint(qt->SW);
        QTreePrint(qt->NE);
        QTreePrint(qt->SE);
    }
}
