/* ddgraphs.h
 * =========================================================================
 * This file is part of the ddgraphs project
 *
 * Copyright (C) 2011 Universiteit Gent
 *
 * Author: Nicolas Van Cleemput
 * In collaboration with Gunnar Brinkmann
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * A copy of the GNU General Public License can be found in the file
 * LICENSE.txt provided with this distribution. This license can also
 * be found on the GNU website at http://www.gnu.org/licenses/gpl.html.
 *
 * If you did not receive a copy of the GNU General Public License along
 * with this program, contact the lead developer, or write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 */

#ifndef _DDGRAPHS_H
#define	_DDGRAPHS_H

/******************Includes**********************/

#include "util.h"
#include "nausparse.h"

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <ctype.h>

/******************Defines*******************************/

#define SEMIEDGE INT_MAX

short endian = LITTLE_ENDIAN; // defines which endian should be used while exporting pregraph code

struct _buildingblock {
    int id; //normally this is the index of this block in the list

    int type; //1, 2, 3 or 4
    int component;
    int parameter;

    int connectorCount;
    struct _buildingblock **connections; //the blocks to which this block is connected
    int *targetConnector; //the connectors to which this block is connected

    int *connectionVertices;
};

typedef struct _buildingblock BBLOCK;

struct _ddgraph {
    /* This is the true order of the Delaney-Dress graph. It may differ from the
     * order of the nauty graph because that graph includes extra vertices to
     * remove the double edges. The vertices are arranged in such a way that
     * the numbering of the vertices are the same for both graph (i.e. all dummy
     * vertices are at the end of the nauty graph).
     */
    int order;

    /* In this graph nv is equal to order + number of double edges, nde is equal
     * to 2*(number of edges (including double edges) + number of double edges).
     * The array e has length 3*order + order (-1) (i.e. +/- 4*order). The first
     * range of 3*order is subdivided in ranges of 3 corresponding to the
     * 3-regular vertices of the graph. The last range is subdivided into ranges
     * of 2, corresponding to the dummy vertices for the double edges.
     */
    sparsegraph *underlyingGraph;

    /* This array is used to count the number of semi-edges incident with each
     * vertex. This is doen, because semi-edges can't be used in nauty, but we
     * do use the number of semi-edges as a color for the vertices to make sure
     * that they are taken into account for the automorphism group.
     */
    int *semiEdges;

    /* This array is used to store a 1-factor (i.e. the complement of 2-factor)
     * for the Delaney-Dress graph. For each vertex the position in the
     * adjacency list is given of the neighbour that is connected to this vertex
     * through the edge of the 1-factor.
     */
    int *oneFactor;

    int dummyVertexCount;
};

typedef struct _ddgraph DDGRAPH;

//================== Component lists =====================

unsigned int componentListsCount;

#define Q1TypeComponentsCount 12

int maximumQ1TypeComponents;
int Q1TypeComponentsConnectors[Q1TypeComponentsCount];
int Q1TypeComponentsSmallestCase[Q1TypeComponentsCount];
int** Q1TypeComponentsComponentCount;

#define Q2TypeComponentsCount 2

int maximumQ2TypeComponents;
int Q2TypeComponentsConnectors[Q2TypeComponentsCount];
int Q2TypeComponentsSmallestCase[Q2TypeComponentsCount];
int** Q2TypeComponentsComponentCount;

#define Q3TypeComponentsCount 2

int maximumQ3TypeComponents;
int Q3TypeComponentsConnectors[Q2TypeComponentsCount];
int Q3TypeComponentsSmallestCase[Q2TypeComponentsCount];
int** Q3TypeComponentsComponentCount;

int Q4ComponentCount;

#define ComponentsTypesCount (Q1TypeComponentsCount + Q2TypeComponentsCount + Q3TypeComponentsCount + 1)

void (*constructBlock[ComponentsTypesCount])(int *, BBLOCK *, DDGRAPH *, int *, int *) = {NULL};

void (*storeBlockAutomorphismGenerators[ComponentsTypesCount])(BBLOCK *, DDGRAPH *) = {NULL};

void (*storeBlocksMapping[ComponentsTypesCount])(BBLOCK *, BBLOCK *, DDGRAPH *) = {NULL};

/******************Global Variables**********************/

int connectionsMade; //the number of connections made at this point
                     //the maximum value is vertexCount/2

int connections[MAXN/2][2];

unsigned int graphsCount;

boolean markedTwoFactors = FALSE;

/* Provide space for the generators at each recursion depth (maximum depth = MAXN/2)
 * There are at most n<=MAXN generators in a graph with n vertices and the length
 * of each generator is n.
 */
permutation automorphismGroupGenerators[MAXN/2][MAXN][MAXN];
int numberOfGenerators[MAXN/2];

int vertexOrbits[MAXN/2][MAXN];
int vertexOrbitsSizes[MAXN/2][MAXN];

/* Variables for nauty */
int nautyLabelling[MAXN], nautyPtn[MAXN];
DEFAULTOPTIONS_SPARSEGRAPH(nautyOptions);
statsblk nautyStats;
setword nautyWorkspace[50 * MAXM];
sparsegraph canonGraph;

/******************Methods*******************************/

void findNextOrbitToConnect(BBLOCK* blocks, int buildingBlockCount, DDGRAPH *ddgraph,
        int *vertexToBlock, int *vertexToConnector, boolean *freeConnectors);

#endif	/* _DDGRAPHS_H */

