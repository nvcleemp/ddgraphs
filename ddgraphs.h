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

#define SEMIEDGE INT_MAX;

struct _buildingblock {
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
    int order;

    //an adjacency list of the graph
    int **adjList;

};

typedef struct _ddgraph DDGRAPH;

//================== Component lists =====================

unsigned int componentListsCount;

#define Q1TypeComponentsCount 11

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

/******************Global Variables**********************/

int connectionsMade; //the number of connections made at this point
                     //the maximum value is vertexCount/2

/* Provide space for the generators at each recursion depth (maximum depth = MAXN/2)
 * There are at most n<=MAXN generators in a graph with n vertices and the length
 * of each generator is n.
 */
permutation automorphismGroupGenerators[MAXN/2][MAXN][MAXN];
int numberOfGenerators[MAXN/2];
/******************Methods*******************************/


#endif	/* _DDGRAPHS_H */

