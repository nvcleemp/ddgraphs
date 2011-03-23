/* ddgraphs.c
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

//#define _DEBUG
//#define _PROFILING

#include "ddgraphs.h"

#include <sys/times.h>

//========================== Utility methods ================================

DDGRAPH *getNewDDGraph(int order){
    int i;
    
    DDGRAPH *ddgraph = (DDGRAPH *)malloc(sizeof(DDGRAPH));
    ddgraph->order = order;
    ddgraph->dummyVertexCount = 0;
    ddgraph->semiEdges = (int *)malloc(sizeof(int)*order);
    for(i = 0; i < order; i++){
        ddgraph->semiEdges[i] = 0;
    }
    ddgraph->underlyingGraph = (sparsegraph *)malloc(sizeof(sparsegraph));
    if(order==2){
        //order 2 is a bit special because we need to store the theta-graph
        ddgraph->underlyingGraph->nv = 0;
        ddgraph->underlyingGraph->nde = 0;
        ddgraph->underlyingGraph->d = (int *)malloc(4*sizeof(int));
        ddgraph->underlyingGraph->v = (int *)malloc(4*sizeof(int));
        ddgraph->underlyingGraph->e = (int *)malloc(10*sizeof(int));
        for(i = 0; i < order; i++){
            ddgraph->underlyingGraph->d[i]=3;
            ddgraph->underlyingGraph->v[i]=3*i;
        }
        for(i = order; i < 4; i++){
            ddgraph->underlyingGraph->d[i]=2;
            ddgraph->underlyingGraph->v[i]=3*order + 2*i;
        }
        ddgraph->underlyingGraph->dlen = 4;
        ddgraph->underlyingGraph->vlen = 4;
        ddgraph->underlyingGraph->elen = 10;
    } else {
        int maxVertices = order + order/2;
        ddgraph->underlyingGraph->nv = 0;
        ddgraph->underlyingGraph->nde = 0;
        ddgraph->underlyingGraph->d = (int *)malloc(maxVertices*sizeof(int));
        ddgraph->underlyingGraph->v = (int *)malloc(maxVertices*sizeof(int));
        ddgraph->underlyingGraph->e = (int *)malloc((3*order + 2*(order/2))*sizeof(int));
        for(i = 0; i < order; i++){
            ddgraph->underlyingGraph->d[i]=3;
            ddgraph->underlyingGraph->v[i]=3*i;
        }
        for(i = order; i < maxVertices; i++){
            ddgraph->underlyingGraph->d[i]=2;
            ddgraph->underlyingGraph->v[i]=3*order + 2*i;
        }
        ddgraph->underlyingGraph->dlen = maxVertices;
        ddgraph->underlyingGraph->vlen = maxVertices;
        ddgraph->underlyingGraph->elen = 3*order + 2*(order/2);
    }

    return ddgraph;
}

void freeDDGraph(DDGRAPH *ddgraph){
    free(ddgraph->underlyingGraph->d);
    free(ddgraph->underlyingGraph->v);
    free(ddgraph->underlyingGraph->e);
    free(ddgraph->underlyingGraph);
    free(ddgraph->semiEdges);
    free(ddgraph);
}

permutation *getIdentity(int order){
    permutation i;
    permutation *identity = (permutation *)malloc(sizeof(permutation)*order);
    for(i = 0; i<order; i++){
        identity[i]=i;
    }
    return identity;
}

/**
 * Configures and returns the given building block with the given details.
 * 
 * @param block
 * @param type
 * @param component
 * @param parameter
 * @return 
 */
BBLOCK *initBuildingBlock(BBLOCK* block, int type, int component, int parameter){
    block->type = type;
    block->component = component;
    block->parameter = parameter;
    if(type==1){
        block->connectorCount = Q1TypeComponentsConnectors[component];
        block->connections = (BBLOCK **)malloc(sizeof(BBLOCK *)*Q1TypeComponentsConnectors[component]);
        block->targetConnector = (int *)malloc(sizeof(int)*Q1TypeComponentsConnectors[component]);
        block->connectionVertices = (int *)malloc(sizeof(int)*Q1TypeComponentsConnectors[component]);
    } else if(type==2){
        block->connectorCount = Q2TypeComponentsConnectors[component];
        block->connections = (BBLOCK **)malloc(sizeof(BBLOCK *)*Q2TypeComponentsConnectors[component]);
        block->targetConnector = (int *)malloc(sizeof(int)*Q2TypeComponentsConnectors[component]);
        block->connectionVertices = (int *)malloc(sizeof(int)*Q2TypeComponentsConnectors[component]);
    } else if(type==3){
        block->connectorCount = Q3TypeComponentsConnectors[component];
        block->connections = (BBLOCK **)malloc(sizeof(BBLOCK *)*Q3TypeComponentsConnectors[component]);
        block->targetConnector = (int *)malloc(sizeof(int)*Q3TypeComponentsConnectors[component]);
        block->connectionVertices = (int *)malloc(sizeof(int)*Q3TypeComponentsConnectors[component]);
    } else if(type==4){
        block->connectorCount = 1;
        block->connections = (BBLOCK **)malloc(sizeof(BBLOCK *));
        block->targetConnector = (int *)malloc(sizeof(int));
        block->connectionVertices = (int *)malloc(sizeof(int));
    } else {
        fprintf(stderr, "Illegal component type: %d (valid types from 1 to 4)\n", type);
        exit(EXIT_FAILURE);
    }
    return block;
}

/**
 * Creates a new building block with the given details.
 *
 * @param type
 * @param component
 * @param parameter
 * @return
 */
BBLOCK *getBuildingBlock(int type, int component, int parameter){
    BBLOCK *block = (BBLOCK *) malloc(sizeof(BBLOCK));
    return initBuildingBlock(block, type, component, parameter);
}

/**
 * Clears the connectors of the blocks in the given list.
 */
void clearConnectors(BBLOCK* blocks, int buildingBlockCount){
    int i, j;

    for (i = 0; i < buildingBlockCount; i++) {
        for(j = 0; j < (blocks+i)->connectorCount; j++){
            (blocks+i)->connections[j]=NULL;
        }
    }
}

inline void increaseConnectionsMade(){
    connectionsMade++;
    numberOfGenerators[connectionsMade] = 0;
}

void storeGenerator(int count, permutation perm[], nvector orbits[],
        int numorbits, int stabvertex, int n) {
    memcpy(automorphismGroupGenerators[connectionsMade] +
            numberOfGenerators[connectionsMade], perm, sizeof(permutation) * n);

    numberOfGenerators[connectionsMade]++;
}

int buildingBlockTypeToNumber(BBLOCK *block){
    if(block->type==1){
        return block->component;
    } else if(block->type==2){
        return Q1TypeComponentsCount + block->component;
    } else if(block->type==3){
        return Q1TypeComponentsCount + Q2TypeComponentsCount + block->component;
    } else {
        return Q1TypeComponentsCount + Q2TypeComponentsCount + Q3TypeComponentsCount;
    }
}

//====================== Building block construction ==========================

/*
 * Constructs a hub H(n). This block has 4 connectors arranged as follows:
 *
 *         0                 2
 *          \               /
 *           o---o-...-o---o
 *           |   |     |   |
 *           |   |     |   |
 *           o---o-...-o---o
 *          /               \
 *         1                 3
 *
 */
void constructHub(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    int i;
    block->connectionVertices[0] = *currentVertex;
    block->connectionVertices[1] = (*currentVertex)+1;
    block->connectionVertices[2] = (*currentVertex)+(block->parameter-1)*4+2;
    block->connectionVertices[3] = (*currentVertex)+(block->parameter-1)*4+3;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)+2;
    degrees[*currentVertex] = 2;
    positions[*currentVertex]++; //connection points are at position 0

    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //connection points are at position 0

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        edges[positions[*currentVertex]+0] = (*currentVertex)-1;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-2;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;

        (*currentVertex)+=4;
    }

    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)-2;
    degrees[*currentVertex] = 2;
    positions[*currentVertex]++; //connection points are at position 0

    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)-1;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //connection points are at position 0

    (*currentVertex)+=2;
}

void storeHubAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    if(block->parameter==1){
        //cyclic symmetry
        permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

        generator[block->connectionVertices[0]] = block->connectionVertices[1];
        generator[block->connectionVertices[1]] = block->connectionVertices[3];
        generator[block->connectionVertices[2]] = block->connectionVertices[0];
        generator[block->connectionVertices[3]] = block->connectionVertices[2];

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    } else {
        //only mirror-symmetry along long axis
        permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

        int vertex = block->connectionVertices[0];

        int i;
        for(i = 0; i < 2*(block->parameter); i++){
            generator[vertex] = vertex+1;
            generator[vertex+1] = vertex;
            vertex+=2;
        }

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    }
}

void storeHubsMapping(BBLOCK *block1, BBLOCK *block2, DDGRAPH *ddgraph){
    int i, vertexBlock1, vertexBlock2;

    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    vertexBlock1 = block1->connectionVertices[0];
    vertexBlock2 = block2->connectionVertices[0];

    for(i=0; i<4*(block1->parameter); i++){
        generator[vertexBlock1+i] = vertexBlock2+i;
        generator[vertexBlock2+i] = vertexBlock1+i;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

/*
 * Constructs a locked hub LH(n). This block has 3 connectors arranged
 * as follows:
 *
 *                           1
 *          \               /
 *           o---o-...-o---o
 *           |   |     |   |
 *           |   |     |   |
 *           o---o-...-o---o
 *          /               \
 *         0                 2
 *
 */
void constructLockedHub(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    int i;
    block->connectionVertices[0] = (*currentVertex)+1;
    block->connectionVertices[1] = (*currentVertex)+(block->parameter-1)*4+2;
    block->connectionVertices[2] = (*currentVertex)+(block->parameter-1)*4+3;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    edges[positions[*currentVertex]+0] = (*currentVertex)+1;
    edges[positions[*currentVertex]+1] = (*currentVertex)+2;
    edges[positions[*currentVertex]+2] = SEMIEDGE;
    degrees[*currentVertex] = 2;
    ddgraph->semiEdges[*currentVertex] = 1;
    //don't change position: semi-edges are at the end

    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //connection points are at position 0

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        edges[positions[*currentVertex]+0] = (*currentVertex)-1;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-2;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;

        (*currentVertex)+=4;
    }

    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)-2;
    degrees[*currentVertex] = 2;
    positions[*currentVertex]++; //connection points are at position 0

    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)-1;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //connection points are at position 0

    (*currentVertex)+=2;
}

void storeLockedHubAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    if(block->parameter==1){
        //mirror symmetry along diagonal through semi-edge
        permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

        generator[block->connectionVertices[0]] = block->connectionVertices[1];
        generator[block->connectionVertices[1]] = block->connectionVertices[0];

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    } // else no symmetry
}

void storeLockedHubsMapping(BBLOCK *block1, BBLOCK *block2, DDGRAPH *ddgraph){
    int i, vertexBlock1, vertexBlock2;

    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    vertexBlock1 = block1->connectionVertices[0]-1;
    vertexBlock2 = block2->connectionVertices[0]-1;

    for(i=0; i<4*(block1->parameter); i++){
        generator[vertexBlock1+i] = vertexBlock2+i;
        generator[vertexBlock2+i] = vertexBlock1+i;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

/*
 * Constructs a diagonal chain DC(n). This block has 2 connectors arranged
 * as follows:
 *              _______________
 *             /               |
 *            /              1 |
 *           /              /  |
 *           o---o-...-o---o  /
 *           |   |     |   | /
 *           |   |     |   |/
 *           o---o-...-o---o
 *          /
 *         0
 *
 */
void constructDiagonalChain(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    int i, start;
    block->connectionVertices[0] = (*currentVertex)+1;
    block->connectionVertices[1] = (*currentVertex)+(block->parameter-1)*4+2;

    start = *currentVertex;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    edges[positions[*currentVertex]+0] = (*currentVertex)+(block->parameter-1)*4+3;
    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)+2;

    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //connection points are at position 0
    
    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        edges[positions[*currentVertex]+0] = (*currentVertex)-1;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-2;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;

        (*currentVertex)+=4;
    }

    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)-2;
    degrees[*currentVertex] = 2;
    positions[*currentVertex]++; //connection points are at position 0

    edges[positions[(*currentVertex)+1]+1] = start;
    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)-1;

    (*currentVertex)+=2;
}

void storeDiagonalChainAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    if(block->parameter==1){
        //mirror symmetry along diagonal
        permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

        generator[block->connectionVertices[0]] = block->connectionVertices[1];
        generator[block->connectionVertices[1]] = block->connectionVertices[0];

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        //mirror symmetry along other diagonal

        generator[block->connectionVertices[0]] = block->connectionVertices[0];
        generator[block->connectionVertices[1]] = block->connectionVertices[1];
        generator[block->connectionVertices[0]-1] = block->connectionVertices[1]+1;
        generator[block->connectionVertices[1]+1] = block->connectionVertices[0]-1;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    } // else no symmetry
}

void storeDiagonalChainsMapping(BBLOCK *block1, BBLOCK *block2, DDGRAPH *ddgraph){
    int i, vertexBlock1, vertexBlock2;

    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    vertexBlock1 = block1->connectionVertices[0]-1;
    vertexBlock2 = block2->connectionVertices[0]-1;

    for(i=0; i<4*(block1->parameter); i++){
        generator[vertexBlock1+i] = vertexBlock2+i;
        generator[vertexBlock2+i] = vertexBlock1+i;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

/*
 * Constructs a double-roof high building DHB(n). This block has 2 connectors
 * arranged as follows:
 *
 *         0
 *          \
 *           o---o-...-o---o
 *           |   |     |   |\
 *           |   |     |   |/
 *           o---o-...-o---o
 *          /
 *         1
 *
 */
void constructDoubleroofHighBuilding(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    int i;
    block->connectionVertices[0] = *currentVertex;
    block->connectionVertices[1] = (*currentVertex)+1;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)+2;
    degrees[*currentVertex] = 2;
    positions[*currentVertex]++; //connection points are at position 0

    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //connection points are at position 0

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        edges[positions[*currentVertex]+0] = (*currentVertex)-1;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-2;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;

        (*currentVertex)+=4;
    }

    int dummyVertex = ddgraph->order + ddgraph->dummyVertexCount;

    edges[positions[*currentVertex]+0] = dummyVertex;
    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)-2;

    edges[positions[(*currentVertex)+1]+0] = dummyVertex;
    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)-1;

    edges[positions[dummyVertex]+0] = (*currentVertex);
    edges[positions[dummyVertex]+1] = (*currentVertex)+1;

    ddgraph->dummyVertexCount++;
    
    (*currentVertex)+=2;
}

void storeDoubleroofHighBuildingAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    //only mirror-symmetry along long axis
    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    int vertex = block->connectionVertices[0];

    int i;
    for(i = 0; i < 2*(block->parameter); i++){
        generator[vertex] = vertex+1;
        generator[vertex+1] = vertex;
        vertex+=2;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

void storeDoubleroofHighBuildingsMapping(BBLOCK *block1, BBLOCK *block2, DDGRAPH *ddgraph){
    int i, vertexBlock1, vertexBlock2;

    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    vertexBlock1 = block1->connectionVertices[0];
    vertexBlock2 = block2->connectionVertices[0];

    for(i=0; i<4*(block1->parameter); i++){
        generator[vertexBlock1+i] = vertexBlock2+i;
        generator[vertexBlock2+i] = vertexBlock1+i;
    }

    int dummyVertex1 = ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[vertexBlock1 + 4*(block1->parameter)]+0];
    int dummyVertex2 = ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[vertexBlock2 + 4*(block2->parameter)]+0];

    generator[dummyVertex1] = dummyVertex2;
    generator[dummyVertex2] = dummyVertex1;

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

/*
 * Constructs a open-roof high building OHB(n). This block has 2 connectors
 * arranged as follows:
 *
 *         0
 *          \               /
 *           o---o-...-o---o
 *           |   |     |   |
 *           |   |     |   |
 *           o---o-...-o---o
 *          /               \
 *         1
 *
 */
void constructOpenroofHighBuilding(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    int i;
    block->connectionVertices[0] = *currentVertex;
    block->connectionVertices[1] = (*currentVertex)+1;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)+2;
    degrees[*currentVertex] = 2;
    positions[*currentVertex]++; //connection points are at position 0

    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //connection points are at position 0

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        edges[positions[*currentVertex]+0] = (*currentVertex)-1;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-2;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;

        (*currentVertex)+=4;
    }


    edges[positions[*currentVertex]+0] = SEMIEDGE;
    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)-2;
    degrees[*currentVertex] = 2;
    positions[*currentVertex]++; //semi-edge is at position 0
    ddgraph->semiEdges[*currentVertex]=1;

    edges[positions[(*currentVertex)+1]+0] = SEMIEDGE;
    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)-1;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //semi-edge is at position 0
    ddgraph->semiEdges[(*currentVertex)+1]=1;

    (*currentVertex)+=2;
}

void storeOpenroofHighBuildingAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    //only mirror-symmetry along long axis
    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    int vertex = block->connectionVertices[0];

    int i;
    for(i = 0; i < 2*(block->parameter); i++){
        generator[vertex] = vertex+1;
        generator[vertex+1] = vertex;
        vertex+=2;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

void storeOpenroofHighBuildingsMapping(BBLOCK *block1, BBLOCK *block2, DDGRAPH *ddgraph){
    int i, vertexBlock1, vertexBlock2;

    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    vertexBlock1 = block1->connectionVertices[0];
    vertexBlock2 = block2->connectionVertices[0];

    for(i=0; i<4*(block1->parameter); i++){
        generator[vertexBlock1+i] = vertexBlock2+i;
        generator[vertexBlock2+i] = vertexBlock1+i;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

/*
 * Constructs a double-roof long building DLB(n). This block has 2 connectors
 * arranged as follows:
 *
 *             ___________
 *            /           \
 *           o---o-...-o---o
 *           |   |     |   |
 *           |   |     |   |
 *           o---o-...-o---o
 *          /               \
 *         0                 1
 *
 * No special measures need to be made for the case where n is 1, because this
 * is not a legal type.
 */
void constructDoubleroofLongBuilding(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    int i, start;
    block->connectionVertices[0] = (*currentVertex)+1;
    block->connectionVertices[1] = (*currentVertex)+(block->parameter-1)*4+3;

    start = *currentVertex;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    edges[positions[*currentVertex]+0] = (*currentVertex)+(block->parameter-1)*4+2;
    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)+2;

    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //connection points are at position 0

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        edges[positions[*currentVertex]+0] = (*currentVertex)-1;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-2;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;

        (*currentVertex)+=4;
    }

    edges[positions[*currentVertex]+0] = start;
    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)-2;

    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)-1;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //connection points are at position 0

    (*currentVertex)+=2;
}

void storeDoubleroofLongBuildingAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    //only mirror-symmetry along short axis
    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    int vertexLeft = block->connectionVertices[0];
    int vertexRight = block->connectionVertices[1];

    int i;
    for(i = 0; i < block->parameter; i++){
        generator[vertexLeft] = vertexRight;
        generator[vertexLeft - 1] = vertexRight-1;
        generator[vertexRight] = vertexLeft;
        generator[vertexRight - 1] = vertexLeft-1;
        vertexLeft+=2;
        vertexRight-=2;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

void storeDoubleroofLongBuildingsMapping(BBLOCK *block1, BBLOCK *block2, DDGRAPH *ddgraph){
    int i, vertexBlock1, vertexBlock2;

    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    vertexBlock1 = block1->connectionVertices[0]-1;
    vertexBlock2 = block2->connectionVertices[0]-1;

    for(i=0; i<4*(block1->parameter); i++){
        generator[vertexBlock1+i] = vertexBlock2+i;
        generator[vertexBlock2+i] = vertexBlock1+i;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

/*
 * Constructs a open-roof long building OLB(n). This block has 2 connectors
 * arranged as follows:
 *
 *
 *          \               /
 *           o---o-...-o---o
 *           |   |     |   |
 *           |   |     |   |
 *           o---o-...-o---o
 *          /               \
 *         0                 1
 *
 */
void constructOpenroofLongBuilding(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    int i;
    block->connectionVertices[0] = (*currentVertex)+1;
    block->connectionVertices[1] = (*currentVertex)+(block->parameter-1)*4+3;
    

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    edges[positions[*currentVertex]+0] = SEMIEDGE;
    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)+2;
    degrees[*currentVertex] = 2;
    positions[*currentVertex]++; //semi-edge is at position 0
    ddgraph->semiEdges[*currentVertex] = 1;

    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //connection points are at position 0

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        edges[positions[*currentVertex]+0] = (*currentVertex)-1;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-2;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;

        (*currentVertex)+=4;
    }

    edges[positions[*currentVertex]+1] = SEMIEDGE;
    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)-2;
    degrees[*currentVertex] = 2;
    positions[*currentVertex]++; //semi-edge is at position 0
    ddgraph->semiEdges[*currentVertex] = 1;

    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)-1;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //connection points are at position 0

    (*currentVertex)+=2;
}

void storeOpenroofLongBuildingAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    //only mirror-symmetry along short axis
    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    int vertexLeft = block->connectionVertices[0];
    int vertexRight = block->connectionVertices[1];

    int i;
    for(i = 0; i < block->parameter; i++){
        generator[vertexLeft] = vertexRight;
        generator[vertexLeft - 1] = vertexRight-1;
        generator[vertexRight] = vertexLeft;
        generator[vertexRight - 1] = vertexLeft-1;
        vertexLeft+=2;
        vertexRight-=2;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

void storeOpenroofLongBuildingsMapping(BBLOCK *block1, BBLOCK *block2, DDGRAPH *ddgraph){
    int i, vertexBlock1, vertexBlock2;

    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    vertexBlock1 = block1->connectionVertices[0]-1;
    vertexBlock2 = block2->connectionVertices[0]-1;

    for(i=0; i<4*(block1->parameter); i++){
        generator[vertexBlock1+i] = vertexBlock2+i;
        generator[vertexBlock2+i] = vertexBlock1+i;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

/*
 * Constructs a locked diagonal chain LDC(n). This block has 1 connector
 * arranged as follows:
 *
 *              _______________
 *             /               |
 *            /                |
 *           /              /  |
 *           o---o-...-o---o  /
 *           |   |     |   | /
 *           |   |     |   |/
 *           o---o-...-o---o
 *          /
 *         0
 *
 */
void constructLockedDiagonalChain(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    int i, start;
    block->connectionVertices[0] = (*currentVertex)+1;

    start = *currentVertex;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    edges[positions[*currentVertex]+0] = (*currentVertex)+(block->parameter-1)*4+3;
    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)+2;

    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //connection points are at position 0

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        edges[positions[*currentVertex]+0] = (*currentVertex)-1;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-2;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;

        (*currentVertex)+=4;
    }

    edges[positions[*currentVertex]+0] = SEMIEDGE;
    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)-2;
    degrees[*currentVertex] = 2;
    positions[*currentVertex]++; //semi-edge is at position 0
    ddgraph->semiEdges[*currentVertex] = 1;

    edges[positions[(*currentVertex)+1]+1] = start;
    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)-1;

    (*currentVertex)+=2;
}

void storeLockedDiagonalChainAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    if(block->parameter==1){
        //mirror symmetry along diagonal
        permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

        generator[block->connectionVertices[0] - 1] = block->connectionVertices[0] + 2;
        generator[block->connectionVertices[0] + 2] = block->connectionVertices[0] - 1;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    } // else no symmetry
}

void storeLockedDiagonalChainsMapping(BBLOCK *block1, BBLOCK *block2, DDGRAPH *ddgraph){
    int i, vertexBlock1, vertexBlock2;

    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    vertexBlock1 = block1->connectionVertices[0]-1;
    vertexBlock2 = block2->connectionVertices[0]-1;

    for(i=0; i<4*(block1->parameter); i++){
        generator[vertexBlock1+i] = vertexBlock2+i;
        generator[vertexBlock2+i] = vertexBlock1+i;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

/*
 * Constructs a locked double-roof high building LDHB(n). This block has 1
 * connector arranged as follows:
 *
 *         0
 *          \
 *           o---o-...-o---o
 *           |   |     |   |\
 *           |   |     |   |/
 *           o---o-...-o---o
 *          /
 *
 * This block has a trivial symmetry group.
 */
void constructLockedDoubleroofHighBuilding(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    int i;
    block->connectionVertices[0] = *currentVertex;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)+2;
    degrees[*currentVertex] = 2;
    positions[*currentVertex]++; //connection points are at position 0

    edges[positions[(*currentVertex)+1]+1] = SEMIEDGE;
    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //semi-edge is at position 0
    ddgraph->semiEdges[(*currentVertex)+1] = 1;

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        edges[positions[*currentVertex]+0] = (*currentVertex)-1;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-2;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;

        (*currentVertex)+=4;
    }

    int dummyVertex = ddgraph->order + ddgraph->dummyVertexCount;

    edges[positions[*currentVertex]+0] = dummyVertex;
    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)-2;

    edges[positions[(*currentVertex)+1]+0] = dummyVertex;
    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)-1;

    edges[positions[dummyVertex]+0] = (*currentVertex);
    edges[positions[dummyVertex]+1] = (*currentVertex)+1;

    ddgraph->dummyVertexCount++;

    (*currentVertex)+=2;
}

void storeLockedDoubleroofHighBuildingsMapping(BBLOCK *block1, BBLOCK *block2, DDGRAPH *ddgraph){
    int i, vertexBlock1, vertexBlock2;

    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    vertexBlock1 = block1->connectionVertices[0];
    vertexBlock2 = block2->connectionVertices[0];

    for(i=0; i<4*(block1->parameter); i++){
        generator[vertexBlock1+i] = vertexBlock2+i;
        generator[vertexBlock2+i] = vertexBlock1+i;
    }

    int dummyVertex1 = ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[vertexBlock1 + 4*(block1->parameter)]+0];
    int dummyVertex2 = ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[vertexBlock2 + 4*(block2->parameter)]+0];

    generator[dummyVertex1] = dummyVertex2;
    generator[dummyVertex2] = dummyVertex1;

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

/*
 * Constructs a locked open-roof high building LOHB(n). This block has 1
 * connector arranged as follows:
 *
 *         0
 *          \               /
 *           o---o-...-o---o
 *           |   |     |   |
 *           |   |     |   |
 *           o---o-...-o---o
 *          /               \
 *         
 *
 */
void constructLockedOpenroofHighBuilding(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    int i;
    block->connectionVertices[0] = *currentVertex;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)+2;
    degrees[*currentVertex] = 2;
    positions[*currentVertex]++; //connection points are at position 0

    edges[positions[(*currentVertex)+1]+0] = SEMIEDGE;
    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //semi-edge is at position 0
    ddgraph->semiEdges[(*currentVertex)+1]=1;

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        edges[positions[*currentVertex]+0] = (*currentVertex)-1;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-2;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;

        (*currentVertex)+=4;
    }


    edges[positions[*currentVertex]+0] = SEMIEDGE;
    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)-2;
    degrees[*currentVertex] = 2;
    positions[*currentVertex]++; //semi-edge is at position 0
    ddgraph->semiEdges[*currentVertex]=1;

    edges[positions[(*currentVertex)+1]+0] = SEMIEDGE;
    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)-1;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //semi-edge is at position 0
    ddgraph->semiEdges[(*currentVertex)+1]=1;

    (*currentVertex)+=2;
}

void storeLockedOpenroofHighBuildingAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    if(block->parameter==1){
        //mirror symmetry along diagonal
        permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

        generator[block->connectionVertices[0] + 1] = block->connectionVertices[0] + 2;
        generator[block->connectionVertices[0] + 2] = block->connectionVertices[0] + 1;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    } // else no symmetry
}

void storeLockedOpenroofHighBuildingsMapping(BBLOCK *block1, BBLOCK *block2, DDGRAPH *ddgraph){
    int i, vertexBlock1, vertexBlock2;

    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    vertexBlock1 = block1->connectionVertices[0];
    vertexBlock2 = block2->connectionVertices[0];

    for(i=0; i<4*(block1->parameter); i++){
        generator[vertexBlock1+i] = vertexBlock2+i;
        generator[vertexBlock2+i] = vertexBlock1+i;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

/*
 * Constructs a locked double-roof long building LDLB(n). This block has 1
 * connector arranged as follows:
 *
 *             ___________
 *            /           \
 *           o---o-...-o---o
 *           |   |     |   |
 *           |   |     |   |
 *           o---o-...-o---o
 *          /               \
 *         0
 * No special measures need to be made for the case where n is 1, because this
 * is not a legal type.
 * 
 * This block has a trivial symmetry group.
 */
void constructLockedDoubleroofLongBuilding(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    int i, start;
    block->connectionVertices[0] = (*currentVertex)+1;

    start = *currentVertex;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    edges[positions[*currentVertex]+0] = (*currentVertex)+(block->parameter-1)*4+2;
    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)+2;

    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //connection points are at position 0

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        edges[positions[*currentVertex]+0] = (*currentVertex)-1;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-2;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;

        (*currentVertex)+=4;
    }

    edges[positions[*currentVertex]+0] = start;
    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)-2;

    edges[positions[(*currentVertex)+1]+0] = SEMIEDGE;
    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)-1;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //semi-edge is at position 0

    (*currentVertex)+=2;
}

void storeLockedDoubleroofLongBuildingsMapping(BBLOCK *block1, BBLOCK *block2, DDGRAPH *ddgraph){
    int i, vertexBlock1, vertexBlock2;

    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    vertexBlock1 = block1->connectionVertices[0]-1;
    vertexBlock2 = block2->connectionVertices[0]-1;

    for(i=0; i<4*(block1->parameter); i++){
        generator[vertexBlock1+i] = vertexBlock2+i;
        generator[vertexBlock2+i] = vertexBlock1+i;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

/*
 * Constructs a pearl chain PC(n). This block has 2 connectors
 * arranged as follows:
 *         ____         ____
 *        /    \       /    \
 *    0--o      o-...-o      o--1
 *        \____/       \____/
 *
 */
void constructPearlChain(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    int i;
    block->connectionVertices[0] = *currentVertex;
    block->connectionVertices[1] = *currentVertex + 1 + 2*(block->parameter);

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    for(i=0; i<block->parameter; i++){
        int dummyVertex = ddgraph->order + ddgraph->dummyVertexCount;
        ddgraph->dummyVertexCount++;

        edges[positions[*currentVertex]+0] = (*currentVertex)-1;
        //for the first vertex this will be overwritten when making the connections
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = dummyVertex;
        edges[positions[dummyVertex]+0] = (*currentVertex);
        (*currentVertex)++;

        edges[positions[*currentVertex]+0] = (*currentVertex)+1;
        //for the last vertex this will be overwritten when making the connections
        edges[positions[*currentVertex]+1] = (*currentVertex)-1;
        edges[positions[*currentVertex]+2] = dummyVertex;
        edges[positions[dummyVertex]+1] = (*currentVertex);
        (*currentVertex)++;
    }

    degrees[block->connectionVertices[0]] = 2;
    degrees[block->connectionVertices[1]] = 2;
    positions[block->connectionVertices[0]]++;
    positions[block->connectionVertices[1]]++;

}

void storePearlChainAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    int i, vertexLeft, vertexRight, dummyLeft, dummyRight, parameter;

    //mirror symmetry
    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    vertexLeft = block->connectionVertices[0];
    vertexRight = block->connectionVertices[1];
    dummyLeft = ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[vertexLeft]+2];
    dummyRight = ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[vertexRight]+2];
    parameter = block->parameter;

    for(i=0; i<parameter; i++){
        generator[vertexLeft] = vertexRight;
        generator[vertexRight] = vertexLeft;
        vertexLeft++;
        vertexRight--;
    }
    for(i=0; i<parameter/2; i++){
        generator[dummyLeft] = dummyRight;
        generator[dummyRight] = dummyLeft;
        dummyLeft++;
        dummyRight--;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

void storePearlChainsMapping(BBLOCK *block1, BBLOCK *block2, DDGRAPH *ddgraph){
    int i, vertexBlock1, vertexBlock2, dummyBlock1, dummyBlock2;

    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    vertexBlock1 = block1->connectionVertices[0];
    dummyBlock1 = ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[vertexBlock1]+2];
    vertexBlock2 = block2->connectionVertices[0];
    dummyBlock2 = ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[vertexBlock2]+2];

    for(i=0; i<(block1->parameter); i++){
        generator[vertexBlock1+2*i] = vertexBlock2+2*i;
        generator[vertexBlock2+2*i] = vertexBlock1+2*i;
        generator[vertexBlock1+2*i+1] = vertexBlock2+2*i+1;
        generator[vertexBlock2+2*i+1] = vertexBlock1+2*i+1;
        generator[dummyBlock1+i] = dummyBlock2+i;
        generator[dummyBlock2+i] = dummyBlock1+i;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

/*
 * Constructs a locked pearl chain LPC(n). This block has 1 connector
 * arranged as follows:
 *         ____         ____
 *        /    \       /    \
 *    0--o      o-...-o      o-
 *        \____/       \____/
 *
 * This block has a trivial symmetry group. (We don't count symmetries of the
 * double edge).
 */
void constructLockedPearlChain(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    int i;
    block->connectionVertices[0] = *currentVertex;
    
    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    for(i=0; i<block->parameter; i++){
        int dummyVertex = ddgraph->order + ddgraph->dummyVertexCount;
        ddgraph->dummyVertexCount++;

        edges[positions[*currentVertex]+0] = (*currentVertex)-1;
        //for the first vertex this will be overwritten when making the connections
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = dummyVertex;
        edges[positions[dummyVertex]+0] = (*currentVertex);
        (*currentVertex)++;

        edges[positions[*currentVertex]+0] = (*currentVertex)+1;
        //for the last vertex this will be overwritten when making the connections
        edges[positions[*currentVertex]+1] = (*currentVertex)-1;
        edges[positions[*currentVertex]+2] = dummyVertex;
        edges[positions[dummyVertex]+1] = (*currentVertex);
        (*currentVertex)++;
    }

    edges[positions[(*currentVertex)-1]+0] = SEMIEDGE;
    ddgraph->semiEdges[(*currentVertex)-1] = 1;

    degrees[block->connectionVertices[0]] = 2;
    degrees[block->connectionVertices[1]] = 2;
    positions[block->connectionVertices[0]]++;
    positions[block->connectionVertices[1]]++;
}

void storeLockedPearlChainsMapping(BBLOCK *block1, BBLOCK *block2, DDGRAPH *ddgraph){
    int i, vertexBlock1, vertexBlock2, dummyBlock1, dummyBlock2;

    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    vertexBlock1 = block1->connectionVertices[0];
    dummyBlock1 = ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[vertexBlock1]+2];
    vertexBlock2 = block2->connectionVertices[0];
    dummyBlock2 = ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[vertexBlock2]+2];

    for(i=0; i<(block1->parameter); i++){
        generator[vertexBlock1+2*i] = vertexBlock2+2*i;
        generator[vertexBlock2+2*i] = vertexBlock1+2*i;
        generator[vertexBlock1+2*i+1] = vertexBlock2+2*i+1;
        generator[vertexBlock2+2*i+1] = vertexBlock1+2*i+1;
        generator[dummyBlock1+i] = dummyBlock2+i;
        generator[dummyBlock2+i] = dummyBlock1+i;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

/*
 * Constructs a barb wire BW(n). This block has 2 connectors
 * arranged as follows:
 *
 *       |    |     |    |
 *    0--o----o-...-o----o--1
 *
 *
 */
void constructBarbWire(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    int i;
    block->connectionVertices[0] = *currentVertex;
    block->connectionVertices[1] = *currentVertex + 1 + 2*(block->parameter);

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    for(i=0; i<block->parameter; i++){

        edges[positions[*currentVertex]+0] = (*currentVertex)-1;
        //for the first vertex this will be overwritten when making the connections
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = SEMIEDGE;
        degrees[*currentVertex] = 2;
        (*currentVertex)++;

        edges[positions[*currentVertex]+0] = (*currentVertex)+1;
        //for the last vertex this will be overwritten when making the connections
        edges[positions[*currentVertex]+1] = (*currentVertex)-1;
        edges[positions[*currentVertex]+2] = SEMIEDGE;
        degrees[*currentVertex] = 2;
        (*currentVertex)++;
    }

    degrees[block->connectionVertices[0]] = 1;
    degrees[block->connectionVertices[1]] = 1;
    positions[block->connectionVertices[0]]++;
    positions[block->connectionVertices[1]]++;
}

void storeBarbWireAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    int i, vertexLeft, vertexRight, parameter;

    //mirror symmetry
    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    vertexLeft = block->connectionVertices[0];
    vertexRight = block->connectionVertices[1];
    parameter = block->parameter;

    for(i=0; i<parameter; i++){
        generator[vertexLeft] = vertexRight;
        generator[vertexRight] = vertexLeft;
        vertexLeft++;
        vertexRight--;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

void storeBarbWiresMapping(BBLOCK *block1, BBLOCK *block2, DDGRAPH *ddgraph){
    int i, vertexBlock1, vertexBlock2;

    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    vertexBlock1 = block1->connectionVertices[0];
    vertexBlock2 = block2->connectionVertices[0];

    for(i=0; i<2*(block1->parameter); i++){
        generator[vertexBlock1+i] = vertexBlock2+i;
        generator[vertexBlock2+i] = vertexBlock1+i;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

/*
 * Constructs a locked barb wire LBW(n). This block has 1 connector
 * arranged as follows:
 *
 *       |    |     |    |
 *    0--o----o-...-o----o-
 *
 * This block has a trivial symmetry group.
 */
void constructLockedBarbWire(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    int i;
    block->connectionVertices[0] = *currentVertex;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    for(i=0; i<block->parameter; i++){

        edges[positions[*currentVertex]+0] = (*currentVertex)-1;
        //for the first vertex this will be overwritten when making the connections
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = SEMIEDGE;
        degrees[*currentVertex] = 2;
        ddgraph->semiEdges[(*currentVertex)] = 1;
        (*currentVertex)++;

        edges[positions[*currentVertex]+0] = (*currentVertex)+1;
        //for the last vertex this will be overwritten when making the connections
        edges[positions[*currentVertex]+1] = (*currentVertex)-1;
        edges[positions[*currentVertex]+2] = SEMIEDGE;
        degrees[*currentVertex] = 2;
        ddgraph->semiEdges[(*currentVertex)] = 1;
        (*currentVertex)++;
    }
    edges[positions[(*currentVertex)-1]+0] = SEMIEDGE;
    degrees[(*currentVertex)-1] = 1;
    ddgraph->semiEdges[(*currentVertex)-1] = 2;
    positions[(*currentVertex)-1]++;

    degrees[block->connectionVertices[0]] = 1;
    positions[block->connectionVertices[0]]++;
}

void storeLockedBarbWiresMapping(BBLOCK *block1, BBLOCK *block2, DDGRAPH *ddgraph){
    int i, vertexBlock1, vertexBlock2;

    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    vertexBlock1 = block1->connectionVertices[0];
    vertexBlock2 = block2->connectionVertices[0];

    for(i=0; i<2*(block1->parameter); i++){
        generator[vertexBlock1+i] = vertexBlock2+i;
        generator[vertexBlock2+i] = vertexBlock1+i;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

/*
 * Constructs a q4 block. This block has 1 connector
 * arranged as follows:
 *
 *    0--o<
 *
 * This block has a trivial symmetry group (We don't consider switching the two
 * semi-edges).
 */
void constructQ4(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    ddgraph->underlyingGraph->e[(*currentVertex) + 1] = ddgraph->underlyingGraph->e[(*currentVertex) + 2] = SEMIEDGE;
    ddgraph->underlyingGraph->d[(*currentVertex)] = 0;
    ddgraph->semiEdges[(*currentVertex)] = 2;

    //it is not necessary to adjust ddgraph->underlyingGraph->v because it will
    //not be used when the degree is 0.
    block->connectionVertices[0] = *currentVertex;
    (*currentVertex)++;
}

void storeQ4sMapping(BBLOCK *block1, BBLOCK *block2, DDGRAPH *ddgraph){
    int vertexBlock1, vertexBlock2;

    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    vertexBlock1 = block1->connectionVertices[0];
    vertexBlock2 = block2->connectionVertices[0];

    generator[vertexBlock1] = vertexBlock2;
    generator[vertexBlock2] = vertexBlock1;

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

void constructBuildingBlockListAsGraph(BBLOCK* blocks, int buildingBlockCount, DDGRAPH *ddgraph){
    int i, currentVertex=0;

    for (i = 0; i < buildingBlockCount; i++) {
        int number = buildingBlockTypeToNumber(blocks + i);
        (*constructBlock[number])(&currentVertex, (blocks+i), ddgraph);
    }

    ddgraph->underlyingGraph->nv = ddgraph->order + ddgraph->dummyVertexCount;

    for(i = 0; i < ddgraph->underlyingGraph->nv; i++){
        ddgraph->underlyingGraph->nde += ddgraph->underlyingGraph->d[i];
    }
}


//=============== PHASE 2: CONNECTING THE BUILDING BLOCKS ===================

BBLOCK *constructComponentList(int *blockListSize){
    int i,j,k;
    int currentBlock;
    int blockCount = Q4ComponentCount;

    for(i = 0; i < Q1TypeComponentsCount; i++){
        for(j = 0; j < maximumQ1TypeComponents; j++){
            blockCount += Q1TypeComponentsComponentCount[i][j];
        }
    }
    for(i = 0; i < Q2TypeComponentsCount; i++){
        for(j = 0; j < maximumQ2TypeComponents; j++){
            blockCount += Q2TypeComponentsComponentCount[i][j];
        }
    }
    for(i = 0; i < Q3TypeComponentsCount; i++){
        for(j = 0; j < maximumQ3TypeComponents; j++){
            blockCount += Q3TypeComponentsComponentCount[i][j];
        }
    }

    *blockListSize = blockCount;
    BBLOCK *block = (BBLOCK *) malloc(sizeof(BBLOCK)*blockCount);
    currentBlock = 0;

    for(i = 0; i < Q1TypeComponentsCount; i++){
        for(j = 0; j < maximumQ1TypeComponents; j++){
            for(k = 0; k < Q1TypeComponentsComponentCount[i][j]; k++){
                initBuildingBlock(block + currentBlock, 1, i, j);
            }
        }
    }
    for(i = 0; i < Q2TypeComponentsCount; i++){
        for(j = 0; j < maximumQ2TypeComponents; j++){
            for(k = 0; k < Q2TypeComponentsComponentCount[i][j]; k++){
                initBuildingBlock(block + currentBlock, 2, i, j);
            }
        }
    }
    for(i = 0; i < Q3TypeComponentsCount; i++){
        for(j = 0; j < maximumQ3TypeComponents; j++){
            for(k = 0; k < Q3TypeComponentsComponentCount[i][j]; k++){
                initBuildingBlock(block + currentBlock, 3, i, j);
            }
        }
    }
    for(k = 0; k < Q4ComponentCount; k++){
        initBuildingBlock(block + currentBlock, 4, 0, 0);
    }

    return block;
}

void storeInitialGenerators(BBLOCK* blocks, int buildingBlockCount, DDGRAPH *ddgraph){
    int firstBlocks[16]; //stores the index of the first block of each type
                         //-1 in case no block of this type was encountered yet
    int i;
    for(i=0;i<16;i++) firstBlocks[i]=-1;

    for (i = 0; i < buildingBlockCount; i++) {
        int number = buildingBlockTypeToNumber(blocks+i);
        if(firstBlocks[number]==-1){
            firstBlocks[number]=i;
            if(storeBlockAutomorphismGenerators[number]!=NULL){
                (*storeBlockAutomorphismGenerators[number])(blocks+i, ddgraph);
            }
        } else {
            (*storeBlocksMapping[number])(blocks+firstBlocks[number], blocks+i, ddgraph);
        }
    }
}

void connectComponentList(int vertexCount){
    int blockCount = 0;
    //create an array of blocks based upon the numbers in the global arrays
    BBLOCK *blocks = constructComponentList(&blockCount);
    DDGRAPH *graph = getNewDDGraph(vertexCount);

    //reset counters
    connectionsMade=0;
    numberOfGenerators[connectionsMade]=0;

    //convert the blocks to a disconnected graph
    constructBuildingBlockListAsGraph(blocks, blockCount, graph);
    //store the generators of the automorphism group of this disconnected graph
    storeInitialGenerators(blocks, blockCount, graph);

    

    //free the memory allocated at the beginning of this method
    freeDDGraph(graph);
    free(blocks);
}

//================ PHASE 1: GENERATION OF COMPONENT LISTS ===================

/**
 * This method checks whether the current list of components is realizable as a
 * metagraph. It does not take into account extra restrictions about which 
 * components may be connected to which other components.
 * 
 * @return FALSE if it is certain that the given list of components is impossible 
 */
boolean isMaybeRealizableComponentList(){
    int degreeSequence[4];
    int i, j;

    for (i = 0; i < 4; i++) {
        degreeSequence[i]=0;
    }

    for(i = 0; i < Q1TypeComponentsCount; i++){
        for(j = 0; j < maximumQ1TypeComponents; j++){
            degreeSequence[Q1TypeComponentsConnectors[i]-1] += Q1TypeComponentsComponentCount[i][j];
        }
    }

    for(i = 0; i < Q2TypeComponentsCount; i++){
        for(j = 0; j < maximumQ2TypeComponents; j++){
            degreeSequence[Q2TypeComponentsConnectors[i]-1] += Q2TypeComponentsComponentCount[i][j];
        }
    }

    for(i = 0; i < Q3TypeComponentsCount; i++){
        for(j = 0; j < maximumQ3TypeComponents; j++){
            degreeSequence[Q3TypeComponentsConnectors[i]-1] += Q3TypeComponentsComponentCount[i][j];
        }
    }

    degreeSequence[0] += Q4ComponentCount;

    int doubleEdgeSum = 0;
    int verticesCount = 0;
    int maxDegree = 0;
    for (i = 0; i < 4; i++) {
        doubleEdgeSum += (i+1)*degreeSequence[i];
        verticesCount += degreeSequence[i];
        if(degreeSequence[i]>0){
            maxDegree = i+1;
        }
    }

    if(maxDegree==0){
        return FALSE;
    } else if(doubleEdgeSum % 2 == 1){
        return FALSE;
    } else if(maxDegree*2 > doubleEdgeSum){
        return FALSE;
    } else if(verticesCount > doubleEdgeSum/2 + 1){
        return FALSE;
    } else {
        return TRUE;
    }

}

/**
 * Tests whether the current component lists passes some simple tests for
 * valid lists.
 * This method takes into account:
 *   - a q2 component may not be connected to a q2 component
 *   - a q3 component may not be connected to a q3 component
 *
 * @return FALSE if it is certain that the given list of components is impossible 
 */
boolean passesSimpleForbiddenConnectionsTest(){
    int i, j;
    int q2Connections=0, q3Connections=0, totalConnections;

    for(i = 0; i < Q2TypeComponentsCount; i++){
        for(j = 0; j < maximumQ2TypeComponents; j++){
            q2Connections += Q2TypeComponentsComponentCount[i][j]*Q2TypeComponentsConnectors[i];
        }
    }

    for(i = 0; i < Q3TypeComponentsCount; i++){
        for(j = 0; j < maximumQ3TypeComponents; j++){
            q3Connections += Q3TypeComponentsComponentCount[i][j]*Q3TypeComponentsConnectors[i];
        }
    }

    totalConnections = q2Connections + q3Connections + Q4ComponentCount;
    for(i = 0; i < Q1TypeComponentsCount; i++){
        for(j = 0; j < maximumQ1TypeComponents; j++){
            totalConnections += Q1TypeComponentsComponentCount[i][j]*Q1TypeComponentsConnectors[i];
        }
    }

    if(q2Connections*2>totalConnections){
        return FALSE;
    } else if(q3Connections*2>totalConnections){
        return FALSE;
    } else {
        return TRUE;
    }

}

void handleComponentList(int vertexCount){
    if(!isMaybeRealizableComponentList()){
        return;
    } else if(!passesSimpleForbiddenConnectionsTest()){
        return;
    } else {
        componentListsCount++;
#ifdef _DEBUG
        int i, j;
        for(i = 0; i < Q1TypeComponentsCount; i++){
            fprintf(stderr, "(");
            for(j = 0; j < maximumQ1TypeComponents; j++){
                fprintf(stderr, "%d ", Q1TypeComponentsComponentCount[i][j]);
            }
            fprintf(stderr, ")");
        }
        fprintf(stderr, "| ");
        for(i = 0; i < Q2TypeComponentsCount; i++){
            fprintf(stderr, "(");
            for(j = 0; j < maximumQ2TypeComponents; j++){
                fprintf(stderr, "%d ", Q2TypeComponentsComponentCount[i][j]);
            }
            fprintf(stderr, ")");
        }
        fprintf(stderr, "| ");
        for(i = 0; i < Q3TypeComponentsCount; i++){
            fprintf(stderr, "(");
            for(j = 0; j < maximumQ3TypeComponents; j++){
                fprintf(stderr, "%d ", Q3TypeComponentsComponentCount[i][j]);
            }
            fprintf(stderr, ")");
        }
        fprintf(stderr, "| %d\n", Q4ComponentCount);
#endif
        connectComponentList(vertexCount);
    }
}

void q4Components(int targetSize, int currentSize){
    Q4ComponentCount = targetSize - currentSize; //each q4 component has 1 vertex
    handleComponentList(targetSize);
}

void q3Components(int currentType, int currentParameter, int targetSize, int currentSize){
    int i;

    int remainingVertices = targetSize - currentSize;

    for(i = 0; i <= remainingVertices/(2*currentParameter); i++){
        Q3TypeComponentsComponentCount[currentType][currentParameter-1]=i;
        int newSize = currentSize + 2*currentParameter*i; //each q2 component has 2n vertices

        if(targetSize - newSize > 2*(currentParameter+1)){
            q3Components(currentType, currentParameter+1, targetSize, newSize);
        } else if(currentType+1==Q3TypeComponentsCount){
            q4Components(targetSize, newSize);
        } else {
            q3Components(currentType+1, Q3TypeComponentsSmallestCase[currentType+1], targetSize, newSize);
        }
    }
    Q3TypeComponentsComponentCount[currentType][currentParameter-1]=0; //reset this type to 0
}

void q2Components(int currentType, int currentParameter, int targetSize, int currentSize){
    int i;

    int remainingVertices = targetSize - currentSize;

    for(i = 0; i <= remainingVertices/(2*currentParameter); i++){
        Q2TypeComponentsComponentCount[currentType][currentParameter-1]=i;
        int newSize = currentSize + 2*currentParameter*i; //each q2 component has 2n vertices

        if(targetSize - newSize > 2*(currentParameter+1)){
            q2Components(currentType, currentParameter+1, targetSize, newSize);
        } else if(currentType+1==Q2TypeComponentsCount){
            q3Components(0, Q3TypeComponentsSmallestCase[0], targetSize, newSize);
        } else {
            q2Components(currentType+1, Q2TypeComponentsSmallestCase[currentType+1], targetSize, newSize);
        }
    }
    Q2TypeComponentsComponentCount[currentType][currentParameter-1]=0; //reset this type to 0
}

void q1Components(int currentType, int currentParameter, int targetSize, int currentSize){
    int i;

    int remainingVertices = targetSize - currentSize;

    for(i = 0; i <= remainingVertices/(4*currentParameter); i++){
        Q1TypeComponentsComponentCount[currentType][currentParameter-1]=i;
        int newSize = currentSize + 4*currentParameter*i; //each q1 component has 4n vertices

        if(targetSize - newSize > 4*(currentParameter+1)){
            q1Components(currentType, currentParameter+1, targetSize, newSize);
        } else if(currentType+1==Q1TypeComponentsCount){
            q2Components(0, Q2TypeComponentsSmallestCase[0], targetSize, newSize);
        } else {
            q1Components(currentType+1, Q1TypeComponentsSmallestCase[currentType+1], targetSize, newSize);
        }
    }
    Q1TypeComponentsComponentCount[currentType][currentParameter-1]=0; //reset this type to 0
}

//====================== INITIALIZATION =======================

void initComponents(int targetSize){
    int i,j;

    maximumQ1TypeComponents = targetSize/4;
    maximumQ2TypeComponents = maximumQ3TypeComponents = targetSize/2;

    Q1TypeComponentsConnectors[0] = 4; //Hub
    Q1TypeComponentsConnectors[1] = 3; //Locked Hub
    Q1TypeComponentsConnectors[2] = 2; //Diagonal Chain
    Q1TypeComponentsConnectors[3] = 2; //Double-roof High Building
    Q1TypeComponentsConnectors[4] = 2; //Open-roof High Building
    Q1TypeComponentsConnectors[5] = 2; //Double-roof Long Building
    Q1TypeComponentsConnectors[6] = 2; //Open-roof Long Building
    Q1TypeComponentsConnectors[7] = 1; //Locked Diagonal Chain
    Q1TypeComponentsConnectors[8] = 1; //Locked Double-roof High Building
    Q1TypeComponentsConnectors[9] = 1; //Locked Open-roof High Building
    Q1TypeComponentsConnectors[10] = 1; //Locked Double-roof Long Building

    Q1TypeComponentsSmallestCase[0] = 1;
    Q1TypeComponentsSmallestCase[1] = 1;
    Q1TypeComponentsSmallestCase[2] = 1;
    Q1TypeComponentsSmallestCase[3] = 1;
    Q1TypeComponentsSmallestCase[4] = 1;
    Q1TypeComponentsSmallestCase[5] = 2;
    Q1TypeComponentsSmallestCase[6] = 2;
    Q1TypeComponentsSmallestCase[7] = 1;
    Q1TypeComponentsSmallestCase[8] = 1;
    Q1TypeComponentsSmallestCase[9] = 1;
    Q1TypeComponentsSmallestCase[10] = 2;

    Q1TypeComponentsComponentCount = (int**)malloc(sizeof(int*)*Q1TypeComponentsCount);

    for(i = 0; i < Q1TypeComponentsCount; i++){
        Q1TypeComponentsComponentCount[i] = (int*)malloc(sizeof(int)*maximumQ1TypeComponents);
        for(j = 0; j < maximumQ1TypeComponents; j++){
            Q1TypeComponentsComponentCount[i][j]=0;
        }
    }

    Q2TypeComponentsConnectors[0] = 2; //Pearl Chain
    Q2TypeComponentsConnectors[1] = 1; //Locked Pearl Chain

    Q2TypeComponentsSmallestCase[0] = 1;
    Q2TypeComponentsSmallestCase[1] = 1;

    Q2TypeComponentsComponentCount = (int**)malloc(sizeof(int*)*Q2TypeComponentsCount);

    for(i = 0; i < Q2TypeComponentsCount; i++){
        Q2TypeComponentsComponentCount[i] = (int*)malloc(sizeof(int)*maximumQ2TypeComponents);
        for(j = 0; j < maximumQ2TypeComponents; j++){
            Q2TypeComponentsComponentCount[i][j]=0;
        }
    }

    Q3TypeComponentsConnectors[0] = 2; //Barb Wire
    Q3TypeComponentsConnectors[1] = 1; //Locked Barb Wire

    Q3TypeComponentsSmallestCase[0] = 1;
    Q3TypeComponentsSmallestCase[1] = 1;

    Q3TypeComponentsComponentCount = (int**)malloc(sizeof(int*)*Q3TypeComponentsCount);

    for(i = 0; i < Q3TypeComponentsCount; i++){
        Q3TypeComponentsComponentCount[i] = (int*)malloc(sizeof(int)*maximumQ3TypeComponents);
        for(j = 0; j < maximumQ3TypeComponents; j++){
            Q3TypeComponentsComponentCount[i][j]=0;
        }
    }

    Q4ComponentCount = 0;

    //set up arrays with function pointers

    constructBlock[0] = &constructHub;
    constructBlock[1] = &constructLockedHub;
    constructBlock[2] = &constructDiagonalChain;
    constructBlock[3] = &constructDoubleroofHighBuilding;
    constructBlock[4] = &constructOpenroofHighBuilding;
    constructBlock[5] = &constructDoubleroofLongBuilding;
    constructBlock[6] = &constructOpenroofLongBuilding;
    constructBlock[7] = &constructLockedDiagonalChain;
    constructBlock[8] = &constructLockedDoubleroofHighBuilding;
    constructBlock[9] = &constructLockedOpenroofHighBuilding;
    constructBlock[10] = &constructLockedDoubleroofLongBuilding;
    constructBlock[11] = &constructPearlChain;
    constructBlock[12] = &constructLockedPearlChain;
    constructBlock[13] = &constructBarbWire;
    constructBlock[14] = &constructLockedBarbWire;
    constructBlock[15] = &constructQ4;

    storeBlockAutomorphismGenerators[0] = &storeHubAutomorphismGenerators;
    storeBlockAutomorphismGenerators[1] = &storeLockedHubAutomorphismGenerators;
    storeBlockAutomorphismGenerators[2] = &storeDiagonalChainAutomorphismGenerators;
    storeBlockAutomorphismGenerators[3] = &storeDoubleroofHighBuildingAutomorphismGenerators;
    storeBlockAutomorphismGenerators[4] = &storeOpenroofHighBuildingAutomorphismGenerators;
    storeBlockAutomorphismGenerators[5] = &storeDoubleroofLongBuildingAutomorphismGenerators;
    storeBlockAutomorphismGenerators[6] = &storeOpenroofLongBuildingAutomorphismGenerators;
    storeBlockAutomorphismGenerators[7] = &storeLockedDiagonalChainAutomorphismGenerators;
    storeBlockAutomorphismGenerators[8] = NULL;
    storeBlockAutomorphismGenerators[9] = &storeLockedOpenroofHighBuildingAutomorphismGenerators;
    storeBlockAutomorphismGenerators[10] = NULL;
    storeBlockAutomorphismGenerators[11] = &storePearlChainAutomorphismGenerators;
    storeBlockAutomorphismGenerators[12] = NULL;
    storeBlockAutomorphismGenerators[13] = &storeBarbWireAutomorphismGenerators;
    storeBlockAutomorphismGenerators[14] = NULL;
    storeBlockAutomorphismGenerators[15] = NULL;

    storeBlocksMapping[0] = &storeHubsMapping;
    storeBlocksMapping[1] = &storeLockedHubsMapping;
    storeBlocksMapping[2] = &storeDiagonalChainsMapping;
    storeBlocksMapping[3] = &storeDoubleroofHighBuildingsMapping;
    storeBlocksMapping[4] = &storeOpenroofHighBuildingsMapping;
    storeBlocksMapping[5] = &storeDoubleroofLongBuildingsMapping;
    storeBlocksMapping[6] = &storeOpenroofLongBuildingsMapping;
    storeBlocksMapping[7] = &storeLockedDiagonalChainsMapping;
    storeBlocksMapping[8] = &storeLockedDoubleroofHighBuildingsMapping;
    storeBlocksMapping[9] = &storeLockedOpenroofHighBuildingsMapping;
    storeBlocksMapping[10] = &storeLockedDoubleroofLongBuildingsMapping;
    storeBlocksMapping[11] = &storePearlChainsMapping;
    storeBlocksMapping[12] = &storeLockedPearlChainsMapping;
    storeBlocksMapping[13] = &storeBarbWiresMapping;
    storeBlocksMapping[14] = &storeLockedBarbWiresMapping;
    storeBlocksMapping[15] = &storeQ4sMapping;
}

void initStatistics(){
    componentListsCount = 0;
}

//====================== START =======================

void startGeneration(int targetSize){

    initComponents(targetSize);
    initStatistics();

    q1Components(0, Q1TypeComponentsSmallestCase[0], targetSize, 0);

    fprintf(stderr, "Found %d component lists.\n", componentListsCount);

}

//====================== USAGE =======================

void help(char *name){

}

void usage(char *name){

}


#ifdef DDGRAPHS_NO_MAIN
    #define DDGRAPHS_MAIN_FUNCTION ddgraphsnomain
#else
    #define DDGRAPHS_MAIN_FUNCTION main
#endif
/*
 *
 */
int DDGRAPHS_MAIN_FUNCTION(int argc, char** argv) {

    /*=========== commandline parsing ===========*/

    int c;
    char *name = argv[0];

    while ((c = getopt(argc, argv, "h")) != -1) {
        switch (c) {
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            default:
                fprintf(stderr, "Illegal option %c.\n", c);
                usage(name);
                return EXIT_FAILURE;
        }
    }

    // check the non-option arguments
    if (argc - optind != 1) {
        usage(name);
        return EXIT_FAILURE;
    }

    //parse the order
    int vertexCount = strtol(argv[optind], NULL, 10);
    DEBUGDUMP(vertexCount, "%d")

    /*=========== initialization ===========*/
    struct tms TMS;
    unsigned int oldtime = 0;

    startGeneration(vertexCount);



    times(&TMS);
    unsigned int savetime = oldtime + (unsigned int) TMS.tms_utime;
    fprintf(stderr, "CPU time: %.1f seconds.\n", (double) savetime / time_factor);

    return EXIT_SUCCESS;
}
