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
    DDGRAPH *ddgraph = (DDGRAPH *)malloc(sizeof(DDGRAPH));
    ddgraph->order=order;
    ddgraph->adjList = (int**)malloc(sizeof(int *)*order);
    int i;
    for(i=0;i<order;i++){
        ddgraph->adjList[i] = (int*)malloc(sizeof(int)*3);
    }
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

/**
 * Translates an array of (connected!) blocks to a Delaney-Dress graph. This
 * method assumes that the different building blocks already have been
 * constructed in the Delaney-Dress graph and that the free neighbour of each
 * connector is the first neighbour (position 0 in the adjacency list).
 *
 * @param blocks
 * @param buildingBlockCount
 */
void metaGraph2DDgraph(BBLOCK* blocks, int buildingBlockCount, DDGRAPH *ddgraph){
    int i, j;

    for (i = 0; i < buildingBlockCount; i++) {
        for(j = 0; j < (blocks+i)->connectorCount; j++){
            ddgraph->adjList[(blocks+i)->connectionVertices[j]][0] =
                    (blocks+i)->connections[j]->connectionVertices[(blocks+i)->targetConnector[j]];
        }
    }
}

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

    ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][2] = (*currentVertex)+2;
    ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)+3;
    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        ddgraph->adjList[*currentVertex][0] = (*currentVertex)-1;
        ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
        ddgraph->adjList[*currentVertex][2] = (*currentVertex)+2;

        ddgraph->adjList[(*currentVertex)+1][0] = (*currentVertex)-2;
        ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
        ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)+3;

        ddgraph->adjList[(*currentVertex)+2][0] = *currentVertex;
        ddgraph->adjList[(*currentVertex)+2][1] = (*currentVertex)+3;
        ddgraph->adjList[(*currentVertex)+2][2] = (*currentVertex)+4;

        ddgraph->adjList[(*currentVertex)+3][0] = (*currentVertex)+1;
        ddgraph->adjList[(*currentVertex)+3][1] = (*currentVertex)+2;
        ddgraph->adjList[(*currentVertex)+3][2] = (*currentVertex)+5;
        (*currentVertex)+=4;
    }

    ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][2] = (*currentVertex)-2;
    ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)-1;
    (*currentVertex)+=2;
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

    ddgraph->adjList[*currentVertex][0] = SEMIEDGE;
    ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][2] = (*currentVertex)+2;
    ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)+3;
    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        ddgraph->adjList[*currentVertex][0] = (*currentVertex)-1;
        ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
        ddgraph->adjList[*currentVertex][2] = (*currentVertex)+2;

        ddgraph->adjList[(*currentVertex)+1][0] = (*currentVertex)-2;
        ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
        ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)+3;

        ddgraph->adjList[(*currentVertex)+2][0] = *currentVertex;
        ddgraph->adjList[(*currentVertex)+2][1] = (*currentVertex)+3;
        ddgraph->adjList[(*currentVertex)+2][2] = (*currentVertex)+4;

        ddgraph->adjList[(*currentVertex)+3][0] = (*currentVertex)+1;
        ddgraph->adjList[(*currentVertex)+3][1] = (*currentVertex)+2;
        ddgraph->adjList[(*currentVertex)+3][2] = (*currentVertex)+5;
        (*currentVertex)+=4;
    }

    ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][2] = (*currentVertex)-2;
    ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)-1;
    (*currentVertex)+=2;
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

    ddgraph->adjList[*currentVertex][0] = (*currentVertex)+(block->parameter-1)*4+3;
    ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][2] = (*currentVertex)+2;
    ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)+3;
    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        ddgraph->adjList[*currentVertex][0] = (*currentVertex)-1;
        ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
        ddgraph->adjList[*currentVertex][2] = (*currentVertex)+2;

        ddgraph->adjList[(*currentVertex)+1][0] = (*currentVertex)-2;
        ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
        ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)+3;

        ddgraph->adjList[(*currentVertex)+2][0] = *currentVertex;
        ddgraph->adjList[(*currentVertex)+2][1] = (*currentVertex)+3;
        ddgraph->adjList[(*currentVertex)+2][2] = (*currentVertex)+4;

        ddgraph->adjList[(*currentVertex)+3][0] = (*currentVertex)+1;
        ddgraph->adjList[(*currentVertex)+3][1] = (*currentVertex)+2;
        ddgraph->adjList[(*currentVertex)+3][2] = (*currentVertex)+5;
        (*currentVertex)+=4;
    }

    ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][2] = (*currentVertex)-2;
    ddgraph->adjList[(*currentVertex)+1][0] = start;
    ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)-1;
    (*currentVertex)+=2;
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

    ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][2] = (*currentVertex)+2;
    ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)+3;
    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        ddgraph->adjList[*currentVertex][0] = (*currentVertex)-1;
        ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
        ddgraph->adjList[*currentVertex][2] = (*currentVertex)+2;

        ddgraph->adjList[(*currentVertex)+1][0] = (*currentVertex)-2;
        ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
        ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)+3;

        ddgraph->adjList[(*currentVertex)+2][0] = *currentVertex;
        ddgraph->adjList[(*currentVertex)+2][1] = (*currentVertex)+3;
        ddgraph->adjList[(*currentVertex)+2][2] = (*currentVertex)+4;

        ddgraph->adjList[(*currentVertex)+3][0] = (*currentVertex)+1;
        ddgraph->adjList[(*currentVertex)+3][1] = (*currentVertex)+2;
        ddgraph->adjList[(*currentVertex)+3][2] = (*currentVertex)+5;
        (*currentVertex)+=4;
    }

    ddgraph->adjList[*currentVertex][0] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][2] = (*currentVertex)-2;
    ddgraph->adjList[(*currentVertex)+1][0] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)-1;
    (*currentVertex)+=2;
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

    ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][2] = (*currentVertex)+2;
    ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)+3;
    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        ddgraph->adjList[*currentVertex][0] = (*currentVertex)-1;
        ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
        ddgraph->adjList[*currentVertex][2] = (*currentVertex)+2;

        ddgraph->adjList[(*currentVertex)+1][0] = (*currentVertex)-2;
        ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
        ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)+3;

        ddgraph->adjList[(*currentVertex)+2][0] = *currentVertex;
        ddgraph->adjList[(*currentVertex)+2][1] = (*currentVertex)+3;
        ddgraph->adjList[(*currentVertex)+2][2] = (*currentVertex)+4;

        ddgraph->adjList[(*currentVertex)+3][0] = (*currentVertex)+1;
        ddgraph->adjList[(*currentVertex)+3][1] = (*currentVertex)+2;
        ddgraph->adjList[(*currentVertex)+3][2] = (*currentVertex)+5;
        (*currentVertex)+=4;
    }

    ddgraph->adjList[*currentVertex][0] = SEMIEDGE;
    ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][2] = (*currentVertex)-2;
    ddgraph->adjList[(*currentVertex)+1][0] = SEMIEDGE;
    ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)-1;
    (*currentVertex)+=2;
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
 */
void constructDoubleroofLongBuilding(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    int i, start;
    block->connectionVertices[0] = (*currentVertex)+1;
    block->connectionVertices[1] = (*currentVertex)+(block->parameter-1)*4+3;

    start = *currentVertex;

    ddgraph->adjList[*currentVertex][0] = (*currentVertex)+(block->parameter-1)*4+2;
    ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][2] = (*currentVertex)+2;
    ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)+3;
    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        ddgraph->adjList[*currentVertex][0] = (*currentVertex)-1;
        ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
        ddgraph->adjList[*currentVertex][2] = (*currentVertex)+2;

        ddgraph->adjList[(*currentVertex)+1][0] = (*currentVertex)-2;
        ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
        ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)+3;

        ddgraph->adjList[(*currentVertex)+2][0] = *currentVertex;
        ddgraph->adjList[(*currentVertex)+2][1] = (*currentVertex)+3;
        ddgraph->adjList[(*currentVertex)+2][2] = (*currentVertex)+4;

        ddgraph->adjList[(*currentVertex)+3][0] = (*currentVertex)+1;
        ddgraph->adjList[(*currentVertex)+3][1] = (*currentVertex)+2;
        ddgraph->adjList[(*currentVertex)+3][2] = (*currentVertex)+5;
        (*currentVertex)+=4;
    }

    ddgraph->adjList[*currentVertex][0] = start;
    ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][2] = (*currentVertex)-2;
    ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)-1;
    (*currentVertex)+=2;
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
    

    ddgraph->adjList[*currentVertex][0] = SEMIEDGE;
    ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][2] = (*currentVertex)+2;
    ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)+3;
    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        ddgraph->adjList[*currentVertex][0] = (*currentVertex)-1;
        ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
        ddgraph->adjList[*currentVertex][2] = (*currentVertex)+2;

        ddgraph->adjList[(*currentVertex)+1][0] = (*currentVertex)-2;
        ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
        ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)+3;

        ddgraph->adjList[(*currentVertex)+2][0] = *currentVertex;
        ddgraph->adjList[(*currentVertex)+2][1] = (*currentVertex)+3;
        ddgraph->adjList[(*currentVertex)+2][2] = (*currentVertex)+4;

        ddgraph->adjList[(*currentVertex)+3][0] = (*currentVertex)+1;
        ddgraph->adjList[(*currentVertex)+3][1] = (*currentVertex)+2;
        ddgraph->adjList[(*currentVertex)+3][2] = (*currentVertex)+5;
        (*currentVertex)+=4;
    }

    ddgraph->adjList[*currentVertex][0] = SEMIEDGE;
    ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][2] = (*currentVertex)-2;
    ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)-1;
    (*currentVertex)+=2;
}

/*
 * Constructs a locked diagonal chain LDC(n). This block has 2 connectors
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

    ddgraph->adjList[*currentVertex][0] = (*currentVertex)+(block->parameter-1)*4+3;
    ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][2] = (*currentVertex)+2;
    ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)+3;
    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        ddgraph->adjList[*currentVertex][0] = (*currentVertex)-1;
        ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
        ddgraph->adjList[*currentVertex][2] = (*currentVertex)+2;

        ddgraph->adjList[(*currentVertex)+1][0] = (*currentVertex)-2;
        ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
        ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)+3;

        ddgraph->adjList[(*currentVertex)+2][0] = *currentVertex;
        ddgraph->adjList[(*currentVertex)+2][1] = (*currentVertex)+3;
        ddgraph->adjList[(*currentVertex)+2][2] = (*currentVertex)+4;

        ddgraph->adjList[(*currentVertex)+3][0] = (*currentVertex)+1;
        ddgraph->adjList[(*currentVertex)+3][1] = (*currentVertex)+2;
        ddgraph->adjList[(*currentVertex)+3][2] = (*currentVertex)+5;
        (*currentVertex)+=4;
    }

    ddgraph->adjList[*currentVertex][0] = SEMIEDGE;
    ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][2] = (*currentVertex)-2;
    ddgraph->adjList[(*currentVertex)+1][0] = start;
    ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)-1;
    (*currentVertex)+=2;
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
 *
 */
void constructLockedDoubleroofHighBuilding(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    int i;
    block->connectionVertices[0] = *currentVertex;

    ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][2] = (*currentVertex)+2;
    ddgraph->adjList[(*currentVertex)+1][0] = SEMIEDGE;
    ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)+3;
    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        ddgraph->adjList[*currentVertex][0] = (*currentVertex)-1;
        ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
        ddgraph->adjList[*currentVertex][2] = (*currentVertex)+2;

        ddgraph->adjList[(*currentVertex)+1][0] = (*currentVertex)-2;
        ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
        ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)+3;

        ddgraph->adjList[(*currentVertex)+2][0] = *currentVertex;
        ddgraph->adjList[(*currentVertex)+2][1] = (*currentVertex)+3;
        ddgraph->adjList[(*currentVertex)+2][2] = (*currentVertex)+4;

        ddgraph->adjList[(*currentVertex)+3][0] = (*currentVertex)+1;
        ddgraph->adjList[(*currentVertex)+3][1] = (*currentVertex)+2;
        ddgraph->adjList[(*currentVertex)+3][2] = (*currentVertex)+5;
        (*currentVertex)+=4;
    }

    ddgraph->adjList[*currentVertex][0] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][2] = (*currentVertex)-2;
    ddgraph->adjList[(*currentVertex)+1][0] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)-1;
    (*currentVertex)+=2;
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

    ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][2] = (*currentVertex)+2;
    ddgraph->adjList[(*currentVertex)+1][0] = SEMIEDGE;
    ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)+3;
    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        ddgraph->adjList[*currentVertex][0] = (*currentVertex)-1;
        ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
        ddgraph->adjList[*currentVertex][2] = (*currentVertex)+2;

        ddgraph->adjList[(*currentVertex)+1][0] = (*currentVertex)-2;
        ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
        ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)+3;

        ddgraph->adjList[(*currentVertex)+2][0] = *currentVertex;
        ddgraph->adjList[(*currentVertex)+2][1] = (*currentVertex)+3;
        ddgraph->adjList[(*currentVertex)+2][2] = (*currentVertex)+4;

        ddgraph->adjList[(*currentVertex)+3][0] = (*currentVertex)+1;
        ddgraph->adjList[(*currentVertex)+3][1] = (*currentVertex)+2;
        ddgraph->adjList[(*currentVertex)+3][2] = (*currentVertex)+5;
        (*currentVertex)+=4;
    }

    ddgraph->adjList[*currentVertex][0] = SEMIEDGE;
    ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][2] = (*currentVertex)-2;
    ddgraph->adjList[(*currentVertex)+1][0] = SEMIEDGE;
    ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)-1;
    (*currentVertex)+=2;
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
 *
 */
void constructLockedDoubleroofLongBuilding(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    int i, start;
    block->connectionVertices[0] = (*currentVertex)+1;

    start = *currentVertex;

    ddgraph->adjList[*currentVertex][0] = (*currentVertex)+(block->parameter-1)*4+2;
    ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][2] = (*currentVertex)+2;
    ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)+3;
    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        ddgraph->adjList[*currentVertex][0] = (*currentVertex)-1;
        ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
        ddgraph->adjList[*currentVertex][2] = (*currentVertex)+2;

        ddgraph->adjList[(*currentVertex)+1][0] = (*currentVertex)-2;
        ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
        ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)+3;

        ddgraph->adjList[(*currentVertex)+2][0] = *currentVertex;
        ddgraph->adjList[(*currentVertex)+2][1] = (*currentVertex)+3;
        ddgraph->adjList[(*currentVertex)+2][2] = (*currentVertex)+4;

        ddgraph->adjList[(*currentVertex)+3][0] = (*currentVertex)+1;
        ddgraph->adjList[(*currentVertex)+3][1] = (*currentVertex)+2;
        ddgraph->adjList[(*currentVertex)+3][2] = (*currentVertex)+5;
        (*currentVertex)+=4;
    }

    ddgraph->adjList[*currentVertex][0] = start;
    ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
    ddgraph->adjList[*currentVertex][2] = (*currentVertex)-2;
    ddgraph->adjList[(*currentVertex)+1][0] = SEMIEDGE;
    ddgraph->adjList[(*currentVertex)+1][1] = *currentVertex;
    ddgraph->adjList[(*currentVertex)+1][2] = (*currentVertex)-1;
    (*currentVertex)+=2;
}

void constructPearlChain(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    int i;
    block->connectionVertices[0] = *currentVertex;
    for(i=0; i<block->parameter; i++){
        ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
        ddgraph->adjList[*currentVertex][2] = (*currentVertex)+1;
        ddgraph->adjList[*currentVertex][0] = (*currentVertex)-1;
        //for the first vertex this will be overwritten when making the connections
        (*currentVertex)++;
        ddgraph->adjList[*currentVertex][1] = (*currentVertex)-1;
        ddgraph->adjList[*currentVertex][2] = (*currentVertex)-1;
        ddgraph->adjList[*currentVertex][0] = (*currentVertex)+1;
        //for the last vertex this will be overwritten when making the connections
        (*currentVertex)++;
    }
}

void constructLockedPearlChain(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    int i;
    block->connectionVertices[0] = *currentVertex;
    for(i=0; i<block->parameter; i++){
        ddgraph->adjList[*currentVertex][1] = (*currentVertex)+1;
        ddgraph->adjList[*currentVertex][2] = (*currentVertex)+1;
        ddgraph->adjList[*currentVertex][0] = (*currentVertex)-1;
        //for the first vertex this will be overwritten when making the connections
        (*currentVertex)++;
        ddgraph->adjList[*currentVertex][1] = (*currentVertex)-1;
        ddgraph->adjList[*currentVertex][2] = (*currentVertex)-1;
        ddgraph->adjList[*currentVertex][0] = (*currentVertex)+1;
        //for the last vertex this will be overwritten when making the connections
        (*currentVertex)++;
    }
    ddgraph->adjList[(*currentVertex-1)][0] = SEMIEDGE;
}

void constructBarbWire(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    int i;
    block->connectionVertices[0] = *currentVertex;
    for(i=0; i<block->parameter; i++){
        ddgraph->adjList[*currentVertex][1] = SEMIEDGE;
        ddgraph->adjList[*currentVertex][2] = (*currentVertex)+1;
        ddgraph->adjList[*currentVertex][0] = (*currentVertex)-1;
        //for the first vertex this will be overwritten when making the connections
        (*currentVertex)++;
        ddgraph->adjList[*currentVertex][1] = SEMIEDGE;
        ddgraph->adjList[*currentVertex][2] = (*currentVertex)-1;
        ddgraph->adjList[*currentVertex][0] = (*currentVertex)+1;
        //for the last vertex this will be overwritten when making the connections
        (*currentVertex)++;
    }
}

void constructLockedBarbWire(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    int i;
    block->connectionVertices[0] = *currentVertex;
    for(i=0; i<2*(block->parameter); i++){
        ddgraph->adjList[*currentVertex][1] = SEMIEDGE;
        ddgraph->adjList[*currentVertex][2] = (*currentVertex)+1;
        ddgraph->adjList[*currentVertex][0] = (*currentVertex)-1;
        //for the first vertex this will be overwritten when making the connections
        (*currentVertex)++;
    }
    ddgraph->adjList[(*currentVertex)-1][2] = SEMIEDGE;
}

void constructQ4(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph){
    ddgraph->adjList[*currentVertex][1] = ddgraph->adjList[*currentVertex][2] = SEMIEDGE;
    block->connectionVertices[0] = *currentVertex;
    (*currentVertex)++;
}

void constructBuildingBlockListAsGraph(BBLOCK* blocks, int buildingBlockCount, DDGRAPH *ddgraph){
    int i, currentVertex=0;

    for (i = 0; i < buildingBlockCount; i++) {
        if((blocks+i)->type==1){
            if((blocks+i)->component==0){
                constructHub(&currentVertex, (blocks+i), ddgraph);
            } else if((blocks+i)->component==1){
                constructLockedHub(&currentVertex, (blocks+i), ddgraph);
            } else if((blocks+i)->component==2){
                constructDiagonalChain(&currentVertex, (blocks+i), ddgraph);
            } else if((blocks+i)->component==3){
                constructDoubleroofHighBuilding(&currentVertex, (blocks+i), ddgraph);
            } else if((blocks+i)->component==4){
                constructOpenroofHighBuilding(&currentVertex, (blocks+i), ddgraph);
            } else if((blocks+i)->component==5){
                constructDoubleroofLongBuilding(&currentVertex, (blocks+i), ddgraph);
            } else if((blocks+i)->component==6){
                constructOpenroofLongBuilding(&currentVertex, (blocks+i), ddgraph);
            } else if((blocks+i)->component==7){
                constructLockedDiagonalChain(&currentVertex, (blocks+i), ddgraph);
            } else if((blocks+i)->component==8){
                constructLockedDoubleroofHighBuilding(&currentVertex, (blocks+i), ddgraph);
            } else if((blocks+i)->component==9){
                constructLockedOpenroofHighBuilding(&currentVertex, (blocks+i), ddgraph);
            } else if((blocks+i)->component==10){
                constructLockedDoubleroofLongBuilding(&currentVertex, (blocks+i), ddgraph);
            } else {
                fprintf(stderr, "Illegal component number for type 1: %d (valid numbers from 0 to 10)\n", (blocks+i)->component);
                exit(EXIT_FAILURE);
            }
        } else if((blocks+i)->type==2){
            if((blocks+i)->component==0){
                constructPearlChain(&currentVertex, (blocks+i), ddgraph);
            } else if((blocks+i)->component==1){
                constructLockedPearlChain(&currentVertex, (blocks+i), ddgraph);
            } else {
                fprintf(stderr, "Illegal component number for type 2: %d (valid numbers are 0 or 1)\n", (blocks+i)->component);
                exit(EXIT_FAILURE);
            }
        } else if((blocks+i)->type==3){
            if((blocks+i)->component==0){
                constructBarbWire(&currentVertex, (blocks+i), ddgraph);
            } else if((blocks+i)->component==1){
                constructLockedBarbWire(&currentVertex, (blocks+i), ddgraph);
            } else {
                fprintf(stderr, "Illegal component number for type 3: %d (valid numbers are 0 or 1)\n", (blocks+i)->component);
                exit(EXIT_FAILURE);
            }
        } else if((blocks+i)->type==4){
            constructQ4(&currentVertex, (blocks+i), ddgraph);
        } else {
            fprintf(stderr, "Illegal component type: %d (valid types from 1 to 4)\n", (blocks+i)->type);
            exit(EXIT_FAILURE);
        }
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

void connectComponentList(int vertexCount){
    int blockCount = 0;
    BBLOCK *blocks = constructComponentList(&blockCount);
    DDGRAPH *graph = getNewDDGraph(vertexCount);

    constructBuildingBlockListAsGraph(blocks, blockCount, graph);
    
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
        int i, j;
        componentListsCount++;
/*        for(i = 0; i < Q1TypeComponentsCount; i++){
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
        fprintf(stderr, "| %d\n", Q4ComponentCount);*/
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

