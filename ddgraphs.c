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

#ifdef _DEBUG
void show_stackframe() {
  void *trace[16];
  char **messages = (char **)NULL;
  int i, trace_size = 0;

  trace_size = backtrace(trace, 16);
  messages = backtrace_symbols(trace, trace_size);
  fprintf(stderr, "[bt] Execution path:\n");
  for (i=0; i<trace_size; ++i){
	fprintf(stderr, "[bt] %s\n", messages[i]);
  }
}
#endif



#ifdef _DEBUGMETHODS
void printGenerators (DDGRAPH *ddgraph, int printDepth){
    int i, j;
    fprintf(stderr, "Generators at level %2d:\n", printDepth);
    for(i=0; i<numberOfGenerators[printDepth]; i++){
        boolean done[ddgraph->underlyingGraph->nv];
        for(j=0; j<ddgraph->underlyingGraph->nv; j++){
            done[j]=FALSE;
        }
        fprintf(stderr, "  Generator %2d: ", i+1);
        for(j=0; j<ddgraph->underlyingGraph->nv; j++){
            if(!done[j] && (*(automorphismGroupGenerators + printDepth))[i][j] != j){
                fprintf(stderr, "(%d", j);
                done[j] = TRUE;
                int next = (*(automorphismGroupGenerators + printDepth))[i][j];
                done[next] = TRUE;
                int teller=0;
                while(next!=j){
                    fprintf(stderr, " %d", next);
                    next = (*(automorphismGroupGenerators + printDepth))[i][next];
                    done[next] = TRUE;
                    teller++;
                    if(teller==MAXN) break;
                }
                fprintf(stderr, ") ");
            }
        }
        fprintf(stderr, "\n");
    }
}

void printOrbits(int n, int printDepth){
    int i;
    fprintf(stderr, "Orbits: [%d", vertexOrbits[printDepth][0]);
    for(i=1; i<n; i++){
        fprintf(stderr, ", %d", vertexOrbits[printDepth][i]);
    }
    fprintf(stderr, "]\n");
}

void printConnectorOrbits(BBLOCK* blocks, int buildingBlockCount, int printDepth){
    int i, j;
    fprintf(stderr, "Orbits of connection vertices:\n");
    for(i=0; i<buildingBlockCount; i++){
        for(j=0; j<(blocks+i)->connectorCount; j++){
            if((blocks+i)->connectionVertices[j] == vertexOrbits[printDepth][(blocks+i)->connectionVertices[j]]){
                fprintf(stderr, "Orbit representative: %2d (size: %2d)\n", (blocks+i)->connectionVertices[j],
                        vertexOrbitsSizes[printDepth][(blocks+i)->connectionVertices[j]]);
            }
        }
    }
}

void printComponentList(){
    int i, j;
    for(i = 0; i < Q1TypeComponentsCount; i++){
        fprintf(stderr, "(%d", Q1TypeComponentsComponentCount[i][0]);
        for(j = 1; j < maximumQ1TypeComponents; j++){
            fprintf(stderr, " %d", Q1TypeComponentsComponentCount[i][j]);
        }
        fprintf(stderr, ")");
    }
    fprintf(stderr, "| ");
    for(i = 0; i < Q2TypeComponentsCount; i++){
        fprintf(stderr, "(%d", Q2TypeComponentsComponentCount[i][0]);
        for(j = 1; j < maximumQ2TypeComponents; j++){
            fprintf(stderr, " %d", Q2TypeComponentsComponentCount[i][j]);
        }
        fprintf(stderr, ")");
    }
    fprintf(stderr, "| ");
    for(i = 0; i < Q3TypeComponentsCount; i++){
        fprintf(stderr, "(%d", Q3TypeComponentsComponentCount[i][0]);
        for(j = 1; j < maximumQ3TypeComponents; j++){
            fprintf(stderr, " %d", Q3TypeComponentsComponentCount[i][j]);
        }
        fprintf(stderr, ")");
    }
    fprintf(stderr, "| %d\n", Q4ComponentCount);
}

char (*blockName[ComponentsTypesCount]);

void printBlockName(int count, int blockNumber, int parameter, boolean first){
    if(!first){
        fprintf(stderr, ", ");
    }
    fprintf(stderr, "%dx", count);
    fprintf(stderr, blockName[blockNumber], parameter+1);
}

void printHumanReadableComponentList(){
    int i, j;
    boolean first = TRUE;
    for(i = 0; i < Q1TypeComponentsCount; i++){
        for(j = 0; j < maximumQ1TypeComponents; j++){
            if(Q1TypeComponentsComponentCount[i][j]){
                printBlockName(Q1TypeComponentsComponentCount[i][j], i, j, first);
                first = FALSE;
            }
        }
    }
    for(i = 0; i < Q2TypeComponentsCount; i++){
        for(j = 0; j < maximumQ2TypeComponents; j++){
            if(Q2TypeComponentsComponentCount[i][j]){
                printBlockName(Q2TypeComponentsComponentCount[i][j], Q1TypeComponentsCount + i, j, first);
                first = FALSE;
            }
        }
    }
    for(i = 0; i < Q3TypeComponentsCount; i++){
        for(j = 0; j < maximumQ3TypeComponents; j++){
            if(Q3TypeComponentsComponentCount[i][j]){
                printBlockName(Q3TypeComponentsComponentCount[i][j], Q1TypeComponentsCount + Q2TypeComponentsCount + i, j, first);
                first = FALSE;
            }
        }
    }
    if(Q4ComponentCount){
        printBlockName(Q4ComponentCount, ComponentsTypesCount - 1, 0, first);
    }
    fprintf(stderr, "\n");
}
#endif

void printDDGraph(DDGRAPH *graph){
    fprintf(stderr, "DDGRAPH %p\n", graph);
    fprintf(stderr, "================\n");
    fprintf(stderr, "order      %d\n", graph->order);
    //fprintf(stderr, "semi-edges %d\n", graph->semiEdges);
    fprintf(stderr, "dummies    %d\n", graph->dummyVertexCount);
    fprintf(stderr, "\n");
    fprintf(stderr, "nv         %d\n", graph->underlyingGraph->nv);
    fprintf(stderr, "nde        %d\n", graph->underlyingGraph->nde);
#ifdef _DEBUG
    fprintf(stderr, "\n");
    fprintf(stderr, "vlen       %d\n", graph->underlyingGraph->vlen);
    fprintf(stderr, "dlen       %d\n", graph->underlyingGraph->dlen);
    fprintf(stderr, "elen       %d\n", graph->underlyingGraph->elen);
#endif

    int i,j;
    fprintf(stderr, "v :");
    for(i=0; i<graph->underlyingGraph->nv; i++){
        fprintf(stderr, "%d ", graph->underlyingGraph->v[i]);
    }
    fprintf(stderr, "\n");

    for(i=0; i<graph->underlyingGraph->nv; i++){
        fprintf(stderr, "%d :", i);
        for(j=graph->underlyingGraph->v[i]; j < graph->underlyingGraph->v[i] + graph->underlyingGraph->d[i]; j++){
            int s = SEMIEDGE;
            if(graph->underlyingGraph->e[j] == s){
                fprintf(stderr, "S ");
            } else {
                fprintf(stderr, "%d ", graph->underlyingGraph->e[j]);
            }
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

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
            ddgraph->underlyingGraph->v[i]=3*order + 2*(i-order);
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
            ddgraph->underlyingGraph->v[i]=3*order + 2*(i-order);
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
BBLOCK *initBuildingBlock(BBLOCK* block, int type, int component, int parameter, int id){
    block->id = id;
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
    int i;
    for(i=0; i<block->connectorCount; i++){
        block->connections[i]=NULL;
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
BBLOCK *getBuildingBlock(int type, int component, int parameter, int id){
    BBLOCK *block = (BBLOCK *) malloc(sizeof(BBLOCK));
    return initBuildingBlock(block, type, component, parameter, id);
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

//union-find

/* Searches the given forest to find the root of the given element and performs
 * path compression at the same time.
 */
int findRootOfElement(int *forest, int element) {
    /* find with path-compression: once the recursion finds the root, the parent
     * for each element it met on its search is set to the root
     */
    if(element!=forest[element]){
        forest[element]=findRootOfElement(forest, forest[element]);
    }
    return forest[element];
}

/* Connects element1 and element2 in the given forest. This is done by making
 * the root of the largest tree the parent of the root of the smallest tree.
 */
void unionElements(int *forest, int *treeSizes, int *numberOfComponents, int element1, int element2){
    int root1 = findRootOfElement(forest, element1);
    int root2 = findRootOfElement(forest, element2);

    DEBUGMSG("Union:")
    DEBUGDUMP(element1, "%d")
    DEBUGDUMP(element2, "%d")
    DEBUGDUMP(root1, "%d")
    DEBUGDUMP(root2, "%d")

    //if these elements are already in the same tree, we can just return
    if(root1==root2) return;

    if(treeSizes[root1]<treeSizes[root2]){
        forest[root1]=root2;
        treeSizes[root2]+=treeSizes[root1];
    } else {
        forest[root2]=root1;
        treeSizes[root1]+=treeSizes[root2];
    }
    (*numberOfComponents)--;
}

/* Determines the vertex orbits and their sizes based upon the provided generators
 */
void determineVertexOrbits(int vertexCount, int *vertexOrbits, int *orbitSizes,
        int *orbitCount, permutation (*generators)[MAXN][MAXN], int generatorCount){
    int i, j;

    //initialize variables
    for(i=0; i<vertexCount; i++){
        vertexOrbits[i] = i;
        orbitSizes[i] = 1;
    }
    *orbitCount = vertexCount;

    //return in case of a trivial automorphism group
    if(generatorCount==0) return;

    for(i=0; i<generatorCount; i++){
        for(j=0; j<vertexCount; j++){
            unionElements(vertexOrbits, orbitSizes, orbitCount, j, (*generators)[i][j]);
        }
    }

    //make sure that each element is connected to its root
    for(i = 0; i < vertexCount; i++){
        findRootOfElement(vertexOrbits, i);
    }
}

void callNautyForGraph(DDGRAPH *ddgraph){
    int i, j;
    //calculate automorphisms of the new graph
    int currentOrbits[ddgraph->underlyingGraph->nv];

    for(i=0; i<ddgraph->underlyingGraph->nv; i++){
        nautyPtn[i] = 1;
    }
    int counter = 0;
    for(j = 2; j>=0; j--){
        for(i=0; i<ddgraph->order; i++){
            if(ddgraph->semiEdges[i]==j){
                nautyLabelling[counter] = i;
                counter++;
            }
        }
        if(counter>0){
            nautyPtn[counter-1]=0;
        }
    }
    for(i=ddgraph->order; i<ddgraph->underlyingGraph->nv; i++){
        nautyLabelling[i] = i;
    }
    nautyPtn[ddgraph->underlyingGraph->nv-1]=0;

#ifdef _DEBUG
    fprintf(stderr, "nautyLab: ");
    for(i=0; i<ddgraph->underlyingGraph->nv; i++){
        fprintf(stderr, "%d ", nautyLabelling[i]);
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "nautyPtn: ");
    for(i=0; i<ddgraph->underlyingGraph->nv; i++){
        fprintf(stderr, "%d ", nautyPtn[i]);
    }
    fprintf(stderr, "\n");
#endif

    nauty((graph*)(ddgraph->underlyingGraph), nautyLabelling, nautyPtn, NULL, currentOrbits,
            &nautyOptions, &nautyStats, nautyWorkspace, 50 * MAXM, MAXM,
            ddgraph->underlyingGraph->nv, (graph*)&canonGraph);

#ifdef _DEBUG
    fprintf(stderr, "nauty Orbits: [%d", currentOrbits[0]);
    for(i=1; i<ddgraph->underlyingGraph->nv; i++){
        fprintf(stderr, ", %d", currentOrbits[i]);
    }
    fprintf(stderr, "]\n");
#endif
}

void determineVertexPairsOrbits(int (*vertexPairList)[2], int vertexPairListSize, int *vertexPairOrbits, int *orbitCount,
        permutation (*currentGenerators)[MAXN][MAXN] , int currentNumberOfGenerators){

    int i, j, k, temp;
    int orbitSize[vertexPairListSize];

    //initialization of the variables
    for(i=0; i<vertexPairListSize; i++){
        vertexPairOrbits[i]=i;
        orbitSize[i]=1;
    }
    *orbitCount=vertexPairListSize;

    if(currentNumberOfGenerators==0){
        return;
    }

    permutation *permutation;
    int pair[2];
    for(i = 0; i < currentNumberOfGenerators; i++) {
        permutation = (*currentGenerators)[i];

        for(j = 0; j<vertexPairListSize; j++){
            //apply permutation to current vertex pair
            pair[0] = permutation[vertexPairList[j][0]];
            pair[1] = permutation[vertexPairList[j][1]];

            //canonical form of vertex pair
            if(pair[0]>pair[1]){
                temp = pair[1];
                pair[1] = pair[0];
                pair[0] = temp;
            }

            //search the pair in the list
            for(k = 0; k<vertexPairListSize; k++){
                if(pair[0] == vertexPairList[k][0] && pair[1] == vertexPairList[k][1]){
                    unionElements(vertexPairOrbits, orbitSize, orbitCount, j, k);
                    break; //the list of pairs doesn't contain any duplicates so we can stop
                }
            }
        }
    }

    //make sure that each element is connected to its root
    for(i = 0; i < vertexPairListSize; i++){
        findRootOfElement(vertexPairOrbits, i);
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
void constructHub(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    int i;
    block->connectionVertices[0] = *currentVertex;
    block->connectionVertices[1] = (*currentVertex)+1;
    block->connectionVertices[2] = (*currentVertex)+(block->parameter-1)*4+2;
    block->connectionVertices[3] = (*currentVertex)+(block->parameter-1)*4+3;
    vertexToConnector[*currentVertex] = 0;
    vertexToConnector[(*currentVertex)+1] = 1;
    vertexToConnector[(*currentVertex)+(block->parameter-1)*4+2] = 2;
    vertexToConnector[(*currentVertex)+(block->parameter-1)*4+3] = 3;
    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;

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
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;

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

    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;
    
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
        //mirror-symmetry along both axis

        //long axis
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

        //short axis
        generator = getIdentity(ddgraph->underlyingGraph->nv);

        int vertexLeft = block->connectionVertices[0];
        int vertexRight = block->connectionVertices[2];

        for(i = 0; i < block->parameter; i++){
            generator[vertexLeft] = vertexRight;
            generator[vertexRight] = vertexLeft;
            generator[vertexLeft+1] = vertexRight+1;
            generator[vertexRight+1] = vertexLeft+1;
            vertexLeft+=2;
            vertexRight-=2;
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
void constructLockedHub(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    int i;
    block->connectionVertices[0] = (*currentVertex)+1;
    block->connectionVertices[1] = (*currentVertex)+(block->parameter-1)*4+2;
    block->connectionVertices[2] = (*currentVertex)+(block->parameter-1)*4+3;
    vertexToConnector[(*currentVertex)+1] = 0;
    vertexToConnector[(*currentVertex)+(block->parameter-1)*4+2] = 1;
    vertexToConnector[(*currentVertex)+(block->parameter-1)*4+3] = 2;
    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;

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
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;

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

    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;

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
 * Constructs a double locked hub DLH(n). This block has 2 connectors arranged
 * as follows:
 *
 *                           1
 *          \               /
 *           o---o-...-o---o
 *           |   |     |   |
 *           |   |     |   |
 *           o---o-...-o---o
 *          /               \
 *         0
 *
 */
void constructDoubleLockedHub(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    int i, start;
    block->connectionVertices[0] = (*currentVertex)+1;
    block->connectionVertices[1] = (*currentVertex)+(block->parameter-1)*4+2;
    vertexToConnector[(*currentVertex)+1] = 0;
    vertexToConnector[(*currentVertex)+(block->parameter-1)*4+2] = 1;
    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;

    start = *currentVertex;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    edges[positions[*currentVertex]+0] = (*currentVertex)+1;
    edges[positions[*currentVertex]+1] = (*currentVertex)+2;
    edges[positions[*currentVertex]+2] = SEMIEDGE;
    degrees[*currentVertex] = 2;
    ddgraph->semiEdges[*currentVertex] = 1;

    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //connection points are at position 0

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;

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

    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;

    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)-2;
    degrees[*currentVertex] = 2;
    positions[*currentVertex]++; //connection points are at position 0

    edges[positions[(*currentVertex)+1]+0] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+1] = (*currentVertex)-1;
    edges[positions[(*currentVertex)+1]+2] = SEMIEDGE;
    degrees[(*currentVertex)+1] = 2;
    ddgraph->semiEdges[(*currentVertex)+1] = 1;

    (*currentVertex)+=2;
}

void storeDoubleLockedHubAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
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
    } else {
        //rotation of 180
        permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

        int vertexBottom = block->connectionVertices[0];
        int vertexTop = block->connectionVertices[1];

        int i;
        for(i = 0; i < 2*(block->parameter); i++){
            generator[vertexBottom] = vertexTop;
            generator[vertexTop] = vertexBottom;
            vertexBottom+=2;
            vertexTop-=2;
        }

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    }
}

void storeDoubleLockedHubsMapping(BBLOCK *block1, BBLOCK *block2, DDGRAPH *ddgraph){
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
void constructDiagonalChain(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    int i, start;
    block->connectionVertices[0] = (*currentVertex)+1;
    block->connectionVertices[1] = (*currentVertex)+(block->parameter-1)*4+2;
    vertexToConnector[(*currentVertex)+1] = 0;
    vertexToConnector[(*currentVertex)+(block->parameter-1)*4+2] = 1;
    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;

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
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;

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

    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;

    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)-2;
    degrees[*currentVertex] = 2;
    positions[*currentVertex]++; //connection points are at position 0

    edges[positions[(*currentVertex)+1]+0] = start;
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
    } else {
        //rotation of 180
        permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

        int vertexBottom = block->connectionVertices[0];
        int vertexTop = block->connectionVertices[1];

        int i;
        for(i = 0; i < 2*(block->parameter); i++){
            generator[vertexBottom] = vertexTop;
            generator[vertexTop] = vertexBottom;
            vertexBottom+=2;
            vertexTop-=2;
        }

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    }
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
void constructDoubleroofHighBuilding(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    int i;
    block->connectionVertices[0] = *currentVertex;
    block->connectionVertices[1] = (*currentVertex)+1;
    vertexToConnector[*currentVertex] = 0;
    vertexToConnector[(*currentVertex)+1] = 1;
    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;

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
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;

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

    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;
    vertexToBlock[dummyVertex] = block->id;

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

    int dummyVertex1 = ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[vertexBlock1 + 4*(block1->parameter)-1]+0];
    int dummyVertex2 = ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[vertexBlock2 + 4*(block2->parameter)-1]+0];

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
void constructOpenroofHighBuilding(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    int i;
    block->connectionVertices[0] = *currentVertex;
    block->connectionVertices[1] = (*currentVertex)+1;
    vertexToConnector[*currentVertex] = 0;
    vertexToConnector[(*currentVertex)+1] = 1;
    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;


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
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;

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

    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;

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
void constructDoubleroofLongBuilding(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    int i, start;
    block->connectionVertices[0] = (*currentVertex)+1;
    block->connectionVertices[1] = (*currentVertex)+(block->parameter-1)*4+3;
    vertexToConnector[(*currentVertex)+1] = 0;
    vertexToConnector[(*currentVertex)+(block->parameter-1)*4+3] = 1;
    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;


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
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;

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

    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;

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
void constructOpenroofLongBuilding(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    int i;
    block->connectionVertices[0] = (*currentVertex)+1;
    block->connectionVertices[1] = (*currentVertex)+(block->parameter-1)*4+3;
    vertexToConnector[(*currentVertex)+1] = 0;
    vertexToConnector[(*currentVertex)+(block->parameter-1)*4+3] = 1;
    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;
    

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
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;


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
    
    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;

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
void constructLockedDiagonalChain(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    int i, start;
    block->connectionVertices[0] = (*currentVertex)+1;
    vertexToConnector[(*currentVertex)+1] = 0;
    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;

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
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;

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
    
    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;

    edges[positions[*currentVertex]+0] = SEMIEDGE;
    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)-2;
    degrees[*currentVertex] = 2;
    positions[*currentVertex]++; //semi-edge is at position 0
    ddgraph->semiEdges[*currentVertex] = 1;

    edges[positions[(*currentVertex)+1]+0] = start;
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
void constructLockedDoubleroofHighBuilding(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    int i;
    block->connectionVertices[0] = *currentVertex;
    vertexToConnector[*currentVertex] = 0;
    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;

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
    ddgraph->semiEdges[(*currentVertex)+1] = 1;

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;

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

    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;
    vertexToBlock[dummyVertex] = block->id;

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

    int dummyVertex1 = ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[vertexBlock1 + 4*(block1->parameter)-1]+0];
    int dummyVertex2 = ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[vertexBlock2 + 4*(block2->parameter)-1]+0];

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
void constructLockedOpenroofHighBuilding(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    int i;
    block->connectionVertices[0] = *currentVertex;
    vertexToConnector[*currentVertex] = 0;
    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;

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
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;
        
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

    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;

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
void constructLockedDoubleroofLongBuilding(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    int i, start;
    block->connectionVertices[0] = (*currentVertex)+1;
    vertexToConnector[(*currentVertex)+1] = 0;
    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;

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
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;

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
    
    vertexToBlock[*currentVertex] = block->id;
    vertexToBlock[(*currentVertex)+1] = block->id;

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
void constructPearlChain(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    int i;
    block->connectionVertices[0] = *currentVertex;
    block->connectionVertices[1] = *currentVertex + 1 + 2*(block->parameter-1);
    vertexToConnector[*currentVertex] = 0;
    vertexToConnector[(*currentVertex)+1+2*(block->parameter-1)] = 1;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    for(i=0; i<block->parameter; i++){
        int dummyVertex = ddgraph->order + ddgraph->dummyVertexCount;
        ddgraph->dummyVertexCount++;

        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[dummyVertex] = block->id;

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
    dummyLeft = ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[vertexLeft]-1+2];
    dummyRight = ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[vertexRight]-1+2];
    //-1 because the start of the edges of the connections points was moved
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
    dummyBlock1 = ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[vertexBlock1]-1+2];
    vertexBlock2 = block2->connectionVertices[0];
    dummyBlock2 = ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[vertexBlock2]-1+2];

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
void constructLockedPearlChain(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    int i;
    block->connectionVertices[0] = *currentVertex;
    vertexToConnector[*currentVertex] = 0;
    
    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    for(i=0; i<block->parameter; i++){
        int dummyVertex = ddgraph->order + ddgraph->dummyVertexCount;
        ddgraph->dummyVertexCount++;

        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[dummyVertex] = block->id;

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
    degrees[(*currentVertex)-1] = 2;
    positions[(*currentVertex)-1]++;

    degrees[block->connectionVertices[0]] = 2;
    positions[block->connectionVertices[0]]++;
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
void constructBarbWire(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    int i;
    block->connectionVertices[0] = *currentVertex;
    block->connectionVertices[1] = *currentVertex + 1 + 2*(block->parameter-1);
    vertexToConnector[*currentVertex] = 0;
    vertexToConnector[(*currentVertex)+1+2*(block->parameter-1)] = 1;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    for(i=0; i<block->parameter; i++){
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;

        edges[positions[*currentVertex]+0] = (*currentVertex)-1;
        //for the first vertex this will be overwritten when making the connections
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = SEMIEDGE;
        ddgraph->semiEdges[(*currentVertex)] = 1;
        degrees[*currentVertex] = 2;
        (*currentVertex)++;

        edges[positions[*currentVertex]+0] = (*currentVertex)+1;
        //for the last vertex this will be overwritten when making the connections
        edges[positions[*currentVertex]+1] = (*currentVertex)-1;
        edges[positions[*currentVertex]+2] = SEMIEDGE;
        ddgraph->semiEdges[(*currentVertex)] = 1;
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
void constructLockedBarbWire(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    int i;
    block->connectionVertices[0] = *currentVertex;
    vertexToConnector[*currentVertex] = 0;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    for(i=0; i<block->parameter; i++){
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;

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
void constructQ4(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[(*currentVertex)] + 1]
            = ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[(*currentVertex)] + 2] = SEMIEDGE;
    ddgraph->underlyingGraph->v[(*currentVertex)]++;
    ddgraph->underlyingGraph->d[(*currentVertex)] = 0;
    ddgraph->semiEdges[(*currentVertex)] = 2;

    //it is not necessary to adjust ddgraph->underlyingGraph->v because it will
    //not be used when the degree is 0.
    block->connectionVertices[0] = *currentVertex;
    vertexToConnector[*currentVertex] = 0;
    vertexToBlock[*currentVertex] = block->id;
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

void constructBuildingBlockListAsGraph(BBLOCK* blocks, int buildingBlockCount, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    int i, currentVertex=0;

    for (i = 0; i < buildingBlockCount; i++) {
        int number = buildingBlockTypeToNumber(blocks + i);
        (*constructBlock[number])(&currentVertex, (blocks+i), ddgraph, vertexToBlock, vertexToConnector);
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
                initBuildingBlock(block + currentBlock, 1, i, j+1, currentBlock);
                currentBlock++;
            }
        }
    }
    for(i = 0; i < Q2TypeComponentsCount; i++){
        for(j = 0; j < maximumQ2TypeComponents; j++){
            for(k = 0; k < Q2TypeComponentsComponentCount[i][j]; k++){
                initBuildingBlock(block + currentBlock, 2, i, j+1, currentBlock);
                currentBlock++;
            }
        }
    }
    for(i = 0; i < Q3TypeComponentsCount; i++){
        for(j = 0; j < maximumQ3TypeComponents; j++){
            for(k = 0; k < Q3TypeComponentsComponentCount[i][j]; k++){
                initBuildingBlock(block + currentBlock, 3, i, j+1, currentBlock);
                currentBlock++;
            }
        }
    }
    for(k = 0; k < Q4ComponentCount; k++){
        initBuildingBlock(block + currentBlock, 4, 0, 1, currentBlock);
        currentBlock++;
    }

    return block;
}

void storeInitialGenerators(BBLOCK* blocks, int buildingBlockCount, DDGRAPH *ddgraph){
    int firstBlocks[ComponentsTypesCount][MAXN]; //stores the index of the first block of each type
                         //-1 in case no block of this type was encountered yet
    int i,j;
    for(i=0;i<ComponentsTypesCount;i++) for(j=0;j<MAXN;j++) firstBlocks[i][j]=-1;

    for (i = 0; i < buildingBlockCount; i++) {
        int number = buildingBlockTypeToNumber(blocks+i);
        if(firstBlocks[number][(blocks+i)->parameter]==-1){
            firstBlocks[number][(blocks+i)->parameter]=i;
            if(storeBlockAutomorphismGenerators[number]!=NULL){
                (*storeBlockAutomorphismGenerators[number])(blocks+i, ddgraph);
            }
        } else {
            (*storeBlocksMapping[number])(blocks+firstBlocks[number][(blocks+i)->parameter], blocks+i, ddgraph);
        }
    }
}

inline void connectConnectors(BBLOCK* blocks, DDGRAPH *ddgraph,
        int *vertexToBlock, int *vertexToConnector, int connector1, int connector2){
    int block1 = vertexToBlock[connector1];
    int block2 = vertexToBlock[connector2];
    int connectorPosition1 = vertexToConnector[connector1];
    int connectorPosition2 = vertexToConnector[connector2];

    (blocks+block1)->connections[connectorPosition1] = blocks+block2;
    (blocks+block2)->connections[connectorPosition2] = blocks+block1;
    (blocks+block1)->targetConnector[connectorPosition1] = connectorPosition2;
    (blocks+block2)->targetConnector[connectorPosition2] = connectorPosition1;
    ddgraph->underlyingGraph->d[connector1]++;
    ddgraph->underlyingGraph->v[connector1]--;
    ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[connector1]] = connector2;
    ddgraph->underlyingGraph->d[connector2]++;
    ddgraph->underlyingGraph->v[connector2]--;
    ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[connector2]] = connector1;
    ddgraph->underlyingGraph->nde+=2;
}

inline void disconnectConnectors(BBLOCK* blocks, DDGRAPH *ddgraph,
        int *vertexToBlock, int *vertexToConnector, int connector1, int connector2){
    int block1 = vertexToBlock[connector1];
    int block2 = vertexToBlock[connector2];
    int connectorPosition1 = vertexToConnector[connector1];
    int connectorPosition2 = vertexToConnector[connector2];

    (blocks+block1)->connections[connectorPosition1] = NULL;
    (blocks+block2)->connections[connectorPosition2] = NULL;
    ddgraph->underlyingGraph->d[connector1]--;
    ddgraph->underlyingGraph->v[connector1]++;
    ddgraph->underlyingGraph->d[connector2]--;
    ddgraph->underlyingGraph->v[connector2]++;
    ddgraph->underlyingGraph->nde-=2;
}

boolean areNeighbouringConnections(int family, int parameter, int connection1, int connection2){
    DEBUGASSERT(Q1TypeComponentsCount==12)
    if(family == 4 || family == 5){
        return TRUE;
    } else if(family == 0){
        if(parameter==1){
            return (connection1 + connection2) != 3;
        } else {
            int min = connection1 < connection2 ? connection1 : connection2;
            int max = connection1 < connection2 ? connection2 : connection1;
            return (min==0 && max==1) || (min==2 && max==3);
        }
    } else if(family == 1){
        if(parameter==1){
            return (connection1 + connection2) != 1;
        } else {
            return !(connection1==0 || connection2==0);
        }
    }
    return FALSE;
}

boolean isLegalConnection(BBLOCK* blocks, int buildingBlockCount, DDGRAPH *ddgraph,
        int *vertexToBlock, int *vertexToConnector, int connector1, int connector2){
    DEBUGASSERT(Q1TypeComponentsCount==12)
    int i;
    int block1 = vertexToBlock[connector1];
    int block2 = vertexToBlock[connector2];
    int connectorPosition1 = vertexToConnector[connector1];
    int connectorPosition2 = vertexToConnector[connector2];
    int type1 = (blocks+block1)->type;
    int type2 = (blocks+block2)->type;

    if(type1 == type2){
        if(type1 == 2 || type1 == 3){
            //type 2 and type 3 can't be connected to blocks of the same type
            return FALSE;
        }

        if(type1 == 4 && ddgraph->order != 2){
            //connecting two type 4 blocks leads to a component without open connections
            return FALSE;
        }

        if(type1 == 1){
            int family1 = (blocks+block1)->component;
            int family2 = (blocks+block2)->component;

            if((family1 == 0 || family1 == 1 || family1 == 4 || family1 == 5) &&
                    (family2 == 0 || family2 == 1 || family2 == 4 || family2 == 5)){
                for(i=0; i<(blocks+block1)->connectorCount; i++){
                    if((blocks+block1)->connections[i]==(blocks+block2)){
                        //there already exists a connection between these blocks
                        //so we need to verify that that connection doesn't conflict with this one
                        if(areNeighbouringConnections(family1,
                                                    (blocks+block1)->parameter,
                                                     connectorPosition1, i) &&
                                areNeighbouringConnections(family2,
                                                    (blocks+block2)->parameter,
                                                     connectorPosition2,
                                                    (blocks+block1)->targetConnector[i])){
                            return FALSE;
                        }
                    }
                }
            }
        }
    }

    int currentBlock;

    int queue[buildingBlockCount];
    int visited[buildingBlockCount];
    for(i=0; i<buildingBlockCount; i++){
        visited[i]=FALSE;
    }
    int head = -1, tail = 1;
    queue[0] = block1;
    queue[1] = block2;
    visited[block1]=TRUE;
    visited[block2]=TRUE;
    
    while(head!=tail){
        head++;
        currentBlock = queue[head];
        for(i=0; i<(blocks+currentBlock)->connectorCount; i++){
            if((currentBlock==block1 && i == connectorPosition1)||
                    (currentBlock==block2 && i == connectorPosition2)){
                //ignore the new connections
                continue;
            }
            if((blocks+currentBlock)->connections[i]==NULL){
                //there still is an open connection in this component after making the new connection
                return TRUE;
            } else {
                int neighbour = (blocks+currentBlock)->connections[i]->id;
                if(!visited[neighbour]){
                    tail++;
                    queue[tail]=neighbour;
                    visited[neighbour]=TRUE;
                }
            }
        }
    }

    i=0;
    while(i<buildingBlockCount && visited[i]) i++;

    return i==buildingBlockCount;
}

boolean isCanonicalConnection(BBLOCK* blocks, int buildingBlockCount, DDGRAPH *ddgraph,
        int *vertexToBlock, int *vertexToConnector, int orbit, int depth, int connector1, int connector2){
    int i, j;

    int smallest = (connector1 < connector2 ? connector1 : connector2);
    int biggest = (connector1 < connector2 ? connector2 : connector1);
    int newConnection = -1;

    //construct list of made connections in this orbit
    int madeConnections[vertexOrbitsSizes[depth][orbit]][2];
    int madeConnectionsCount = 0;
    for(i=0; i<ddgraph->order; i++){ //dummy vertices can never be connectors
        if(vertexOrbits[depth][i] == orbit && (blocks+vertexToBlock[i])->connections[vertexToConnector[i]] != NULL){
            int opposite = (blocks+vertexToBlock[i])->connections[vertexToConnector[i]]->connectionVertices[(blocks+vertexToBlock[i])->targetConnector[vertexToConnector[i]]];
            if(vertexOrbits[depth][opposite] != orbit || i<opposite){
                int localSmallest = (i < opposite ? i : opposite);
                int localBiggest = (i < opposite ? opposite : i);
                if(localSmallest==smallest && localBiggest==biggest){
                    newConnection = madeConnectionsCount;
                }
                madeConnections[madeConnectionsCount][0]=localSmallest;
                madeConnections[madeConnectionsCount][1]=localBiggest;
                madeConnectionsCount++;
            }
        }
    }
    DEBUGASSERT(newConnection!=-1)

    //calculate automorphisms of the new graph
    int currentOrbits[ddgraph->underlyingGraph->nv];

    for(i=0; i<ddgraph->underlyingGraph->nv; i++){
        nautyPtn[i] = 1;
    }
    int counter = 0;
    for(j = 2; j>=0; j--){
        for(i=0; i<ddgraph->order; i++){
            if(ddgraph->semiEdges[i]==j){
                nautyLabelling[counter] = i;
                counter++;
            }
        }
        if(counter>0){
            nautyPtn[counter-1]=0;
        }
    }
    for(i=ddgraph->order; i<ddgraph->underlyingGraph->nv; i++){
        nautyLabelling[i] = i;
    }
    nautyPtn[ddgraph->underlyingGraph->nv-1]=0;

#ifdef _DEBUG
    printComponentList();
    fprintf(stderr, "nautyLab: ");
    for(i=0; i<ddgraph->underlyingGraph->nv; i++){
        fprintf(stderr, "%d ", nautyLabelling[i]);
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "nautyPtn: ");
    for(i=0; i<ddgraph->underlyingGraph->nv; i++){
        fprintf(stderr, "%d ", nautyPtn[i]);
    }
    fprintf(stderr, "\n");
#endif

    nauty((graph*)(ddgraph->underlyingGraph), nautyLabelling, nautyPtn, NULL, currentOrbits,
            &nautyOptions, &nautyStats, nautyWorkspace, 50 * MAXM, MAXM,
            ddgraph->underlyingGraph->nv, (graph*)&canonGraph);

#ifdef _DEBUG
    fprintf(stderr, "nauty Orbits: [%d", currentOrbits[0]);
    for(i=1; i<ddgraph->underlyingGraph->nv; i++){
        fprintf(stderr, ", %d", currentOrbits[i]);
    }
    fprintf(stderr, "]\n");
#endif


    int madeConnectionOrbits[madeConnectionsCount];
    int madeConnectionOrbitsCount;
    determineVertexPairsOrbits(madeConnections, madeConnectionsCount, madeConnectionOrbits,
            &madeConnectionOrbitsCount, automorphismGroupGenerators + connectionsMade, numberOfGenerators[connectionsMade]);

    int smallestConnection = 0;
    int temp1 = currentOrbits[madeConnections[smallestConnection][0]];
    int temp2 = currentOrbits[madeConnections[smallestConnection][1]];
    
    int smallestConnectionSmallest = temp1 < temp2 ? temp1 : temp2;
    int smallestConnectionBiggest = temp1 < temp2 ? temp2 : temp1;
    for(i=1; i<madeConnectionsCount; i++){
        if(madeConnectionOrbits[i]==i){
            int localTemp1 = currentOrbits[madeConnections[i][0]];
            int localTemp2 = currentOrbits[madeConnections[i][1]];

            int localSmallest = localTemp1 < localTemp2 ? localTemp1 : localTemp2;
            int localBiggest = localTemp1 < localTemp2 ? localTemp2 : localTemp1;

            if(localSmallest < smallestConnectionSmallest){
                smallestConnectionSmallest = localSmallest;
                smallestConnectionBiggest = localBiggest;
                smallestConnection = i;
            } else if(localSmallest == smallestConnectionSmallest && localBiggest < smallestConnectionBiggest){
                smallestConnectionSmallest = localSmallest;
                smallestConnectionBiggest = localBiggest;
                smallestConnection = i;
            }
        }
    }
    return madeConnectionOrbits[smallestConnection] == madeConnectionOrbits[newConnection];
}

void connectCompleteOrbit(BBLOCK* blocks, int buildingBlockCount, DDGRAPH *ddgraph,
        int *vertexToBlock, int *vertexToConnector, int orbit, int depth,
        int openConnectionsLeftInOrbit, int totalConnectionsLeft, boolean *freeConnectors){
    DEBUGTRACE_ENTER
    //find all pairs of possible connections with the current orbit
    int i, j;
    int possibleConnections[openConnectionsLeftInOrbit*totalConnectionsLeft][2];
    int possibleConnectionsCount = 0;
    for(i=0; i<ddgraph->order; i++){ //dummy vertices can never be connectors
        if(vertexOrbits[depth][i] == orbit && (blocks+vertexToBlock[i])->connections[vertexToConnector[i]] == NULL){
            for(j=0; j<ddgraph->order; j++){
                if(vertexOrbits[depth][j] == orbit && j<i) continue;
                if(i!=j && freeConnectors[j] && vertexToBlock[i]!=vertexToBlock[j]){
                    possibleConnections[possibleConnectionsCount][0] = (i<j ? i : j);
                    possibleConnections[possibleConnectionsCount][1] = (i<j ? j : i);
                    possibleConnectionsCount++;
                    DEBUGASSERT(possibleConnectionsCount <= openConnectionsLeftInOrbit*totalConnectionsLeft)
                }
            }
        }
    }
    int connectorOrbits[possibleConnectionsCount];
    int connectorOrbitCount;
    determineVertexPairsOrbits(possibleConnections, possibleConnectionsCount, connectorOrbits,
            &connectorOrbitCount, automorphismGroupGenerators + connectionsMade, numberOfGenerators[connectionsMade]);

    for(i=0; i<possibleConnectionsCount; i++){
        if(connectorOrbits[i]==i){
            //try connection
            int v1 = possibleConnections[i][0];
            int v2 = possibleConnections[i][1];

            //check if connecting v1 to v2 is legal
            if(isLegalConnection(blocks, buildingBlockCount, ddgraph, vertexToBlock, vertexToConnector, v1, v2)){
                
                //connect v1 to v2
                connectConnectors(blocks, ddgraph, vertexToBlock, vertexToConnector, v1, v2);
                connectionsMade++;
                numberOfGenerators[connectionsMade]=0;
                freeConnectors[v1]=FALSE;
                freeConnectors[v2]=FALSE;

                //check canonicity of operation and recurse
                if(isCanonicalConnection(blocks, buildingBlockCount, ddgraph, vertexToBlock, vertexToConnector, orbit, depth, v1, v2)){
                    int inCurrentOrbit = 0;
                    if(vertexOrbits[depth][v1] == orbit) inCurrentOrbit++;
                    if(vertexOrbits[depth][v2] == orbit) inCurrentOrbit++;

                    //
                    if(totalConnectionsLeft - 2 == 0){
                        graphsCount++;
                    } else if(openConnectionsLeftInOrbit - inCurrentOrbit == 0){
                        //we've just made the final connection for the orbit currently under consideration
                        findNextOrbitToConnect(blocks, buildingBlockCount, ddgraph, vertexToBlock, vertexToConnector, freeConnectors);
                    } else {
                        connectCompleteOrbit(blocks, buildingBlockCount, ddgraph,
                                vertexToBlock, vertexToConnector, orbit, depth,
                                openConnectionsLeftInOrbit - inCurrentOrbit, totalConnectionsLeft-2, freeConnectors);
                    }
                }

                //disconnect v1 from v2
                freeConnectors[v1]=TRUE;
                freeConnectors[v2]=TRUE;
                connectionsMade--;
                disconnectConnectors(blocks, ddgraph, vertexToBlock, vertexToConnector, v1, v2);
            }
        }
    }
    DEBUGTRACE_EXIT

}

void findNextOrbitToConnect(BBLOCK* blocks, int buildingBlockCount, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector, boolean *freeConnectors){
    DEBUGTRACE_ENTER
    //first we need the vertex orbits

    int orbitCount = 0;
    int i,j;
    
#ifdef _DEBUGINTERMEDIATE
    printGenerators(ddgraph, connectionsMade);
#endif
    determineVertexOrbits(
            ddgraph->underlyingGraph->nv,
            vertexOrbits[connectionsMade],
            vertexOrbitsSizes[connectionsMade],
            &orbitCount,
            automorphismGroupGenerators + connectionsMade,
            numberOfGenerators[connectionsMade]);
#ifdef _DEBUGINTERMEDIATE
    printConnectorOrbits(blocks, buildingBlockCount, connectionsMade);
#endif

    int minimumOrbitSize = ddgraph->order + 1;
    int smallestOrbit = INT_MAX;
    int connectionCount = 0;

    for(i=0; i<buildingBlockCount; i++){
        for(j=0; j<(blocks+i)->connectorCount; j++){
            if((blocks+i)->connections[j]==NULL){
                connectionCount++;
                //only look at connections that haven't been made
                if((blocks+i)->connectionVertices[j] == vertexOrbits[connectionsMade][(blocks+i)->connectionVertices[j]]){
                    if(vertexOrbitsSizes[connectionsMade][(blocks+i)->connectionVertices[j]] < minimumOrbitSize){
                        minimumOrbitSize = vertexOrbitsSizes[connectionsMade][(blocks+i)->connectionVertices[j]];
                        smallestOrbit = (blocks+i)->connectionVertices[j];
                    }
                }
            }
        }
    }
    connectCompleteOrbit(blocks, buildingBlockCount, ddgraph, vertexToBlock,
            vertexToConnector, smallestOrbit, connectionsMade, minimumOrbitSize,
            connectionCount, freeConnectors);
    DEBUGTRACE_EXIT
}

void connectComponentList(int vertexCount){
    DEBUGTRACE_ENTER
    int blockCount = 0;
    //create an array of blocks based upon the numbers in the global arrays
    BBLOCK *blocks = constructComponentList(&blockCount);
    DDGRAPH *ddgraph = getNewDDGraph(vertexCount);
    int vertexToBlock[ddgraph->underlyingGraph->vlen];
    int vertexToConnector[ddgraph->underlyingGraph->vlen];

    //reset counters
    connectionsMade=0;
    numberOfGenerators[connectionsMade]=0;

    //convert the blocks to a disconnected graph
    constructBuildingBlockListAsGraph(blocks, blockCount, ddgraph, vertexToBlock, vertexToConnector);
    //store the generators of the automorphism group of this disconnected graph
    storeInitialGenerators(blocks, blockCount, ddgraph);

    boolean freeConnectors[vertexCount];
    int i, j;
    for(i=0; i<vertexCount; i++){
        freeConnectors[i] = FALSE;
    }
    for(i=0; i<blockCount; i++){
        for(j=0; j<(blocks+i)->connectorCount; j++){
            freeConnectors[(blocks+i)->connectionVertices[j]] = TRUE;
        }
    }


    findNextOrbitToConnect(blocks, blockCount, ddgraph, vertexToBlock, vertexToConnector, freeConnectors);

    //free the memory allocated at the beginning of this method
    freeDDGraph(ddgraph);
    free(blocks);
    DEBUGTRACE_EXIT
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
    DEBUGTRACE_ENTER
    if(!isMaybeRealizableComponentList()){
        DEBUGTRACE_EXIT
        return;
    } else if(!passesSimpleForbiddenConnectionsTest()){
        DEBUGTRACE_EXIT
        return;
    } else {
        componentListsCount++;
#ifdef _DEBUGINTERMEDIATE
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
    DEBUGTRACE_EXIT
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

    for(i = 0; i < Q1TypeComponentsCount; i++){
        Q1TypeComponentsComponentCount[i] = (int*)malloc(sizeof(int)*maximumQ1TypeComponents);
        for(j = 0; j < maximumQ1TypeComponents; j++){
            Q1TypeComponentsComponentCount[i][j]=0;
        }
    }

    for(i = 0; i < Q2TypeComponentsCount; i++){
        Q2TypeComponentsComponentCount[i] = (int*)malloc(sizeof(int)*maximumQ2TypeComponents);
        for(j = 0; j < maximumQ2TypeComponents; j++){
            Q2TypeComponentsComponentCount[i][j]=0;
        }
    }

    for(i = 0; i < Q3TypeComponentsCount; i++){
        Q3TypeComponentsComponentCount[i] = (int*)malloc(sizeof(int)*maximumQ3TypeComponents);
        for(j = 0; j < maximumQ3TypeComponents; j++){
            Q3TypeComponentsComponentCount[i][j]=0;
        }
    }

    Q4ComponentCount = 0;
}

void freeComponents(){
    int i;

    for(i = 0; i < Q1TypeComponentsCount; i++){
        free(Q1TypeComponentsComponentCount[i]);
    }

    for(i = 0; i < Q2TypeComponentsCount; i++){
        free(Q2TypeComponentsComponentCount[i]);
    }

    for(i = 0; i < Q3TypeComponentsCount; i++){
        free(Q3TypeComponentsComponentCount[i]);
    }
}

void initComponentsStatic(){
    
    Q1TypeComponentsConnectors[0] = 4; //Hub
    Q1TypeComponentsConnectors[1] = 3; //Locked Hub
    Q1TypeComponentsConnectors[2] = 2; //Double Locked Hub
    Q1TypeComponentsConnectors[3] = 2; //Diagonal Chain
    Q1TypeComponentsConnectors[4] = 2; //Double-roof High Building
    Q1TypeComponentsConnectors[5] = 2; //Open-roof High Building
    Q1TypeComponentsConnectors[6] = 2; //Double-roof Long Building
    Q1TypeComponentsConnectors[7] = 2; //Open-roof Long Building
    Q1TypeComponentsConnectors[8] = 1; //Locked Diagonal Chain
    Q1TypeComponentsConnectors[9] = 1; //Locked Double-roof High Building
    Q1TypeComponentsConnectors[10] = 1; //Locked Open-roof High Building
    Q1TypeComponentsConnectors[11] = 1; //Locked Double-roof Long Building

    Q1TypeComponentsSmallestCase[0] = 1;
    Q1TypeComponentsSmallestCase[1] = 1;
    Q1TypeComponentsSmallestCase[2] = 1;
    Q1TypeComponentsSmallestCase[3] = 1;
    Q1TypeComponentsSmallestCase[4] = 1;
    Q1TypeComponentsSmallestCase[5] = 1;
    Q1TypeComponentsSmallestCase[6] = 2;
    Q1TypeComponentsSmallestCase[7] = 2;
    Q1TypeComponentsSmallestCase[8] = 1;
    Q1TypeComponentsSmallestCase[9] = 1;
    Q1TypeComponentsSmallestCase[10] = 1;
    Q1TypeComponentsSmallestCase[11] = 2;

    Q1TypeComponentsComponentCount = (int**)malloc(sizeof(int*)*Q1TypeComponentsCount);

    Q2TypeComponentsConnectors[0] = 2; //Pearl Chain
    Q2TypeComponentsConnectors[1] = 1; //Locked Pearl Chain

    Q2TypeComponentsSmallestCase[0] = 1;
    Q2TypeComponentsSmallestCase[1] = 1;

    Q2TypeComponentsComponentCount = (int**)malloc(sizeof(int*)*Q2TypeComponentsCount);

    Q3TypeComponentsConnectors[0] = 2; //Barb Wire
    Q3TypeComponentsConnectors[1] = 1; //Locked Barb Wire

    Q3TypeComponentsSmallestCase[0] = 1;
    Q3TypeComponentsSmallestCase[1] = 1;

    Q3TypeComponentsComponentCount = (int**)malloc(sizeof(int*)*Q3TypeComponentsCount);


    //set up arrays with function pointers

    constructBlock[0] = &constructHub;
    constructBlock[1] = &constructLockedHub;
    constructBlock[2] = &constructDoubleLockedHub;
    constructBlock[3] = &constructDiagonalChain;
    constructBlock[4] = &constructDoubleroofHighBuilding;
    constructBlock[5] = &constructOpenroofHighBuilding;
    constructBlock[6] = &constructDoubleroofLongBuilding;
    constructBlock[7] = &constructOpenroofLongBuilding;
    constructBlock[8] = &constructLockedDiagonalChain;
    constructBlock[9] = &constructLockedDoubleroofHighBuilding;
    constructBlock[10] = &constructLockedOpenroofHighBuilding;
    constructBlock[11] = &constructLockedDoubleroofLongBuilding;
    constructBlock[12] = &constructPearlChain;
    constructBlock[13] = &constructLockedPearlChain;
    constructBlock[14] = &constructBarbWire;
    constructBlock[15] = &constructLockedBarbWire;
    constructBlock[16] = &constructQ4;

    storeBlockAutomorphismGenerators[0] = &storeHubAutomorphismGenerators;
    storeBlockAutomorphismGenerators[1] = &storeLockedHubAutomorphismGenerators;
    storeBlockAutomorphismGenerators[2] = &storeDoubleLockedHubAutomorphismGenerators;
    storeBlockAutomorphismGenerators[3] = &storeDiagonalChainAutomorphismGenerators;
    storeBlockAutomorphismGenerators[4] = &storeDoubleroofHighBuildingAutomorphismGenerators;
    storeBlockAutomorphismGenerators[5] = &storeOpenroofHighBuildingAutomorphismGenerators;
    storeBlockAutomorphismGenerators[6] = &storeDoubleroofLongBuildingAutomorphismGenerators;
    storeBlockAutomorphismGenerators[7] = &storeOpenroofLongBuildingAutomorphismGenerators;
    storeBlockAutomorphismGenerators[8] = &storeLockedDiagonalChainAutomorphismGenerators;
    storeBlockAutomorphismGenerators[9] = NULL;
    storeBlockAutomorphismGenerators[10] = &storeLockedOpenroofHighBuildingAutomorphismGenerators;
    storeBlockAutomorphismGenerators[11] = NULL;
    storeBlockAutomorphismGenerators[12] = &storePearlChainAutomorphismGenerators;
    storeBlockAutomorphismGenerators[13] = NULL;
    storeBlockAutomorphismGenerators[14] = &storeBarbWireAutomorphismGenerators;
    storeBlockAutomorphismGenerators[15] = NULL;
    storeBlockAutomorphismGenerators[16] = NULL;

    storeBlocksMapping[0] = &storeHubsMapping;
    storeBlocksMapping[1] = &storeLockedHubsMapping;
    storeBlocksMapping[2] = &storeDoubleLockedHubsMapping;
    storeBlocksMapping[3] = &storeDiagonalChainsMapping;
    storeBlocksMapping[4] = &storeDoubleroofHighBuildingsMapping;
    storeBlocksMapping[5] = &storeOpenroofHighBuildingsMapping;
    storeBlocksMapping[6] = &storeDoubleroofLongBuildingsMapping;
    storeBlocksMapping[7] = &storeOpenroofLongBuildingsMapping;
    storeBlocksMapping[8] = &storeLockedDiagonalChainsMapping;
    storeBlocksMapping[9] = &storeLockedDoubleroofHighBuildingsMapping;
    storeBlocksMapping[10] = &storeLockedOpenroofHighBuildingsMapping;
    storeBlocksMapping[11] = &storeLockedDoubleroofLongBuildingsMapping;
    storeBlocksMapping[12] = &storePearlChainsMapping;
    storeBlocksMapping[13] = &storeLockedPearlChainsMapping;
    storeBlocksMapping[14] = &storeBarbWiresMapping;
    storeBlocksMapping[15] = &storeLockedBarbWiresMapping;
    storeBlocksMapping[16] = &storeQ4sMapping;

#ifdef _DEBUGMETHODS
    blockName[0] = "H(%d)";
    blockName[1] = "LH(%d)";
    blockName[2] = "DLH(%d)";
    blockName[3] = "DC(%d)";
    blockName[4] = "DHB(%d)";
    blockName[5] = "OHB(%d)";
    blockName[6] = "DLB(%d)";
    blockName[7] = "OLB(%d)";
    blockName[8] = "LDC(%d)";
    blockName[9] = "LDHB(%d)";
    blockName[10] = "LOHB(%d)";
    blockName[11] = "LDLB(%d)";
    blockName[12] = "PC(%d)";
    blockName[13] = "LPC(%d)";
    blockName[14] = "BW(%d)";
    blockName[15] = "LBW(%d)";
    blockName[16] = "Q4";

#endif
}

void initStatistics(){
    componentListsCount = 0;
    graphsCount = 0;
}

void initNautyOptions(int order) {
    nautyOptions.defaultptn = FALSE;
    nautyOptions.getcanon = TRUE;
    nautyOptions.userautomproc = storeGenerator;
    #ifdef _DEBUGINTERMEDIATE
    nautyOptions.writeautoms = TRUE;
    nautyOptions.writemarkers = TRUE;
    nautyOptions.outfile = stderr;
    #endif

    int maxVertices = order + order/2 + 1;
    canonGraph.nv = 0;
    canonGraph.nde = 0;
    canonGraph.d = (int *)malloc(maxVertices*sizeof(int));
    canonGraph.v = (int *)malloc(maxVertices*sizeof(int));
    canonGraph.e = (int *)malloc((3*order + 2*(order/2) + 2)*sizeof(int));
    canonGraph.dlen = maxVertices;
    canonGraph.vlen = maxVertices;
    canonGraph.elen = 3*order + 2*(order/2) + 2;
}

//====================== START =======================

void startGeneration(int targetSize){

    initComponentsStatic();
    initComponents(targetSize);
    initStatistics();
    initNautyOptions(targetSize);

    q1Components(0, Q1TypeComponentsSmallestCase[0], targetSize, 0);

    fprintf(stderr, "Found %d component lists.\n", componentListsCount);
    fprintf(stderr, "Found %d Delaney-Dress graphs.\n", graphsCount);

}

void startFromListFile(char *filename){
    //read a list of components from a file
    FILE *f = fopen(filename, "r");
    if(f==NULL)
        ERRORMSG("Could not open file for input.")

    initComponentsStatic();
    initStatistics();

    int vertexCount = 0;
    while(fscanf(f, "%d", &vertexCount)!=-1){
        if(vertexCount<=0){
            ERRORMSG("Error while parsing file: illegal graph order.")
        }
        int realVertexCount = 0;

        initComponents(vertexCount);
        initNautyOptions(vertexCount);

        int type = 0;
        int family = 0;
        int parameter = 0;
        int count = 0;
        while(fscanf(f, "%d", &type)!=-1){
            if(type==0){
                break;
                //end of graph
            } else if(type==1){
                if(!fscanf(f, "%d", &family)){
                    ERRORMSG("Error while parsing file.")
                }
                if(!fscanf(f, "%d", &parameter)){
                    ERRORMSG("Error while parsing file.")
                }
                if(!fscanf(f, "%d", &count)){
                    ERRORMSG("Error while parsing file.")
                }
                if(family >= Q1TypeComponentsCount){
                    ERRORMSG("Error while parsing file: illegal family for type 1.")
                }
                if(parameter < Q1TypeComponentsSmallestCase[family] ||
                        parameter >maximumQ1TypeComponents){
                    ERRORMSG("Error while parsing file: illegal parameter.")
                }
                if(Q1TypeComponentsComponentCount[family][parameter-1]!=0){
                    ERRORMSG("Error while parsing file: this block type was already set.")
                }
                if(count <= 0){
                    ERRORMSG("Error while parsing file: illegal count.")
                }
                Q1TypeComponentsComponentCount[family][parameter-1] = count;
                realVertexCount += 4*count;
            } else if(type==2){
                if(!fscanf(f, "%d", &family)){
                    ERRORMSG("Error while parsing file.")
                }
                if(!fscanf(f, "%d", &parameter)){
                    ERRORMSG("Error while parsing file.")
                }
                if(!fscanf(f, "%d", &count)){
                    ERRORMSG("Error while parsing file.")
                }
                if(family >= Q2TypeComponentsCount){
                    ERRORMSG("Error while parsing file: illegal family for type 2.")
                }
                if(parameter < Q2TypeComponentsSmallestCase[family] ||
                        parameter >maximumQ2TypeComponents){
                    ERRORMSG("Error while parsing file: illegal parameter.")
                }
                if(Q2TypeComponentsComponentCount[family][parameter-1]!=0){
                    ERRORMSG("Error while parsing file: this block type was already set.")
                }
                if(count <= 0){
                    ERRORMSG("Error while parsing file: illegal count.")
                }
                Q2TypeComponentsComponentCount[family][parameter-1] = count;
                realVertexCount += 2*count;
            } else if(type==3){
                if(!fscanf(f, "%d", &family)){
                    ERRORMSG("Error while parsing file.")
                }
                if(!fscanf(f, "%d", &parameter)){
                    ERRORMSG("Error while parsing file.")
                }
                if(!fscanf(f, "%d", &count)){
                    ERRORMSG("Error while parsing file.")
                }
                if(family >= Q3TypeComponentsCount){
                    ERRORMSG("Error while parsing file: illegal family for type 3.")
                }
                if(parameter < Q3TypeComponentsSmallestCase[family] ||
                        parameter >maximumQ3TypeComponents){
                    ERRORMSG("Error while parsing file: illegal parameter.")
                }
                if(Q3TypeComponentsComponentCount[family][parameter-1]!=0){
                    ERRORMSG("Error while parsing file: this block type was already set.")
                }
                if(count <= 0){
                    ERRORMSG("Error while parsing file: illegal count.")
                }
                Q3TypeComponentsComponentCount[family][parameter-1] = count;
                realVertexCount += 2*count;
            } else if(type==4){
                if(fscanf(f, "%d", &count)){
                    if(Q4ComponentCount!=0){
                        ERRORMSG("Error while parsing file: this block type was already set.")
                    }
                    if(count <= 0){
                        ERRORMSG("Error while parsing file: illegal count.")
                    }
                    Q4ComponentCount = count;
                    realVertexCount += count;
                } else {
                    ERRORMSG("Error while parsing file.")
                }
            } else {
                ERRORMSG("Error while parsing file: illegal type.")
            }
        }
        
        if(realVertexCount!=vertexCount){
            ERRORMSG("Error while parsing file: incorrect vertex count.")
        }

        handleComponentList(vertexCount);
        freeComponents();
    }

    fprintf(stderr, "Read %d component lists.\n", componentListsCount);
    fprintf(stderr, "Found %d Delaney-Dress graphs.\n", graphsCount);
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
    char *listFilename = NULL;

    while ((c = getopt(argc, argv, "hl:")) != -1) {
        switch (c) {
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            case 'l': //(defaults to stdin)
                listFilename = optarg;
                break;
            default:
                fprintf(stderr, "Illegal option %c.\n", c);
                usage(name);
                return EXIT_FAILURE;
        }
    }

    // check the non-option arguments
    if (argc - optind != 1 && listFilename==NULL) {
        usage(name);
        return EXIT_FAILURE;
    }

    /*=========== initialization ===========*/
    struct tms TMS;
    unsigned int oldtime = 0;

    if(listFilename!=NULL){
        startFromListFile(listFilename);
    } else {

        //parse the order
        int vertexCount = strtol(argv[optind], NULL, 10);
        DEBUGDUMP(vertexCount, "%d")

        startGeneration(vertexCount);
    }



    times(&TMS);
    unsigned int savetime = oldtime + (unsigned int) TMS.tms_utime;
    fprintf(stderr, "CPU time: %.1f seconds.\n", (double) savetime / time_factor);

    return EXIT_SUCCESS;
}
