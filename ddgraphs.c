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

void printCanonicalLabelling(DDGRAPH *ddgraph){
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

    void (*userautomproc) (int,permutation*,int*,int,int,int) = nautyOptions.userautomproc;
    nautyOptions.userautomproc = NULL;
    nauty((graph*)(ddgraph->underlyingGraph), nautyLabelling, nautyPtn, NULL, currentOrbits,
            &nautyOptions, &nautyStats, nautyWorkspace, 50 * MAXM, MAXM,
            ddgraph->underlyingGraph->nv, (graph*)&canonGraph);
    nautyOptions.userautomproc = userautomproc;

    fprintf(stderr, "Canonical order of vertices: [%d", nautyLabelling[0]);
    for(i=1; i<ddgraph->underlyingGraph->nv; i++){
        fprintf(stderr, ", %d", nautyLabelling[i]);
    }
    fprintf(stderr, "]\n");
}

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

void printConnectionsMade(){
    int i;
    fprintf(stderr, "Connections:\n");
    for(i=0; i<connectionsMade; i++){
        fprintf(stderr, "  %2d) %2d - %2d\n", i, connections[i][0], connections[i][1]);
    }
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

#endif

void printBlockName(FILE *file, int count, int blockNumber, int parameter, boolean first){
    if(!first){
        fprintf(file, ", ");
    }
    fprintf(file, "%dx", count);
    fprintf(file, blockName[blockNumber], parameter);
}

void printHumanReadableComponentList(FILE *file){
    int i, j;
    boolean first = TRUE;
    for(i = 0; i < Q1TypeComponentsCount; i++){
        for(j = 0; j < maximumQ1TypeComponents; j++){
            if(Q1TypeComponentsComponentCount[i][j]){
                printBlockName(file, Q1TypeComponentsComponentCount[i][j], i, j+1, first);
                first = FALSE;
            }
        }
    }
    for(i = 0; i < Q2TypeComponentsCount; i++){
        for(j = 0; j < maximumQ2TypeComponents; j++){
            if(Q2TypeComponentsComponentCount[i][j]){
                printBlockName(file, Q2TypeComponentsComponentCount[i][j], Q1TypeComponentsCount + i, j+1, first);
                first = FALSE;
            }
        }
    }
    for(i = 0; i < Q3TypeComponentsCount; i++){
        for(j = 0; j < maximumQ3TypeComponents; j++){
            if(Q3TypeComponentsComponentCount[i][j]){
                printBlockName(file, Q3TypeComponentsComponentCount[i][j], Q1TypeComponentsCount + Q2TypeComponentsCount + i, j+1, first);
                first = FALSE;
            }
        }
    }
    if(Q4ComponentCount){
        printBlockName(file, Q4ComponentCount, Q1TypeComponentsCount + Q2TypeComponentsCount + Q3TypeComponentsCount, 0, first);
    }
    fprintf(file, "\n");
}

void writeListCodeMultipleBlockList(FILE *f, int vertexCount){
    int i, j;

    fprintf(f, "%d  ", vertexCount);
    for(i = 0; i < Q1TypeComponentsCount; i++){
        for(j = 0; j < maximumQ1TypeComponents; j++){
            if(Q1TypeComponentsComponentCount[i][j]){
                fprintf(f, "1 %d %d %d  ", i, j+1, Q1TypeComponentsComponentCount[i][j]);
            }
        }
    }
    for(i = 0; i < Q2TypeComponentsCount; i++){
        for(j = 0; j < maximumQ2TypeComponents; j++){
            if(Q2TypeComponentsComponentCount[i][j]){
                fprintf(f, "2 %d %d %d  ", i, j+1, Q2TypeComponentsComponentCount[i][j]);
            }
        }
    }
    for(i = 0; i < Q3TypeComponentsCount; i++){
        for(j = 0; j < maximumQ3TypeComponents; j++){
            if(Q3TypeComponentsComponentCount[i][j]){
                fprintf(f, "3 %d %d %d  ", i, j+1, Q3TypeComponentsComponentCount[i][j]);
            }
        }
    }
    if(Q4ComponentCount){
        fprintf(f, "4 %d  ", Q4ComponentCount);
    }
    fprintf(f, "0\n"); //end of list
}

void writeListCodeSingleBlockList(FILE *f, int vertexCount, BBLOCK *block){
    fprintf(f, "%d  ", vertexCount);
    fprintf(f, "%d %d %d  ", block->type, block->component, block->parameter);
    fprintf(f, "0\n"); //end of list
}

char write_2byte_number(FILE *f, unsigned short n, short writeEndian) {
    if (writeEndian == BIG_ENDIAN) {
        fprintf(f, "%c%c", n / 256, n % 256);
    } else {
        fprintf(f, "%c%c", n % 256, n / 256);
    }
    return (ferror(f) ? 2 : 1);
}

char writePregraphCode(FILE *f, DDGRAPH *ddgraph, boolean firstInFile) {
    unsigned short i, j;
    unsigned short semiEdge = ddgraph->order + 1;
    if (firstInFile) { //if first graph
        fprintf(f, ">>pregraph_code %s<<", (endian == LITTLE_ENDIAN ? "le" : "be"));
    }
    if (ddgraph->order + 1 <= UCHAR_MAX) {
        fprintf(f, "%c", (unsigned char) ddgraph->order);
    } else {
        fprintf(f, "%c", 0);
        /* big graph */
        if (write_2byte_number(f, (unsigned short) ddgraph->order, endian) == 2) {
            return (2);
        }
    }

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    for (i = 0; i < ddgraph->order; i++) {
        for (j = 0; j < degrees[i]; j++) {
            int neighbour = edges[positions[i]+j];
            if(neighbour >= ddgraph->order){
                //neighbour is a dummy vertex
                neighbour = edges[positions[neighbour]+0] + edges[positions[neighbour]+1]-i;
            }
            if(neighbour>i){
                //only include adjacency information for vertices with a larger index
                if (ddgraph->order + 1 <= UCHAR_MAX) {
                    fprintf(f, "%c", (unsigned char) (neighbour + 1));
                } else {
                    if (write_2byte_number(f, neighbour + 1, endian) == 2) {
                        return (2);
                    }
                }
            }
        }
        //add semi-edges
        for(j = 0; j < ddgraph->semiEdges[i]; j++){
            if (ddgraph->order + 1 <= UCHAR_MAX) {
                fprintf(f, "%c", (unsigned char)semiEdge);
            } else {
                if (write_2byte_number(f, semiEdge, endian) == 2) {
                    return (2);
                }
            }
        }
        //closing 0
        if (ddgraph->order + 1 <= UCHAR_MAX) {
            fprintf(f, "%c", 0);
        } else {
            if (write_2byte_number(f, 0, endian) == 2) {
                return (2);
            }
        }
    }
    return (ferror(f) ? 2 : 1);
}


char writePregraphColorCode2Factor(FILE *f, DDGRAPH *ddgraph, boolean firstInFile) {
    unsigned short i, j;
    unsigned short semiEdge = ddgraph->order + 1;
    if (firstInFile) { //if first graph
        fprintf(f, ">>pregraphcolor_code %s<<", (endian == LITTLE_ENDIAN ? "le" : "be"));
    }
    if (ddgraph->order + 1 <= UCHAR_MAX) {
        fprintf(f, "%c", (unsigned char) ddgraph->order);
    } else {
        fprintf(f, "%c", 0);
        /* big graph */
        if (write_2byte_number(f, (unsigned short) ddgraph->order, endian) == 2) {
            return (2);
        }
    }

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;

    int colours[3];
    int adjacencyListSize;

    for (i = 0; i < ddgraph->order; i++) {
        adjacencyListSize = 0;
        for (j = 0; j < 3; j++) {
            int neighbour = edges[3*i+j]; //don't use the current positions, but the initial ones!!
            if(neighbour==SEMIEDGE){
                if (ddgraph->order + 1 <= UCHAR_MAX) {
                    fprintf(f, "%c", (unsigned char)semiEdge);
                } else {
                    if (write_2byte_number(f, semiEdge, endian) == 2) {
                        return (2);
                    }
                }
                if(ddgraph->oneFactor[i]==j){
                    colours[adjacencyListSize] = 2;
                } else {
                    colours[adjacencyListSize] = 1;
                }
                adjacencyListSize++;
            } else {
                if(neighbour >= ddgraph->order){
                    //neighbour is a dummy vertex
                    neighbour = edges[positions[neighbour]+0] + edges[positions[neighbour]+1]-i;
                }
                if(neighbour>i){
                    //only include adjacency information for vertices with a larger index
                    if (ddgraph->order + 1 <= UCHAR_MAX) {
                        fprintf(f, "%c", (unsigned char) (neighbour + 1));
                    } else {
                        if (write_2byte_number(f, neighbour + 1, endian) == 2) {
                            return (2);
                        }
                    }
                    if(ddgraph->oneFactor[i]==j){
                        colours[adjacencyListSize] = 2;
                    } else {
                        colours[adjacencyListSize] = 1;
                    }
                    adjacencyListSize++;
                }
            }
        }
        //closing 0
        if (ddgraph->order + 1 <= UCHAR_MAX) {
            fprintf(f, "%c", 0);
        } else {
            if (write_2byte_number(f, 0, endian) == 2) {
                return (2);
            }
        }
        //colour list
        for(j=0; j<adjacencyListSize; j++){
            if (ddgraph->order + 1 <= UCHAR_MAX) {
                fprintf(f, "%c", colours[j]);
            } else {
                if (write_2byte_number(f, colours[j], endian) == 2) {
                    return (2);
                }
            }
        }
        //closing 0
        if (ddgraph->order + 1 <= UCHAR_MAX) {
            fprintf(f, "%c", 0);
        } else {
            if (write_2byte_number(f, 0, endian) == 2) {
                return (2);
            }
        }
    }
    return (ferror(f) ? 2 : 1);
}


char writePregraphColorCodeEdgeColouring(FILE *f, DDGRAPH *ddgraph, boolean firstInFile) {
    unsigned short i, j;
    unsigned short semiEdge = ddgraph->order + 1;
    if (firstInFile) { //if first graph
        fprintf(f, ">>pregraphcolor_code %s<<", (endian == LITTLE_ENDIAN ? "le" : "be"));
    }
    if (ddgraph->order + 1 <= UCHAR_MAX) {
        fprintf(f, "%c", (unsigned char) ddgraph->order);
    } else {
        fprintf(f, "%c", 0);
        /* big graph */
        if (write_2byte_number(f, (unsigned short) ddgraph->order, endian) == 2) {
            return (2);
        }
    }

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;

    int colours[3];
    int adjacencyListSize;

    for (i = 0; i < ddgraph->order; i++) {
        adjacencyListSize = 0;
        for (j = 0; j < 3; j++) {
            int neighbour = edges[3*i+j]; //don't use the current positions, but the initial ones!!
            if(neighbour==SEMIEDGE){
                if (ddgraph->order + 1 <= UCHAR_MAX) {
                    fprintf(f, "%c", (unsigned char)semiEdge);
                } else {
                    if (write_2byte_number(f, semiEdge, endian) == 2) {
                        return (2);
                    }
                }
                colours[adjacencyListSize] = ddgraph->colours[4*i+j] == INT_MAX ? 4 : ddgraph->colours[4*i+j] + 1;
                adjacencyListSize++;
            } else {
                if(neighbour >= ddgraph->order){
                    //neighbour is a dummy vertex
                    neighbour = edges[positions[neighbour]+0] + edges[positions[neighbour]+1]-i;
                }
                if(neighbour>i){
                    //only include adjacency information for vertices with a larger index
                    if (ddgraph->order + 1 <= UCHAR_MAX) {
                        fprintf(f, "%c", (unsigned char) (neighbour + 1));
                    } else {
                        if (write_2byte_number(f, neighbour + 1, endian) == 2) {
                            return (2);
                        }
                    }
                    colours[adjacencyListSize] = ddgraph->colours[4*i+j] == INT_MAX ? 4 : ddgraph->colours[4*i+j] + 1;
                    adjacencyListSize++;
                }
            }
        }
        //closing 0
        if (ddgraph->order + 1 <= UCHAR_MAX) {
            fprintf(f, "%c", 0);
        } else {
            if (write_2byte_number(f, 0, endian) == 2) {
                return (2);
            }
        }
        //colour list
        for(j=0; j<adjacencyListSize; j++){
            if (ddgraph->order + 1 <= UCHAR_MAX) {
                fprintf(f, "%c", colours[j]);
            } else {
                if (write_2byte_number(f, colours[j], endian) == 2) {
                    return (2);
                }
            }
        }
        //closing 0
        if (ddgraph->order + 1 <= UCHAR_MAX) {
            fprintf(f, "%c", 0);
        } else {
            if (write_2byte_number(f, 0, endian) == 2) {
                return (2);
            }
        }
    }
    return (ferror(f) ? 2 : 1);
}

void printRealDDGraph(DDGRAPH *graph){
    fprintf(stderr, "DDGRAPH %p\n", graph);
    fprintf(stderr, "================\n");
    fprintf(stderr, "order      %d\n", graph->order);

    int i,j;

    for(i=0; i<graph->order; i++){
        fprintf(stderr, "%d :", i);
        for(j=0; j < 3; j++){
            int s = SEMIEDGE;
            int neighbour = graph->underlyingGraph->e[3*i+j];
            if(neighbour == s){
                fprintf(stderr, "S ");
            } else if(neighbour >= graph->order){
                int dummy = neighbour;
                neighbour = graph->underlyingGraph->e[3*dummy+0] + graph->underlyingGraph->e[3*dummy+1] - i;
                fprintf(stderr, "%d ", neighbour);
            } else {
                fprintf(stderr, "%d ", neighbour);
            }
        }
        fprintf(stderr, "|");
        for(j=0; j<3; j++){
            fprintf(stderr, "%d ", graph->colours[4*i+j]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

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
    ddgraph->oneFactor = (int *)malloc(sizeof(int)*order);
    for(i = 0; i < order; i++){
        ddgraph->oneFactor[i] = INT_MAX;
    }
    ddgraph->colours = (int *)malloc(sizeof(int)*order*4);
    for(i = 0; i < 4*order; i++){
        ddgraph->colours[i] = INT_MAX;
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

    ddgraph->vertex2UncolouredFactor = (int *)malloc(sizeof(int)*order);
    ddgraph->vertex2FactorType = (int *)malloc(sizeof(int)*order);
    for(i = 0; i < order; i++){
        ddgraph->vertex2UncolouredFactor[i] = INT_MAX;
        ddgraph->vertex2FactorType[i] = INT_MAX;
    }

    ddgraph->uncolouredFactorCount = 0;

    ddgraph->uncolouredFactor2Vertex = (int *)malloc(sizeof(int)*(order/2));
    for(i = 0; i < (order/2); i++){
        ddgraph->uncolouredFactor2Vertex[i] = INT_MAX;
    }

    return ddgraph;
}

void cleanDDGraph(DDGRAPH * ddgraph){
    int i, order;

    order = ddgraph->order;
    ddgraph->dummyVertexCount = 0;
    for(i = 0; i < order; i++){
        ddgraph->semiEdges[i] = 0;
    }
    for(i = 0; i < order; i++){
        ddgraph->oneFactor[i] = INT_MAX;
    }
    for(i = 0; i < 4*order; i++){
        ddgraph->colours[i] = INT_MAX;
    }
    if(order==2){
        //order 2 is a bit special because we need to store the theta-graph
        ddgraph->underlyingGraph->nv = 0;
        ddgraph->underlyingGraph->nde = 0;
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

    //free old memory before reallocating new
    free(ddgraph->vertex2UncolouredFactor);
    free(ddgraph->vertex2FactorType);
    free(ddgraph->uncolouredFactor2Vertex);

    ddgraph->vertex2UncolouredFactor = (int *)malloc(sizeof(int)*order);
    ddgraph->vertex2FactorType = (int *)malloc(sizeof(int)*order);
    for(i = 0; i < order; i++){
        ddgraph->vertex2UncolouredFactor[i] = INT_MAX;
        ddgraph->vertex2FactorType[i] = INT_MAX;
    }

    ddgraph->uncolouredFactorCount = 0;

    ddgraph->uncolouredFactor2Vertex = (int *)malloc(sizeof(int)*(order/2));
    for(i = 0; i < (order/2); i++){
        ddgraph->uncolouredFactor2Vertex[i] = INT_MAX;
    }
}

void freeDDGraph(DDGRAPH *ddgraph){
    free(ddgraph->underlyingGraph->d);
    free(ddgraph->underlyingGraph->v);
    free(ddgraph->underlyingGraph->e);
    free(ddgraph->underlyingGraph);
    free(ddgraph->semiEdges);
    free(ddgraph->oneFactor);
    free(ddgraph->colours);
    free(ddgraph->uncolouredFactor2Vertex);
    free(ddgraph->vertex2UncolouredFactor);
    free(ddgraph->vertex2FactorType);
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
    } else if(type==5){
        block->connectorCount = 0;
        block->connections = NULL;
        block->targetConnector = NULL;
        block->connectionVertices = NULL;
    } else if(type==6){
        block->connectorCount = 0;
        block->connections = NULL;
        block->targetConnector = NULL;
        block->connectionVertices = NULL;
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

void freeBuildingBlock(BBLOCK *block){
    if(block->type<5){
        free(block->connections);
        free(block->targetConnector);
        free(block->connectionVertices);
    }
    free(block);
}

void freeBuildingBlocks(BBLOCK *block, int blockCount){
    int i;
    for(i=0; i<blockCount; i++){
        if((block+i)->type<5){
            free((block+i)->connections);
            free((block+i)->targetConnector);
            free((block+i)->connectionVertices);
        }
    }
    free(block);
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
    } else if(block->type==4){
        return Q1TypeComponentsCount + Q2TypeComponentsCount + Q3TypeComponentsCount;
    } else if(block->type==5){
        return Q1TypeComponentsCount + Q2TypeComponentsCount
                + Q3TypeComponentsCount + 1 + block->component;
    } else {
        return Q1TypeComponentsCount + Q2TypeComponentsCount 
                + Q3TypeComponentsCount + 1 + NoConnectorsFixedColouringComponentsCount
                + block->component;
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

    //the connections correspond to the one-factor
    ddgraph->oneFactor[*currentVertex] = 0;
    ddgraph->oneFactor[(*currentVertex)+1] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+2] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+3] = 0;
    //colours
    ddgraph->colours[4*(*currentVertex)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+1)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+2)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+3)+0] = 1;

    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)+2;
    degrees[*currentVertex] = 2;
    positions[*currentVertex]++; //connection points are at position 0

    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //connection points are at position 0

    //store q1-factor for colouring
    ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex);
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+0] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+1] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+2] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+3] = ddgraph->uncolouredFactorCount;
    ddgraph->uncolouredFactorCount++;
    ddgraph->vertex2FactorType[(*currentVertex)+0] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+1] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+2] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+3] = 1;

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;

        edges[positions[*currentVertex]+0] = (*currentVertex)-2;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-1;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;
        
        //set the one factor
        ddgraph->oneFactor[*currentVertex] = 2;
        ddgraph->oneFactor[(*currentVertex)+1] = 2;
        ddgraph->oneFactor[(*currentVertex)+2] = 0;
        ddgraph->oneFactor[(*currentVertex)+3] = 0;
        //colours
        ddgraph->colours[4*(*currentVertex)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+1)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+2)+0] = 1;
        ddgraph->colours[4*((*currentVertex)+3)+0] = 1;

        //store q1-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex)+2;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+0] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+1] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+2] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+3] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[(*currentVertex)+2+0] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+1] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+2] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+3] = 1;

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

    //the connections correspond to the one-factor
    ddgraph->oneFactor[*currentVertex] = 2;
    ddgraph->oneFactor[(*currentVertex)+1] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+2] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+3] = 0;
    //colours
    ddgraph->colours[4*(*currentVertex)+2] = 1;
    ddgraph->colours[4*((*currentVertex)+1)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+2)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+3)+0] = 1;

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

    //store q1-factor for colouring
    ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex);
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+0] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+1] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+2] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+3] = ddgraph->uncolouredFactorCount;
    ddgraph->uncolouredFactorCount++;
    ddgraph->vertex2FactorType[(*currentVertex)+0] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+1] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+2] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+3] = 1;

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;

        edges[positions[*currentVertex]+0] = (*currentVertex)-2;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-1;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;

        //set the one factor
        ddgraph->oneFactor[*currentVertex] = 2;
        ddgraph->oneFactor[(*currentVertex)+1] = 2;
        ddgraph->oneFactor[(*currentVertex)+2] = 0;
        ddgraph->oneFactor[(*currentVertex)+3] = 0;
        //colours
        ddgraph->colours[4*(*currentVertex)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+1)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+2)+0] = 1;
        ddgraph->colours[4*((*currentVertex)+3)+0] = 1;

        //store q1-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex)+2;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+0] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+1] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+2] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+3] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[(*currentVertex)+2+0] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+1] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+2] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+3] = 1;

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

    //the connections correspond to the one-factor
    ddgraph->oneFactor[*currentVertex] = 2;
    ddgraph->oneFactor[(*currentVertex)+1] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+2] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+3] = 2;
    //colours
    ddgraph->colours[4*(*currentVertex)+2] = 1;
    ddgraph->colours[4*((*currentVertex)+1)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+2)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+3)+2] = 1;

    edges[positions[*currentVertex]+0] = (*currentVertex)+1;
    edges[positions[*currentVertex]+1] = (*currentVertex)+2;
    edges[positions[*currentVertex]+2] = SEMIEDGE;
    degrees[*currentVertex] = 2;
    ddgraph->semiEdges[*currentVertex] = 1;

    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //connection points are at position 0

    //store q1-factor for colouring
    ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex);
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+0] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+1] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+2] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+3] = ddgraph->uncolouredFactorCount;
    ddgraph->uncolouredFactorCount++;
    ddgraph->vertex2FactorType[(*currentVertex)+0] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+1] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+2] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+3] = 1;

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;

        edges[positions[*currentVertex]+0] = (*currentVertex)-2;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-1;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;

        //set the one factor
        ddgraph->oneFactor[*currentVertex] = 2;
        ddgraph->oneFactor[(*currentVertex)+1] = 2;
        ddgraph->oneFactor[(*currentVertex)+2] = 0;
        ddgraph->oneFactor[(*currentVertex)+3] = 0;
        //colours
        ddgraph->colours[4*(*currentVertex)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+1)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+2)+0] = 1;
        ddgraph->colours[4*((*currentVertex)+3)+0] = 1;

        //store q1-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex)+2;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+0] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+1] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+2] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+3] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[(*currentVertex)+2+0] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+1] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+2] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+3] = 1;

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

    //the connections correspond to the one-factor
    ddgraph->oneFactor[*currentVertex] = 0;
    ddgraph->oneFactor[(*currentVertex)+1] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+2] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+3] = 0;
    //colours
    ddgraph->colours[4*(*currentVertex)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+1)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+2)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+3)+0] = 1;

    edges[positions[*currentVertex]+0] = (*currentVertex)+(block->parameter-1)*4+3;
    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)+2;

    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //connection points are at position 0

    //store q1-factor for colouring
    ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex);
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+0] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+1] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+2] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+3] = ddgraph->uncolouredFactorCount;
    ddgraph->uncolouredFactorCount++;
    ddgraph->vertex2FactorType[(*currentVertex)+0] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+1] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+2] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+3] = 1;

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;

        edges[positions[*currentVertex]+0] = (*currentVertex)-2;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-1;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;

        //set the one factor
        ddgraph->oneFactor[*currentVertex] = 2;
        ddgraph->oneFactor[(*currentVertex)+1] = 2;
        ddgraph->oneFactor[(*currentVertex)+2] = 0;
        ddgraph->oneFactor[(*currentVertex)+3] = 0;
        //colours
        ddgraph->colours[4*(*currentVertex)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+1)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+2)+0] = 1;
        ddgraph->colours[4*((*currentVertex)+3)+0] = 1;

        //store q1-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex)+2;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+0] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+1] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+2] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+3] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[(*currentVertex)+2+0] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+1] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+2] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+3] = 1;

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

    //the connections correspond to the one-factor
    ddgraph->oneFactor[*currentVertex] = 0;
    ddgraph->oneFactor[(*currentVertex)+1] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+2] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+3] = 0;
    //colours
    ddgraph->colours[4*(*currentVertex)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+1)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+2)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+3)+0] = 1;

    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)+2;
    degrees[*currentVertex] = 2;
    positions[*currentVertex]++; //connection points are at position 0

    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //connection points are at position 0

    //store q1-factor for colouring
    ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex);
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+0] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+1] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+2] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+3] = ddgraph->uncolouredFactorCount;
    ddgraph->uncolouredFactorCount++;
    ddgraph->vertex2FactorType[(*currentVertex)+0] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+1] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+2] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+3] = 1;

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;

        edges[positions[*currentVertex]+0] = (*currentVertex)-2;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-1;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;

        //set the one factor
        ddgraph->oneFactor[*currentVertex] = 2;
        ddgraph->oneFactor[(*currentVertex)+1] = 2;
        ddgraph->oneFactor[(*currentVertex)+2] = 0;
        ddgraph->oneFactor[(*currentVertex)+3] = 0;
        //colours
        ddgraph->colours[4*(*currentVertex)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+1)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+2)+0] = 1;
        ddgraph->colours[4*((*currentVertex)+3)+0] = 1;

        //store q1-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex)+2;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+0] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+1] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+2] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+3] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[(*currentVertex)+2+0] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+1] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+2] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+3] = 1;

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

    //the connections correspond to the one-factor
    ddgraph->oneFactor[*currentVertex] = 0;
    ddgraph->oneFactor[(*currentVertex)+1] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+2] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+3] = 0;
    //colours
    ddgraph->colours[4*(*currentVertex)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+1)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+2)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+3)+0] = 1;

    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)+2;
    degrees[*currentVertex] = 2;
    positions[*currentVertex]++; //connection points are at position 0

    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //connection points are at position 0

    //store q1-factor for colouring
    ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex);
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+0] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+1] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+2] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+3] = ddgraph->uncolouredFactorCount;
    ddgraph->uncolouredFactorCount++;
    ddgraph->vertex2FactorType[(*currentVertex)+0] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+1] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+2] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+3] = 1;

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;

        edges[positions[*currentVertex]+0] = (*currentVertex)-2;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-1;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;

        //set the one factor
        ddgraph->oneFactor[*currentVertex] = 2;
        ddgraph->oneFactor[(*currentVertex)+1] = 2;
        ddgraph->oneFactor[(*currentVertex)+2] = 0;
        ddgraph->oneFactor[(*currentVertex)+3] = 0;
        //colours
        ddgraph->colours[4*(*currentVertex)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+1)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+2)+0] = 1;
        ddgraph->colours[4*((*currentVertex)+3)+0] = 1;

        //store q1-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex)+2;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+0] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+1] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+2] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+3] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[(*currentVertex)+2+0] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+1] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+2] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+3] = 1;

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

    //the connections correspond to the one-factor
    ddgraph->oneFactor[*currentVertex] = 0;
    ddgraph->oneFactor[(*currentVertex)+1] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+2] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+3] = 0;
    //colours
    ddgraph->colours[4*(*currentVertex)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+1)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+2)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+3)+0] = 1;

    edges[positions[*currentVertex]+0] = (*currentVertex)+(block->parameter-1)*4+2;
    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)+2;

    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //connection points are at position 0

    //store q1-factor for colouring
    ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex);
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+0] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+1] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+2] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+3] = ddgraph->uncolouredFactorCount;
    ddgraph->uncolouredFactorCount++;
    ddgraph->vertex2FactorType[(*currentVertex)+0] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+1] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+2] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+3] = 1;

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;

        edges[positions[*currentVertex]+0] = (*currentVertex)-2;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-1;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;

        //set the one factor
        ddgraph->oneFactor[*currentVertex] = 2;
        ddgraph->oneFactor[(*currentVertex)+1] = 2;
        ddgraph->oneFactor[(*currentVertex)+2] = 0;
        ddgraph->oneFactor[(*currentVertex)+3] = 0;
        //colours
        ddgraph->colours[4*(*currentVertex)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+1)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+2)+0] = 1;
        ddgraph->colours[4*((*currentVertex)+3)+0] = 1;

        //store q1-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex)+2;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+0] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+1] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+2] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+3] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[(*currentVertex)+2+0] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+1] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+2] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+3] = 1;

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

    if(block->parameter == 2){
        generator = getIdentity(ddgraph->underlyingGraph->nv);

        int firstVertex = block->connectionVertices[0]-1;

        generator[firstVertex] = firstVertex + 3;
        generator[firstVertex + 3] = firstVertex;
        generator[firstVertex + 5] = firstVertex + 6;
        generator[firstVertex + 6] = firstVertex + 5;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);
    }

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

    //the connections correspond to the one-factor
    ddgraph->oneFactor[*currentVertex] = 0;
    ddgraph->oneFactor[(*currentVertex)+1] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+2] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+3] = 0;
    //colours
    ddgraph->colours[4*(*currentVertex)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+1)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+2)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+3)+0] = 1;

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

    //store q1-factor for colouring
    ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex);
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+0] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+1] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+2] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+3] = ddgraph->uncolouredFactorCount;
    ddgraph->uncolouredFactorCount++;
    ddgraph->vertex2FactorType[(*currentVertex)+0] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+1] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+2] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+3] = 1;

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;


        edges[positions[*currentVertex]+0] = (*currentVertex)-2;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-1;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;

        //set the one factor
        ddgraph->oneFactor[*currentVertex] = 2;
        ddgraph->oneFactor[(*currentVertex)+1] = 2;
        ddgraph->oneFactor[(*currentVertex)+2] = 0;
        ddgraph->oneFactor[(*currentVertex)+3] = 0;
        //colours
        ddgraph->colours[4*(*currentVertex)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+1)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+2)+0] = 1;
        ddgraph->colours[4*((*currentVertex)+3)+0] = 1;

        //store q1-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex)+2;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+0] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+1] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+2] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+3] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[(*currentVertex)+2+0] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+1] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+2] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+3] = 1;

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

    //the connections correspond to the one-factor
    ddgraph->oneFactor[*currentVertex] = 0;
    ddgraph->oneFactor[(*currentVertex)+1] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+2] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+3] = 0;
    //colours
    ddgraph->colours[4*(*currentVertex)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+1)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+2)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+3)+0] = 1;

    edges[positions[*currentVertex]+0] = (*currentVertex)+(block->parameter-1)*4+3;
    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)+2;

    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //connection points are at position 0

    //store q1-factor for colouring
    ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex);
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+0] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+1] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+2] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+3] = ddgraph->uncolouredFactorCount;
    ddgraph->uncolouredFactorCount++;
    ddgraph->vertex2FactorType[(*currentVertex)+0] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+1] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+2] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+3] = 1;

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;

        edges[positions[*currentVertex]+0] = (*currentVertex)-2;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-1;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;

        //set the one factor
        ddgraph->oneFactor[*currentVertex] = 2;
        ddgraph->oneFactor[(*currentVertex)+1] = 2;
        ddgraph->oneFactor[(*currentVertex)+2] = 0;
        ddgraph->oneFactor[(*currentVertex)+3] = 0;
        //colours
        ddgraph->colours[4*(*currentVertex)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+1)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+2)+0] = 1;
        ddgraph->colours[4*((*currentVertex)+3)+0] = 1;

        //store q1-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex)+2;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+0] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+1] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+2] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+3] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[(*currentVertex)+2+0] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+1] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+2] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+3] = 1;

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

    //the connections correspond to the one-factor
    ddgraph->oneFactor[*currentVertex] = 0;
    ddgraph->oneFactor[(*currentVertex)+1] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+2] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+3] = 0;
    //colours
    ddgraph->colours[4*(*currentVertex)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+1)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+2)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+3)+0] = 1;

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

    //store q1-factor for colouring
    ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex);
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+0] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+1] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+2] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+3] = ddgraph->uncolouredFactorCount;
    ddgraph->uncolouredFactorCount++;
    ddgraph->vertex2FactorType[(*currentVertex)+0] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+1] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+2] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+3] = 1;

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;

        edges[positions[*currentVertex]+0] = (*currentVertex)-2;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-1;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;

        //set the one factor
        ddgraph->oneFactor[*currentVertex] = 2;
        ddgraph->oneFactor[(*currentVertex)+1] = 2;
        ddgraph->oneFactor[(*currentVertex)+2] = 0;
        ddgraph->oneFactor[(*currentVertex)+3] = 0;
        //colours
        ddgraph->colours[4*(*currentVertex)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+1)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+2)+0] = 1;
        ddgraph->colours[4*((*currentVertex)+3)+0] = 1;

        //store q1-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex)+2;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+0] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+1] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+2] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+3] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[(*currentVertex)+2+0] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+1] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+2] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+3] = 1;

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

    //the connections correspond to the one-factor
    ddgraph->oneFactor[*currentVertex] = 0;
    ddgraph->oneFactor[(*currentVertex)+1] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+2] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+3] = 0;
    //colours
    ddgraph->colours[4*(*currentVertex)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+1)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+2)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+3)+0] = 1;

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

    //store q1-factor for colouring
    ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex);
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+0] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+1] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+2] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+3] = ddgraph->uncolouredFactorCount;
    ddgraph->uncolouredFactorCount++;
    ddgraph->vertex2FactorType[(*currentVertex)+0] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+1] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+2] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+3] = 1;

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;
        
        edges[positions[*currentVertex]+0] = (*currentVertex)-2;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-1;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;

        //set the one factor
        ddgraph->oneFactor[*currentVertex] = 2;
        ddgraph->oneFactor[(*currentVertex)+1] = 2;
        ddgraph->oneFactor[(*currentVertex)+2] = 0;
        ddgraph->oneFactor[(*currentVertex)+3] = 0;
        //colours
        ddgraph->colours[4*(*currentVertex)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+1)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+2)+0] = 1;
        ddgraph->colours[4*((*currentVertex)+3)+0] = 1;

        //store q1-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex)+2;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+0] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+1] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+2] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+3] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[(*currentVertex)+2+0] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+1] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+2] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+3] = 1;

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

    //the connections correspond to the one-factor
    ddgraph->oneFactor[*currentVertex] = 0;
    ddgraph->oneFactor[(*currentVertex)+1] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+2] = 0;
    ddgraph->oneFactor[(*currentVertex)+(block->parameter-1)*4+3] = 0;
    //colours
    ddgraph->colours[4*(*currentVertex)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+1)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+2)+0] = 1;
    ddgraph->colours[4*((*currentVertex)+(block->parameter-1)*4+3)+0] = 1;

    edges[positions[*currentVertex]+0] = (*currentVertex)+(block->parameter-1)*4+2;
    edges[positions[*currentVertex]+1] = (*currentVertex)+1;
    edges[positions[*currentVertex]+2] = (*currentVertex)+2;

    edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
    edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;
    degrees[(*currentVertex)+1] = 2;
    positions[(*currentVertex)+1]++; //connection points are at position 0

    //store q1-factor for colouring
    ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex);
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+0] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+1] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+2] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[(*currentVertex)+3] = ddgraph->uncolouredFactorCount;
    ddgraph->uncolouredFactorCount++;
    ddgraph->vertex2FactorType[(*currentVertex)+0] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+1] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+2] = 1;
    ddgraph->vertex2FactorType[(*currentVertex)+3] = 1;

    (*currentVertex)+=2;

    for(i=0; i<block->parameter-1; i++){
        vertexToBlock[(*currentVertex)+0] = block->id;
        vertexToBlock[(*currentVertex)+1] = block->id;
        vertexToBlock[(*currentVertex)+2] = block->id;
        vertexToBlock[(*currentVertex)+3] = block->id;

        edges[positions[*currentVertex]+0] = (*currentVertex)-2;
        edges[positions[*currentVertex]+1] = (*currentVertex)+1;
        edges[positions[*currentVertex]+2] = (*currentVertex)+2;

        edges[positions[(*currentVertex)+1]+0] = (*currentVertex)-1;
        edges[positions[(*currentVertex)+1]+1] = (*currentVertex);
        edges[positions[(*currentVertex)+1]+2] = (*currentVertex)+3;

        edges[positions[(*currentVertex)+2]+0] = (*currentVertex);
        edges[positions[(*currentVertex)+2]+1] = (*currentVertex)+3;
        edges[positions[(*currentVertex)+2]+2] = (*currentVertex)+4;

        edges[positions[(*currentVertex)+3]+0] = (*currentVertex)+1;
        edges[positions[(*currentVertex)+3]+1] = (*currentVertex)+2;
        edges[positions[(*currentVertex)+3]+2] = (*currentVertex)+5;

        //set the one factor
        ddgraph->oneFactor[*currentVertex] = 2;
        ddgraph->oneFactor[(*currentVertex)+1] = 2;
        ddgraph->oneFactor[(*currentVertex)+2] = 0;
        ddgraph->oneFactor[(*currentVertex)+3] = 0;
        //colours
        ddgraph->colours[4*(*currentVertex)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+1)+2] = 1;
        ddgraph->colours[4*((*currentVertex)+2)+0] = 1;
        ddgraph->colours[4*((*currentVertex)+3)+0] = 1;

        //store q1-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex)+2;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+0] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+1] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+2] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[(*currentVertex)+2+3] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[(*currentVertex)+2+0] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+1] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+2] = 1;
        ddgraph->vertex2FactorType[(*currentVertex)+2+3] = 1;

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
    ddgraph->semiEdges[(*currentVertex)+1] = 1;

    (*currentVertex)+=2;
}

void storeLockedDoubleroofLongBuildingAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    if(block->parameter == 2){
        permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

        int firstVertex = block->connectionVertices[0]-1;

        generator[firstVertex] = firstVertex + 3;
        generator[firstVertex + 3] = firstVertex;
        generator[firstVertex + 5] = firstVertex + 6;
        generator[firstVertex + 6] = firstVertex + 5;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    }
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
        ddgraph->oneFactor[*currentVertex] = 0;
        //colours
        ddgraph->colours[4*(*currentVertex)+0] = 1;
        ddgraph->colours[4*(*currentVertex)+1] = 0;
        ddgraph->colours[4*(*currentVertex)+2] = 2;
        ddgraph->vertex2FactorType[(*currentVertex)] = 2;
        (*currentVertex)++;

        edges[positions[*currentVertex]+0] = (*currentVertex)+1;
        //for the last vertex this will be overwritten when making the connections
        edges[positions[*currentVertex]+1] = (*currentVertex)-1;
        edges[positions[*currentVertex]+2] = dummyVertex;
        edges[positions[dummyVertex]+1] = (*currentVertex);
        ddgraph->oneFactor[*currentVertex] = 0;
        //colours
        ddgraph->colours[4*(*currentVertex)+0] = 1;
        ddgraph->colours[4*(*currentVertex)+1] = 0;
        ddgraph->colours[4*(*currentVertex)+2] = 2;
        ddgraph->vertex2FactorType[(*currentVertex)] = 2;
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
        ddgraph->oneFactor[*currentVertex] = 0;
        //colours
        ddgraph->colours[4*(*currentVertex)+0] = 1;
        ddgraph->colours[4*(*currentVertex)+1] = 0;
        ddgraph->colours[4*(*currentVertex)+2] = 2;
        ddgraph->vertex2FactorType[(*currentVertex)] = 2;
        (*currentVertex)++;

        edges[positions[*currentVertex]+0] = (*currentVertex)+1;
        //for the last vertex this will be overwritten when making the connections
        edges[positions[*currentVertex]+1] = (*currentVertex)-1;
        edges[positions[*currentVertex]+2] = dummyVertex;
        edges[positions[dummyVertex]+1] = (*currentVertex);
        ddgraph->oneFactor[*currentVertex] = 0;
        //colours
        ddgraph->colours[4*(*currentVertex)+0] = 1;
        ddgraph->colours[4*(*currentVertex)+1] = 0;
        ddgraph->colours[4*(*currentVertex)+2] = 2;
        ddgraph->vertex2FactorType[(*currentVertex)] = 2;
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
    dummyBlock1 = ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[vertexBlock1]+1];
    vertexBlock2 = block2->connectionVertices[0];
    dummyBlock2 = ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[vertexBlock2]+1];
    //+1 instead of +2 because the positions have been adjusted

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
 * Constructs a barbed wire BW(n). This block has 2 connectors
 * arranged as follows:
 *
 *       |    |     |    |
 *    0--o----o-...-o----o--1
 *
 *
 */
void constructBarbedWire(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
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
        ddgraph->oneFactor[*currentVertex] = 0;
        //colours
        ddgraph->colours[4*(*currentVertex)+0] = 1;
        //store q3-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex);
        ddgraph->vertex2UncolouredFactor[(*currentVertex)] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2FactorType[(*currentVertex)] = 3;
        (*currentVertex)++;

        edges[positions[*currentVertex]+0] = (*currentVertex)+1;
        //for the last vertex this will be overwritten when making the connections
        edges[positions[*currentVertex]+1] = (*currentVertex)-1;
        edges[positions[*currentVertex]+2] = SEMIEDGE;
        ddgraph->semiEdges[(*currentVertex)] = 1;
        degrees[*currentVertex] = 2;
        ddgraph->oneFactor[*currentVertex] = 0;
        //colours
        ddgraph->colours[4*(*currentVertex)+0] = 1;

        //store q3-factor for colouring
        ddgraph->vertex2UncolouredFactor[(*currentVertex)] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[(*currentVertex)] = 3;

        (*currentVertex)++;
    }

    degrees[block->connectionVertices[0]] = 1;
    degrees[block->connectionVertices[1]] = 1;
    positions[block->connectionVertices[0]]++;
    positions[block->connectionVertices[1]]++;
}

void storeBarbedWireAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
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

void storeBarbedWiresMapping(BBLOCK *block1, BBLOCK *block2, DDGRAPH *ddgraph){
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
 * Constructs a locked barbed wire LBW(n). This block has 1 connector
 * arranged as follows:
 *
 *       |    |     |    |
 *    0--o----o-...-o----o-
 *
 * This block has a trivial symmetry group.
 */
void constructLockedBarbedWire(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
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
        ddgraph->oneFactor[*currentVertex] = 0;
        //colours
        ddgraph->colours[4*(*currentVertex)+0] = 1;
        //store q3-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = (*currentVertex);
        ddgraph->vertex2UncolouredFactor[(*currentVertex)] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2FactorType[(*currentVertex)] = 3;
        (*currentVertex)++;

        edges[positions[*currentVertex]+0] = (*currentVertex)+1;
        //for the last vertex this will be overwritten when making the connections
        edges[positions[*currentVertex]+1] = (*currentVertex)-1;
        edges[positions[*currentVertex]+2] = SEMIEDGE;
        degrees[*currentVertex] = 2;
        ddgraph->semiEdges[(*currentVertex)] = 1;
        ddgraph->oneFactor[*currentVertex] = 0;
        //colours
        ddgraph->colours[4*(*currentVertex)+0] = 1;

        //store q3-factor for colouring
        ddgraph->vertex2UncolouredFactor[(*currentVertex)] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[(*currentVertex)] = 3;

        (*currentVertex)++;
    }
    edges[positions[(*currentVertex)-1]+0] = SEMIEDGE;
    degrees[(*currentVertex)-1] = 1;
    ddgraph->semiEdges[(*currentVertex)-1] = 2;
    positions[(*currentVertex)-1]++;

    degrees[block->connectionVertices[0]] = 1;
    positions[block->connectionVertices[0]]++;
}

void storeLockedBarbedWiresMapping(BBLOCK *block1, BBLOCK *block2, DDGRAPH *ddgraph){
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
    ddgraph->oneFactor[*currentVertex] = 0;
    //colours
    ddgraph->colours[4*(*currentVertex)+0] = 1;
    ddgraph->colours[4*(*currentVertex)+1] = 0;
    ddgraph->colours[4*(*currentVertex)+2] = 2;

    //it is not necessary to adjust ddgraph->underlyingGraph->v because it will
    //not be used when the degree is 0.
    block->connectionVertices[0] = *currentVertex;
    vertexToConnector[*currentVertex] = 0;
    vertexToBlock[*currentVertex] = block->id;
    ddgraph->vertex2FactorType[*currentVertex] = 4;
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

/*
 * Constructs a tristar. This block has 0 connectors.
 *
 *       /
 *    --o
 *       \
 *
 */
void constructTristar(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    DEBUGASSERT((*currentVertex) == 0)
    ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[0]+0] = SEMIEDGE;
    ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[0]+1] = SEMIEDGE;
    ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[0]+2] = SEMIEDGE;
    ddgraph->semiEdges[0] = 3;
    ddgraph->oneFactor[0] = 0;
    //colours
    ddgraph->colours[0+0] = 1;
    ddgraph->colours[0+1] = 0;
    ddgraph->colours[0+2] = 2;
    ddgraph->underlyingGraph->d[0] = 0;
    ddgraph->vertex2FactorType[0] = 4;
}

/*
 * Constructs a double locked pearl chain DLPC(n). This block has 0 connectors.
 *        ____         ____
 *       /    \       /    \
 *    --o      o-...-o      o--
 *       \____/       \____/
 *
 */
void constructDoubleLockedPearlChain(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    DEBUGASSERT((*currentVertex) == 0)
    int parameter = block->parameter;
    int i;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    for(i=0; i<parameter; i++){
        int dummyVertex = ddgraph->order + ddgraph->dummyVertexCount;
        ddgraph->dummyVertexCount++;

        edges[positions[2*i]+0] = 2*i-1;
        //for the first vertex this will be overwritten when making the connections
        edges[positions[2*i]+1] = 2*i+1;
        edges[positions[2*i]+2] = dummyVertex;
        edges[positions[dummyVertex]+0] = 2*i;
        ddgraph->oneFactor[2*i] = 0;
        //colours
        ddgraph->colours[4*(2*i)+0] = 1;
        ddgraph->colours[4*(2*i)+1] = 0;
        ddgraph->colours[4*(2*i)+2] = 2;

        edges[positions[2*i+1]+0] = 2*i+2;
        //for the last vertex this will be overwritten when making the connections
        edges[positions[2*i+1]+1] = 2*i;
        edges[positions[2*i+1]+2] = dummyVertex;
        edges[positions[dummyVertex]+1] = 2*i+1;
        ddgraph->oneFactor[2*i+1] = 0;
        //colours
        ddgraph->colours[4*(2*i+1)+0] = 1;
        ddgraph->colours[4*(2*i+1)+1] = 0;
        ddgraph->colours[4*(2*i+1)+2] = 2;
        ddgraph->vertex2FactorType[2*i] = 2;
        ddgraph->vertex2FactorType[2*i+1] = 2;
    }

    edges[positions[0]+0] = SEMIEDGE;
    ddgraph->semiEdges[0] = 1;
    degrees[0] = 2;
    positions[0]++;

    edges[positions[2*parameter-1]+0] = SEMIEDGE;
    ddgraph->semiEdges[2*parameter-1] = 1;
    degrees[2*parameter-1] = 2;
    positions[2*parameter-1]++;
}

void storeDoubleLockedPearlChainAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    int parameter = block->parameter;
    int i, vertexLeft, vertexRight, dummyLeft, dummyRight;

    //mirror symmetry
    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    vertexLeft = 0;
    vertexRight = 2*parameter-1;
    dummyLeft = 2*parameter;
    dummyRight = 3*parameter-1;

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

/*
 * Constructs a pearl necklace DLPC(n). This block has 0 connectors.
 *     _______________________
 *    |   ____         ____   |
 *    |  /    \       /    \  |
 *    `-o      o-...-o      o-'
 *       \____/       \____/
 *
 */
void constructPearlNecklace(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    DEBUGASSERT((*currentVertex) == 0)
    int parameter = block->parameter;
    int i;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;

    for(i=0; i<parameter; i++){
        int dummyVertex = ddgraph->order + ddgraph->dummyVertexCount;
        ddgraph->dummyVertexCount++;

        edges[positions[2*i]+0] = 2*i-1;
        //for the first vertex this will be overwritten when making the connections
        edges[positions[2*i]+1] = 2*i+1;
        edges[positions[2*i]+2] = dummyVertex;
        edges[positions[dummyVertex]+0] = 2*i;
        ddgraph->oneFactor[2*i] = 0;
        //colours
        ddgraph->colours[4*(2*i)+0] = 1;
        ddgraph->colours[4*(2*i)+1] = 0;
        ddgraph->colours[4*(2*i)+2] = 2;

        edges[positions[2*i+1]+0] = 2*i+2;
        //for the last vertex this will be overwritten when making the connections
        edges[positions[2*i+1]+1] = 2*i;
        edges[positions[2*i+1]+2] = dummyVertex;
        edges[positions[dummyVertex]+1] = 2*i+1;
        ddgraph->oneFactor[2*i+1] = 0;
        //colours
        ddgraph->colours[4*(2*i+1)+0] = 1;
        ddgraph->colours[4*(2*i+1)+1] = 0;
        ddgraph->colours[4*(2*i+1)+2] = 2;

        ddgraph->vertex2FactorType[2*i] = 2;
        ddgraph->vertex2FactorType[2*i+1] = 2;
    }

    edges[positions[0]+0] = 2*parameter-1;
    edges[positions[2*parameter-1]+0] = 0;
}

void storePearlNecklaceAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    int parameter = block->parameter;
    int i, vertexLeft, vertexRight, dummyLeft, dummyRight;

    //mirror symmetry
    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    vertexLeft = 0;
    vertexRight = 2*parameter-1;
    dummyLeft = 2*parameter;
    dummyRight = 3*parameter-1;

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

    if(parameter>1){
        //rotation symmetry
        generator = getIdentity(ddgraph->underlyingGraph->nv);

        for(i=0; i<parameter-1; i++){
            generator[2*i] = 2*(i+1);
            generator[2*i+1] = 2*(i+1)+1;
            generator[2*parameter+i] = 2*parameter+i+1;
        }

        generator[2*(parameter-1)] = 0;
        generator[2*(parameter-1)+1] = 1;
        generator[3*parameter-1] = 2*parameter;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    }
}

/*
 * Constructs a double locked barbed wire DLBW(n). This block has 0 connectors.
 *
 *      |    |     |    |
 *    --o----o-...-o----o--
 *
 *
 */
void constructDoubleLockedBarbedWire(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    DEBUGASSERT((*currentVertex) == 0)
    int parameter = block->parameter;
    int i;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    for(i=0; i<parameter; i++){

        edges[positions[2*i]+0] = 2*i-1;
        //for the first vertex this will be overwritten when making the connections
        edges[positions[2*i]+1] = 2*i+1;
        edges[positions[2*i]+2] = SEMIEDGE;
        ddgraph->semiEdges[2*i] = 1;
        degrees[2*i] = 2;
        ddgraph->oneFactor[2*i] = 0;
        //colours
        ddgraph->colours[4*(2*i)+0] = 1;

        edges[positions[2*i+1]+0] = 2*i+2;
        //for the last vertex this will be overwritten when making the connections
        edges[positions[2*i+1]+1] = 2*i;
        edges[positions[2*i+1]+2] = SEMIEDGE;
        ddgraph->semiEdges[2*i+1] = 1;
        degrees[2*i+1] = 2;
        ddgraph->oneFactor[2*i+1] = 0;
        //colours
        ddgraph->colours[4*(2*i+1)+0] = 1;

        //store q3-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = 2*i;
        ddgraph->vertex2UncolouredFactor[2*i+0] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[2*i+1] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[2*i+0] = 3;
        ddgraph->vertex2FactorType[2*i+1] = 3;
    }
    edges[positions[0]+0] = SEMIEDGE;
    degrees[0] = 1;
    ddgraph->semiEdges[0] = 2;
    positions[0]++;
    edges[positions[2*parameter-1]+0] = SEMIEDGE;
    degrees[2*parameter-1] = 1;
    ddgraph->semiEdges[2*parameter-1] = 2;
    positions[2*parameter-1]++;
}

void storeDoubleLockedBarbedWireAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    int parameter = block->parameter;
    int i, vertexLeft, vertexRight;

    //mirror symmetry
    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    vertexLeft = 0;
    vertexRight = 2*parameter-1;

    for(i=0; i<parameter; i++){
        generator[vertexLeft] = vertexRight;
        generator[vertexRight] = vertexLeft;
        vertexLeft++;
        vertexRight--;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);
}

/*
 * Constructs a barbed wire necklace BWN(n). This block has 0 connectors.
 *     _____________________
 *    |                     |
 *    |  |    |     |    |  |
 *    `--o----o-...-o----o--'
 *
 *
 */
void constructBarbedWireNecklace(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    DEBUGASSERT((*currentVertex) == 0)
    int parameter = block->parameter;
    int i;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    for(i=0; i<parameter; i++){

        edges[positions[2*i]+0] = 2*i-1;
        //for the first vertex this will be overwritten when making the connections
        edges[positions[2*i]+1] = 2*i+1;
        edges[positions[2*i]+2] = SEMIEDGE;
        ddgraph->semiEdges[2*i] = 1;
        degrees[2*i] = 2;
        ddgraph->oneFactor[2*i] = 0;
        //colours
        ddgraph->colours[4*(2*i)+0] = 1;

        edges[positions[2*i+1]+0] = 2*i+2;
        //for the last vertex this will be overwritten when making the connections
        edges[positions[2*i+1]+1] = 2*i;
        edges[positions[2*i+1]+2] = SEMIEDGE;
        ddgraph->semiEdges[2*i+1] = 1;
        degrees[2*i+1] = 2;
        ddgraph->oneFactor[2*i+1] = 0;
        //colours
        ddgraph->colours[4*(2*i+1)+0] = 1;

        //store q3-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = 2*i;
        ddgraph->vertex2UncolouredFactor[2*i+0] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[2*i+1] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[2*i+0] = 3;
        ddgraph->vertex2FactorType[2*i+1] = 3;
    }
    edges[positions[0]+0] = 2*parameter-1;
    edges[positions[2*parameter-1]+0] = 0;
}

void storeBarbedWireNecklaceAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    int parameter = block->parameter;
    int i, vertexLeft, vertexRight;

    //mirror symmetry
    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    vertexLeft = 0;
    vertexRight = 2*parameter-1;

    for(i=0; i<parameter; i++){
        generator[vertexLeft] = vertexRight;
        generator[vertexRight] = vertexLeft;
        vertexLeft++;
        vertexRight--;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);

    if(parameter>1){
        //rotation symmetry
        generator = getIdentity(ddgraph->underlyingGraph->nv);

        for(i=0; i<parameter-1; i++){
            generator[2*i] = 2*(i+1);
            generator[2*i+1] = 2*(i+1)+1;
        }

        generator[2*(parameter-1)] = 0;
        generator[2*(parameter-1)+1] = 1;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    }
}

/*
 * Constructs a double locked diagonal chain DLDC(n). This block has 0 connectors.
 *              _______________
 *             /               |
 *            /                |
 *           /              /  |
 *           o---o-...-o---o  /
 *           |   |     |   | /
 *           |   |     |   |/
 *           o---o-...-o---o
 *          /
 *
 */
void constructDoubleLockedDiagonalChain(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    DEBUGASSERT((*currentVertex) == 0)
    int parameter = block->parameter;
    int i;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    //the connections correspond to the one-factor
    ddgraph->oneFactor[0] = 0;
    ddgraph->oneFactor[1] = 0;
    ddgraph->oneFactor[parameter*4-2] = 0;
    ddgraph->oneFactor[parameter*4-1] = 0;
    //colours
    ddgraph->colours[4*(0)+0] = 1;
    ddgraph->colours[4*(1)+0] = 1;
    ddgraph->colours[4*(4*parameter-2)+0] = 1;
    ddgraph->colours[4*(4*parameter-1)+0] = 1;

    edges[positions[0]+0] = parameter*4-1;
    edges[positions[0]+1] = 1;
    edges[positions[0]+2] = 2;

    edges[positions[1]+0] = SEMIEDGE;
    edges[positions[1]+1] = 0;
    edges[positions[1]+2] = 3;
    degrees[1] = 2;
    positions[1]++;
    ddgraph->semiEdges[1] = 1;

    //store q1-factor for colouring
    ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = 0;
    ddgraph->vertex2UncolouredFactor[0] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[1] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[2] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[3] = ddgraph->uncolouredFactorCount;
    ddgraph->uncolouredFactorCount++;
    ddgraph->vertex2FactorType[0] = 1;
    ddgraph->vertex2FactorType[1] = 1;
    ddgraph->vertex2FactorType[2] = 1;
    ddgraph->vertex2FactorType[3] = 1;

    for(i=1; i<parameter; i++){

        edges[positions[4*i-2]+0] = 4*i-4;
        edges[positions[4*i-2]+1] = 4*i-1;
        edges[positions[4*i-2]+2] = 4*i;

        edges[positions[4*i-1]+0] = 4*i-3;
        edges[positions[4*i-1]+1] = 4*i-2;
        edges[positions[4*i-1]+2] = 4*i+1;

        edges[positions[4*i]+0] = 4*i-2;
        edges[positions[4*i]+1] = 4*i+1;
        edges[positions[4*i]+2] = 4*i+2;

        edges[positions[4*i+1]+0] = 4*i-1;
        edges[positions[4*i+1]+1] = 4*i;
        edges[positions[4*i+1]+2] = 4*i+3;

        //set the one factor
        ddgraph->oneFactor[4*i-2] = 2;
        ddgraph->oneFactor[4*i-1] = 2;
        ddgraph->oneFactor[4*i] = 0;
        ddgraph->oneFactor[4*i+1] = 0;
        //colours
        ddgraph->colours[4*(4*i-2)+2] = 1;
        ddgraph->colours[4*(4*i-1)+2] = 1;
        ddgraph->colours[4*(4*i)+0] = 1;
        ddgraph->colours[4*(4*i+1)+0] = 1;

        //store q1-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = 4*i;
        ddgraph->vertex2UncolouredFactor[4*i+0] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[4*i+1] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[4*i+2] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[4*i+3] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[4*i+0] = 1;
        ddgraph->vertex2FactorType[4*i+1] = 1;
        ddgraph->vertex2FactorType[4*i+2] = 1;
        ddgraph->vertex2FactorType[4*i+3] = 1;
    }

    edges[positions[parameter*4-2]+0] = SEMIEDGE;
    edges[positions[parameter*4-2]+1] = parameter*4-1;
    edges[positions[parameter*4-2]+2] = parameter*4-4;
    degrees[parameter*4-2] = 2;
    positions[parameter*4-2]++;
    ddgraph->semiEdges[parameter*4-2] = 1;

    edges[positions[parameter*4-1]+0] = 0;
    edges[positions[parameter*4-1]+1] = parameter*4-2;
    edges[positions[parameter*4-1]+2] = parameter*4-3;
}

void storeDoubleLockedDiagonalChainAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    int parameter = block->parameter;
    if(parameter==1){
        //mirror symmetry along diagonal
        permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

        generator[0] = 3;
        generator[3] = 0;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        //mirror symmetry along other diagonal

        generator[0] = 0;
        generator[3] = 3;
        generator[1] = 2;
        generator[2] = 1;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    } else {
        //rotation of 180
        permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

        int vertexBottom = 1;
        int vertexTop = parameter*4-2;

        int i;
        for(i = 0; i < 2*parameter; i++){
            generator[vertexBottom] = vertexTop;
            generator[vertexTop] = vertexBottom;
            vertexBottom+=2;
            vertexTop-=2;
        }

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    }
}

/*
 * Constructs a Mobius ladder ML(n). This block has 0 connectors.
 *              _______________
 *             /               |
 *            /              __|_
 *           /              /  | |
 *           o---o-...-o---o  / /
 *           |   |     |   | / /
 *           |   |     |   |/ /
 *           o---o-...-o---o /
 *            \_____________/
 *
 */
void constructMobiusLadder(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    DEBUGASSERT((*currentVertex) == 0)
    int parameter = block->parameter;
    int i;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;

    //the connections correspond to the one-factor
    ddgraph->oneFactor[0] = 0;
    ddgraph->oneFactor[1] = 0;
    ddgraph->oneFactor[parameter*4-2] = 0;
    ddgraph->oneFactor[parameter*4-1] = 0;
    //colours
    ddgraph->colours[4*(0)+0] = 1;
    ddgraph->colours[4*(1)+0] = 1;
    ddgraph->colours[4*(4*parameter-2)+0] = 1;
    ddgraph->colours[4*(4*parameter-1)+0] = 1;

    edges[positions[0]+0] = parameter*4-1;
    edges[positions[0]+1] = 1;
    edges[positions[0]+2] = 2;

    edges[positions[1]+0] = parameter*4-2;
    edges[positions[1]+1] = 0;
    edges[positions[1]+2] = 3;

    //store q1-factor for colouring
    ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = 0;
    ddgraph->vertex2UncolouredFactor[0] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[1] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[2] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[3] = ddgraph->uncolouredFactorCount;
    ddgraph->uncolouredFactorCount++;
    ddgraph->vertex2FactorType[0] = 1;
    ddgraph->vertex2FactorType[1] = 1;
    ddgraph->vertex2FactorType[2] = 1;
    ddgraph->vertex2FactorType[3] = 1;

    for(i=1; i<parameter; i++){

        edges[positions[4*i-2]+0] = 4*i-4;
        edges[positions[4*i-2]+1] = 4*i-1;
        edges[positions[4*i-2]+2] = 4*i;

        edges[positions[4*i-1]+0] = 4*i-3;
        edges[positions[4*i-1]+1] = 4*i-2;
        edges[positions[4*i-1]+2] = 4*i+1;

        edges[positions[4*i]+0] = 4*i-2;
        edges[positions[4*i]+1] = 4*i+1;
        edges[positions[4*i]+2] = 4*i+2;

        edges[positions[4*i+1]+0] = 4*i-1;
        edges[positions[4*i+1]+1] = 4*i;
        edges[positions[4*i+1]+2] = 4*i+3;

        //set the one factor
        ddgraph->oneFactor[4*i-2] = 2;
        ddgraph->oneFactor[4*i-1] = 2;
        ddgraph->oneFactor[4*i] = 0;
        ddgraph->oneFactor[4*i+1] = 0;
        //colours
        ddgraph->colours[4*(4*i-2)+2] = 1;
        ddgraph->colours[4*(4*i-1)+2] = 1;
        ddgraph->colours[4*(4*i)+0] = 1;
        ddgraph->colours[4*(4*i+1)+0] = 1;

        //store q1-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = 4*i;
        ddgraph->vertex2UncolouredFactor[4*i+0] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[4*i+1] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[4*i+2] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[4*i+3] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[4*i+0] = 1;
        ddgraph->vertex2FactorType[4*i+1] = 1;
        ddgraph->vertex2FactorType[4*i+2] = 1;
        ddgraph->vertex2FactorType[4*i+3] = 1;
    }

    edges[positions[parameter*4-2]+0] = 1;
    edges[positions[parameter*4-2]+1] = parameter*4-1;
    edges[positions[parameter*4-2]+2] = parameter*4-4;

    edges[positions[parameter*4-1]+0] = 0;
    edges[positions[parameter*4-1]+1] = parameter*4-2;
    edges[positions[parameter*4-1]+2] = parameter*4-3;
}

void storeMobiusLadderAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    int parameter = block->parameter;
    if(parameter==1){
        //mirror symmetry along diagonal
        permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

        generator[0] = 3;
        generator[3] = 0;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        //mirror symmetry along other diagonal

        generator[0] = 0;
        generator[3] = 3;
        generator[1] = 2;
        generator[2] = 1;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        //mirror symmetry along vertical axis

        generator[0] = 2;
        generator[1] = 3;
        generator[2] = 0;
        generator[3] = 1;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        //mirror symmetry along horizontal axis

        generator[0] = 1;
        generator[1] = 0;
        generator[2] = 3;
        generator[3] = 2;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    } else {
        //rotation of 180
        permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

        int vertexBottom = 1;
        int vertexTop = parameter*4-2;

        int i;
        for(i = 0; i < 2*parameter; i++){
            generator[vertexBottom] = vertexTop;
            generator[vertexTop] = vertexBottom;
            vertexBottom+=2;
            vertexTop-=2;
        }

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);

        generator = getIdentity(ddgraph->underlyingGraph->nv);

        for(i = 0; i < parameter-1; i++){
            generator[4*i] = 4*(i+1);
            generator[4*i+1] = 4*(i+1)+1;
            generator[4*i+2] = 4*(i+1)+2;
            generator[4*i+3] = 4*(i+1)+3;
        }
        generator[4*(parameter-1)] = 1;
        generator[4*(parameter-1)+1] = 0;
        generator[4*(parameter-1)+2] = 3;
        generator[4*(parameter-1)+3] = 2;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);

        generator = getIdentity(ddgraph->underlyingGraph->nv);

        for(i = 0; i < parameter; i++){
            generator[4*i] = 4*i+1;
            generator[4*i+1] = 4*i;
            generator[4*i+2] = 4*i+3;
            generator[4*i+3] = 4*i+2;
        }

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);

        generator = getIdentity(ddgraph->underlyingGraph->nv);

        for(i = 0; i < parameter; i++){
            generator[4*i] = 4*(parameter-i) - 2;
            generator[4*i+1] = 4*(parameter-i) - 1;
            generator[4*i+2] = 4*(parameter-i) - 4;
            generator[4*i+3] = 4*(parameter-i) - 3;
        }

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    }
}

/*
 * Constructs a prism P(n). This block has 0 connectors.
 *             ___________
 *            /           \
 *           o---o-...-o---o
 *           |   |     |   |
 *           |   |     |   |
 *           o---o-...-o---o
 *            \___________/
 *
 */
void constructPrism(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    DEBUGASSERT((*currentVertex) == 0)
    int parameter = block->parameter;
    int i;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;

    //the connections correspond to the one-factor
    ddgraph->oneFactor[0] = 0;
    ddgraph->oneFactor[1] = 0;
    ddgraph->oneFactor[parameter*4-2] = 0;
    ddgraph->oneFactor[parameter*4-1] = 0;
    //colours
    ddgraph->colours[4*(0)+0] = 1;
    ddgraph->colours[4*(1)+0] = 1;
    ddgraph->colours[4*(4*parameter-2)+0] = 1;
    ddgraph->colours[4*(4*parameter-1)+0] = 1;

    edges[positions[0]+0] = parameter*4-2;
    edges[positions[0]+1] = 1;
    edges[positions[0]+2] = 2;

    edges[positions[1]+0] = parameter*4-1;
    edges[positions[1]+1] = 0;
    edges[positions[1]+2] = 3;

    //store q1-factor for colouring
    ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = 0;
    ddgraph->vertex2UncolouredFactor[0] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[1] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[2] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[3] = ddgraph->uncolouredFactorCount;
    ddgraph->uncolouredFactorCount++;
    ddgraph->vertex2FactorType[0] = 1;
    ddgraph->vertex2FactorType[1] = 1;
    ddgraph->vertex2FactorType[2] = 1;
    ddgraph->vertex2FactorType[3] = 1;

    for(i=1; i<parameter; i++){

        edges[positions[4*i-2]+0] = 4*i-4;
        edges[positions[4*i-2]+1] = 4*i-1;
        edges[positions[4*i-2]+2] = 4*i;

        edges[positions[4*i-1]+0] = 4*i-3;
        edges[positions[4*i-1]+1] = 4*i-2;
        edges[positions[4*i-1]+2] = 4*i+1;

        edges[positions[4*i]+0] = 4*i-2;
        edges[positions[4*i]+1] = 4*i+1;
        edges[positions[4*i]+2] = 4*i+2;

        edges[positions[4*i+1]+0] = 4*i-1;
        edges[positions[4*i+1]+1] = 4*i;
        edges[positions[4*i+1]+2] = 4*i+3;

        //set the one factor
        ddgraph->oneFactor[4*i-2] = 2;
        ddgraph->oneFactor[4*i-1] = 2;
        ddgraph->oneFactor[4*i] = 0;
        ddgraph->oneFactor[4*i+1] = 0;
        //colours
        ddgraph->colours[4*(4*i-2)+2] = 1;
        ddgraph->colours[4*(4*i-1)+2] = 1;
        ddgraph->colours[4*(4*i)+0] = 1;
        ddgraph->colours[4*(4*i+1)+0] = 1;

        //store q1-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = 4*i;
        ddgraph->vertex2UncolouredFactor[4*i+0] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[4*i+1] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[4*i+2] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[4*i+3] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[4*i+0] = 1;
        ddgraph->vertex2FactorType[4*i+1] = 1;
        ddgraph->vertex2FactorType[4*i+2] = 1;
        ddgraph->vertex2FactorType[4*i+3] = 1;
    }

    edges[positions[parameter*4-2]+0] = 0;
    edges[positions[parameter*4-2]+1] = parameter*4-1;
    edges[positions[parameter*4-2]+2] = parameter*4-4;

    edges[positions[parameter*4-1]+0] = 1;
    edges[positions[parameter*4-1]+1] = parameter*4-2;
    edges[positions[parameter*4-1]+2] = parameter*4-3;
}

void storePrismAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    int parameter = block->parameter;
    DEBUGASSERTMSG(parameter>1, "Use DDHB instead.")
    //rotation of 180
    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    int vertexBottom = 1;
    int vertexTop = parameter*4-2;

    int i;
    for(i = 0; i < 2*parameter; i++){
        generator[vertexBottom] = vertexTop;
        generator[vertexTop] = vertexBottom;
        vertexBottom+=2;
        vertexTop-=2;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);

    generator = getIdentity(ddgraph->underlyingGraph->nv);

    for(i = 0; i < parameter-1; i++){
        generator[4*i] = 4*(i+1);
        generator[4*i+1] = 4*(i+1)+1;
        generator[4*i+2] = 4*(i+1)+2;
        generator[4*i+3] = 4*(i+1)+3;
    }
    generator[4*(parameter-1)] = 0;
    generator[4*(parameter-1)+1] = 1;
    generator[4*(parameter-1)+2] = 2;
    generator[4*(parameter-1)+3] = 3;

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);

    generator = getIdentity(ddgraph->underlyingGraph->nv);

    for(i = 0; i < parameter; i++){
        generator[4*i] = 4*i+1;
        generator[4*i+1] = 4*i;
        generator[4*i+2] = 4*i+3;
        generator[4*i+3] = 4*i+2;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);

    generator = getIdentity(ddgraph->underlyingGraph->nv);

    for(i = 0; i < parameter; i++){
        generator[4*i] = 4*(parameter-i) - 2;
        generator[4*i+1] = 4*(parameter-i) - 1;
        generator[4*i+2] = 4*(parameter-i) - 4;
        generator[4*i+3] = 4*(parameter-i) - 3;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

    free(generator);

    if(parameter==2){
        generator = getIdentity(ddgraph->underlyingGraph->nv);

        generator[0] = 1;
        generator[1] = 3;
        generator[2] = 0;
        generator[3] = 2;
        generator[4] = 6;
        generator[5] = 4;
        generator[6] = 7;
        generator[7] = 5;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    }
}

/*
 * Constructs a double locked double roof long building DLDLB(n). This block has 0 connectors.
 *             ___________
 *            /           \
 *           o---o-...-o---o
 *           |   |     |   |
 *           |   |     |   |
 *           o---o-...-o---o
 *          /               \
 *
 */
void constructDoubleLockedDoubleroofLongBuilding(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    DEBUGASSERT((*currentVertex) == 0)
    int parameter = block->parameter;
    int i;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    //the connections correspond to the one-factor
    ddgraph->oneFactor[0] = 0;
    ddgraph->oneFactor[1] = 0;
    ddgraph->oneFactor[parameter*4-2] = 0;
    ddgraph->oneFactor[parameter*4-1] = 0;
    //colours
    ddgraph->colours[4*(0)+0] = 1;
    ddgraph->colours[4*(1)+0] = 1;
    ddgraph->colours[4*(4*parameter-2)+0] = 1;
    ddgraph->colours[4*(4*parameter-1)+0] = 1;

    edges[positions[0]+0] = parameter*4-2;
    edges[positions[0]+1] = 1;
    edges[positions[0]+2] = 2;

    edges[positions[1]+0] = SEMIEDGE;
    edges[positions[1]+1] = 0;
    edges[positions[1]+2] = 3;
    degrees[1] = 2;
    positions[1]++;
    ddgraph->semiEdges[1] = 1;

    //store q1-factor for colouring
    ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = 0;
    ddgraph->vertex2UncolouredFactor[0] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[1] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[2] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[3] = ddgraph->uncolouredFactorCount;
    ddgraph->uncolouredFactorCount++;
    ddgraph->vertex2FactorType[0] = 1;
    ddgraph->vertex2FactorType[1] = 1;
    ddgraph->vertex2FactorType[2] = 1;
    ddgraph->vertex2FactorType[3] = 1;

    for(i=1; i<parameter; i++){

        edges[positions[4*i-2]+0] = 4*i-4;
        edges[positions[4*i-2]+1] = 4*i-1;
        edges[positions[4*i-2]+2] = 4*i;

        edges[positions[4*i-1]+0] = 4*i-3;
        edges[positions[4*i-1]+1] = 4*i-2;
        edges[positions[4*i-1]+2] = 4*i+1;

        edges[positions[4*i]+0] = 4*i-2;
        edges[positions[4*i]+1] = 4*i+1;
        edges[positions[4*i]+2] = 4*i+2;

        edges[positions[4*i+1]+0] = 4*i-1;
        edges[positions[4*i+1]+1] = 4*i;
        edges[positions[4*i+1]+2] = 4*i+3;

        //set the one factor
        ddgraph->oneFactor[4*i-2] = 2;
        ddgraph->oneFactor[4*i-1] = 2;
        ddgraph->oneFactor[4*i] = 0;
        ddgraph->oneFactor[4*i+1] = 0;
        //colours
        ddgraph->colours[4*(4*i-2)+2] = 1;
        ddgraph->colours[4*(4*i-1)+2] = 1;
        ddgraph->colours[4*(4*i)+0] = 1;
        ddgraph->colours[4*(4*i+1)+0] = 1;

        //store q1-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = 4*i;
        ddgraph->vertex2UncolouredFactor[4*i+0] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[4*i+1] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[4*i+2] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[4*i+3] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[4*i+0] = 1;
        ddgraph->vertex2FactorType[4*i+1] = 1;
        ddgraph->vertex2FactorType[4*i+2] = 1;
        ddgraph->vertex2FactorType[4*i+3] = 1;
    }

    edges[positions[parameter*4-2]+0] = 0;
    edges[positions[parameter*4-2]+1] = parameter*4-1;
    edges[positions[parameter*4-2]+2] = parameter*4-4;

    edges[positions[parameter*4-1]+0] = SEMIEDGE;
    edges[positions[parameter*4-1]+1] = parameter*4-2;
    edges[positions[parameter*4-1]+2] = parameter*4-3;
    degrees[parameter*4-1] = 2;
    positions[parameter*4-1]++;
    ddgraph->semiEdges[parameter*4-1] = 1;
}

void storeDoubleLockedDoubleroofLongBuildingAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    int parameter = block->parameter;
    if(parameter==1){
        permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

        //mirror symmetry along vertical axis

        generator[0] = 2;
        generator[1] = 3;
        generator[2] = 0;
        generator[3] = 1;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    } else if(parameter==2){
        int i;
        permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

        generator[0] = 3;
        generator[3] = 0;
        generator[5] = 6;
        generator[6] = 5;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        generator = getIdentity(ddgraph->underlyingGraph->nv);

        for(i = 0; i < parameter; i++){
            generator[4*i] = 4*(parameter-i) - 2;
            generator[4*i+1] = 4*(parameter-i) - 1;
            generator[4*i+2] = 4*(parameter-i) - 4;
            generator[4*i+3] = 4*(parameter-i) - 3;
        }

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    } else {
        int i;
        permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

        for(i = 0; i < parameter; i++){
            generator[4*i] = 4*(parameter-i) - 2;
            generator[4*i+1] = 4*(parameter-i) - 1;
            generator[4*i+2] = 4*(parameter-i) - 4;
            generator[4*i+3] = 4*(parameter-i) - 3;
        }

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    }
}

/*
 * Constructs a completely locked hub CLH(n). This block has 0 connectors.
 *
 *          \               /
 *           o---o-...-o---o
 *           |   |     |   |
 *           |   |     |   |
 *           o---o-...-o---o
 *          /               \
 *
 */
void constructCompletelyLockedHub(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    DEBUGASSERT((*currentVertex) == 0)
    int parameter = block->parameter;
    int i;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    //the connections correspond to the one-factor
    ddgraph->oneFactor[0] = 0;
    ddgraph->oneFactor[1] = 0;
    ddgraph->oneFactor[parameter*4-2] = 0;
    ddgraph->oneFactor[parameter*4-1] = 0;
    //colours
    ddgraph->colours[4*(0)+0] = 1;
    ddgraph->colours[4*(1)+0] = 1;
    ddgraph->colours[4*(4*parameter-2)+0] = 1;
    ddgraph->colours[4*(4*parameter-1)+0] = 1;

    edges[positions[0]+0] = SEMIEDGE;
    edges[positions[0]+1] = 1;
    edges[positions[0]+2] = 2;
    degrees[0] = 2;
    positions[0]++;
    ddgraph->semiEdges[0] = 1;

    edges[positions[1]+0] = SEMIEDGE;
    edges[positions[1]+1] = 0;
    edges[positions[1]+2] = 3;
    degrees[1] = 2;
    positions[1]++;
    ddgraph->semiEdges[1] = 1;

    //store q1-factor for colouring
    ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = 0;
    ddgraph->vertex2UncolouredFactor[0] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[1] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[2] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[3] = ddgraph->uncolouredFactorCount;
    ddgraph->uncolouredFactorCount++;
    ddgraph->vertex2FactorType[0] = 1;
    ddgraph->vertex2FactorType[1] = 1;
    ddgraph->vertex2FactorType[2] = 1;
    ddgraph->vertex2FactorType[3] = 1;

    for(i=1; i<parameter; i++){

        edges[positions[4*i-2]+0] = 4*i-4;
        edges[positions[4*i-2]+1] = 4*i-1;
        edges[positions[4*i-2]+2] = 4*i;

        edges[positions[4*i-1]+0] = 4*i-3;
        edges[positions[4*i-1]+1] = 4*i-2;
        edges[positions[4*i-1]+2] = 4*i+1;

        edges[positions[4*i]+0] = 4*i-2;
        edges[positions[4*i]+1] = 4*i+1;
        edges[positions[4*i]+2] = 4*i+2;

        edges[positions[4*i+1]+0] = 4*i-1;
        edges[positions[4*i+1]+1] = 4*i;
        edges[positions[4*i+1]+2] = 4*i+3;

        //set the one factor
        ddgraph->oneFactor[4*i-2] = 2;
        ddgraph->oneFactor[4*i-1] = 2;
        ddgraph->oneFactor[4*i] = 0;
        ddgraph->oneFactor[4*i+1] = 0;
        //colours
        ddgraph->colours[4*(4*i-2)+2] = 1;
        ddgraph->colours[4*(4*i-1)+2] = 1;
        ddgraph->colours[4*(4*i)+0] = 1;
        ddgraph->colours[4*(4*i+1)+0] = 1;

        //store q1-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = 4*i;
        ddgraph->vertex2UncolouredFactor[4*i+0] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[4*i+1] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[4*i+2] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[4*i+3] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[4*i+0] = 1;
        ddgraph->vertex2FactorType[4*i+1] = 1;
        ddgraph->vertex2FactorType[4*i+2] = 1;
        ddgraph->vertex2FactorType[4*i+3] = 1;
    }

    edges[positions[parameter*4-2]+0] = SEMIEDGE;
    edges[positions[parameter*4-2]+1] = parameter*4-1;
    edges[positions[parameter*4-2]+2] = parameter*4-4;
    degrees[parameter*4-2] = 2;
    positions[parameter*4-2]++;
    ddgraph->semiEdges[parameter*4-2] = 1;

    edges[positions[parameter*4-1]+0] = SEMIEDGE;
    edges[positions[parameter*4-1]+1] = parameter*4-2;
    edges[positions[parameter*4-1]+2] = parameter*4-3;
    degrees[parameter*4-1] = 2;
    positions[parameter*4-1]++;
    ddgraph->semiEdges[parameter*4-1] = 1;
}

void storeCompletelyLockedHubAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    int parameter = block->parameter;
    if(parameter==1){
        //mirror symmetry along diagonal
        permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

        generator[0] = 3;
        generator[3] = 0;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        //mirror symmetry along other diagonal

        generator[0] = 0;
        generator[3] = 3;
        generator[1] = 2;
        generator[2] = 1;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        //mirror symmetry along vertical axis

        generator[0] = 2;
        generator[1] = 3;
        generator[2] = 0;
        generator[3] = 1;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        //mirror symmetry along horizontal axis

        generator[0] = 1;
        generator[1] = 0;
        generator[2] = 3;
        generator[3] = 2;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    } else {
        int i;
        permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

        for(i = 0; i < parameter; i++){
            generator[4*i] = 4*i+1;
            generator[4*i+1] = 4*i;
            generator[4*i+2] = 4*i+3;
            generator[4*i+3] = 4*i+2;
        }

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);

        generator = getIdentity(ddgraph->underlyingGraph->nv);

        for(i = 0; i < parameter; i++){
            generator[4*i] = 4*(parameter-i) - 2;
            generator[4*i+1] = 4*(parameter-i) - 1;
            generator[4*i+2] = 4*(parameter-i) - 4;
            generator[4*i+3] = 4*(parameter-i) - 3;
        }

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    }
}

/*
 * Constructs a double roof double floor high building DDHB(n). This block has 0 connectors.
 *
 *           o---o-...-o---o
 *          /|   |     |   |\
 *          \|   |     |   |/
 *           o---o-...-o---o
 *
 */
void constructDoubleroofDoublefloorHighBuilding(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    DEBUGASSERT((*currentVertex) == 0)
    int parameter = block->parameter;
    int i;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;

    //the connections correspond to the one-factor
    ddgraph->oneFactor[0] = 0;
    ddgraph->oneFactor[1] = 0;
    ddgraph->oneFactor[parameter*4-2] = 0;
    ddgraph->oneFactor[parameter*4-1] = 0;
    //colours
    ddgraph->colours[4*(0)+0] = 1;
    ddgraph->colours[4*(1)+0] = 1;
    ddgraph->colours[4*(4*parameter-2)+0] = 1;
    ddgraph->colours[4*(4*parameter-1)+0] = 1;

    edges[positions[0]+0] = 4*parameter;
    edges[positions[0]+1] = 1;
    edges[positions[0]+2] = 2;

    edges[positions[1]+0] = 4*parameter;
    edges[positions[1]+1] = 0;
    edges[positions[1]+2] = 3;

    edges[positions[4*parameter]+0] = 0;
    edges[positions[4*parameter]+1] = 1;

    //store q1-factor for colouring
    ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = 0;
    ddgraph->vertex2UncolouredFactor[0] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[1] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[2] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[3] = ddgraph->uncolouredFactorCount;
    ddgraph->uncolouredFactorCount++;
    ddgraph->vertex2FactorType[0] = 1;
    ddgraph->vertex2FactorType[1] = 1;
    ddgraph->vertex2FactorType[2] = 1;
    ddgraph->vertex2FactorType[3] = 1;

    for(i=1; i<parameter; i++){

        edges[positions[4*i-2]+0] = 4*i-4;
        edges[positions[4*i-2]+1] = 4*i-1;
        edges[positions[4*i-2]+2] = 4*i;

        edges[positions[4*i-1]+0] = 4*i-3;
        edges[positions[4*i-1]+1] = 4*i-2;
        edges[positions[4*i-1]+2] = 4*i+1;

        edges[positions[4*i]+0] = 4*i-2;
        edges[positions[4*i]+1] = 4*i+1;
        edges[positions[4*i]+2] = 4*i+2;

        edges[positions[4*i+1]+0] = 4*i-1;
        edges[positions[4*i+1]+1] = 4*i;
        edges[positions[4*i+1]+2] = 4*i+3;

        //set the one factor
        ddgraph->oneFactor[4*i-2] = 2;
        ddgraph->oneFactor[4*i-1] = 2;
        ddgraph->oneFactor[4*i] = 0;
        ddgraph->oneFactor[4*i+1] = 0;
        //colours
        ddgraph->colours[4*(4*i-2)+2] = 1;
        ddgraph->colours[4*(4*i-1)+2] = 1;
        ddgraph->colours[4*(4*i)+0] = 1;
        ddgraph->colours[4*(4*i+1)+0] = 1;

        //store q1-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = 4*i;
        ddgraph->vertex2UncolouredFactor[4*i+0] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[4*i+1] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[4*i+2] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[4*i+3] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[4*i+0] = 1;
        ddgraph->vertex2FactorType[4*i+1] = 1;
        ddgraph->vertex2FactorType[4*i+2] = 1;
        ddgraph->vertex2FactorType[4*i+3] = 1;
    }

    edges[positions[parameter*4-2]+0] = parameter*4+1;
    edges[positions[parameter*4-2]+1] = parameter*4-1;
    edges[positions[parameter*4-2]+2] = parameter*4-4;

    edges[positions[parameter*4-1]+0] = parameter*4+1;
    edges[positions[parameter*4-1]+1] = parameter*4-2;
    edges[positions[parameter*4-1]+2] = parameter*4-3;

    edges[positions[4*parameter+1]+0] = parameter*4-2;
    edges[positions[4*parameter+1]+1] = parameter*4-1;
}

void storeDoubleroofDoublefloorHighBuildingAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    int parameter = block->parameter;
    if(parameter==1){
        permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

        //mirror symmetry along vertical axis

        generator[0] = 2;
        generator[1] = 3;
        generator[2] = 0;
        generator[3] = 1;
        generator[4] = 5;
        generator[5] = 4;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        //mirror symmetry along horizontal axis

        generator[0] = 1;
        generator[1] = 0;
        generator[2] = 3;
        generator[3] = 2;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    } else {
        permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

        int i;

        for(i = 0; i < parameter; i++){
            generator[4*i] = 4*i+1;
            generator[4*i+1] = 4*i;
            generator[4*i+2] = 4*i+3;
            generator[4*i+3] = 4*i+2;
        }

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);

        generator = getIdentity(ddgraph->underlyingGraph->nv);

        for(i = 0; i < parameter; i++){
            generator[4*i] = 4*(parameter-i) - 2;
            generator[4*i+1] = 4*(parameter-i) - 1;
            generator[4*i+2] = 4*(parameter-i) - 4;
            generator[4*i+3] = 4*(parameter-i) - 3;
        }
        generator[4*parameter] = 4*parameter+1;
        generator[4*parameter+1] = 4*parameter;

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

        free(generator);
    }
}

/*
 * Constructs a double locked double roof high building DLDHB(n). This block has 0 connectors.
 *
 *          \
 *           o---o-...-o---o
 *           |   |     |   |\
 *           |   |     |   |/
 *           o---o-...-o---o
 *          /
 *
 */
void constructDoubleLockedDoubleroofHighBuilding(int *currentVertex, BBLOCK *block, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector){
    DEBUGASSERT((*currentVertex) == 0)
    int parameter = block->parameter;
    int i;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    int *degrees = ddgraph->underlyingGraph->d;

    //the connections correspond to the one-factor
    ddgraph->oneFactor[0] = 0;
    ddgraph->oneFactor[1] = 0;
    ddgraph->oneFactor[parameter*4-2] = 0;
    ddgraph->oneFactor[parameter*4-1] = 0;
    //colours
    ddgraph->colours[4*(0)+0] = 1;
    ddgraph->colours[4*(1)+0] = 1;
    ddgraph->colours[4*(4*parameter-2)+0] = 1;
    ddgraph->colours[4*(4*parameter-1)+0] = 1;

    edges[positions[0]+0] = SEMIEDGE;
    edges[positions[0]+1] = 1;
    edges[positions[0]+2] = 2;
    degrees[0] = 2;
    positions[0]++;
    ddgraph->semiEdges[0] = 1;

    edges[positions[1]+0] = SEMIEDGE;
    edges[positions[1]+1] = 0;
    edges[positions[1]+2] = 3;
    degrees[1] = 2;
    positions[1]++;
    ddgraph->semiEdges[1] = 1;

    //store q1-factor for colouring
    ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = 0;
    ddgraph->vertex2UncolouredFactor[0] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[1] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[2] = ddgraph->uncolouredFactorCount;
    ddgraph->vertex2UncolouredFactor[3] = ddgraph->uncolouredFactorCount;
    ddgraph->uncolouredFactorCount++;
    ddgraph->vertex2FactorType[0] = 1;
    ddgraph->vertex2FactorType[1] = 1;
    ddgraph->vertex2FactorType[2] = 1;
    ddgraph->vertex2FactorType[3] = 1;

    for(i=1; i<parameter; i++){

        edges[positions[4*i-2]+0] = 4*i-4;
        edges[positions[4*i-2]+1] = 4*i-1;
        edges[positions[4*i-2]+2] = 4*i;

        edges[positions[4*i-1]+0] = 4*i-3;
        edges[positions[4*i-1]+1] = 4*i-2;
        edges[positions[4*i-1]+2] = 4*i+1;

        edges[positions[4*i]+0] = 4*i-2;
        edges[positions[4*i]+1] = 4*i+1;
        edges[positions[4*i]+2] = 4*i+2;

        edges[positions[4*i+1]+0] = 4*i-1;
        edges[positions[4*i+1]+1] = 4*i;
        edges[positions[4*i+1]+2] = 4*i+3;

        //set the one factor
        ddgraph->oneFactor[4*i-2] = 2;
        ddgraph->oneFactor[4*i-1] = 2;
        ddgraph->oneFactor[4*i] = 0;
        ddgraph->oneFactor[4*i+1] = 0;
        //colours
        ddgraph->colours[4*(4*i-2)+2] = 1;
        ddgraph->colours[4*(4*i-1)+2] = 1;
        ddgraph->colours[4*(4*i)+0] = 1;
        ddgraph->colours[4*(4*i+1)+0] = 1;

        //store q1-factor for colouring
        ddgraph->uncolouredFactor2Vertex[ddgraph->uncolouredFactorCount] = 4*i;
        ddgraph->vertex2UncolouredFactor[4*i+0] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[4*i+1] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[4*i+2] = ddgraph->uncolouredFactorCount;
        ddgraph->vertex2UncolouredFactor[4*i+3] = ddgraph->uncolouredFactorCount;
        ddgraph->uncolouredFactorCount++;
        ddgraph->vertex2FactorType[4*i+0] = 1;
        ddgraph->vertex2FactorType[4*i+1] = 1;
        ddgraph->vertex2FactorType[4*i+2] = 1;
        ddgraph->vertex2FactorType[4*i+3] = 1;
    }

    edges[positions[parameter*4-2]+0] = parameter*4;
    edges[positions[parameter*4-2]+1] = parameter*4-1;
    edges[positions[parameter*4-2]+2] = parameter*4-4;

    edges[positions[parameter*4-1]+0] = parameter*4;
    edges[positions[parameter*4-1]+1] = parameter*4-2;
    edges[positions[parameter*4-1]+2] = parameter*4-3;

    edges[positions[4*parameter]+0] = parameter*4-2;
    edges[positions[4*parameter]+1] = parameter*4-1;
}

void storeDoubleLockedDoubleroofHighBuildingAutomorphismGenerators(BBLOCK *block, DDGRAPH *ddgraph){
    int parameter = block->parameter;
    permutation *generator = getIdentity(ddgraph->underlyingGraph->nv);

    int i;

    for(i = 0; i < parameter; i++){
        generator[4*i] = 4*i+1;
        generator[4*i+1] = 4*i;
        generator[4*i+2] = 4*i+3;
        generator[4*i+3] = 4*i+2;
    }

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

//========= PHASE 4: ENUMERATION OF DELANEY-DRESS SYMBOLS ===================
void findComponent(DDGRAPH *ddgraph, int *components, int*componentsSize, int *count, int colour1, int colour2){
    int i, j, k, v, nextV, c, size;
    boolean visited[MAXN];
    for(i=0; i<MAXN; i++){
        visited[i]=FALSE;
    }
    int colours[2];
    colours[0]=colour1;
    colours[1]=colour2;
    
    *count = 0;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    
    for(i=0; i<ddgraph->order; i++){
        if(!visited[i]){
            size = 0;
            for(k=0; k<2; k++){
                c = k;
                v = i;
                visited[v]=FALSE;
                while(v!=SEMIEDGE && !visited[v]){
                    size++;
                    visited[v]=TRUE; //TODO: this is a bit messy
                    j = 0;
                    while(ddgraph->colours[4*v+j]!=colours[c]) j++;
                    nextV = edges[3*v+j];
                    if(nextV >= ddgraph->order && nextV!=SEMIEDGE){
                        //dummy edge, so we go directly to the neighbour
                        nextV = edges[positions[nextV]+0] + edges[positions[nextV]+1] - v;
                    }
                    c = (c+1)%2;
                    v = nextV;
                }
            }
            components[*count] = i;
            componentsSize[*count] = size-1; //TODO: this is a bit messy
            (*count)++;
        }
    }
}

void assignComponentLabels(DDGRAPH *ddgraph){
    int s0s1Components[MAXN];
    int s0s1ComponentsSize[MAXN];
    int s0s1ComponentCount = 0;
    int s1s2Components[MAXN];
    int s1s2ComponentsSize[MAXN];
    int s1s2ComponentCount = 0;

    findComponent(ddgraph, s0s1Components, s0s1ComponentsSize, &s0s1ComponentCount, 0, 1);
    findComponent(ddgraph, s1s2Components, s1s2ComponentsSize, &s1s2ComponentCount, 1, 2);
    
}

//========= PHASE 3: HANDLING THE GENERATED DELANEY-DRESS GRAPHS ============
boolean first = TRUE;

void handleColouredDelaneyDressGraph(DDGRAPH *ddgraph){
    edgeColouredGraphsCount++;
    if(symbols){
        assignComponentLabels(ddgraph);
    } else {
        if(outputType=='c'){
            writePregraphColorCodeEdgeColouring(stdout, ddgraph, first);
            first = FALSE;
        }
    }
}

void colourEdge(DDGRAPH *ddgraph, int v1, int v2, int colour){
    int i;
    for(i=0; i<3; i++){
        if(ddgraph->underlyingGraph->e[3*v1+i]==v2){
            ddgraph->colours[4*v1+i]=colour;
        }
        if(ddgraph->underlyingGraph->e[3*v2+i]==v1){
            ddgraph->colours[4*v2+i]=colour;
        }
    }
}

int getImageOfColouring(int colourAssignment, permutation *automorphism, DDGRAPH *ddgraph){
    int i, image = 0;

    for(i = 0; i<ddgraph->uncolouredFactorCount; i++){
        int v = ddgraph->uncolouredFactor2Vertex[i];
        int vImage = automorphism[v];
        int currentState = colourAssignment & (1<<i);
        if(ddgraph->vertex2FactorType[v]==3){
            if(currentState){
                image |= 1 << (ddgraph->vertex2UncolouredFactor[vImage]);
            }
        } else if (ddgraph->vertex2FactorType[v]==1){
            int u = v + 1;
            int uImage = automorphism[u];
            if((vImage-uImage)%2 && currentState){
                image |= 1 << (ddgraph->vertex2UncolouredFactor[vImage]);
            } else if(!((vImage-uImage)%2) && !currentState){
                image |= 1 << (ddgraph->vertex2UncolouredFactor[vImage]);
            }
        } else {
            ERRORMSG("Illegal vertex type for colour assignment.")
        }
    }

    return image;
}

void assignEdgeColours(DDGRAPH *ddgraph){
    int i, colourAssignment;

    int possibleEdgeColouringsCount = (1 << (ddgraph->uncolouredFactorCount));
    int colourOrbits[possibleEdgeColouringsCount];
    int orbitSizes[possibleEdgeColouringsCount];
    int edgeColouringsCount = possibleEdgeColouringsCount;

    for(colourAssignment=0; colourAssignment<possibleEdgeColouringsCount; colourAssignment++){
        colourOrbits[colourAssignment] = colourAssignment;
        orbitSizes[colourAssignment] = 1;
    }

    for(i=0; i<numberOfGenerators[connectionsMade]; i++){
        for(colourAssignment=0; colourAssignment<possibleEdgeColouringsCount; colourAssignment++){
            int colourAssignmentImage = getImageOfColouring(colourAssignment, (*(automorphismGroupGenerators + connectionsMade))[i], ddgraph);
            unionElements(colourOrbits, orbitSizes, &edgeColouringsCount, colourAssignment, colourAssignmentImage);
        }
    }

    //make sure that each element is connected to its root
    for(colourAssignment=0; colourAssignment<possibleEdgeColouringsCount; colourAssignment++){
        findRootOfElement(colourOrbits, colourAssignment);
    }

    for(colourAssignment=0; colourAssignment<possibleEdgeColouringsCount; colourAssignment++){
        if(colourOrbits[colourAssignment] == colourAssignment){
            //assign colours
            for(i=0; i<ddgraph->uncolouredFactorCount; i++){
                int representingVertex = ddgraph->uncolouredFactor2Vertex[i];
                if(ddgraph->vertex2FactorType[representingVertex]==3){
                    if(colourAssignment & (1<<i)){
                        //colour the semi-edges with 0
                        ddgraph->colours[4*representingVertex+2] = 0;
                        ddgraph->colours[4*(representingVertex+1)+2] = 0;
                        ddgraph->colours[4*representingVertex+1] = 2;
                        ddgraph->colours[4*(representingVertex+1)+1] = 2;
                    } else {
                        //colour the semi-edges with 2
                        ddgraph->colours[4*representingVertex+2] = 2;
                        ddgraph->colours[4*(representingVertex+1)+2] = 2;
                        ddgraph->colours[4*representingVertex+1] = 0;
                        ddgraph->colours[4*(representingVertex+1)+1] = 0;
                    }
                } else if(ddgraph->vertex2FactorType[representingVertex]==1){
                    if(colourAssignment & (1<<i)){
                        //colour edge(rv, rv+1) with 2
                        colourEdge(ddgraph, representingVertex, representingVertex+1, 2);
                        colourEdge(ddgraph, representingVertex+2, representingVertex+3, 2);
                        colourEdge(ddgraph, representingVertex, representingVertex+2, 0);
                        colourEdge(ddgraph, representingVertex+1, representingVertex+3, 0);
                    } else {
                        //colour edge(rv, rv+1) with 0
                        colourEdge(ddgraph, representingVertex, representingVertex+1, 0);
                        colourEdge(ddgraph, representingVertex+2, representingVertex+3, 0);
                        colourEdge(ddgraph, representingVertex, representingVertex+2, 2);
                        colourEdge(ddgraph, representingVertex+1, representingVertex+3, 2);
                    }
                } else {
                    ERRORMSG("Incorrect factor type.")
                }
            }
            handleColouredDelaneyDressGraph(ddgraph);
        }
    }

}

void handleDelaneyDressGraph(DDGRAPH *ddgraph){
    graphsCount++;
    if(colouredEdges){
        assignEdgeColours(ddgraph);
    } else {
        if(outputType=='c'){
            if(markedTwoFactors){
                writePregraphColorCode2Factor(stdout, ddgraph, first);
            } else {
                writePregraphCode(stdout, ddgraph, first);
            }
        } else if(outputType=='h'){
            //TODO
            //printDDGraph(ddgraph);
        }
        first = FALSE;
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

    if(!markedTwoFactors){
        int family1 = (blocks+block1)->component;
        int family2 = (blocks+block2)->component;
        int param1 = (blocks+block1)->parameter;
        int param2 = (blocks+block2)->parameter;
        if((type1 == 1 && family1 == 4 && type2 == 3 && family2 == 0 && param2 == 1) ||
                (type1 == 3 && family1 == 0 && type2 == 1 && family2 == 4 && param1 == 1)){
            for(i=0; i<(blocks+block1)->connectorCount; i++){
                if((blocks+block1)->connections[i]==(blocks+block2)){
                    return FALSE;
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
    int i, j, type;

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
    for(type=1; type<=4; type++){
        for(j = 2; j>=0; j--){
            for(i=0; i<ddgraph->order; i++){
                if(ddgraph->semiEdges[i]==j && ddgraph->vertex2FactorType[i]==type){
                    nautyLabelling[counter] = i;
                    counter++;
                }
            }
            if(counter>0){
                nautyPtn[counter-1]=0;
            }
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

    int reverseLabelling[ddgraph->underlyingGraph->nv];
    for(i=0; i<ddgraph->underlyingGraph->nv; i++){
        reverseLabelling[nautyLabelling[i]] = i;
    }

    int smallestConnection = 0;
    int temp1 = reverseLabelling[madeConnections[smallestConnection][0]];
    int temp2 = reverseLabelling[madeConnections[smallestConnection][1]];
    
    int smallestConnectionSmallest = temp1 < temp2 ? temp1 : temp2;
    int smallestConnectionBiggest = temp1 < temp2 ? temp2 : temp1;
    for(i=1; i<madeConnectionsCount; i++){
        int localTemp1 = reverseLabelling[madeConnections[i][0]];
        int localTemp2 = reverseLabelling[madeConnections[i][1]];

        int localSmallest = localTemp1 < localTemp2 ? localTemp1 : localTemp2;
        int localBiggest = localTemp1 < localTemp2 ? localTemp2 : localTemp1;

        if(localSmallest < smallestConnectionSmallest){
            smallestConnectionSmallest = localSmallest;
            smallestConnectionBiggest = localBiggest;
            smallestConnection = madeConnectionOrbits[i];
        } else if(localSmallest == smallestConnectionSmallest && localBiggest < smallestConnectionBiggest){
            smallestConnectionSmallest = localSmallest;
            smallestConnectionBiggest = localBiggest;
            smallestConnection = madeConnectionOrbits[i];
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
                connections[connectionsMade][0] = v1;
                connections[connectionsMade][1] = v2;
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
                        //printConnectionsMade();
                        //printGenerators(ddgraph, 6);
                        //fprintf(stderr, "Graph %2d: ", graphsCount);
                        //printCanonicalLabelling(ddgraph);
                        //fprintf(stderr, "Found graph based on: ");
                        //printHumanReadableComponentList();
                        handleDelaneyDressGraph(ddgraph);
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

void connectComponentList(int vertexCount, DDGRAPH *ddgraph){
    DEBUGTRACE_ENTER
    int blockCount = 0;
    //create an array of blocks based upon the numbers in the global arrays
    BBLOCK *blocks = constructComponentList(&blockCount);
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

    cleanDDGraph(ddgraph);
    //free the memory allocated at the beginning of this method
    freeBuildingBlocks(blocks, blockCount);
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

void handleComponentList(int vertexCount, DDGRAPH *ddgraph){
    DEBUGTRACE_ENTER
    if(!isMaybeRealizableComponentList()){
        DEBUGTRACE_EXIT
        return;
    } else if(!passesSimpleForbiddenConnectionsTest()){
        DEBUGTRACE_EXIT
        return;
    } else {
        if(!moduloEnabled || (splitPointCount%moduloMod == moduloRest)){
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
            if(onlyLists){
                if(outputType=='c'){
                    writeListCodeMultipleBlockList(stdout, vertexCount);
                } else if(outputType=='h'){
                    printHumanReadableComponentList(stdout);
                }
            } else {
                connectComponentList(vertexCount, ddgraph);
            }
        }
        splitPointCount++;
    }
    DEBUGTRACE_EXIT
}

void q4Components(int targetSize, int currentSize, DDGRAPH *ddgraph){
    Q4ComponentCount = targetSize - currentSize; //each q4 component has 1 vertex
    handleComponentList(targetSize, ddgraph);
}

void q3Components(int currentType, int currentParameter, int targetSize, int currentSize, DDGRAPH *ddgraph){
    int i;

    int remainingVertices = targetSize - currentSize;

    for(i = 0; i <= remainingVertices/(2*currentParameter); i++){
        Q3TypeComponentsComponentCount[currentType][currentParameter-1]=i;
        int newSize = currentSize + 2*currentParameter*i; //each q2 component has 2n vertices

        if(targetSize - newSize >= 2*(currentParameter+1)){
            q3Components(currentType, currentParameter+1, targetSize, newSize, ddgraph);
        } else if(currentType+1==Q3TypeComponentsCount){
            q4Components(targetSize, newSize, ddgraph);
        } else {
            q3Components(currentType+1, Q3TypeComponentsSmallestCase[currentType+1], targetSize, newSize, ddgraph);
        }
    }
    Q3TypeComponentsComponentCount[currentType][currentParameter-1]=0; //reset this type to 0
}

void q2Components(int currentType, int currentParameter, int targetSize, int currentSize, DDGRAPH *ddgraph){
    int i;

    int remainingVertices = targetSize - currentSize;

    for(i = 0; i <= remainingVertices/(2*currentParameter); i++){
        Q2TypeComponentsComponentCount[currentType][currentParameter-1]=i;
        int newSize = currentSize + 2*currentParameter*i; //each q2 component has 2n vertices

        if(targetSize - newSize >= 2*(currentParameter+1)){
            q2Components(currentType, currentParameter+1, targetSize, newSize, ddgraph);
        } else if(currentType+1==Q2TypeComponentsCount){
            q3Components(0, Q3TypeComponentsSmallestCase[0], targetSize, newSize, ddgraph);
        } else {
            q2Components(currentType+1, Q2TypeComponentsSmallestCase[currentType+1], targetSize, newSize, ddgraph);
        }
    }
    Q2TypeComponentsComponentCount[currentType][currentParameter-1]=0; //reset this type to 0
}

void q1Components(int currentType, int currentParameter, int targetSize, int currentSize, DDGRAPH *ddgraph){
    int i;

    int remainingVertices = targetSize - currentSize;

    for(i = 0; i <= remainingVertices/(4*currentParameter); i++){
        Q1TypeComponentsComponentCount[currentType][currentParameter-1]=i;
        int newSize = currentSize + 4*currentParameter*i; //each q1 component has 4n vertices

        if(targetSize - newSize >= 4*(currentParameter+1)){
            q1Components(currentType, currentParameter+1, targetSize, newSize, ddgraph);
        } else if(currentType+1==Q1TypeComponentsCount){
            q2Components(0, Q2TypeComponentsSmallestCase[0], targetSize, newSize, ddgraph);
        } else {
            q1Components(currentType+1, Q1TypeComponentsSmallestCase[currentType+1], targetSize, newSize, ddgraph);
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

    Q3TypeComponentsConnectors[0] = 2; //Barbed Wire
    Q3TypeComponentsConnectors[1] = 1; //Locked Barbed Wire

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
    constructBlock[14] = &constructBarbedWire;
    constructBlock[15] = &constructLockedBarbedWire;
    constructBlock[16] = &constructQ4;
    constructBlock[17] = &constructTristar;
    constructBlock[18] = &constructDoubleLockedPearlChain;
    constructBlock[19] = &constructPearlNecklace;
    constructBlock[20] = &constructDoubleLockedBarbedWire;
    constructBlock[21] = &constructBarbedWireNecklace;
    constructBlock[22] = &constructDoubleLockedDiagonalChain;
    constructBlock[23] = &constructMobiusLadder;
    constructBlock[24] = &constructPrism;
    constructBlock[25] = &constructDoubleLockedDoubleroofLongBuilding;
    constructBlock[26] = &constructCompletelyLockedHub;
    constructBlock[27] = &constructDoubleroofDoublefloorHighBuilding;
    constructBlock[28] = &constructDoubleLockedDoubleroofHighBuilding;

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
    storeBlockAutomorphismGenerators[11] = &storeLockedDoubleroofLongBuildingAutomorphismGenerators;
    storeBlockAutomorphismGenerators[12] = &storePearlChainAutomorphismGenerators;
    storeBlockAutomorphismGenerators[13] = NULL;
    storeBlockAutomorphismGenerators[14] = &storeBarbedWireAutomorphismGenerators;
    storeBlockAutomorphismGenerators[15] = NULL;
    storeBlockAutomorphismGenerators[16] = NULL;
    storeBlockAutomorphismGenerators[17] = NULL;
    storeBlockAutomorphismGenerators[18] = &storeDoubleLockedPearlChainAutomorphismGenerators;
    storeBlockAutomorphismGenerators[19] = &storePearlNecklaceAutomorphismGenerators;
    storeBlockAutomorphismGenerators[20] = &storeDoubleLockedBarbedWireAutomorphismGenerators;
    storeBlockAutomorphismGenerators[21] = &storeBarbedWireNecklaceAutomorphismGenerators;
    storeBlockAutomorphismGenerators[22] = &storeDoubleLockedDiagonalChainAutomorphismGenerators;
    storeBlockAutomorphismGenerators[23] = &storeMobiusLadderAutomorphismGenerators;
    storeBlockAutomorphismGenerators[24] = &storePrismAutomorphismGenerators;
    storeBlockAutomorphismGenerators[25] = &storeDoubleLockedDoubleroofLongBuildingAutomorphismGenerators;
    storeBlockAutomorphismGenerators[26] = &storeCompletelyLockedHubAutomorphismGenerators;
    storeBlockAutomorphismGenerators[27] = &storeDoubleroofDoublefloorHighBuildingAutomorphismGenerators;
    storeBlockAutomorphismGenerators[28] = &storeDoubleLockedDoubleroofHighBuildingAutomorphismGenerators;

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
    storeBlocksMapping[14] = &storeBarbedWiresMapping;
    storeBlocksMapping[15] = &storeLockedBarbedWiresMapping;
    storeBlocksMapping[16] = &storeQ4sMapping;
    storeBlocksMapping[17] = NULL;
    storeBlocksMapping[18] = NULL;
    storeBlocksMapping[19] = NULL;
    storeBlocksMapping[20] = NULL;
    storeBlocksMapping[21] = NULL;
    storeBlocksMapping[22] = NULL;
    storeBlocksMapping[23] = NULL;
    storeBlocksMapping[24] = NULL;
    storeBlocksMapping[25] = NULL;
    storeBlocksMapping[26] = NULL;
    storeBlocksMapping[27] = NULL;
    storeBlocksMapping[28] = NULL;

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
    blockName[17] = "T";
    blockName[18] = "DLPC(%d)";
    blockName[19] = "PN(%d)";
    blockName[20] = "DLBW(%d)";
    blockName[21] = "BWN(%d)";
    blockName[22] = "DLDC(%d)";
    blockName[23] = "ML(%d)";
    blockName[24] = "P(%d)";
    blockName[25] = "DLDLB(%d)";
    blockName[26] = "CLH(%d)";
    blockName[27] = "DDHB(%d)";
    blockName[28] = "DLDHB(%d)";

}

void cleanComponentStatistics(){
    free(Q1TypeComponentsComponentCount);
    free(Q2TypeComponentsComponentCount);
    free(Q3TypeComponentsComponentCount);
}

void initStatistics(){
    componentListsCount = 0;
    graphsCount = 0;
    edgeColouredGraphsCount = 0;
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

void cleanNautyOptions(){
    free(canonGraph.d);
    free(canonGraph.v);
    free(canonGraph.e);
}

//====================== START =======================
void finishGraph(DDGRAPH *ddgraph){
    int i;

    ddgraph->underlyingGraph->nv = ddgraph->order + ddgraph->dummyVertexCount;

    for(i = 0; i < ddgraph->underlyingGraph->nv; i++){
        ddgraph->underlyingGraph->nde += ddgraph->underlyingGraph->d[i];
    }
}

void handleSingleBlockComponentList(BBLOCK * bblock, int order, DDGRAPH * ddgraph){
    if(!moduloEnabled || (splitPointCount%moduloMod == moduloRest)){
        componentListsCount++;
        if(onlyLists){
            if(outputType=='c'){
                writeListCodeSingleBlockList(stdout, order, bblock);
            } else if(outputType=='h'){
                printBlockName(stdout, 1, buildingBlockTypeToNumber(bblock), bblock->parameter, TRUE);
                fprintf(stdout, "\n");
            }
        } else {
            int vertexToBlock[ddgraph->underlyingGraph->vlen];
            int vertexToConnector[ddgraph->underlyingGraph->vlen];
            constructBuildingBlockListAsGraph(bblock, 1, ddgraph, vertexToBlock, vertexToConnector);
            numberOfGenerators[0] = 0;
            storeInitialGenerators(bblock, 1, ddgraph);
            handleDelaneyDressGraph(ddgraph);
            cleanDDGraph(ddgraph);
        }
    }
    splitPointCount++;
}

void addTristar(DDGRAPH * ddgraph){
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 5, 0, 0, 0);
    handleSingleBlockComponentList(bblock, 1, ddgraph);
    free(bblock);
}

void addDoubleLockedPearlChain(DDGRAPH * ddgraph, int parameter){
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 5, 1, parameter, 0);
    handleSingleBlockComponentList(bblock, parameter*2, ddgraph);
    free(bblock);
}

void addPearlNecklace(DDGRAPH * ddgraph, int parameter){
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 5, 2, parameter, 0);
    handleSingleBlockComponentList(bblock, parameter*2, ddgraph);
    free(bblock);
}

void addDoubleLockedBarbedWire(DDGRAPH * ddgraph, int parameter){
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 6, 0, parameter, 0);
    handleSingleBlockComponentList(bblock, parameter*2, ddgraph);
    free(bblock);
}

void addBarbedWireNecklace(DDGRAPH * ddgraph, int parameter){
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 6, 1, parameter, 0);
    handleSingleBlockComponentList(bblock, parameter*2, ddgraph);
    free(bblock);
}

void addDoubleLockedDiagonalChain(DDGRAPH * ddgraph, int parameter){
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 6, 2, parameter, 0);
    handleSingleBlockComponentList(bblock, parameter*4, ddgraph);
    free(bblock);
}

void addMobiusLadder(DDGRAPH * ddgraph, int parameter){
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 6, 3, parameter, 0);
    handleSingleBlockComponentList(bblock, parameter*4, ddgraph);
    free(bblock);
}

void addPrism(DDGRAPH * ddgraph, int parameter){
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 6, 4, parameter, 0);
    handleSingleBlockComponentList(bblock, parameter*4, ddgraph);
    free(bblock);
}

void addDoubleLockedDoubleroofLongBuilding(DDGRAPH * ddgraph, int parameter){
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 6, 5, parameter, 0);
    handleSingleBlockComponentList(bblock, parameter*4, ddgraph);
    free(bblock);
}

void addCompletelyLockedHub(DDGRAPH * ddgraph, int parameter){
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 6, 6, parameter, 0);
    handleSingleBlockComponentList(bblock, parameter*4, ddgraph);
    free(bblock);
}

void addDoubleroofDoublefloorHighBuilding(DDGRAPH * ddgraph, int parameter){
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 6, 7, parameter, 0);
    handleSingleBlockComponentList(bblock, parameter*4, ddgraph);
    free(bblock);
}

void addDoubleLockedDoubleroofHighBuilding(DDGRAPH * ddgraph, int parameter){
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 6, 8, parameter, 0);
    handleSingleBlockComponentList(bblock, parameter*4, ddgraph);
    free(bblock);
}

void extraUnconstructableGraphs_1(DDGRAPH * ddgraph){
    addTristar(ddgraph);
}

void extraUnconstructableGraphs_2(DDGRAPH * ddgraph){
    addDoubleLockedPearlChain(ddgraph, 1);
    addPearlNecklace(ddgraph, 1);

    if(markedTwoFactors){
        addDoubleLockedBarbedWire(ddgraph, 1);

        //isomorph to DLPC(1) in unmarked case
        addBarbedWireNecklace(ddgraph, 1);
    }
}

void extraUnconstructableGraphs_4(DDGRAPH * ddgraph){
    addDoubleLockedPearlChain(ddgraph, 2);
    addPearlNecklace(ddgraph, 2);
    addBarbedWireNecklace(ddgraph, 2);
    addDoubleLockedDiagonalChain(ddgraph, 1);
    addMobiusLadder(ddgraph, 1);

    if(markedTwoFactors){
        addDoubleLockedBarbedWire(ddgraph, 2);
        addCompletelyLockedHub(ddgraph, 1);
        addDoubleroofDoublefloorHighBuilding(ddgraph, 1);
        addDoubleLockedDoubleroofHighBuilding(ddgraph, 1);
    }
}

void extraUnconstructableGraphs_4n(DDGRAPH * ddgraph, int targetSize){
    addDoubleLockedPearlChain(ddgraph, targetSize/2);
    addPearlNecklace(ddgraph, targetSize/2);
    addBarbedWireNecklace(ddgraph, targetSize/2);
    addDoubleLockedDiagonalChain(ddgraph, targetSize/4);
    addMobiusLadder(ddgraph, targetSize/4);
    addPrism(ddgraph, targetSize/4);
    addDoubleLockedDoubleroofLongBuilding(ddgraph, targetSize/4);

    if(markedTwoFactors){
        addDoubleLockedBarbedWire(ddgraph, targetSize/2);
        addCompletelyLockedHub(ddgraph, targetSize/4);
        addDoubleroofDoublefloorHighBuilding(ddgraph, targetSize/4);
        addDoubleLockedDoubleroofHighBuilding(ddgraph, targetSize/4);
    }
}

void extraUnconstructableGraphs_4n2(DDGRAPH * ddgraph, int targetSize){
    addDoubleLockedPearlChain(ddgraph, targetSize/2);
    addPearlNecklace(ddgraph, targetSize/2);
    addBarbedWireNecklace(ddgraph, targetSize/2);

    if(markedTwoFactors){
        addDoubleLockedBarbedWire(ddgraph, targetSize/2);
    }
    
}

void extraUnconstructableGraphs(DDGRAPH * ddgraph, int targetSize){
    if(targetSize == 1){
        extraUnconstructableGraphs_1(ddgraph);
    } else if(targetSize == 2){
        extraUnconstructableGraphs_2(ddgraph);
    } else if(targetSize == 4){
        extraUnconstructableGraphs_4(ddgraph);
    } else if(targetSize % 4 == 2){
        extraUnconstructableGraphs_4n2(ddgraph, targetSize);
    } else if(targetSize % 4 == 0){
        extraUnconstructableGraphs_4n(ddgraph, targetSize);
    }
}

void startGeneration(int targetSize){

    initComponentsStatic();
    initComponents(targetSize);
    initStatistics();
    initNautyOptions(targetSize);

    DDGRAPH * ddgraph = getNewDDGraph(targetSize);

    q1Components(0, Q1TypeComponentsSmallestCase[0], targetSize, 0, ddgraph);

    extraUnconstructableGraphs(ddgraph, targetSize);

    freeDDGraph(ddgraph);

    fprintf(stderr, "Found %llu component list%s.\n", componentListsCount, componentListsCount==1 ? (char *)"" : (char *)"s");
    if(!onlyLists){
        fprintf(stderr, "Found %llu Delaney-Dress graph%s%s.\n",
                graphsCount,
                graphsCount==1 ? (char *)"" : (char *)"s",
                markedTwoFactors ? (char *)" with marked 2-factors" : (char *)"");
        if(colouredEdges){
            fprintf(stderr, "Found %llu edge-coloured Delaney-Dress graph%s.\n",
                edgeColouredGraphsCount,
                edgeColouredGraphsCount==1 ? (char *)"" : (char *)"s");
        }
    }
    if(moduloEnabled){
        fprintf(stderr, "Generated only part %llu of %llu.\n", moduloRest+1, moduloMod);
    }
    
    cleanNautyOptions();
    cleanComponentStatistics();
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
        boolean singleComponentList = FALSE;
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
                realVertexCount += 4*count*parameter;
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
                realVertexCount += 2*count*parameter;
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
                realVertexCount += 2*count*parameter;
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
            } else if(type==5){
                if(!fscanf(f, "%d", &family)){
                    ERRORMSG("Error while parsing file.")
                }
                if(!fscanf(f, "%d", &parameter)){
                    ERRORMSG("Error while parsing file.")
                }
                if(family==0){
                    realVertexCount=1;
                } else if(family==1 || family==2){
                    realVertexCount=2*parameter;
                } else {
                    ERRORMSG("Error while parsing file: illegal family for type 5.")
                }
                int last = -1;
                if(fscanf(f, "%d", &last)!=-1 && last!=0){
                    ERRORMSG("Error while parsing file: Block of type 5 can only be in list of size 1.")
                }
                singleComponentList = TRUE;
                break;
            } else if(type==6){
                if(!fscanf(f, "%d", &family)){
                    ERRORMSG("Error while parsing file.")
                }
                if(!fscanf(f, "%d", &parameter)){
                    ERRORMSG("Error while parsing file.")
                }
                if(family==0 || family==1){
                    realVertexCount=2*parameter;
                } else if(family==2 || family==3 || family==4 || family==5 || family==6 || family==7 || family==8){
                    realVertexCount=4*parameter;
                } else {
                    ERRORMSG("Error while parsing file: illegal family for type 6.")
                }
                int last = -1;
                if(fscanf(f, "%d", &last)!=-1 && last!=0){
                    ERRORMSG("Error while parsing file: Block of type 6 can only be in list of size 1.")
                }
                singleComponentList = TRUE;
                break;
            } else {
                ERRORMSG("Error while parsing file: illegal type.")
            }
        }

        if(realVertexCount!=vertexCount){
            ERRORMSG("Error while parsing file: incorrect vertex count.")
        }

        DDGRAPH * ddgraph = getNewDDGraph(vertexCount);

        if(singleComponentList){
            BBLOCK buildingBlock;
            initBuildingBlock(&buildingBlock, type, family, parameter, 0);
            handleSingleBlockComponentList(&buildingBlock, realVertexCount, ddgraph);
        } else {
            handleComponentList(vertexCount, ddgraph);
            freeComponents();
        }
        freeDDGraph(ddgraph);

        cleanNautyOptions();
    }

    fprintf(stderr, "Read %llu component list%s.\n", componentListsCount, componentListsCount==1 ? (char *)"" : (char *)"s");
    if(!onlyLists){
        fprintf(stderr, "Found %llu Delaney-Dress graph%s%s.\n",
                graphsCount,
                graphsCount==1 ? (char *)"" : (char *)"s",
                markedTwoFactors ? (char *)" with marked 2-factors" : (char *)"");
        if(colouredEdges){
            fprintf(stderr, "Found %llu edge-coloured Delaney-Dress graph%s.\n",
                edgeColouredGraphsCount,
                edgeColouredGraphsCount==1 ? (char *)"" : (char *)"s");
        }
    }
    if(moduloEnabled){
        fprintf(stderr, "Generated only part %llu of %llu.\n", moduloRest+1, moduloMod);
    }

    cleanComponentStatistics();
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
    char *moduloString;

    while ((c = getopt(argc, argv, "hl:Ltcso:m:")) != -1) {
        switch (c) {
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            case 'l':
                listFilename = optarg;
                break;
            case 'L':
                onlyLists = TRUE;
                break;
            case 't':
                markedTwoFactors = TRUE;
                break;
            case 'c':
                markedTwoFactors = TRUE;
                colouredEdges = TRUE;
                break;
            case 's':
                markedTwoFactors = TRUE;
                colouredEdges = TRUE;
                symbols = TRUE;
                break;
            case 'o':
                outputType = optarg[0];
                switch (outputType) {
                    case 'n': //no output (default)
                    case 'c': //listcode, pregraphcode, pregraphcolorcode
                    case 'h': //human-readable
                        break;
                    default:
                        fprintf(stderr, "Illegal output format %c.\n", c);
                        usage(name);
                        return 1;
                }
                break;
            case 'm':
                //modulo
                moduloEnabled = TRUE;
                moduloString = optarg;
                moduloRest = atoi(moduloString);
                moduloString = strchr(moduloString, ':');
                if(moduloString==NULL){
                    fprintf(stderr, "Illegal format for modulo.\n");
                    usage(name);
                    return EXIT_FAILURE;
                }
                moduloMod = atoi(moduloString+1);
                if (moduloRest >= moduloMod) {
                    fprintf(stderr, "Illegal format for modulo: rest must be smaller than mod.\n");
                    usage(name);
                    return EXIT_FAILURE;
                }
                if (moduloRest < 0) {
                    fprintf(stderr, "Illegal format for modulo: rest must be positive.\n");
                    usage(name);
                    return EXIT_FAILURE;
                }
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
