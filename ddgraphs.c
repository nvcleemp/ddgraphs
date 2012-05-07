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

//======================== Mathematical methods =============================

int gcd(int a, int b){
    if(b==0){
       return a;
    } else {
       return gcd(b, a%b);
    }
}

int lcm(int a, int b){
    return a/gcd(a,b)*b;
}

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

void writeDDSymbol(FILE *f, DDGRAPH *ddgraph, COLOURCOMPONENTS *s0s1Components, COLOURCOMPONENTS *s1s2Components, unsigned long long int n1, unsigned long long int n2){
    int i, j, v;
    int *edges = ddgraph->underlyingGraph->e;
    int *positions = ddgraph->underlyingGraph->v;
    
    fprintf(f, "<%llu.%llu:%d 2:", n1, n2, ddgraph->order);
    
    //s_i actions
    for(i=0; i<2; i++){
        j = 0;
        while(ddgraph->colours[j]!=i) j++;
        if(edges[j]==SEMIEDGE){
            fprintf(f, "%d", 1);
        } else if(edges[j]>=ddgraph->order) {
            fprintf(f, "%d", edges[positions[edges[j]]+0] + edges[positions[edges[j]]+1]+1);
        } else {
            fprintf(f, "%d", edges[j]+1);
        }
        for(v=1; v<ddgraph->order; v++){
            j = 0;
            while(ddgraph->colours[4*v+j]!=i) j++;
            if(edges[3*v+j]==SEMIEDGE){
                fprintf(f, " %d", v+1);
            } else {
                int n = edges[3*v+j];
                if(n>=ddgraph->order) {
                    n = edges[positions[edges[3*v+j]]+0] + edges[positions[edges[3*v+j]]+1]-v;
                }
                if(v<=n){
                    fprintf(f, " %d", n+1);
                }
            }
        }
        fprintf(f, ",");
    }
    j = 0;
    while(ddgraph->colours[j]!=2) j++;
    if(edges[j]==SEMIEDGE){
        fprintf(f, "%d", 1);
    } else if(edges[j]>=ddgraph->order) {
        fprintf(f, "%d", edges[positions[edges[j]]+0] + edges[positions[edges[j]]+1]+1);
    } else {
        fprintf(f, "%d", edges[j]+1);
    }
    for(v=1; v<ddgraph->order; v++){
        j = 0;
        while(ddgraph->colours[4*v+j]!=2) j++;
        if(edges[3*v+j]==SEMIEDGE){
            fprintf(f, " %d", v+1);
        } else {
            int n = edges[3*v+j];
            if(n>=ddgraph->order) {
                n = edges[positions[edges[3*v+j]]+0] + edges[positions[edges[3*v+j]]+1]-v;
            }
            if(v<=n){
                fprintf(f, " %d", n+1);
            }
        }
    }
    fprintf(f, ":");
    
    //m_i,i+1 values
    fprintf(f, "%d", s0s1Components->componentLabels[0]);
    for(i=1; i<s0s1Components->componentCount; i++){
        fprintf(f, " %d", s0s1Components->componentLabels[i]);
    }
    fprintf(f, ",");
    fprintf(f, "%d", s1s2Components->componentLabels[0]);
    for(i=1; i<s1s2Components->componentCount; i++){
        fprintf(f, " %d", s1s2Components->componentLabels[i]);
    }
    
    fprintf(f, ">\n");
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



void printColouredDDGraph(CDDGRAPH *graph){
    fprintf(stderr, "CDDGRAPH %p\n", graph);
    fprintf(stderr, "===================\n");
    
    int i, j;
    for(i=0; i<graph->order; i++){
        fprintf(stderr, "%2d) ", i+1);
        for(j=0; j<3; j++){
            if(graph->adjacency[i][j]==SEMIEDGE){
                fprintf(stderr, " S ");
            } else {
                fprintf(stderr, "%2d ", graph->adjacency[i][j]+1);
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

permutation globalIdentity[MAXN];

permutation *getIdentity(int order){
    if(order>MAXN){
        ERRORMSG("Order to large: recompile with a larger MAXN.")
    }
    permutation i;
    for(i = 0; i<order; i++){
        globalIdentity[i]=i;
    }
    return globalIdentity;
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

int buildingBlockParametersToNumber(int type, int component){
    if(type==1){
        return component;
    } else if(type==2){
        return Q1TypeComponentsCount + component;
    } else if(type==3){
        return Q1TypeComponentsCount + Q2TypeComponentsCount + component;
    } else if(type==4){
        return Q1TypeComponentsCount + Q2TypeComponentsCount + Q3TypeComponentsCount;
    } else if(type==5){
        return Q1TypeComponentsCount + Q2TypeComponentsCount
                + Q3TypeComponentsCount + 1 + component;
    } else {
        return Q1TypeComponentsCount + Q2TypeComponentsCount 
                + Q3TypeComponentsCount + 1 + NoConnectorsFixedColouringComponentsCount
                + component;
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
    int smallestVertex = *currentVertex;
    for(i=smallestVertex; i<smallestVertex + 4*(block->parameter); i++){
        ddgraph->vertexComponents[0][i] = smallestVertex;
    }
    ddgraph->vertexPartition[0][*currentVertex] = 0;
    ddgraph->vertexPartition[0][(*currentVertex)+1] = 1;
    ddgraph->vertexPartition[0][(*currentVertex)+(block->parameter-1)*4+2] = 1;
    ddgraph->vertexPartition[0][(*currentVertex)+(block->parameter-1)*4+3] = 0;

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
    int smallestVertex = *currentVertex;
    for(i=smallestVertex; i<smallestVertex + 4*(block->parameter); i++){
        ddgraph->vertexComponents[0][i] = smallestVertex;
    }
    ddgraph->vertexPartition[0][(*currentVertex)+1] = 1;
    ddgraph->vertexPartition[0][(*currentVertex)+(block->parameter-1)*4+2] = 1;
    ddgraph->vertexPartition[0][(*currentVertex)+(block->parameter-1)*4+3] = 0;

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
    int smallestVertex = *currentVertex;
    for(i=smallestVertex; i<smallestVertex + 4*(block->parameter); i++){
        ddgraph->vertexComponents[0][i] = smallestVertex;
    }
    ddgraph->vertexPartition[0][(*currentVertex)+1] = 1;
    ddgraph->vertexPartition[0][(*currentVertex)+(block->parameter-1)*4+2] = 1;

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
    int smallestVertex = *currentVertex;
    for(i=smallestVertex; i<smallestVertex + 4*(block->parameter); i++){
        ddgraph->vertexComponents[0][i] = smallestVertex;
    }

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
    int smallestVertex = *currentVertex;
    for(i=smallestVertex; i<smallestVertex + 4*(block->parameter); i++){
        ddgraph->vertexComponents[0][i] = smallestVertex;
    }
    ddgraph->vertexPartition[0][*currentVertex] = 0;
    ddgraph->vertexPartition[0][(*currentVertex)+1] = 1;

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
    int smallestVertex = *currentVertex;
    for(i=smallestVertex; i<smallestVertex + 4*(block->parameter); i++){
        ddgraph->vertexComponents[0][i] = smallestVertex;
    }
    ddgraph->vertexPartition[0][*currentVertex] = 0;
    ddgraph->vertexPartition[0][(*currentVertex)+1] = 1;


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
    int smallestVertex = *currentVertex;
    for(i=smallestVertex; i<smallestVertex + 4*(block->parameter); i++){
        ddgraph->vertexComponents[0][i] = smallestVertex;
    }
    ddgraph->vertexPartition[0][(*currentVertex)+1] = 1;
    ddgraph->vertexPartition[0][(*currentVertex)+(block->parameter-1)*4+3] = 0;


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
    int smallestVertex = *currentVertex;
    for(i=smallestVertex; i<smallestVertex + 4*(block->parameter); i++){
        ddgraph->vertexComponents[0][i] = smallestVertex;
    }
    ddgraph->vertexPartition[0][(*currentVertex)+1] = 1;
    ddgraph->vertexPartition[0][(*currentVertex)+(block->parameter-1)*4+3] = 0;
    

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
    int smallestVertex = *currentVertex;
    for(i=smallestVertex; i<smallestVertex + 4*(block->parameter); i++){
        ddgraph->vertexComponents[0][i] = smallestVertex;
    }
    ddgraph->vertexPartition[0][(*currentVertex)+1] = 1;

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
    int smallestVertex = *currentVertex;
    for(i=smallestVertex; i<smallestVertex + 4*(block->parameter); i++){
        ddgraph->vertexComponents[0][i] = smallestVertex;
    }
    ddgraph->vertexPartition[0][*currentVertex] = 0;

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
    int smallestVertex = *currentVertex;
    for(i=smallestVertex; i<smallestVertex + 4*(block->parameter); i++){
        ddgraph->vertexComponents[0][i] = smallestVertex;
    }
    ddgraph->vertexPartition[0][*currentVertex] = 0;

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
    int smallestVertex = *currentVertex;
    for(i=smallestVertex; i<smallestVertex + 4*(block->parameter); i++){
        ddgraph->vertexComponents[0][i] = smallestVertex;
    }
    ddgraph->vertexPartition[0][(*currentVertex)+1] = 1;

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
    int smallestVertex = *currentVertex;
    for(i=smallestVertex; i<smallestVertex + 2*(block->parameter); i++){
        ddgraph->vertexComponents[0][i] = smallestVertex;
    }
    ddgraph->vertexPartition[0][*currentVertex] = 0;
    ddgraph->vertexPartition[0][*currentVertex + 1 + 2*(block->parameter-1)] = 1;

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
    int smallestVertex = *currentVertex;
    for(i=smallestVertex; i<smallestVertex + 2*(block->parameter); i++){
        ddgraph->vertexComponents[0][i] = smallestVertex;
    }
    ddgraph->vertexPartition[0][*currentVertex] = 0;
    
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
    int smallestVertex = *currentVertex;
    for(i=smallestVertex; i<smallestVertex + 2*(block->parameter); i++){
        ddgraph->vertexComponents[0][i] = smallestVertex;
    }
    ddgraph->vertexPartition[0][*currentVertex] = 0;
    ddgraph->vertexPartition[0][*currentVertex + 1 + 2*(block->parameter-1)] = 1;

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
    int smallestVertex = *currentVertex;
    for(i=smallestVertex; i<smallestVertex + 2*(block->parameter); i++){
        ddgraph->vertexComponents[0][i] = smallestVertex;
    }
    ddgraph->vertexPartition[0][*currentVertex] = 0;

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
    ddgraph->vertexComponents[0][*currentVertex] = *currentVertex;
    ddgraph->vertexPartition[0][*currentVertex] = 0;
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


        generator = getIdentity(ddgraph->underlyingGraph->nv);

        for(i = 0; i < parameter; i++){
            generator[4*i] = 4*i+1;
            generator[4*i+1] = 4*i;
            generator[4*i+2] = 4*i+3;
            generator[4*i+3] = 4*i+2;
        }

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);


        generator = getIdentity(ddgraph->underlyingGraph->nv);

        for(i = 0; i < parameter; i++){
            generator[4*i] = 4*(parameter-i) - 2;
            generator[4*i+1] = 4*(parameter-i) - 1;
            generator[4*i+2] = 4*(parameter-i) - 4;
            generator[4*i+3] = 4*(parameter-i) - 3;
        }

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

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


    generator = getIdentity(ddgraph->underlyingGraph->nv);

    for(i = 0; i < parameter; i++){
        generator[4*i] = 4*i+1;
        generator[4*i+1] = 4*i;
        generator[4*i+2] = 4*i+3;
        generator[4*i+3] = 4*i+2;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);


    generator = getIdentity(ddgraph->underlyingGraph->nv);

    for(i = 0; i < parameter; i++){
        generator[4*i] = 4*(parameter-i) - 2;
        generator[4*i+1] = 4*(parameter-i) - 1;
        generator[4*i+2] = 4*(parameter-i) - 4;
        generator[4*i+3] = 4*(parameter-i) - 3;
    }

    storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);


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


        generator = getIdentity(ddgraph->underlyingGraph->nv);

        for(i = 0; i < parameter; i++){
            generator[4*i] = 4*(parameter-i) - 2;
            generator[4*i+1] = 4*(parameter-i) - 1;
            generator[4*i+2] = 4*(parameter-i) - 4;
            generator[4*i+3] = 4*(parameter-i) - 3;
        }

        storeGenerator(0, generator, NULL, 0, 0, ddgraph->underlyingGraph->nv);

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

/* Translates a DDGRAPH (which has its adjacency lists ordered in construction
 * order) to a CDDGRAPH (which has its adjacency lists ordered according to
 * colour of the corresponding edges).
 */
void translateDdgraphToCddgraph(DDGRAPH *ddgraph, CDDGRAPH *cddgraph){
    int i, j;
    cddgraph->order = ddgraph->order;
    
    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;

    for (i = 0; i < ddgraph->order; i++) {
        for (j = 0; j < 3; j++) {
            int neighbour = edges[3*i+j]; //don't use the current positions, but the initial ones!!
            if(neighbour >= ddgraph->order && neighbour != SEMIEDGE){
                //neighbour is a dummy vertex
                neighbour = edges[positions[neighbour]+0] + edges[positions[neighbour]+1]-i;
            }
            cddgraph->adjacency[i][ddgraph->colours[4*i+j]] = neighbour;
        }
    }
}



/* Calculates a certificate and a canonical labelling for a given coloured 
 * Delaney-Dress graph and a given vertex that needs to receive label 1.
 */
void getCertificateFromVertex(CDDGRAPH *cddgraph, COLOURCOMPONENTS *s0s1Components,
        COLOURCOMPONENTS *s1s2Components, int *vertex2s0s1components, int *vertex2s1s2components,
        int *certificate, int *old2new, int *new2old, int vertex){
    int i,j;
    int enqueued[MAXN];
    for(i=0;i<MAXN;i++) enqueued[i]=FALSE;
            
    int queue[MAXN];
    int head = 0, tail = 0;
    int currentLabel = 0;
    
    queue[tail++] = vertex;
    enqueued[vertex]=TRUE;

    while(head<tail){
        int v = queue[head++];
        new2old[currentLabel]=v;
        old2new[v]=currentLabel++;
        for(i=0; i<3; i++){
            if(cddgraph->adjacency[v][i]!=SEMIEDGE && !enqueued[cddgraph->adjacency[v][i]]){
                queue[tail++] = cddgraph->adjacency[v][i];
                enqueued[cddgraph->adjacency[v][i]] = TRUE;
            }
        }
    }

    for(i=0; i<cddgraph->order; i++){
        int old = new2old[i];
        for(j=0; j<3; j++){
            int neighbour = cddgraph->adjacency[old][j];
            if(neighbour == SEMIEDGE){
                certificate[i*3 + j] = SEMIEDGE;
            } else {
                certificate[i*3 + j] = old2new[neighbour];
            }
        }
        certificate[(cddgraph->order)*3 + i*2 + 0] =
                s0s1Components->componentLabels[vertex2s0s1components[old]];
        certificate[(cddgraph->order)*3 + i*2 + 1] =
                s1s2Components->componentLabels[vertex2s1s2components[old]];
    }

}

/* Returns a certificate and a canonical labelling for a given coloured 
 * Delaney-Dress graph and a given vertex that needs to receive label 1.
 */
void getSmallestCertificateForComponent(CDDGRAPH *cddgraph, COLOURCOMPONENTS *s0s1Components,
        COLOURCOMPONENTS *s1s2Components, int *vertex2s0s1components, int *vertex2s1s2components,
        int *certificate, int currentComponent, int *v2c){
    int old2new[MAXN];
    int new2old[MAXN];
    int tempold2new[MAXN];
    int tempnew2old[MAXN];
    int tempCertificate[MAXN*5];
    int currentVertex = 0;

    //first we calculate the certificate if we label vertex 0 with 0 and look
    //for the canonical labeling
    
    while(currentVertex < cddgraph->order && v2c[currentVertex]!=currentComponent){
        currentVertex++;
    }
    
    if(currentVertex == cddgraph->order){
        ERRORMSG("Error while calculating canonical form of coloured graph.")
    }
    
    getCertificateFromVertex(cddgraph, s0s1Components, s1s2Components,
            vertex2s0s1components, vertex2s1s2components, certificate, 
            old2new, new2old, currentVertex);
    
    for(currentVertex++ ; currentVertex<cddgraph->order; currentVertex++){
        if(v2c[currentVertex] == currentComponent){
            getCertificateFromVertex(cddgraph, s0s1Components, s1s2Components,
                        vertex2s0s1components, vertex2s1s2components, tempCertificate, 
                        tempold2new, tempnew2old, currentVertex);
            
            int i=0, j;
            while(i<5*cddgraph->order && certificate[i]==tempCertificate[i]) i++;
            if(i<5*(cddgraph->order)){
                if(certificate[i]>tempCertificate[i]){
                    //copy tempCertificate to certificate
                    //TODO: use memcpy
                    for(j=i; j<5*cddgraph->order; j++){
                        certificate[j] = tempCertificate[j];
                    }
                    for(j=0; j<cddgraph->order; j++){
                        old2new[j] = tempold2new[j];
                        new2old[j] = tempnew2old[j];
                    }
                }
            }
        }
    }
}

/*
 * This method is used by the quickSortComponentsList method to partition the
 * given list between left and right in elements larger than the pivot and
 * elements smaller than the pivot. A colour component a is larger then a colour
 * component b if a contains more vertices than b, or, in case they have equal 
 * number of vertices, a has semi-edges and b doesn't.
 */
int partitionComponentsList(int *indices, COLOURCOMPONENTS *componentsList, int left, int right, int pivot){
    int temp, i, newPivotPosition;
    
    int pivotValue = componentsList->componentSizes[indices[pivot]];
    boolean pivotHasSemiEdges = componentsList->containsSemiEdge[indices[pivot]];
    
    //switch the right element and the pivot
    temp = indices[pivot];
    indices[pivot] = indices[right];
    indices[right] = temp;
   
    newPivotPosition = left;
    if(pivotHasSemiEdges){
        for(i = left; i<right; i++){
            if(componentsList->componentSizes[indices[i]]>pivotValue){
                temp = indices[newPivotPosition];
                indices[newPivotPosition] = indices[i];
                indices[i] = temp;
                newPivotPosition++;
            }
        }
    } else {
        for(i = left; i<right; i++){
            if(componentsList->componentSizes[indices[i]]>pivotValue || 
                    (componentsList->componentSizes[indices[i]] == pivotValue &&
                     componentsList->containsSemiEdge[indices[i]])){
                temp = indices[newPivotPosition];
                indices[newPivotPosition] = indices[i];
                indices[i] = temp;
                newPivotPosition++;
            }
        }
    }
    
    //put pivot in its right place
    temp = indices[newPivotPosition];
    indices[newPivotPosition] = indices[right];
    indices[right] = temp;
    
    return newPivotPosition;
}

/*
 * Sorts the given list of indices such that the components they refer to are
 * sorted from large to small.
 */
void quickSortComponentsList(int *indices, COLOURCOMPONENTS *componentsList, int left, int right){
    if(left<right){
        //take the center element as pivot
        int pivot = (left + right)/2;
        
        //partition the list
        int newPivotPosition = partitionComponentsList(indices, componentsList, left, right, pivot);
        
        //recursively sort the two parts
        quickSortComponentsList(indices, componentsList, left, newPivotPosition-1);
        quickSortComponentsList(indices, componentsList, newPivotPosition+1, right);
    }
}

void handleDelaneyDressSymbol(DDGRAPH *ddgraph, COLOURCOMPONENTS *s0s1Components, COLOURCOMPONENTS *s1s2Components){
        symbolsCount++;
        if(giveStatistics){
            //Delaney-Dress graph
            if(ddgraph->order > maxOrderStatistic){
                maxOrderStatistic = ddgraph->order;
            }
            if(ddgraph->order < minOrderStatistic){
                minOrderStatistic = ddgraph->order;
            }
            
            int i;
            //faces
            if(s0s1Components->componentCount > maxFaceOrbitCountStatistic){
                maxFaceOrbitCountStatistic = s0s1Components->componentCount;
            }
            if(s0s1Components->componentCount < minFaceOrbitCountStatistic){
                minFaceOrbitCountStatistic = s0s1Components->componentCount;
            }
            for(i=0; i<s0s1Components->componentCount; i++){
                if(s0s1Components->componentLabels[i] > maxFaceSizeStatistic){
                    maxFaceSizeStatistic = s0s1Components->componentLabels[i];
                }
                if(s0s1Components->componentLabels[i] < minFaceSizeStatistic){
                    minFaceSizeStatistic = s0s1Components->componentLabels[i];
                }
            }
            
            //vertices
            if(s1s2Components->componentCount > maxVertexOrbitCountStatistic){
                maxVertexOrbitCountStatistic = s1s2Components->componentCount;
            }
            if(s1s2Components->componentCount < minVertexOrbitCountStatistic){
                minVertexOrbitCountStatistic = s1s2Components->componentCount;
            }
            for(i=0; i<s1s2Components->componentCount; i++){
                if(s1s2Components->componentLabels[i] > maxVertexDegreeStatistic){
                    maxVertexDegreeStatistic = s1s2Components->componentLabels[i];
                }
                if(s1s2Components->componentLabels[i] < minVertexDegreeStatistic){
                    minVertexDegreeStatistic = s1s2Components->componentLabels[i];
                }
            }
            
            //edges
            if(edgeOrbitCount > maxEdgeOrbitCountStatistic){
                maxEdgeOrbitCountStatistic = edgeOrbitCount;
            }
            if(edgeOrbitCount < minEdgeOrbitCountStatistic){
                minEdgeOrbitCountStatistic = edgeOrbitCount;
            }
        }
        if(outputType=='c'){
            writeDDSymbol(stdout, ddgraph, s0s1Components, s1s2Components, 1, symbolsCount);
        }
}

void findComponents(DDGRAPH *ddgraph, COLOURCOMPONENTS *colourComponents, int *vertex2component){
    int i, j, v, nextV, c, size;
    boolean visited[MAXN];
    for(i=0; i<MAXN; i++){
        visited[i]=FALSE;
    }
    int colours[2];
    colours[0]=colourComponents->colour1;
    colours[1]=colourComponents->colour2;

    colourComponents->componentCount = 0;

    //store some pointers to limit the amount of typing in the next lines
    int *positions = ddgraph->underlyingGraph->v;
    int *edges = ddgraph->underlyingGraph->e;
    
    for(i=0; i<ddgraph->order; i++){
        if(!visited[i]){
            size = 0;
            //visit the component at the side of colour1
            c = 0;
            v = i;
            while(v!=SEMIEDGE && !visited[v]){
                size++;
                visited[v] = TRUE;
                vertex2component[v] = colourComponents->componentCount;
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
            //visit the component at the side of colour2
            c = 1;
            j = 0;
            while(ddgraph->colours[4*i+j]!=colours[c]) j++;
            v = edges[3*i+j];
            if(v >= ddgraph->order && v!=SEMIEDGE){
                //dummy edge, so we go directly to the neighbour
                v = edges[positions[v]+0] + edges[positions[v]+1] - i;
            }
            c = 0;
            while(v!=SEMIEDGE && !visited[v]){
                size++;
                visited[v]=TRUE;
                vertex2component[v] = colourComponents->componentCount;
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
            //component has been completely visited
            colourComponents->components[colourComponents->componentCount] = i;
            colourComponents->componentSizes[colourComponents->componentCount] = size;
            if(v==SEMIEDGE){
                colourComponents->containsSemiEdge[colourComponents->componentCount] = TRUE;
            } else {
                colourComponents->containsSemiEdge[colourComponents->componentCount] = FALSE;
            }
            (colourComponents->componentCount)++;
        }
    }
}

void validateDelaneyDressSymbol(DDGRAPH *ddgraph, COLOURCOMPONENTS *s0s1Components, COLOURCOMPONENTS *s1s2Components){
    int i;
    
    PROFILINGINCREMENT(possibleAssignments)
    int denominator = s0s1Components->componentLabels[0];
    
    for(i=1; i<s0s1Components->componentCount; i++){
        denominator = lcm(denominator, s0s1Components->componentLabels[i]);
    }
    
    for(i=0; i<s1s2Components->componentCount; i++){
        denominator = lcm(denominator, s1s2Components->componentLabels[i]);
    }
    
    int numerator = 0;
    for(i=0; i<s0s1Components->componentCount; i++){
        numerator += 2*denominator/s0s1Components->componentLabels[i]*s0s1Components->componentSizes[i];
    }
    for(i=0; i<s1s2Components->componentCount; i++){
        numerator += 2*denominator/s1s2Components->componentLabels[i]*s1s2Components->componentSizes[i];
    }
    
    if(numerator==(ddgraph->order)*denominator){
        PROFILINGINCREMENT(validAssignments)
        handleDelaneyDressSymbol(ddgraph, s0s1Components, s1s2Components);
    }
}

/*
 * Old method used to assign the labels to the s1s2 components.
 * This method just tries all the combinations and let the user
 * filter out the isomorphic once himself
 */
void bruteForces1s2LabelAssignment(DDGRAPH *ddgraph, COLOURCOMPONENTS *s0s1Components, COLOURCOMPONENTS *s1s2Components, int current){
    int i;
    for(i=minVertexDegree; i<=maxVertexDegree; i++){
        int r = s1s2Components->containsSemiEdge[current] ? 
                        s1s2Components->componentSizes[current] : 
                        s1s2Components->componentSizes[current]/2;
        if(i%r==0){
            s1s2Components->componentLabels[current]=i;
            if(current+1==s1s2Components->componentCount){
                validateDelaneyDressSymbol(ddgraph, s0s1Components, s1s2Components);
            } else {
                bruteForces1s2LabelAssignment(ddgraph, s0s1Components, s1s2Components, current+1);
            }
        }
    }
}

/*
 * Old method used to assign the labels to the s0s1 components.
 * This method just tries all the combinations and let the user
 * filter out the isomorphic once himself
 */
void bruteForces0s1LabelAssignment(DDGRAPH *ddgraph, COLOURCOMPONENTS *s0s1Components, COLOURCOMPONENTS *s1s2Components, int current){
    int i;
    for(i=minFaceSize; i<=maxFaceSize; i++){
        int r = s0s1Components->containsSemiEdge[current] ? 
                        s0s1Components->componentSizes[current] : 
                        s0s1Components->componentSizes[current]/2;
        if(i%r==0){
            s0s1Components->componentLabels[current]=i;
            if(current+1==s0s1Components->componentCount){
                bruteForces1s2LabelAssignment(ddgraph, s0s1Components, s1s2Components, 0);
            } else {
                bruteForces0s1LabelAssignment(ddgraph, s0s1Components, s1s2Components, current+1);
            }
        }
    }
}

//this array is used to store the different certificates to compare
int refineCCPWorkspace[MAXN][MAXN*5];
int refineCCPComparisons[MAXN][MAXN];

void refineColourComponentPartition(CDDGRAPH *cddgraph, COLOURCOMPONENTS *s0s1Components, COLOURCOMPONENTS *s1s2Components,
                                       int *vertex2s0s1components, int *vertex2s1s2components, int depth, int currentPartitionStart,
                                       int partitioning[][MAXN], int partitioningSize[][MAXN], int labelling[][MAXN],
                                       COLOURCOMPONENTS *components, int *v2c){
    int i, j, k;
    int newOrder[partitioningSize[depth][currentPartitionStart]];
    int newLabelling[partitioningSize[depth][currentPartitionStart]]; //only stores the part that needs to be modified
    for(i = 0; i<partitioningSize[depth][currentPartitionStart]; i++){ //for each component in this part
        getSmallestCertificateForComponent(cddgraph, s0s1Components, s1s2Components, vertex2s0s1components, vertex2s1s2components,
                                        refineCCPWorkspace[i], labelling[depth][currentPartitionStart + i], v2c);
        newOrder[i] = i;
    }
    
    int certificateLength = (cddgraph->order)*5;
    for(i = 0; i<partitioningSize[depth][currentPartitionStart]-1; i++){
        for(j = i+1; j<partitioningSize[depth][currentPartitionStart]; j++){ //for pair of components in this part
            //compare certificate of i with certificate of j
            k = 0;
            while(k < certificateLength && refineCCPWorkspace[i][k] == refineCCPWorkspace[j][k]){
                k++;
            }
            if(k == certificateLength){
                refineCCPComparisons[i][j] = 0;
                refineCCPComparisons[j][i] = 0;
            } else if(refineCCPWorkspace[i][k] < refineCCPWorkspace[j][k]){
                refineCCPComparisons[i][j] = -1;
                refineCCPComparisons[j][i] = 1;
            } else {
                refineCCPComparisons[i][j] = 1;
                refineCCPComparisons[j][i] = -1;
            }
        }
    }
    
    //rearrange newOrder to reflect the new order of the current partition
    //we use insertion sort because this will usually be very small lists
    for(i = 1; i < partitioningSize[depth][currentPartitionStart]; i++){
        j = i-1;
        while(j >= 0 && refineCCPComparisons[i][newOrder[j]] < 0){
            newOrder[j+1] = newOrder[j];
            j--;
        }
        newOrder[j+1] = i;
    }
    
    //store the modified part of the new labelling
    for(i = 0; i < partitioningSize[depth][currentPartitionStart]; i++){
        newLabelling[i] = labelling[depth][currentPartitionStart + newOrder[i]];
    }
    
    //modify the labelling and the partitioning
    labelling[depth][currentPartitionStart] = newLabelling[0];
    int currentPartitionSize = 1;
    int lastPartitionStart = currentPartitionStart;
    int originalPartitionSize = partitioningSize[depth][currentPartitionStart];
    for(i = 1; i < originalPartitionSize; i++){
        labelling[depth][currentPartitionStart + i] = newLabelling[i];
        if(refineCCPComparisons[newOrder[i]][newOrder[i-1]] != 0){
            //start new partition
            partitioning[depth][currentPartitionStart + i] = 0;
            for(j = lastPartitionStart; j < currentPartitionStart + i; j++){
                partitioningSize[depth][j] = currentPartitionSize;
            }
            currentPartitionSize = 1;
            lastPartitionStart = currentPartitionStart + i;
        } else {
            currentPartitionSize++;
        }
    }
    for(j = lastPartitionStart; j < currentPartitionStart + originalPartitionSize; j++){
        partitioningSize[depth][j] = currentPartitionSize;
    }
}

void assignLabelsToNexts1s2Component(DDGRAPH *ddgraph, CDDGRAPH *cddgraph, COLOURCOMPONENTS *s0s1Components, COLOURCOMPONENTS *s1s2Components,
                                       int *vertex2s0s1components, int *vertex2s1s2components, int currentComponent,
                                       int depth){
    
    //first determine maximum and minimum for the label
    int localMaximumVertexDegree, localMinimumVertexDegree;
    if(s1s2Components->containsSemiEdge[s1s2labelling[depth][currentComponent]]){
        localMaximumVertexDegree = 6*s1s2Components->componentSizes[s1s2labelling[depth][currentComponent]];
        localMinimumVertexDegree = s1s2Components->componentSizes[s1s2labelling[depth][currentComponent]];
    } else {
        DEBUGASSERTMSG(s1s2Components->componentSizes[s1s2labelling[depth][currentComponent]]%2==0,
                "A component without semi-edges should have an even size.")
        localMaximumVertexDegree = 3*s1s2Components->componentSizes[s1s2labelling[depth][currentComponent]];
        localMinimumVertexDegree = s1s2Components->componentSizes[s1s2labelling[depth][currentComponent]]/2;
    }
    if(localMaximumVertexDegree > maxVertexDegree){
        localMaximumVertexDegree = maxVertexDegree;
    }
    if(localMinimumVertexDegree < minVertexDegree){
        localMinimumVertexDegree = minVertexDegree;
    }
    
    if(localMinimumVertexDegree < minimumPartitionVertexDegree[depth][currentComponent]){
        localMinimumVertexDegree = minimumPartitionVertexDegree[depth][currentComponent];
    }
    
    int i, j;
    int current = s1s2labelling[depth][currentComponent];
    for(i=localMinimumVertexDegree; i<=localMaximumVertexDegree; i++){
        if(!forbiddenFaceSizesTable[i]){
            int r = s1s2Components->containsSemiEdge[current] ? 
                        s1s2Components->componentSizes[current] : 
                        s1s2Components->componentSizes[current]/2;
            if(i%r==0){
                s1s2Components->componentLabels[current]=i;
                
                //update denominator and numerator
                if(depth > 0){
                    partialCurvatureDenominator2[depth] = lcm(partialCurvatureDenominator2[depth-1], i);
                    partialCurvatureNumerator2[depth] = 
                            partialCurvatureNumerator2[depth-1]*partialCurvatureDenominator2[depth]/partialCurvatureDenominator2[depth-1];
                } else {
                    partialCurvatureDenominator2[depth] = partialCurvatureDenominatorLasts0s1;
                    partialCurvatureNumerator2[depth] = partialCurvatureNumeratorLasts0s1;
                }
                partialCurvatureNumerator2[depth] += 2*partialCurvatureDenominator2[depth]/i*s1s2Components->componentSizes[current];
                
                if(partialCurvatureNumerator2[depth] <= partialCurvatureDenominator2[depth] * ddgraph->order) {
                    //proceed to the next component
                    if(currentComponent==s1s2Components->componentCount-1){
                        validateDelaneyDressSymbol(ddgraph, s0s1Components, s1s2Components);
                    } else {
                        //copy necessary array values to the next level
                        for(j = currentComponent+1; j < s1s2Components->componentCount; j++){
                            minimumPartitionVertexDegree[depth+1][j] = minimumPartitionVertexDegree[depth][j];
                            s1s2labelling[depth+1][j] = s1s2labelling[depth][j];
                        }
                        
                        
                        //copy partitioning information to next level
                        //in case the current component belongs to a partition with size > 1 we split this partition
                        int partitioningPosition = currentComponent + 1;
                        if(s1s2componentsPartitionSize[depth][currentComponent]>1){
                            /* The component we just labelled no longer belongs to this partition
                             * because it is labelled and the remaining components in this partition
                             * aren't. So the second component now marks the start of the partition
                             * and the size of the partition is one smaller.
                             * We however do have to set the minimum label for all the components
                             * in this partition to i to prevent isomorphic copies from being generated.
                             */
                            s1s2componentsPartitioning[depth+1][currentComponent+1] = 0;
                            s1s2componentsPartitionSize[depth+1][currentComponent+1] = 
                                    s1s2componentsPartitionSize[depth][currentComponent] - 1;
                            minimumPartitionVertexDegree[depth+1][currentComponent+1] = i;
                            
                            for(j = 1; j < s1s2componentsPartitionSize[depth][currentComponent] - 1; j++){
                                s1s2componentsPartitioning[depth+1][currentComponent+1+j] = 
                                        s1s2componentsPartitioning[depth][currentComponent+1+j];
                                s1s2componentsPartitionSize[depth+1][currentComponent+1+j] = 
                                        s1s2componentsPartitionSize[depth][currentComponent] - 1;
                                
                                //in case this partition has size > 1: set minima for remaining components in this partition
                                minimumPartitionVertexDegree[depth+1][currentComponent+1+j] = i;
                            }
                            
                            partitioningPosition = currentComponent + s1s2componentsPartitionSize[depth][currentComponent];
                        }
                        
                        
                        for(j = partitioningPosition; j < s1s2Components->componentCount; j++){
                            s1s2componentsPartitioning[depth+1][j] = 
                                    s1s2componentsPartitioning[depth+1][j];
                            s1s2componentsPartitionSize[depth+1][j] = 
                                    s1s2componentsPartitionSize[depth][j];
                        }
                        
                        //if the partition still has a size > 1 we try to refine this partition
                        if(s1s2componentsPartitionSize[depth][currentComponent]>2){
                            //refine
                            refineColourComponentPartition(cddgraph, s0s1Components, s1s2Components,
                                    vertex2s0s1components, vertex2s1s2components, depth + 1,
                                    currentComponent + 1, s1s2componentsPartitioning, s1s2componentsPartitionSize,
                                    s1s2labelling, s1s2Components, vertex2s1s2components);
                            /* If the partition still has size > 1, then the components
                             * in this partition really are isomorphic so we don't need to check
                             * if this is the case: this will be handled by the next recursion.
                             */
                        }
                        
                        assignLabelsToNexts1s2Component(ddgraph, cddgraph, s0s1Components,
                                s1s2Components, vertex2s0s1components, 
                                vertex2s1s2components, currentComponent+1, depth + 1);
                    }
                }
            }
        }
    }
}

void startAssignings1s2ComponentLabels(DDGRAPH *ddgraph, CDDGRAPH *cddgraph, COLOURCOMPONENTS *s0s1Components,
                        COLOURCOMPONENTS *s1s2Components, int *vertex2s0s1components, 
                        int *vertex2s1s2components, int s0s1depth){
    
    int i, j;
    
    //initialize some arrays
    for(i=0; i<MAXN; i++){
        s1s2labelling[0][i] = i;
        minimumPartitionVertexDegree[0][i] = 0;
    }
    
    //sort the components from large to small
    quickSortComponentsList(s1s2labelling[0], s1s2Components, 0, s1s2Components->componentCount-1);
    
    //determine the partitioning: currently only for s1s2-components
    s0s1componentsPartitioning[0][0]=0;
    int lastStart = 0;
    int currentSize = 1;
    for(i = 1; i < s1s2Components->componentCount; i++){
        int prevLabel = s1s2labelling[0][i-1];
        int currLabel = s1s2labelling[0][i];
        if((s1s2Components->containsSemiEdge[prevLabel] != s1s2Components->containsSemiEdge[currLabel]) ||
                (s1s2Components->componentSizes[prevLabel] != s1s2Components->componentSizes[currLabel])){
            s1s2componentsPartitioning[0][i] = 0;
            for(j = lastStart; j < i; j++){
                s1s2componentsPartitionSize[0][j] = currentSize;
            }
            currentSize = 1;
            lastStart = i;
        } else {
            s1s2componentsPartitioning[0][i] = i;
            currentSize++;
        }
    }
    for(j = lastStart; j < s1s2Components->componentCount; j++){
        s1s2componentsPartitionSize[0][j] = currentSize;
    }
    if(s1s2componentsPartitionSize[0][0]>1){
        //try to refine partition
        refineColourComponentPartition(cddgraph, s0s1Components, s1s2Components, vertex2s0s1components, vertex2s1s2components,
                                        0, 0, s1s2componentsPartitioning, s1s2componentsPartitionSize,  s1s2labelling,
                                        s1s2Components, vertex2s1s2components);
    }
    
    partialCurvatureDenominatorLasts0s1 = partialCurvatureDenominator[s0s1depth];
    partialCurvatureNumeratorLasts0s1 = partialCurvatureNumerator[s0s1depth];
    
    assignLabelsToNexts1s2Component(ddgraph, cddgraph, s0s1Components, s1s2Components, vertex2s0s1components, vertex2s1s2components, 0, 0);
}

void assignLabelsToNexts0s1Component(DDGRAPH *ddgraph, CDDGRAPH *cddgraph, COLOURCOMPONENTS *s0s1Components, COLOURCOMPONENTS *s1s2Components,
                                       int *vertex2s0s1components, int *vertex2s1s2components, int currentComponent,
                                       int depth){
    
    //first determine maximum and minimum for the label
    int localMaximumFaceSize, localMinimumFaceSize;
    if(s0s1Components->containsSemiEdge[s0s1labelling[depth][currentComponent]]){
        localMaximumFaceSize = 6*s0s1Components->componentSizes[s0s1labelling[depth][currentComponent]];
        localMinimumFaceSize = s0s1Components->componentSizes[s0s1labelling[depth][currentComponent]];
    } else {
        DEBUGASSERTMSG(s0s1Components->componentSizes[s0s1labelling[depth][currentComponent]]%2==0,
                "A component without semi-edges should have an even size.")
        localMaximumFaceSize = 3*s0s1Components->componentSizes[s0s1labelling[depth][currentComponent]];
        localMinimumFaceSize = s0s1Components->componentSizes[s0s1labelling[depth][currentComponent]]/2;
    }
    if(localMaximumFaceSize > maxFaceSize){
        localMaximumFaceSize = maxFaceSize;
    }
    if(localMinimumFaceSize < minFaceSize){
        localMinimumFaceSize = minFaceSize;
    }
    
    if(localMinimumFaceSize < minimumPartitionFaceSize[depth][currentComponent]){
        localMinimumFaceSize = minimumPartitionFaceSize[depth][currentComponent];
    }
    
    int i, j;
    int current = s0s1labelling[depth][currentComponent];
    for(i=localMinimumFaceSize; i<=localMaximumFaceSize; i++){
        if(!forbiddenFaceSizesTable[i]){
            int r = s0s1Components->containsSemiEdge[current] ? 
                        s0s1Components->componentSizes[current] : 
                        s0s1Components->componentSizes[current]/2;
            if(i%r==0){
                s0s1Components->componentLabels[current]=i;
                
                //update denominator and numerator
                if(depth > 0){
                    partialCurvatureDenominator[depth] = lcm(partialCurvatureDenominator[depth-1], i);
                    partialCurvatureNumerator[depth] = 
                            partialCurvatureNumerator[depth-1]*partialCurvatureDenominator[depth]/partialCurvatureDenominator[depth-1];
                } else {
                    partialCurvatureDenominator[depth] = i;
                    partialCurvatureNumerator[depth] = 0;
                }
                partialCurvatureNumerator[depth] += 2*partialCurvatureDenominator[depth]/i*s0s1Components->componentSizes[current];
                
                if(partialCurvatureNumerator[depth] < partialCurvatureDenominator[depth] * ddgraph->order) {
                    //proceed to the next component
                    if(currentComponent==s0s1Components->componentCount-1){
                        startAssignings1s2ComponentLabels(ddgraph, cddgraph, s0s1Components,
                                s1s2Components, vertex2s0s1components, vertex2s1s2components,
                                depth);
                    } else {
                        //copy necessary array values to the next level
                        for(j = currentComponent+1; j < s0s1Components->componentCount; j++){
                            minimumPartitionFaceSize[depth+1][j] = minimumPartitionFaceSize[depth][j];
                            s0s1labelling[depth+1][j] = s0s1labelling[depth][j];
                        }
                        
                        
                        //copy partitioning information to next level
                        //in case the current component belongs to a partition with size > 1 we split this partition
                        int partitioningPosition = currentComponent + 1;
                        if(s0s1componentsPartitionSize[depth][currentComponent]>1){
                            /* The component we just labelled no longer belongs to this partition
                             * because it is labelled and the remaining components in this partition
                             * aren't. So the second component now marks the start of the partition
                             * and the size of the partition is one smaller.
                             * We however do have to set the minimum label for all the components
                             * in this partition to i to prevent isomorphic copies from being generated.
                             */
                            s0s1componentsPartitioning[depth+1][currentComponent+1] = 0;
                            s0s1componentsPartitionSize[depth+1][currentComponent+1] = 
                                    s0s1componentsPartitionSize[depth][currentComponent] - 1;
                            minimumPartitionFaceSize[depth+1][currentComponent+1] = i;
                            
                            for(j = 1; j < s0s1componentsPartitionSize[depth][currentComponent] - 1; j++){
                                s0s1componentsPartitioning[depth+1][currentComponent+1+j] = 
                                        s0s1componentsPartitioning[depth+1][currentComponent+1+j];
                                s0s1componentsPartitionSize[depth+1][currentComponent+1+j] = 
                                        s0s1componentsPartitionSize[depth][currentComponent] - 1;
                                
                                //in case this partition has size > 1: set minima for remaining components in this partition
                                minimumPartitionFaceSize[depth+1][currentComponent+1+j] = i;
                            }
                            
                            partitioningPosition = currentComponent + s0s1componentsPartitionSize[depth][currentComponent];
                        }
                        
                        
                        for(j = partitioningPosition; j < s0s1Components->componentCount; j++){
                            s0s1componentsPartitioning[depth+1][j] = 
                                    s0s1componentsPartitioning[depth+1][j];
                            s0s1componentsPartitionSize[depth+1][j] = 
                                    s0s1componentsPartitionSize[depth][j];
                        }
                        
                        //if the partition still has a size > 1 we try to refine this partition
                        if(s0s1componentsPartitionSize[depth][currentComponent]>2){
                            //refine
                            refineColourComponentPartition(cddgraph, s0s1Components, s1s2Components,
                                    vertex2s0s1components, vertex2s1s2components, depth + 1,
                                    currentComponent + 1, s0s1componentsPartitioning, s0s1componentsPartitionSize,
                                    s0s1labelling, s0s1Components, vertex2s0s1components);
                            /* If the partition still has size > 1, then the components
                             * in this partition really are isomorphic so we don't need to check
                             * if this is the case: this will be handled by the next recursion.
                             */
                        }
                        
                        assignLabelsToNexts0s1Component(ddgraph, cddgraph, s0s1Components,
                                s1s2Components, vertex2s0s1components, 
                                vertex2s1s2components, currentComponent+1, depth + 1);
                    }
                }
            }
        }
    }
}

void startAssignings0s1ComponentLabels(DDGRAPH *ddgraph, COLOURCOMPONENTS *s0s1Components,
                        COLOURCOMPONENTS *s1s2Components, int *vertex2s0s1components, 
                        int *vertex2s1s2components){
    CDDGRAPH cddgraph;
    
    //we first translate the DDGRAPH to a CDDGRAPH
    translateDdgraphToCddgraph(ddgraph, &cddgraph);
    
    int i, j;
    
    //initialize some arrays
    for(i=0; i<MAXN; i++){
        s0s1labelling[0][i] = i;
        minimumPartitionFaceSize[0][i] = 0;
    }
    
    //sort the components from large to small
    quickSortComponentsList(s0s1labelling[0], s0s1Components, 0, s0s1Components->componentCount-1);
    
    //determine the partitioning: currently only for s0s1-components
    s0s1componentsPartitioning[0][0]=0;
    int lastStart = 0;
    int currentSize = 1;
    for(i = 1; i < s0s1Components->componentCount; i++){
        int prevLabel = s0s1labelling[0][i-1];
        int currLabel = s0s1labelling[0][i];
        if((s0s1Components->containsSemiEdge[prevLabel] != s0s1Components->containsSemiEdge[currLabel]) ||
                (s0s1Components->componentSizes[prevLabel] != s0s1Components->componentSizes[currLabel])){
            s0s1componentsPartitioning[0][i] = 0;
            for(j = lastStart; j < i; j++){
                s0s1componentsPartitionSize[0][j] = currentSize;
            }
            currentSize = 1;
            lastStart = i;
        } else {
            s0s1componentsPartitioning[0][i] = i;
            currentSize++;
        }
    }
    for(j = lastStart; j < s0s1Components->componentCount; j++){
        s0s1componentsPartitionSize[0][j] = currentSize;
    }
    if(s0s1componentsPartitionSize[0][0]>1){
        //try to refine partition
        refineColourComponentPartition(&cddgraph, s0s1Components, s1s2Components, vertex2s0s1components, vertex2s1s2components,
                                        0, 0, s0s1componentsPartitioning, s0s1componentsPartitionSize,  s0s1labelling,
                                        s0s1Components, vertex2s0s1components);
    }
    
    assignLabelsToNexts0s1Component(ddgraph, &cddgraph, s0s1Components, s1s2Components, vertex2s0s1components, vertex2s1s2components, 0, 0);
}

void assignComponentLabels(DDGRAPH *ddgraph){
    COLOURCOMPONENTS s0s1Components;
    s0s1Components.colour1 = 0;
    s0s1Components.colour2 = 1;
    int vertex2s0s1Component[MAXN]; //maps vertices to components
    COLOURCOMPONENTS s1s2Components;
    s1s2Components.colour1 = 1;
    s1s2Components.colour2 = 2;
    int vertex2s1s2Component[MAXN]; //maps vertices to components

    findComponents(ddgraph, &s0s1Components, vertex2s0s1Component);
    findComponents(ddgraph, &s1s2Components, vertex2s1s2Component);
    
    if(s0s1Components.componentCount < minFaceOrbitCount || s0s1Components.componentCount > maxFaceOrbitCount){
        PROFILINGINCREMENT(rejectedColouredGraphBecauseWrongNumberFaceOrbits)
        return;
    }
    
    if(s1s2Components.componentCount < minVertexOrbitCount || s1s2Components.componentCount > maxVertexOrbitCount){
        PROFILINGINCREMENT(rejectedColouredGraphBecauseWrongNumberVertexOrbits)
        return;
    }
    
    int i;
    for(i=0; i<s0s1Components.componentCount; i++){
        if(s0s1Components.containsSemiEdge[i] && s0s1Components.componentSizes[i] > maxFaceSize){
            PROFILINGINCREMENT(rejectedColouredGraphBecauseFaceOrbitTooBig)
            return;
        } else if(!s0s1Components.containsSemiEdge[i] && s0s1Components.componentSizes[i]/2 > maxFaceSize){
            PROFILINGINCREMENT(rejectedColouredGraphBecauseFaceOrbitTooBig)
            return;
        }
        //clear component labels
        s0s1Components.componentLabels[i]=0;
    }
    for(i=0; i<s1s2Components.componentCount; i++){
        if(s1s2Components.containsSemiEdge[i] && s1s2Components.componentSizes[i] > maxVertexDegree){
            PROFILINGINCREMENT(rejectedColouredGraphBecauseVertexOrbitTooBig)
            return;
        } else if(!s1s2Components.containsSemiEdge[i] && s1s2Components.componentSizes[i]/2 > maxVertexDegree){
            PROFILINGINCREMENT(rejectedColouredGraphBecauseVertexOrbitTooBig)
            return;
        }
        //clear component labels
        s1s2Components.componentLabels[i]=0;
    }
    
    if(minVertexDegree * (ddgraph->order - 4 * s0s1Components.componentCount) >
            2*ddgraph->order){
        PROFILINGINCREMENT(rejectedColouredGraphBecauseIncompatibleParameters)
        return;
    }
    
    if(minFaceSize * (ddgraph->order - 4 * s1s2Components.componentCount) >
            2*ddgraph->order){
        PROFILINGINCREMENT(rejectedColouredGraphBecauseIncompatibleParameters)
        return;
    }
    
    if(ddgraph->order > 
            4*(s0s1Components.componentCount + s1s2Components.componentCount)){
        PROFILINGINCREMENT(rejectedColouredGraphBecauseTooManyOrbits)
        return;
    }
    
    PROFILINGINCREMENT(acceptedColouredGraphs)
    
    //bruteForces0s1LabelAssignment(ddgraph, &s0s1Components, &s1s2Components, 0);
    startAssignings0s1ComponentLabels(ddgraph, &s0s1Components, &s1s2Components, vertex2s0s1Component, vertex2s1s2Component);
}

//========= PHASE 3: HANDLING THE GENERATED DELANEY-DRESS GRAPHS ============
boolean first = TRUE;

boolean applyColouredDelaneyDressGraphFilter(DDGRAPH *ddgraph){
    COLOURCOMPONENTS s0s1Components;
    s0s1Components.colour1 = 0;
    s0s1Components.colour2 = 1;
    int vertex2s0s1Component[MAXN]; //maps vertices to components
    COLOURCOMPONENTS s1s2Components;
    s1s2Components.colour1 = 1;
    s1s2Components.colour2 = 2;
    int vertex2s1s2Component[MAXN]; //maps vertices to components
    
    int i;

    findComponents(ddgraph, &s0s1Components, vertex2s0s1Component);
    findComponents(ddgraph, &s1s2Components, vertex2s1s2Component);
    
    if(s0s1Components.componentCount < minFaceOrbitCount || s0s1Components.componentCount > maxFaceOrbitCount){
        PROFILINGINCREMENT(rejectedColouredGraphBecauseWrongNumberFaceOrbits)
        return FALSE;
    }
    
    for(i = 0; i < s0s1Components.componentCount; i++){
        int currentR = s0s1Components.componentSizes[i];
        if(!s0s1Components.containsSemiEdge[i]){
            if(currentR%2){
                ERRORMSG("Illegal size of component")
            }
            currentR /= 2;
        }
        if(currentR < minFaceR || currentR > maxFaceR){
            //TODO: profiling
            return FALSE;
        }
    }
    
    if(s1s2Components.componentCount < minVertexOrbitCount || s1s2Components.componentCount > maxVertexOrbitCount){
        PROFILINGINCREMENT(rejectedColouredGraphBecauseWrongNumberVertexOrbits)
        return FALSE;
    }
    
    for(i = 0; i < s1s2Components.componentCount; i++){
        int currentR = s1s2Components.componentSizes[i];
        if(!s1s2Components.containsSemiEdge[i]){
            if(currentR%2){
                ERRORMSG("Illegal size of component")
            }
            currentR /= 2;
        }
        if(currentR < minVertexR || currentR > maxVertexR){
            //TODO: profiling
            return FALSE;
        }
    }
    
    return TRUE;
}

void handleColouredDelaneyDressGraph(DDGRAPH *ddgraph){
    edgeColouredGraphsCount++;
    if(symbols){
        #ifdef _PROFILING
        unsigned long long int oldSymbolCount = symbolsCount;    
        #endif
        assignComponentLabels(ddgraph);
        #ifdef _PROFILING
        if(oldSymbolCount==symbolsCount){
            colouredDelaneyDressGraphsWithoutSymbol++;
        } else if (oldSymbolCount==symbolsCount-1){
            colouredDelaneyDressGraphsWithOneSymbol++;
        } else {
            colouredDelaneyDressGraphsWithMultipleSymbols++;
        }
        #endif
    } else {
        if(filterDelaneyDressGraphs && !applyColouredDelaneyDressGraphFilter(ddgraph)){
            edgeColouredGraphsCount--;
            return;
        }
        if(outputType=='c'){
            writePregraphColorCodeEdgeColouring(stdout, ddgraph, first);
            first = FALSE;
        } else if(outputType=='h'){
            CDDGRAPH cddgraph;
            translateDdgraphToCddgraph(ddgraph, &cddgraph);
            printColouredDDGraph(&cddgraph);
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
    if(bipartite){
        if(ddgraph->vertexPartition[connectionsMade][connector1]
                == ddgraph->vertexPartition[connectionsMade][connector2] &&
                ddgraph->vertexComponents[connectionsMade][connector1] ==
                ddgraph->vertexComponents[connectionsMade][connector2]){
            return FALSE;
        }
    }
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

void calculateAutomorphisms(DDGRAPH* ddgraph, int lastClosedGraphDepth){
    memcpy(nautyPtn, ddgraph->partitioning[lastClosedGraphDepth], sizeof(int)*(ddgraph->underlyingGraph->nv));
    memcpy(nautyLabelling, ddgraph->labelling[lastClosedGraphDepth], sizeof(int)*(ddgraph->underlyingGraph->nv));

#ifdef _DEBUG
    int i, j;
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

    nauty((graph*)(ddgraph->underlyingGraph), nautyLabelling, nautyPtn, NULL, nautyOrbits,
            &nautyOptions, &nautyStats, nautyWorkspace, 50 * MAXM, MAXM,
            ddgraph->underlyingGraph->nv, (graph*)&canonGraph);

#ifdef _DEBUG
    fprintf(stderr, "nauty Orbits: [%d", nautyOrbits[0]);
    for(i=1; i<ddgraph->underlyingGraph->nv; i++){
        fprintf(stderr, ", %d", nautyOrbits[i]);
    }
    fprintf(stderr, "]\n");
#endif
}

boolean isCanonicalConnection(BBLOCK* blocks, int buildingBlockCount, DDGRAPH *ddgraph,
        int *vertexToBlock, int *vertexToConnector, int orbit, int depth, int connector1, int connector2,
        boolean needSymmetryGroup){
    int i;

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
    if(madeConnectionsCount==1){
        PROFILINGINCREMENT(acceptedBecauseOnlyOne)
        //guaranteed to be canonical
        if(needSymmetryGroup){
            calculateAutomorphisms(ddgraph, depth);
        } else {
            PROFILINGINCREMENT(skippedNautyBecauseOnlyOne)
        }
        return TRUE;
    }
            
    //calculate colour for the connection and try to reject based on this colour
    int colour1 = ddgraph->vertexColours[depth][connector1];
    int colour2 = ddgraph->vertexColours[depth][connector2];
    int smallestColour = (colour1 < colour2 ? colour1 : colour2);
    int smallestCount = 0;

    for(i=0; i<madeConnectionsCount; i++){
        int localColour1 = ddgraph->vertexColours[depth][madeConnections[i][0]];
        int localColour2 = ddgraph->vertexColours[depth][madeConnections[i][1]];

        int localSmallest = localColour1 < localColour2 ? localColour1 : localColour2;

        if(localSmallest < smallestColour){
            PROFILINGINCREMENT(rejectedByColour)
            return FALSE;
        } else if(localSmallest == smallestColour){
            smallestCount++;
        }
    }
    
    if(smallestCount==1){   
        PROFILINGINCREMENT(acceptedBecauseOnlyOneMinimalColour)
        //guaranteed to be canonical
        if(needSymmetryGroup){
            calculateAutomorphisms(ddgraph, depth);
        } else {
            PROFILINGINCREMENT(skippedNautyBecauseOnlyOneMinimalColour)
        }
        return TRUE;
    }
    
    //calculate automorphisms of the new graph

    memcpy(nautyPtn, ddgraph->partitioning[depth], sizeof(int)*(ddgraph->underlyingGraph->nv));
    memcpy(nautyLabelling, ddgraph->labelling[depth], sizeof(int)*(ddgraph->underlyingGraph->nv));

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

    nauty((graph*)(ddgraph->underlyingGraph), nautyLabelling, nautyPtn, NULL, nautyOrbits,
            &nautyOptions, &nautyStats, nautyWorkspace, 50 * MAXM, MAXM,
            ddgraph->underlyingGraph->nv, (graph*)&canonGraph);

#ifdef _DEBUG
    fprintf(stderr, "nauty Orbits: [%d", nautyOrbits[0]);
    for(i=1; i<ddgraph->underlyingGraph->nv; i++){
        fprintf(stderr, ", %d", nautyOrbits[i]);
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

#ifdef _PROFILING

    if(madeConnectionOrbits[smallestConnection] != madeConnectionOrbits[newConnection]){
        rejectedByNauty++;
    } else {
        connectionsAccepted++;
    }

#endif 

    return madeConnectionOrbits[smallestConnection] == madeConnectionOrbits[newConnection];
}

boolean hasTrivialSymmetryForFreeConnections(DDGRAPH *ddgraph, boolean *freeConnectors){
    if(numberOfGenerators[connectionsMade]==0)
        return TRUE;

    int i,j;
    for(i=0; i<ddgraph->order; i++){ //dummy vertices can never be connectors
        if(freeConnectors[i]){
            for(j=0; j<numberOfGenerators[connectionsMade]; j++){
                if(automorphismGroupGenerators[connectionsMade][j][i]!=i){
                    return FALSE;
                }
            }
        }
    }
    return TRUE;
}

void connectCompleteOrbit(BBLOCK* blocks, int buildingBlockCount, DDGRAPH *ddgraph,
        int *vertexToBlock, int *vertexToConnector, int orbit, int depth,
        int openConnectionsLeftInOrbit, int totalConnectionsLeft, boolean *freeConnectors){
    //depth corresponds to the depth at which the last closed graph was
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
#ifdef _DEBUGMETHODS
                connections[connectionsMade][0] = v1;
                connections[connectionsMade][1] = v2;
#endif
                connectionsMade++;
                numberOfGenerators[connectionsMade]=0;
                freeConnectors[v1]=FALSE;
                freeConnectors[v2]=FALSE;
                
                //check canonicity of operation and recurse
                if(isCanonicalConnection(blocks, buildingBlockCount, ddgraph,
                        vertexToBlock, vertexToConnector, orbit, depth, v1,
                        v2, (totalConnectionsLeft - 2 != 0) || colouredEdges)){
#ifdef _PROFILING
                    if(numberOfGenerators[connectionsMade]==0){
                        graphsWithTrivialSymmetry[connectionsMade]++;
                    } else {
                        graphsWithNonTrivialSymmetry[connectionsMade]++;
                    }
                    if(hasTrivialSymmetryForFreeConnections(ddgraph, freeConnectors))
                        graphsWithTrivialSymmetryForRemainingConnections[connectionsMade]++;
#endif
                    //copy the component and partition information
                    if(bipartite){
                        int smallestComponent = MIN(ddgraph->vertexComponents[connectionsMade-1][v1],
                                                    ddgraph->vertexComponents[connectionsMade-1][v2]);
                        int largestComponent = MAX(ddgraph->vertexComponents[connectionsMade-1][v1],
                                                    ddgraph->vertexComponents[connectionsMade-1][v2]);
                        boolean switchColours = (ddgraph->vertexPartition[connectionsMade-1][v1] ==
                                                 ddgraph->vertexPartition[connectionsMade-1][v2]);
                        if(switchColours){
                            for(j=0; j<ddgraph->order; j++){
                                if(freeConnectors[j]){
                                    if(ddgraph->vertexComponents[connectionsMade-1][j]==largestComponent){
                                        ddgraph->vertexComponents[connectionsMade][j] = smallestComponent;
                                        ddgraph->vertexPartition[connectionsMade][j] = 
                                                !(ddgraph->vertexPartition[connectionsMade-1][j]);
                                    } else {
                                        ddgraph->vertexComponents[connectionsMade][j] = 
                                                ddgraph->vertexComponents[connectionsMade-1][j];
                                        ddgraph->vertexPartition[connectionsMade][j] = 
                                                ddgraph->vertexPartition[connectionsMade-1][j];
                                    }
                                }
                            }
                        } else {
                            for(j=0; j<ddgraph->order; j++){
                                if(freeConnectors[j]){
                                    if(ddgraph->vertexComponents[connectionsMade-1][j]==largestComponent){
                                        ddgraph->vertexComponents[connectionsMade][j] = smallestComponent;
                                    } else {
                                        ddgraph->vertexComponents[connectionsMade][j] = 
                                                ddgraph->vertexComponents[connectionsMade-1][j];
                                    }
                                    ddgraph->vertexPartition[connectionsMade][j] = ddgraph->vertexPartition[connectionsMade-1][j];
                                }
                            }
                        }
                    }
                    
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

void finishWithoutCanonicityCheck(BBLOCK* blocks, int buildingBlockCount, DDGRAPH *ddgraph,
        int *vertexToBlock, int *vertexToConnector, int totalConnectionsLeft, boolean *freeConnectors,
        boolean hadTrivialSymmetry, int lastClosedDepth){
    int i, j;

    i = 0;
    while(i<ddgraph->order-1 && !freeConnectors[i]) i++;
    DEBUGASSERT(i<ddgraph->order-1)

    for(j=i+1; j<ddgraph->order; j++){
        if(freeConnectors[j] && vertexToBlock[i]!=vertexToBlock[j]){
            if(isLegalConnection(blocks, buildingBlockCount, ddgraph, vertexToBlock, vertexToConnector, i, j)){

                //connect v1 to v2
                connectConnectors(blocks, ddgraph, vertexToBlock, vertexToConnector, i, j);
#ifdef _DEBUGMETHODS
                connections[connectionsMade][0] = i;
                connections[connectionsMade][1] = j;
#endif
                connectionsMade++;
                numberOfGenerators[connectionsMade]=0;
                freeConnectors[i]=FALSE;
                freeConnectors[j]=FALSE;

                PROFILINGINCREMENT(graphsFromClosedGraphsWithTrivialSymmetry[connectionsMade])

                //copy the component and partition information
                if(bipartite){
                    int k;
                    int smallestComponent = MIN(ddgraph->vertexComponents[connectionsMade-1][i],
                                                ddgraph->vertexComponents[connectionsMade-1][j]);
                    int largestComponent = MAX(ddgraph->vertexComponents[connectionsMade-1][i],
                                                ddgraph->vertexComponents[connectionsMade-1][j]);
                    boolean switchColours = (ddgraph->vertexPartition[connectionsMade-1][i] ==
                                             ddgraph->vertexPartition[connectionsMade-1][j]);
                    if(switchColours){
                        for(k=0; k<ddgraph->order; k++){
                            if(freeConnectors[k]){
                                if(ddgraph->vertexComponents[connectionsMade-1][k]==largestComponent){
                                    ddgraph->vertexComponents[connectionsMade][k] = smallestComponent;
                                    ddgraph->vertexPartition[connectionsMade][k] = 
                                            !(ddgraph->vertexPartition[connectionsMade-1][k]);
                                } else {
                                    ddgraph->vertexComponents[connectionsMade][k] = 
                                            ddgraph->vertexComponents[connectionsMade-1][k];
                                    ddgraph->vertexPartition[connectionsMade][k] = 
                                            ddgraph->vertexPartition[connectionsMade-1][k];
                                }
                            }
                        }
                    } else {
                        for(k=0; k<ddgraph->order; k++){
                            if(freeConnectors[k]){
                                if(ddgraph->vertexComponents[connectionsMade-1][k]==largestComponent){
                                    ddgraph->vertexComponents[connectionsMade][k] = smallestComponent;
                                } else {
                                    ddgraph->vertexComponents[connectionsMade][k] = 
                                            ddgraph->vertexComponents[connectionsMade-1][k];
                                }
                                ddgraph->vertexPartition[connectionsMade][k] = ddgraph->vertexPartition[connectionsMade-1][k];
                            }
                        }
                    }
                }
                
                    //
                if(totalConnectionsLeft - 2 == 0){
                    //printConnectionsMade();
                    //printGenerators(ddgraph, 6);
                    //fprintf(stderr, "Graph %2d: ", graphsCount);
                    //printCanonicalLabelling(ddgraph);
                    //fprintf(stderr, "Found graph based on: ");
                    //printHumanReadableComponentList();
                    if(!hadTrivialSymmetry){
                        //calculate automorphism group
                        calculateAutomorphisms(ddgraph, lastClosedDepth);
                    }
                    handleDelaneyDressGraph(ddgraph);
                } else {
                    finishWithoutCanonicityCheck(blocks, buildingBlockCount, ddgraph,
                            vertexToBlock, vertexToConnector, totalConnectionsLeft-2,
                            freeConnectors, hadTrivialSymmetry, lastClosedDepth);
                }

                //disconnect i from j
                freeConnectors[i]=TRUE;
                freeConnectors[j]=TRUE;
                connectionsMade--;
                disconnectConnectors(blocks, ddgraph, vertexToBlock, vertexToConnector, i, j);
            }
        }
    }
}

void findNextOrbitToConnect(BBLOCK* blocks, int buildingBlockCount, DDGRAPH *ddgraph, int *vertexToBlock, int *vertexToConnector, boolean *freeConnectors){
    DEBUGTRACE_ENTER
#ifdef _PROFILING
    if(numberOfGenerators[connectionsMade]==0){
        closedGraphsWithTrivialSymmetry[connectionsMade]++;
    } else {
        closedGraphsWithNonTrivialSymmetry[connectionsMade]++;
    }
    if(hasTrivialSymmetryForFreeConnections(ddgraph, freeConnectors))
        closedGraphsWithTrivialSymmetryForRemainingConnections[connectionsMade]++;
#endif
    if(hasTrivialSymmetryForFreeConnections(ddgraph, freeConnectors)){
        int i, count = 0;
        for(i=0; i<ddgraph->order; i++){
            if(freeConnectors[i]) count++;
        }
        if(numberOfGenerators[connectionsMade]!=0){
            int j;
            //store the vertex orbits as colours
            for(i=0; i<ddgraph->underlyingGraph->nv; i++){
                ddgraph->vertexColours[connectionsMade][i] = nautyOrbits[i];
            }
            //create the partioning and labeling that will be given to nauty while connecting this orbit
            int counter = 0;
            for(i=0; i<ddgraph->underlyingGraph->nv+14; i++){
                for(j=0; j<ddgraph->underlyingGraph->nv; j++){
                    if(ddgraph->vertexColours[connectionsMade][j]==i){
                        ddgraph->labelling[connectionsMade][counter] = j;
                        ddgraph->partitioning[connectionsMade][counter] = 1;
                        counter++;
                    }
                }
                if(counter>0){
                    ddgraph->partitioning[connectionsMade][counter-1]=0;
                }
                if(counter==ddgraph->underlyingGraph->nv) break;
            }
        }
        finishWithoutCanonicityCheck(blocks, buildingBlockCount, ddgraph, vertexToBlock,
                vertexToConnector, count, freeConnectors, numberOfGenerators[connectionsMade]==0,
                connectionsMade);
    } else {
        int i,j;
        //store the vertex orbits as colours
        for(i=0; i<ddgraph->underlyingGraph->nv; i++){
            ddgraph->vertexColours[connectionsMade][i] = nautyOrbits[i];
        }
        //create the partioning and labeling that will be given to nauty while connecting this orbit
        int counter = 0;
        for(i=0; i<ddgraph->underlyingGraph->nv+14; i++){
            for(j=0; j<ddgraph->underlyingGraph->nv; j++){
                if(ddgraph->vertexColours[connectionsMade][j]==i){
                    ddgraph->labelling[connectionsMade][counter] = j;
                    ddgraph->partitioning[connectionsMade][counter] = 1;
                    counter++;
                }
            }
            if(counter>0){
                ddgraph->partitioning[connectionsMade][counter-1]=0;
            }
            if(counter==ddgraph->underlyingGraph->nv) break;
        }
        
        //first we need the vertex orbits

        int orbitCount = 0;
    
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
    }
    DEBUGTRACE_EXIT
}

void connectComponentList(int vertexCount, DDGRAPH *ddgraph){
    DEBUGTRACE_ENTER
    
    //check to see if we cannot discard this list based upon the constraints
    if(symbols){
        int i;
        int maxPearlChainParameter = MIN(maxFaceSize-1, maxVertexDegree-1);
        for(i=maximumQ2TypeComponents; i>maxPearlChainParameter; i--){
            if(Q2TypeComponentsComponentCount[0][i-1]>0){
                //check for pearl chains that are too long
                PROFILINGINCREMENT(rejectedListsBecausePearlChainTooLong)
                DEBUGTRACE_EXIT
                return;
            }
        }
        int maxLockedPearlChainParameter = MIN((maxFaceSize-1)/2, (maxVertexDegree-1)/2);
        for(i=maximumQ2TypeComponents; i>maxLockedPearlChainParameter; i--){
            if(Q2TypeComponentsComponentCount[1][i-1]>0){
                //check for pearl chains that are too long
                PROFILINGINCREMENT(rejectedListsBecauseLockedPearlChainTooLong)
                DEBUGTRACE_EXIT
                return;
            }
        }
        
        /* find lower bound to number of face orbits and vertex orbits any graph
         * based on this list may yield. This lower bound is given by
         *  _                                              _
         * | #(semi-edges with colour 1) + #(q4-components) |
         * | ---------------------------------------------- |
         * |                       2                        |
         * 
         * We don't take q3-components (accept for the locked barb wire) into
         * account here, because these semi-edges can close either a s0s1 or a
         * s1s2 orbit depending on the colour they receive and this is not
         * decided until the colouring phase. 
         */
        
        int minEndingSemiEdgesCount = 0, j;
        for(i=0; i<Q1TypeComponentsCount; i++){
            if(Q1TypeComponentsSemiEdgeCount[i]){
                for(j=0; j<maximumQ1TypeComponents; j++){
                    minEndingSemiEdgesCount += Q1TypeComponentsComponentCount[i][j]*Q1TypeComponentsSemiEdgeCount[i];
                }
            }
        }
        for(i=0; i<maximumQ2TypeComponents; i++){
            minEndingSemiEdgesCount += Q2TypeComponentsComponentCount[1][i];
        }
        for(i=0; i<maximumQ3TypeComponents; i++){
            minEndingSemiEdgesCount += Q3TypeComponentsComponentCount[1][i];
        }
        int s1Count = minEndingSemiEdgesCount;
        minEndingSemiEdgesCount += Q4ComponentCount;
        
        if(minEndingSemiEdgesCount%2){
            minEndingSemiEdgesCount += 1;
            //plus 1 because the ceil of n/2 is equal to (n + 1)/2 for n odd
        }
        
        if(minEndingSemiEdgesCount/2 > MIN(maxFaceOrbitCount, maxVertexOrbitCount)){
            PROFILINGINCREMENT(rejectedListsBecauseTooManySemiEdgesForVertexOrFaceOrbitCount)
            DEBUGTRACE_EXIT
            return;
        }
        
        int q3TwoFactorCount = 0;
        for(i=0; i<maximumQ3TypeComponents; i++){
            q3TwoFactorCount += Q3TypeComponentsComponentCount[0][i] + 
                                Q3TypeComponentsComponentCount[1][i];
        }
        if(minEndingSemiEdgesCount + q3TwoFactorCount > maxFaceOrbitCount + maxVertexOrbitCount){
            PROFILINGINCREMENT(rejectedListsBecauseTooManySemiEdgesForCombinedOrbitCount)
            DEBUGTRACE_EXIT
            return;
        }
        
        //store improved bounds for this list of components
        listLevelMinFaceOrbitCount = MAX(minFaceOrbitCount, minEndingSemiEdgesCount);
        listLevelMinVertexOrbitCount = MAX(minVertexOrbitCount, minEndingSemiEdgesCount);
        listLevelMaxFaceOrbitCount = maxFaceOrbitCount;
        listLevelMaxVertexOrbitCount = MIN(maxVertexOrbitCount, (vertexCount+s1Count)/4);
        listLevelMaxFaceSize = maxFaceSize;
        listLevelMaxVertexDegree = maxVertexDegree;
        listLevelMinFaceSize = minFaceSize;
        listLevelMinVertexDegree = minVertexDegree;

        PROFILINGINCREMENT(numberOfListsAcceptedForSymbols)
    }
    
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

#ifdef _PROFILING
    if(numberOfGenerators[connectionsMade]==0){
        graphsWithTrivialSymmetry[connectionsMade]++;
    } else {
        graphsWithNonTrivialSymmetry[connectionsMade]++;
    }
#endif
    
    //store colours in nautyOrbits because findNextOrbitToConnect uses this to colour the vertices
    for(i=0; i<vertexCount; i++){
        nautyOrbits[i]=ddgraph->vertex2FactorType[i]*3 - ddgraph->semiEdges[i] - 1;
        DEBUGASSERT(nautyOrbits[i]>=0 && nautyOrbits[i]<12)
    }
    for(i=vertexCount; i<ddgraph->underlyingGraph->nv; i++){
        nautyOrbits[i]=12 + ddgraph->vertex2FactorType[ddgraph->underlyingGraph->e[ddgraph->underlyingGraph->v[i]]] - 1;
        DEBUGASSERT(nautyOrbits[i]>=12 && nautyOrbits[i]<14)
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
                if(intermediateStructures){
                    if(outputType=='c'){
                        writeListCodeMultipleBlockList(stdout, vertexCount);
                    } else if(outputType=='h'){
                        printHumanReadableComponentList(stdout);
                    }
                }
//                unsigned long long int oldSymbolsCount = symbolsCount;
//                unsigned long long int oldGraphsCount = graphsCount;
                connectComponentList(vertexCount, ddgraph);
//                if(oldSymbolsCount==symbolsCount){
//                    if(oldGraphsCount==graphsCount){
//                        fprintf(stdout, "nograph: ");
//                        printHumanReadableComponentList(stdout);
//                    } else {
//                        fprintf(stdout, "unused: ");
//                        printHumanReadableComponentList(stdout);
//                    }
//                } else {
//                    fprintf(stdout, "used: ");
//                    printHumanReadableComponentList(stdout);
//                }
            }
        }
        splitPointCount++;
    }
    DEBUGTRACE_EXIT
}

void q4Components(int targetSize, int currentSize, DDGRAPH *ddgraph){
    if(edgeOrbitCount + targetSize - currentSize > maxEdgeOrbitCount) return;
    if(edgeOrbitCount + targetSize - currentSize < minEdgeOrbitCount) return;
    
    if(!orientable || targetSize == currentSize){
        //Q4 component is not orientable
        edgeOrbitCount += targetSize - currentSize;
        Q4ComponentCount = targetSize - currentSize; //each q4 component has 1 vertex
        edgeOrbitCount -= targetSize - currentSize;
        handleComponentList(targetSize, ddgraph);
    }
}

void q3Components(int currentType, int currentParameter, int targetSize, int currentSize, DDGRAPH *ddgraph){
    if(edgeOrbitCount > maxEdgeOrbitCount) return;
    
    boolean skip = FALSE;
    int blockTypeNumber = buildingBlockParametersToNumber(3, currentType);
    if(bipartite && !isBipartiteBlock[blockTypeNumber]) skip = TRUE;
    if(orientable && !isOrientableBlock[blockTypeNumber]) skip = TRUE;
    if(skip){
        if(currentType+1==Q3TypeComponentsCount){
            q4Components(targetSize, currentSize, ddgraph);
        } else {
            q3Components(currentType+1, Q3TypeComponentsSmallestCase[currentType+1], targetSize, currentSize, ddgraph);
        }
        return;
    }
    
    int i;

    int remainingVertices = targetSize - currentSize;

    for(i = 0; i <= remainingVertices/(2*currentParameter); i++){
        Q3TypeComponentsComponentCount[currentType][currentParameter-1]=i;
        int newSize = currentSize + 2*currentParameter*i; //each q2 component has 2n vertices

        if(targetSize - newSize >= 2*(currentParameter+1)){
            edgeOrbitCount += currentParameter*i;
            q3Components(currentType, currentParameter+1, targetSize, newSize, ddgraph);
            edgeOrbitCount -= currentParameter*i;
        } else if(currentType+1==Q3TypeComponentsCount){
            edgeOrbitCount += currentParameter*i;
            q4Components(targetSize, newSize, ddgraph);
            edgeOrbitCount -= currentParameter*i;
        } else {
            edgeOrbitCount += currentParameter*i;
            q3Components(currentType+1, Q3TypeComponentsSmallestCase[currentType+1], targetSize, newSize, ddgraph);
            edgeOrbitCount -= currentParameter*i;
        }
    }
    Q3TypeComponentsComponentCount[currentType][currentParameter-1]=0; //reset this type to 0
}

void q2Components(int currentType, int currentParameter, int targetSize, int currentSize, DDGRAPH *ddgraph){
    if(edgeOrbitCount > maxEdgeOrbitCount) return;
    
    boolean skip = FALSE;
    int blockTypeNumber = buildingBlockParametersToNumber(2, currentType);
    if(bipartite && !isBipartiteBlock[blockTypeNumber]) skip = TRUE;
    if(orientable && !isOrientableBlock[blockTypeNumber]) skip = TRUE;
    if(skip){
        if(currentType+1==Q2TypeComponentsCount){
            q3Components(0, Q3TypeComponentsSmallestCase[0], targetSize, currentSize, ddgraph);
        } else {
            q2Components(currentType+1, Q2TypeComponentsSmallestCase[currentType+1], targetSize, currentSize, ddgraph);
        }
        return;
    }
    
    int i;

    int remainingVertices = targetSize - currentSize;

    for(i = 0; i <= remainingVertices/(2*currentParameter); i++){
        Q2TypeComponentsComponentCount[currentType][currentParameter-1]=i;
        int newSize = currentSize + 2*currentParameter*i; //each q2 component has 2n vertices

        if(targetSize - newSize >= 2*(currentParameter+1)){
            edgeOrbitCount += currentParameter*i;
            q2Components(currentType, currentParameter+1, targetSize, newSize, ddgraph);
            edgeOrbitCount -= currentParameter*i;
        } else if(currentType+1==Q2TypeComponentsCount){
            edgeOrbitCount += currentParameter*i;
            q3Components(0, Q3TypeComponentsSmallestCase[0], targetSize, newSize, ddgraph);
            edgeOrbitCount -= currentParameter*i;
        } else {
            edgeOrbitCount += currentParameter*i;
            q2Components(currentType+1, Q2TypeComponentsSmallestCase[currentType+1], targetSize, newSize, ddgraph);
            edgeOrbitCount -= currentParameter*i;
        }
    }
    Q2TypeComponentsComponentCount[currentType][currentParameter-1]=0; //reset this type to 0
}

void q1Components(int currentType, int currentParameter, int targetSize, int currentSize, DDGRAPH *ddgraph){
    if(edgeOrbitCount > maxEdgeOrbitCount) return;
    
    boolean skip = FALSE;
    int blockTypeNumber = buildingBlockParametersToNumber(1, currentType);
    if(bipartite && !isBipartiteBlock[blockTypeNumber]) skip = TRUE;
    if(orientable && !isOrientableBlock[blockTypeNumber]) skip = TRUE;
    if(skip){
        if(currentType+1==Q1TypeComponentsCount){
            q2Components(0, Q2TypeComponentsSmallestCase[0], targetSize, currentSize, ddgraph);
        } else {
            q1Components(currentType+1, Q1TypeComponentsSmallestCase[currentType+1], targetSize, currentSize, ddgraph);
        }
        return;
    }
    
    int i;

    int remainingVertices = targetSize - currentSize;

    for(i = 0; i <= remainingVertices/(4*currentParameter); i++){
        Q1TypeComponentsComponentCount[currentType][currentParameter-1]=i;
        int newSize = currentSize + 4*currentParameter*i; //each q1 component has 4n vertices

        if(targetSize - newSize >= 4*(currentParameter+1)){
            edgeOrbitCount += currentParameter*i;
            q1Components(currentType, currentParameter+1, targetSize, newSize, ddgraph);
            edgeOrbitCount -= currentParameter*i;
        } else if(currentType+1==Q1TypeComponentsCount){
            edgeOrbitCount += currentParameter*i;
            q2Components(0, Q2TypeComponentsSmallestCase[0], targetSize, newSize, ddgraph);
            edgeOrbitCount -= currentParameter*i;
        } else {
            edgeOrbitCount += currentParameter*i;
            q1Components(currentType+1, Q1TypeComponentsSmallestCase[currentType+1], targetSize, newSize, ddgraph);
            edgeOrbitCount -= currentParameter*i;
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
    
    Q1TypeComponentsSemiEdgeCount[0] = 0;
    Q1TypeComponentsSemiEdgeCount[1] = 1;
    Q1TypeComponentsSemiEdgeCount[2] = 1;
    Q1TypeComponentsSemiEdgeCount[3] = 0;
    Q1TypeComponentsSemiEdgeCount[4] = 0;
    Q1TypeComponentsSemiEdgeCount[5] = 2;
    Q1TypeComponentsSemiEdgeCount[6] = 0;
    Q1TypeComponentsSemiEdgeCount[7] = 2;
    Q1TypeComponentsSemiEdgeCount[8] = 1;
    Q1TypeComponentsSemiEdgeCount[9] = 1;
    Q1TypeComponentsSemiEdgeCount[10] = 3;
    Q1TypeComponentsSemiEdgeCount[11] = 3;

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
    
    isBipartiteBlock[0] = TRUE;
    isBipartiteBlock[1] = TRUE;
    isBipartiteBlock[2] = TRUE;
    isBipartiteBlock[3] = FALSE;
    isBipartiteBlock[4] = TRUE;
    isBipartiteBlock[5] = TRUE;
    isBipartiteBlock[6] = TRUE;
    isBipartiteBlock[7] = TRUE;
    isBipartiteBlock[8] = FALSE;
    isBipartiteBlock[9] = TRUE;
    isBipartiteBlock[10] = TRUE;
    isBipartiteBlock[11] = TRUE;
    isBipartiteBlock[12] = TRUE;
    isBipartiteBlock[13] = TRUE;
    isBipartiteBlock[14] = TRUE;
    isBipartiteBlock[15] = TRUE;
    isBipartiteBlock[16] = TRUE;
    isBipartiteBlock[17] = TRUE;
    isBipartiteBlock[18] = TRUE;
    isBipartiteBlock[19] = TRUE;
    isBipartiteBlock[20] = TRUE;
    isBipartiteBlock[21] = TRUE;
    isBipartiteBlock[22] = FALSE;
    isBipartiteBlock[23] = FALSE;
    isBipartiteBlock[24] = TRUE;
    isBipartiteBlock[25] = TRUE;
    isBipartiteBlock[26] = TRUE;
    isBipartiteBlock[27] = TRUE;
    isBipartiteBlock[28] = TRUE;
    
    isOrientableBlock[0] = TRUE;
    isOrientableBlock[1] = FALSE;
    isOrientableBlock[2] = FALSE;
    isOrientableBlock[3] = FALSE;
    isOrientableBlock[4] = TRUE;
    isOrientableBlock[5] = FALSE;
    isOrientableBlock[6] = TRUE;
    isOrientableBlock[7] = FALSE;
    isOrientableBlock[8] = FALSE;
    isOrientableBlock[9] = FALSE;
    isOrientableBlock[10] = FALSE;
    isOrientableBlock[11] = FALSE;
    isOrientableBlock[12] = TRUE;
    isOrientableBlock[13] = FALSE;
    isOrientableBlock[14] = FALSE;
    isOrientableBlock[15] = FALSE;
    isOrientableBlock[16] = FALSE;
    isOrientableBlock[17] = FALSE;
    isOrientableBlock[18] = FALSE;
    isOrientableBlock[19] = TRUE;
    isOrientableBlock[20] = FALSE;
    isOrientableBlock[21] = FALSE;
    isOrientableBlock[22] = FALSE;
    isOrientableBlock[23] = FALSE;
    isOrientableBlock[24] = TRUE;
    isOrientableBlock[25] = FALSE;
    isOrientableBlock[26] = FALSE;
    isOrientableBlock[27] = TRUE;
    isOrientableBlock[28] = FALSE;

    isEdgeColourableBlock[0] = TRUE;
    isEdgeColourableBlock[1] = FALSE;
    isEdgeColourableBlock[2] = FALSE;
    isEdgeColourableBlock[3] = TRUE;
    isEdgeColourableBlock[4] = FALSE;
    isEdgeColourableBlock[5] = FALSE;
    isEdgeColourableBlock[6] = TRUE;
    isEdgeColourableBlock[7] = FALSE;
    isEdgeColourableBlock[8] = FALSE;
    isEdgeColourableBlock[9] = FALSE;
    isEdgeColourableBlock[10] = FALSE;
    isEdgeColourableBlock[11] = FALSE;
    isEdgeColourableBlock[12] = TRUE;
    isEdgeColourableBlock[13] = FALSE;
    isEdgeColourableBlock[14] = TRUE;
    isEdgeColourableBlock[15] = FALSE;
    isEdgeColourableBlock[16] = TRUE;
    isEdgeColourableBlock[17] = FALSE;
    isEdgeColourableBlock[18] = FALSE;
    isEdgeColourableBlock[19] = TRUE;
    isEdgeColourableBlock[20] = FALSE;
    isEdgeColourableBlock[21] = TRUE;
    isEdgeColourableBlock[22] = FALSE;
    isEdgeColourableBlock[23] = TRUE;
    isEdgeColourableBlock[24] = TRUE;
    isEdgeColourableBlock[25] = FALSE;
    isEdgeColourableBlock[26] = FALSE;
    isEdgeColourableBlock[27] = FALSE;
    isEdgeColourableBlock[28] = FALSE;
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
    symbolsCount = 0;
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
            if(intermediateStructures){
                if(outputType=='c'){
                    writeListCodeSingleBlockList(stdout, order, bblock);
                } else if(outputType=='h'){
                    printBlockName(stdout, 1, buildingBlockTypeToNumber(bblock), bblock->parameter, TRUE);
                    fprintf(stdout, "\n");
                }
            }
            int vertexToBlock[ddgraph->underlyingGraph->vlen];
            int vertexToConnector[ddgraph->underlyingGraph->vlen];
            constructBuildingBlockListAsGraph(bblock, 1, ddgraph, vertexToBlock, vertexToConnector);
            numberOfGenerators[0] = 0;
            storeInitialGenerators(bblock, 1, ddgraph);
#ifdef _PROFILING
    if(numberOfGenerators[connectionsMade]==0){
        graphsWithTrivialSymmetry[connectionsMade]++;
    } else {
        graphsWithNonTrivialSymmetry[connectionsMade]++;
    }
#endif
            handleDelaneyDressGraph(ddgraph);
            cleanDDGraph(ddgraph);
        }
    }
    splitPointCount++;
}

void addTristar(DDGRAPH * ddgraph){
    if(minEdgeOrbitCount>1){
        return;
    }
    edgeOrbitCount = 1;
    
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 5, 0, 0, 0);
    handleSingleBlockComponentList(bblock, 1, ddgraph);
    free(bblock);
}

void addDoubleLockedPearlChain(DDGRAPH * ddgraph, int parameter){
    if(minEdgeOrbitCount > parameter || maxEdgeOrbitCount < parameter){
        return;
    }
    edgeOrbitCount = parameter;
    
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 5, 1, parameter, 0);
    handleSingleBlockComponentList(bblock, parameter*2, ddgraph);
    free(bblock);
}

void addPearlNecklace(DDGRAPH * ddgraph, int parameter){
    if(minEdgeOrbitCount > parameter || maxEdgeOrbitCount < parameter){
        return;
    }
    edgeOrbitCount = parameter;
    
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 5, 2, parameter, 0);
    handleSingleBlockComponentList(bblock, parameter*2, ddgraph);
    free(bblock);
}

void addDoubleLockedBarbedWire(DDGRAPH * ddgraph, int parameter){
    if(minEdgeOrbitCount > parameter || maxEdgeOrbitCount < parameter){
        return;
    }
    edgeOrbitCount = parameter;
    
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 6, 0, parameter, 0);
    handleSingleBlockComponentList(bblock, parameter*2, ddgraph);
    free(bblock);
}

void addBarbedWireNecklace(DDGRAPH * ddgraph, int parameter){
    if(minEdgeOrbitCount > parameter || maxEdgeOrbitCount < parameter){
        return;
    }
    edgeOrbitCount = parameter;
    
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 6, 1, parameter, 0);
    handleSingleBlockComponentList(bblock, parameter*2, ddgraph);
    free(bblock);
}

void addDoubleLockedDiagonalChain(DDGRAPH * ddgraph, int parameter){
    if(minEdgeOrbitCount > parameter || maxEdgeOrbitCount < parameter){
        return;
    }
    edgeOrbitCount = parameter;
    
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 6, 2, parameter, 0);
    handleSingleBlockComponentList(bblock, parameter*4, ddgraph);
    free(bblock);
}

void addMobiusLadder(DDGRAPH * ddgraph, int parameter){
    if(minEdgeOrbitCount > parameter || maxEdgeOrbitCount < parameter){
        return;
    }
    edgeOrbitCount = parameter;
    
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 6, 3, parameter, 0);
    handleSingleBlockComponentList(bblock, parameter*4, ddgraph);
    free(bblock);
}

void addPrism(DDGRAPH * ddgraph, int parameter){
    if(minEdgeOrbitCount > parameter || maxEdgeOrbitCount < parameter){
        return;
    }
    edgeOrbitCount = parameter;
    
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 6, 4, parameter, 0);
    handleSingleBlockComponentList(bblock, parameter*4, ddgraph);
    free(bblock);
}

void addDoubleLockedDoubleroofLongBuilding(DDGRAPH * ddgraph, int parameter){
    if(minEdgeOrbitCount > parameter || maxEdgeOrbitCount < parameter){
        return;
    }
    edgeOrbitCount = parameter;
    
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 6, 5, parameter, 0);
    handleSingleBlockComponentList(bblock, parameter*4, ddgraph);
    free(bblock);
}

void addCompletelyLockedHub(DDGRAPH * ddgraph, int parameter){
    if(minEdgeOrbitCount > parameter || maxEdgeOrbitCount < parameter){
        return;
    }
    edgeOrbitCount = parameter;
    
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 6, 6, parameter, 0);
    handleSingleBlockComponentList(bblock, parameter*4, ddgraph);
    free(bblock);
}

void addDoubleroofDoublefloorHighBuilding(DDGRAPH * ddgraph, int parameter){
    if(minEdgeOrbitCount > parameter || maxEdgeOrbitCount < parameter){
        return;
    }
    edgeOrbitCount = parameter;
    
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 6, 7, parameter, 0);
    handleSingleBlockComponentList(bblock, parameter*4, ddgraph);
    free(bblock);
}

void addDoubleLockedDoubleroofHighBuilding(DDGRAPH * ddgraph, int parameter){
    if(minEdgeOrbitCount > parameter || maxEdgeOrbitCount < parameter){
        return;
    }
    edgeOrbitCount = parameter;
    
    BBLOCK * bblock = (BBLOCK *)malloc(sizeof(BBLOCK));
    initBuildingBlock(bblock, 6, 8, parameter, 0);
    handleSingleBlockComponentList(bblock, parameter*4, ddgraph);
    free(bblock);
}

void extraUnconstructableGraphs_1(DDGRAPH * ddgraph){
    if(!orientable){
        addTristar(ddgraph);
    }
}

void extraUnconstructableGraphs_2(DDGRAPH * ddgraph){
    if(!orientable){
        addDoubleLockedPearlChain(ddgraph, 1);
    }
    addPearlNecklace(ddgraph, 1);

    if(markedTwoFactors && !orientable){
        addDoubleLockedBarbedWire(ddgraph, 1);

        //isomorph to DLPC(1) in unmarked case
        addBarbedWireNecklace(ddgraph, 1);
    }
}

void extraUnconstructableGraphs_4(DDGRAPH * ddgraph){
    addPearlNecklace(ddgraph, 2);
    if(!orientable){
        addDoubleLockedPearlChain(ddgraph, 2);
        addBarbedWireNecklace(ddgraph, 2);
        if(!bipartite){
            addDoubleLockedDiagonalChain(ddgraph, 1);
            addMobiusLadder(ddgraph, 1);
        }
    }

    if(markedTwoFactors){
        if(!orientable){
            addDoubleLockedBarbedWire(ddgraph, 2);
            addCompletelyLockedHub(ddgraph, 1);
            addDoubleLockedDoubleroofHighBuilding(ddgraph, 1);
        }
        addDoubleroofDoublefloorHighBuilding(ddgraph, 1);
    }
}

void extraUnconstructableGraphs_4n(DDGRAPH * ddgraph, int targetSize){
    addPearlNecklace(ddgraph, targetSize/2);
    if(!orientable){
        addDoubleLockedPearlChain(ddgraph, targetSize/2);
        addBarbedWireNecklace(ddgraph, targetSize/2);
        addDoubleLockedDoubleroofLongBuilding(ddgraph, targetSize/4);
        if(!bipartite){
            addDoubleLockedDiagonalChain(ddgraph, targetSize/4);
            addMobiusLadder(ddgraph, targetSize/4);
        }
    }
    addPrism(ddgraph, targetSize/4);

    if(markedTwoFactors){
        if(!orientable){
            addDoubleLockedBarbedWire(ddgraph, targetSize/2);
            addCompletelyLockedHub(ddgraph, targetSize/4);
            addDoubleLockedDoubleroofHighBuilding(ddgraph, targetSize/4);
        }
        addDoubleroofDoublefloorHighBuilding(ddgraph, targetSize/4);
    }
}

void extraUnconstructableGraphs_4n2(DDGRAPH * ddgraph, int targetSize){
    addPearlNecklace(ddgraph, targetSize/2);
    if(!orientable){
        addDoubleLockedPearlChain(ddgraph, targetSize/2);
        addBarbedWireNecklace(ddgraph, targetSize/2);
    }
        
    if(markedTwoFactors && !orientable){
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
    if(orientable && targetSize%2){
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
        return;
    }

    initComponentsStatic();
    initComponents(targetSize);
    initStatistics();
    initNautyOptions(targetSize);

    DDGRAPH * ddgraph = getNewDDGraph(targetSize);

    q1Components(0, Q1TypeComponentsSmallestCase[0], targetSize, 0, ddgraph);

    extraUnconstructableGraphs(ddgraph, targetSize);

    freeDDGraph(ddgraph);

    if(!symbols || verbose){
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
    }
    if(symbols){
        fprintf(stderr, "Found %llu Delaney-Dress symbol%s.\n",
                    symbolsCount,
                    symbolsCount==1 ? (char *)"" : (char *)"s");
    }
    if(moduloEnabled){
        fprintf(stderr, "Generated only part %llu of %llu.\n", moduloRest+1, moduloMod);
    }

    freeComponents();
    cleanNautyOptions();
    cleanComponentStatistics();
}

void startFromListFile(char *filename, int minimumVertexCountAllowed, int maximumVertexCountAllowed){
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
        
        //skip lists with an incorrect number of vertices
        if(vertexCount < minimumVertexCountAllowed || vertexCount > maximumVertexCountAllowed) continue;
        
        int realVertexCount = 0;

        initComponents(vertexCount);
        initNautyOptions(vertexCount);

        int type = 0;
        int family = 0;
        int parameter = 0;
        int count = 0;
        boolean singleComponentList = FALSE;
        boolean canBeBipartite = TRUE;
        boolean canBeOrientable = TRUE;
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
                int blockNumber = buildingBlockParametersToNumber(1, family);
                if(!isBipartiteBlock[blockNumber]) canBeBipartite = FALSE;
                if(!isOrientableBlock[blockNumber]) canBeOrientable = FALSE;
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
                int blockNumber = buildingBlockParametersToNumber(2, family);
                if(!isBipartiteBlock[blockNumber]) canBeBipartite = FALSE;
                if(!isOrientableBlock[blockNumber]) canBeOrientable = FALSE;
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
                int blockNumber = buildingBlockParametersToNumber(3, family);
                if(!isBipartiteBlock[blockNumber]) canBeBipartite = FALSE;
                if(!isOrientableBlock[blockNumber]) canBeOrientable = FALSE;
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
                    int blockNumber = buildingBlockParametersToNumber(4, 0);
                    if(!isBipartiteBlock[blockNumber]) canBeBipartite = FALSE;
                    if(!isOrientableBlock[blockNumber]) canBeOrientable = FALSE;
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
                int blockNumber = buildingBlockParametersToNumber(5, family);
                if(!isBipartiteBlock[blockNumber]) canBeBipartite = FALSE;
                if(!isOrientableBlock[blockNumber]) canBeOrientable = FALSE;
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
                int blockNumber = buildingBlockParametersToNumber(6, family);
                if(!isBipartiteBlock[blockNumber]) canBeBipartite = FALSE;
                if(!isOrientableBlock[blockNumber]) canBeOrientable = FALSE;
                singleComponentList = TRUE;
                break;
            } else {
                ERRORMSG("Error while parsing file: illegal type.")
            }
        }

        if(realVertexCount!=vertexCount){
            ERRORMSG("Error while parsing file: incorrect vertex count.")
        }
        
        //skip this list in case it conflicts with any constraints
        if(orientable && !canBeOrientable) continue;
        if(bipartite && !canBeBipartite) continue;

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
    if(!symbols || verbose){
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
    }
    if(symbols){
        fprintf(stderr, "Found %llu Delaney-Dress symbol%s.\n",
                    symbolsCount,
                    symbolsCount==1 ? (char *)"" : (char *)"s");
    }
    if(moduloEnabled){
        fprintf(stderr, "Generated only part %llu of %llu.\n", moduloRest+1, moduloMod);
    }

    cleanComponentStatistics();
}

void startMultipleGenerations(int startSize, int endSize){
    int targetSize;
    
    initComponentsStatic();
    initStatistics();
    
    int step = 1;
    if(orientable){
        if(startSize%2){
            startSize++;
        }
        step = 2;
    }
    
    for(targetSize = startSize; targetSize<=endSize; targetSize+=step){
        initComponents(targetSize);
        initNautyOptions(targetSize);

        DDGRAPH * ddgraph = getNewDDGraph(targetSize);
        
        edgeOrbitCount = 0;
        //this could be set to a non-zero value in one of the single block lists
        
        q1Components(0, Q1TypeComponentsSmallestCase[0], targetSize, 0, ddgraph);
        
        extraUnconstructableGraphs(ddgraph, targetSize);

        freeDDGraph(ddgraph);

        freeComponents();
        cleanNautyOptions();
    }
    
    if(!symbols || verbose){
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
    }
    if(symbols){
        fprintf(stderr, "Found %llu Delaney-Dress symbol%s.\n",
                    symbolsCount,
                    symbolsCount==1 ? (char *)"" : (char *)"s");
    }
    if(moduloEnabled){
        fprintf(stderr, "Generated only part %llu of %llu.\n", moduloRest+1, moduloMod);
    }
    
    cleanComponentStatistics();
}

//====================== USAGE =======================

void help(char *name){
    fprintf(stderr, "The program %s calculates Delaney-Dress symbols or related structures.\n\n", name);
    fprintf(stderr, "Usage\n=====\n");
    fprintf(stderr, " %s [options] n\n", name);
    fprintf(stderr, "       Generate Delaney-Dress graphs with n vertices.\n");
    fprintf(stderr, " %s -l file [options] \n", name);
    fprintf(stderr, "       Generate Delaney-Dress graphs based on the component lists in file.\n");
    fprintf(stderr, " %s -s [options]\n", name);
    fprintf(stderr, "       Generate Delaney-Dress symbols.\n");
    fprintf(stderr, "\nThis program can handle graphs up to %d vertices. Recompile if you need larger\n", MAXN);
    fprintf(stderr, "graphs.\n\n");
    fprintf(stderr, "Valid options\n=============\n");
    fprintf(stderr, "* Various options\n");
    fprintf(stderr, "    -h, --help\n");
    fprintf(stderr, "       Print this help and return.\n");
    fprintf(stderr, "    -l, --listfile file\n");
    fprintf(stderr, "       Only start the generation based on the component lists in file.\n");
    fprintf(stderr, "    -o, --output type\n");
    fprintf(stderr, "       Specifies the export format where type is one of\n");
    fprintf(stderr, "           c, code    code depends on the generated type\n");
    fprintf(stderr, "           h, human   human-readable output\n");
    fprintf(stderr, "           n, none    no output: only count (default)\n");
    fprintf(stderr, "    -m, --modulo r:n\n");
    fprintf(stderr, "       Split the generation in multiple parts. The generation is split\n");
    fprintf(stderr, "       into n parts and only part r is generated. The number n needs to\n");
    fprintf(stderr, "       be an integer larger than 0 and r should be a positive integer\n");
    fprintf(stderr, "       smaller than n.\n");
    fprintf(stderr, "    --restrictionsonly\n");
    fprintf(stderr, "       calculate the restrictions for the Delaney-Dress symbol, but don't\n");
    fprintf(stderr, "       generate any structures.\n");
    fprintf(stderr, "\n* Generated types\n");
    fprintf(stderr, "    -L, --lists\n");
    fprintf(stderr, "       Generate component lists.\n");
    fprintf(stderr, "    -t, --marked\n");
    fprintf(stderr, "       Generate Delaney-Dress graphs with a marked 2-factor, i.e. graphs that\n");
    fprintf(stderr, "       are the underlying graphs of Delaney-Dress symbols and in which the\n");
    fprintf(stderr, "       s0s2 orbits are marked.\n");
    fprintf(stderr, "    -c, --coloured\n");
    fprintf(stderr, "       Generate coloured Delaney-Dress graphs, i.e. graphs that are the\n");
    fprintf(stderr, "       underlying graphs of Delaney-Dress symbols and in which the edges\n");
    fprintf(stderr, "       are coloured with the colours 0, 1 or 2.\n");
    fprintf(stderr, "    -s, --symbols\n");
    fprintf(stderr, "       Generate Delaney-Dress symbols.\n");
    fprintf(stderr, "\n* Specify constraints\n");
    fprintf(stderr, "    -b, --bipartite\n");
    fprintf(stderr, "       Only generate Delaney-Dress graphs that are bipartite.\n");
    fprintf(stderr, "    -O, --orientable\n");
    fprintf(stderr, "       Only generate Delaney-Dress graphs that (may) correspond to orientable tilings.\n");
    fprintf(stderr, "    -R, --requiredface\n");
    fprintf(stderr, "       Add a face size to the list of required faces.\n");
    fprintf(stderr, "    -F, --forbiddenface\n");
    fprintf(stderr, "       Add a face size to the list of forbidden faces.\n");
    fprintf(stderr, "    -r, --requiredvertex\n");
    fprintf(stderr, "       Add a vertex degree to the list of required vertices.\n");
    fprintf(stderr, "    -f, --forbiddenvertex\n");
    fprintf(stderr, "       Add a vertex degree to the list of forbidden vertices.\n");
    fprintf(stderr, "    --maxfacecount\n");
    fprintf(stderr, "       Specify the maximum number of face orbits in the tiling.\n");
    fprintf(stderr, "    --minfacecount\n");
    fprintf(stderr, "       Specify the minimum number of face orbits in the tiling.\n");
    fprintf(stderr, "    --maxvertexcount\n");
    fprintf(stderr, "       Specify the maximum number of vertex orbits in the tiling.\n");
    fprintf(stderr, "    --minvertexcount\n");
    fprintf(stderr, "       Specify the minimum number of vertex orbits in the tiling.\n");
    fprintf(stderr, "    --minfacesize\n");
    fprintf(stderr, "       Specify the minimum size of a face in the tiling.\n");
    fprintf(stderr, "    --maxfacesize\n");
    fprintf(stderr, "       Specify the maximum size of a face in the tiling.\n");
    fprintf(stderr, "    --minvertexdegree\n");
    fprintf(stderr, "       Specify the minimum degree of a vertex in the tiling.\n");
    fprintf(stderr, "    --maxvertexdegree\n");
    fprintf(stderr, "       Specify the maximum degree of a vertex in the tiling.\n");
    fprintf(stderr, "    -n, --minvertices\n");
    fprintf(stderr, "       Specify the minimum number of vertices in the Delaney-Dress graph.\n");
    fprintf(stderr, "    -N, --maxvertices\n");
    fprintf(stderr, "       Specify the maximum number of vertices in the Delaney-Dress graph.\n");
}

void usage(char *name){
    fprintf(stderr, "Usage: %s [options] n\n", name);
    fprintf(stderr, "       %s [options]\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

boolean checkIntegerValue(const char *paramName, int value, int minimum, int maximum){
    if(value < minimum){
        fprintf(stderr, "Illegal value: %d. The parameter %s should be at least %d.\n",
                value, paramName, minimum);
        return FALSE;
    } else if(value > maximum){
        fprintf(stderr, "Illegal value: %d. The parameter %s should be at most %d.\n",
                value, paramName, maximum);
        return FALSE;
    } else {
        return TRUE;
    }
}

boolean constraintsChanged(){
    //returns true if the constraints have changed since the last time this function was called
    //this function deals with the global constraints and not with the constraints that
    //are being updated during the generation
    
    static int oldMaxFaceOrbitCount = 0;
    static int oldMinFaceOrbitCount = 0;
    static int oldMaxVertexOrbitCount = 0;
    static int oldMinVertexOrbitCount = 0;

    static int oldMinFaceSize = 0;
    static int oldMaxFaceSize = 0;
    static int oldMinVertexDegree = 0;
    static int oldMaxVertexDegree = 0;

    static int oldRequestedFaceSizesCount = 0;
    static int oldRequestedVertexDegreesCount = 0;

    static int oldForbiddenFaceSizesCount = 0;
    static int oldForbiddenVertexDegreesCount = 0;
    
    boolean constraintsChanged = FALSE;
    
    if(oldMaxFaceOrbitCount!=maxFaceOrbitCount){
        constraintsChanged = TRUE;
    }
    if(oldMinFaceOrbitCount!=minFaceOrbitCount){
        constraintsChanged = TRUE;
    }
    if(oldMaxVertexOrbitCount!=maxVertexOrbitCount){
        constraintsChanged = TRUE;
    }
    if(oldMinVertexOrbitCount!=minVertexOrbitCount){
        constraintsChanged = TRUE;
    }
    if(oldMaxFaceSize!=maxFaceSize){
        constraintsChanged = TRUE;
    }
    if(oldMinFaceSize!=minFaceSize){
        constraintsChanged = TRUE;
    }
    if(oldMaxVertexDegree!=maxVertexDegree){
        constraintsChanged = TRUE;
    }
    if(oldMinVertexDegree!=minVertexDegree){
        constraintsChanged = TRUE;
    }
    if(oldRequestedFaceSizesCount!=requestedFaceSizesCount){
        constraintsChanged = TRUE;
    }
    if(oldRequestedVertexDegreesCount!=requestedVertexDegreesCount){
        constraintsChanged = TRUE;
    }
    if(oldForbiddenFaceSizesCount!=forbiddenFaceSizesCount){
        constraintsChanged = TRUE;
    }
    if(oldForbiddenVertexDegreesCount!=forbiddenVertexDegreesCount){
        constraintsChanged = TRUE;
    }
    
    if(constraintsChanged){
        oldMaxFaceOrbitCount = maxFaceOrbitCount;
        oldMinFaceOrbitCount = minFaceOrbitCount;
        oldMaxVertexOrbitCount = maxVertexOrbitCount;
        oldMinVertexOrbitCount = minVertexOrbitCount;
        oldMinFaceSize = minFaceSize;
        oldMaxFaceSize = maxFaceSize;
        oldMinVertexDegree = minVertexDegree;
        oldMaxVertexDegree = maxVertexDegree;
        oldRequestedFaceSizesCount = requestedFaceSizesCount;
        oldRequestedVertexDegreesCount = requestedVertexDegreesCount;
        oldForbiddenFaceSizesCount = forbiddenFaceSizesCount;
        oldForbiddenVertexDegreesCount = forbiddenVertexDegreesCount;
    }
    
    return constraintsChanged;
}

void adjustSymbolConstraints(){
    minFaceOrbitCount = MAX(minFaceOrbitCount, requestedFaceSizesCount);
    minVertexOrbitCount = MAX(minVertexOrbitCount, requestedVertexDegreesCount);
}

boolean validateSymbolConstraints(){
    int i,j;
    boolean validConstraints = TRUE;
    
    //check that the minima are smaller than the respective maxima
    if(maxFaceOrbitCount < minFaceOrbitCount){
        fprintf(stderr, "Maximum number of face orbits needs to be at least the minimum number of face orbits.\n");
        validConstraints = FALSE;
    }
    if(maxVertexOrbitCount < minVertexOrbitCount){
        fprintf(stderr, "Maximum number of vertex orbits needs to be at least the minimum number of vertex orbits.\n");
        validConstraints = FALSE;
    }
    if(maxEdgeOrbitCount < minEdgeOrbitCount){
        fprintf(stderr, "Maximum number of edge orbits needs to be at least the minimum number of edge orbits.\n");
        validConstraints = FALSE;
    }
    if(maxFaceSize < minFaceSize){
        fprintf(stderr, "Maximum face size needs to be at least the minimum face size.\n");
        validConstraints = FALSE;
    }
    if(maxVertexDegree < minVertexDegree){
        fprintf(stderr, "Maximum vertex degree needs to be at least the minimum vertex degree.\n");
        validConstraints = FALSE;
    }
    if(maxVertexCount < minVertexCount){
        fprintf(stderr, "Maximum number of vertices needs to be at least the minimum number of vertices.\n");
        validConstraints = FALSE;
    }
    
    //check that no face/vertex is both required and forbidden
    for(i=0; i<requestedFaceSizesCount; i++){
        for(j=0; j<forbiddenFaceSizesCount; j++){
            if(requestedFaceSizes[i]==forbiddenFaceSizes[j]){
                fprintf(stderr, "Face size %d is both required and forbidden.\n", requestedFaceSizes[i]);
                validConstraints = FALSE;
            }
        }
    }
    for(i=0; i<requestedVertexDegreesCount; i++){
        for(j=0; j<forbiddenVertexDegreesCount; j++){
            if(requestedVertexDegrees[i]==forbiddenVertexDegrees[j]){
                fprintf(stderr, "Vertex degree %d is both required and forbidden.\n", requestedVertexDegrees[i]);
                validConstraints = FALSE;
            }
        }
    }
    
    //calculate the minimum size based upon the required faces and vertices
    int orientableFactor = 1;
    if(orientable) orientableFactor = 2;
    int verticesNeeded = 0;
    for(i=0; i<requestedFaceSizesCount; i++){
        if(requestedFaceSizes[i]%6==0){
            verticesNeeded += requestedFaceSizes[i]/6*orientableFactor;
            continue;
        }
        if(requestedFaceSizes[i]%4==0){
            verticesNeeded += requestedFaceSizes[i]/4*orientableFactor;
            continue;
        }
        if(requestedFaceSizes[i]%3==0){
            verticesNeeded += requestedFaceSizes[i]/3*orientableFactor;
            continue;
        }
        if(requestedFaceSizes[i]%2==0){
            verticesNeeded += requestedFaceSizes[i]/2*orientableFactor;
            continue;
        }
        verticesNeeded += requestedFaceSizes[i]*orientableFactor;
    }
    int minimumVerticesNeededPerFaceForRemainder = maxFaceSize;
    for(i=minFaceSize; i<=maxFaceSize; i++){
        if(i%6==0){
            if(minimumVerticesNeededPerFaceForRemainder > i/6*orientableFactor){
                minimumVerticesNeededPerFaceForRemainder = i/6*orientableFactor;
            }
            continue;
        }
        if(i%4==0){
            if(minimumVerticesNeededPerFaceForRemainder > i/4*orientableFactor){
                minimumVerticesNeededPerFaceForRemainder = i/4*orientableFactor;
            }
            continue;
        }
        if(i%3==0){
            if(minimumVerticesNeededPerFaceForRemainder > i/3*orientableFactor){
                minimumVerticesNeededPerFaceForRemainder = i/3*orientableFactor;
            }
            continue;
        }
        if(i%2==0){
            if(minimumVerticesNeededPerFaceForRemainder > i/2*orientableFactor){
                minimumVerticesNeededPerFaceForRemainder = i/2*orientableFactor;
            }
            continue;
        }
        if(minimumVerticesNeededPerFaceForRemainder > i*orientableFactor){
            minimumVerticesNeededPerFaceForRemainder = i*orientableFactor;
        }
    }
    verticesNeeded += minimumVerticesNeededPerFaceForRemainder*(minFaceOrbitCount-requestedFaceSizesCount);
    if(verticesNeeded > MAXN || verticesNeeded > maxVertexCount){
        fprintf(stderr, "For this list of required face sizes the Delaney-Dress graph\nneeds at least %d vertices.\n", verticesNeeded);
        fprintf(stderr, "This version can only handle %d vertices. Recompile the program\nto be able to handle larger graphs.\n", MAXN);
        validConstraints = FALSE;
    } else if (verticesNeeded > minVertexCount){
        minVertexCount = verticesNeeded;
    }
    
    verticesNeeded = 0;
    for(i=0; i<requestedVertexDegreesCount; i++){
        if(requestedVertexDegrees[i]%6==0){
            verticesNeeded += requestedVertexDegrees[i]/6*orientableFactor;
            continue;
        }
        if(requestedVertexDegrees[i]%4==0){
            verticesNeeded += requestedVertexDegrees[i]/4*orientableFactor;
            continue;
        }
        if(requestedVertexDegrees[i]%3==0){
            verticesNeeded += requestedVertexDegrees[i]/3*orientableFactor;
            continue;
        }
        if(requestedVertexDegrees[i]%2==0){
            verticesNeeded += requestedVertexDegrees[i]/2*orientableFactor;
            continue;
        }
        verticesNeeded += requestedVertexDegrees[i]*orientableFactor;
    }
    int minimumVerticesNeededPerVertexForRemainder = maxVertexDegree;
    for(i=minVertexDegree; i<=maxVertexDegree; i++){
        if(i%6==0){
            if(minimumVerticesNeededPerVertexForRemainder > i/6*orientableFactor){
                minimumVerticesNeededPerVertexForRemainder = i/6*orientableFactor;
            }
            continue;
        }
        if(i%4==0){
            if(minimumVerticesNeededPerVertexForRemainder > i/4*orientableFactor){
                minimumVerticesNeededPerVertexForRemainder = i/4*orientableFactor;
            }
            continue;
        }
        if(i%3==0){
            if(minimumVerticesNeededPerVertexForRemainder > i/3*orientableFactor){
                minimumVerticesNeededPerVertexForRemainder = i/3*orientableFactor;
            }
            continue;
        }
        if(i%2==0){
            if(minimumVerticesNeededPerVertexForRemainder > i/2*orientableFactor){
                minimumVerticesNeededPerVertexForRemainder = i/2*orientableFactor;
            }
            continue;
        }
        if(minimumVerticesNeededPerVertexForRemainder > i*orientableFactor){
            minimumVerticesNeededPerVertexForRemainder = i*orientableFactor;
        }
    }
    verticesNeeded += minimumVerticesNeededPerVertexForRemainder*(minVertexOrbitCount-requestedVertexDegreesCount);
    if(verticesNeeded > MAXN || verticesNeeded > maxVertexCount){
        fprintf(stderr, "For this list of required vertex degrees the Delaney-Dress graph\nneeds at least %d vertices.\n", verticesNeeded);
        fprintf(stderr, "This version can only handle %d vertices. Recompile the program\nto be able to handle larger graphs.\n", MAXN);
        validConstraints = FALSE;
    } else if (verticesNeeded > minVertexCount){
        minVertexCount = verticesNeeded;
    }
    
    //calculate the maximum size based upon the vertex degree, face size and orbit counts
    int verticesAllowed = 4*minVertexDegree*maxFaceOrbitCount/(minVertexDegree-2);
    if(verticesAllowed < minVertexCount){
        fprintf(stderr, "For this minimum vertex degree and maximum number of face orbits the Delaney-Dress graph\ncan have at most %d vertices.\n", verticesAllowed);
        validConstraints = FALSE;
    } else if (verticesAllowed < maxVertexCount){
        maxVertexCount = verticesAllowed;
    }
    verticesAllowed = 4*minFaceSize*maxVertexOrbitCount/(minFaceSize-2);
    if(verticesAllowed < minVertexCount){
        fprintf(stderr, "For this minimum face size and maximum number of vertex orbits the Delaney-Dress graph\ncan have at most %d vertices.\n", verticesAllowed);
        validConstraints = FALSE;
    } else if (verticesAllowed < maxVertexCount){
        maxVertexCount = verticesAllowed;
    }
    verticesAllowed = 4*(maxVertexOrbitCount + maxFaceOrbitCount);
    if(verticesAllowed < minVertexCount){
        fprintf(stderr, "For these maxima for the number of vertex orbits and number of face orbits the Delaney-Dress graph\ncan have at most %d vertices.\n", verticesAllowed);
        validConstraints = FALSE;
    } else if (verticesAllowed < maxVertexCount){
        maxVertexCount = verticesAllowed;
    }
    verticesAllowed = 2*(maxFaceOrbitCount - requestedFaceSizesCount)*maxFaceSize;
    for(i=0; i<requestedFaceSizesCount; i++){
        verticesAllowed += 2*requestedFaceSizes[i];
    }
    if(verticesAllowed < minVertexCount){
        fprintf(stderr, "For these face sizes the Delaney-Dress graph can have at most %d vertices.\n", verticesAllowed);
        validConstraints = FALSE;
    } else if (verticesAllowed < maxVertexCount){
        maxVertexCount = verticesAllowed;
    }
    verticesAllowed = 2*(maxVertexOrbitCount - requestedVertexDegreesCount)*maxVertexDegree;
    for(i=0; i<requestedVertexDegreesCount; i++){
        verticesAllowed += 2*requestedVertexDegrees[i];
    }
    if(verticesAllowed < minVertexCount){
        fprintf(stderr, "For these vertex degrees the Delaney-Dress graph can have at most %d vertices.\n", verticesAllowed);
        validConstraints = FALSE;
    } else if (verticesAllowed < maxVertexCount){
        maxVertexCount = verticesAllowed;
    }
    
    if(2*minVertexDegree + 2*minFaceSize - minVertexDegree*minFaceSize<0){
        fprintf(stderr, "This minimum vertex size and minimum face size can't yield any tiling.\n");
        validConstraints = FALSE;
    }
    
    if(2*maxVertexDegree + 2*maxFaceSize - maxVertexDegree*maxFaceSize>0){
        fprintf(stderr, "This maximum vertex size and maximum face size can't yield any tiling.\n");
        validConstraints = FALSE;
    }
    
//    These bounds only apply for convex tilings!
//
//    int maximumFaceSizeAllowed = 2+maxFaceOrbitCount*4;
//    /* see Two finiteness Theorems for Periodic Tilings of d-Dimensional Euclidean
//     * Space, Dolbilin, Dress, Huson; Discrete Comput GEOM 20:143-153 (1998) for
//     * this bound. The formula given in lemma 4.1 is incorrect: it should be <= 
//     * instead of <.
//     */ 
//    if(minFaceSize>maximumFaceSizeAllowed){
//        fprintf(stderr, "For this number of face orbits the largest face possible\nhas size %d.\n", maximumFaceSizeAllowed);
//        validConstraints = FALSE;
//    } else if (maxFaceSize > maximumFaceSizeAllowed){
//        maxFaceSize = maximumFaceSizeAllowed;
//    }
//    for(i=0; i<requestedFaceSizesCount; i++){
//        if(requestedFaceSizes[i]>maximumFaceSizeAllowed){
//            fprintf(stderr, "For this number of face orbits the largest face possible\nhas size %d.\n", maximumFaceSizeAllowed);
//            validConstraints = FALSE;
//        }
//    }
//    int maximumVertexDegreeAllowed = 2+maxVertexOrbitCount*4;
//    /* the dual of the bond mentioned above
//     */ 
//    if(minVertexDegree>maximumVertexDegreeAllowed){
//        fprintf(stderr, "For this number of vertex orbits the largest vertex degree possible\nhas degree %d.\n", maximumVertexDegreeAllowed);
//        validConstraints = FALSE;
//    } else if (maxVertexDegree > maximumVertexDegreeAllowed){
//        maxVertexDegree = maximumVertexDegreeAllowed;
//    }
//    for(i=0; i<requestedVertexDegreesCount; i++){
//        if(requestedVertexDegrees[i]>maximumVertexDegreeAllowed){
//        fprintf(stderr, "For this number of vertex orbits the largest vertex degree possible\nhas degree %d.\n", maximumVertexDegreeAllowed);
//            validConstraints = FALSE;
//        }
//    }
    
    int maximumVertexOrbitsAllowed = (maxFaceOrbitCount - requestedFaceSizesCount)*maxFaceSize;
    for(i=0; i<requestedFaceSizesCount; i++){
        maximumVertexOrbitsAllowed += requestedFaceSizes[i];
    }
    if(minVertexOrbitCount>maximumVertexOrbitsAllowed){
        fprintf(stderr, "For these faces and minimum vertex degree there can be at most %d vertex orbits.\n", maximumVertexOrbitsAllowed);
        validConstraints = FALSE;
    } else if (maxVertexOrbitCount>maximumVertexOrbitsAllowed){
        maxVertexOrbitCount = maximumVertexOrbitsAllowed;
    }
    
    int maximumFaceOrbitsAllowed = (maxVertexOrbitCount - requestedVertexDegreesCount)*maxVertexDegree;
    for(i=0; i<requestedVertexDegreesCount; i++){
        maximumFaceOrbitsAllowed += requestedVertexDegrees[i];
    }
    if(minFaceOrbitCount>maximumFaceOrbitsAllowed){
        fprintf(stderr, "For these vertices and minimum face size there can be at most %d face orbits.\n", maximumFaceOrbitsAllowed);
        validConstraints = FALSE;
    } else if (maxFaceOrbitCount>maximumFaceOrbitsAllowed){
        maxFaceOrbitCount = maximumFaceOrbitsAllowed;
    }
    /*
     * TODO: these bounds are incorrect: can they be fixed?
    int muV = (maxFaceOrbitCount-requestedFaceSizesCount)*maxFaceSize;
    int muF = (maxVertexOrbitCount-requestedVertexDegreesCount)*maxVertexDegree;
    int maximumVertexDegreeAllowed = 0;
    maximumFaceSizeAllowed = 0;
    for(i=0; i<requestedFaceSizesCount; i++){
        muV += requestedFaceSizes[i];
        muF -= requestedFaceSizes[i];
        if(requestedFaceSizes[i]>maximumFaceSizeAllowed){
            maximumFaceSizeAllowed = requestedFaceSizes[i];
        }
    }
    for(i=0; i<requestedVertexDegreesCount; i++){
        muV -= requestedVertexDegrees[i];
        muF += requestedVertexDegrees[i];
        if(requestedVertexDegrees[i]>maximumVertexDegreeAllowed){
            maximumVertexDegreeAllowed = requestedVertexDegrees[i];
        }
    }
    if(requestedVertexDegreesCount<minVertexOrbitCount){
        muV -= (minVertexOrbitCount-requestedVertexDegreesCount-1)*minVertexDegree;
    }
    if(requestedFaceSizesCount<minFaceOrbitCount){
        muF -= (minFaceOrbitCount-requestedFaceSizesCount-1)*minFaceSize;
    }
    maximumVertexDegreeAllowed = MAX(maximumVertexDegreeAllowed, muV);
    maximumFaceSizeAllowed = MAX(maximumFaceSizeAllowed, muF);
    if(minVertexDegree>maximumVertexDegreeAllowed){
        fprintf(stderr, "For these parameters the vertex degree can be at most %d.\n", maximumVertexDegreeAllowed);
        validConstraints = FALSE;
    } else if(maxVertexDegree>maximumVertexDegreeAllowed){
        maxVertexDegree = maximumVertexDegreeAllowed;
    }
    if(minFaceSize>maximumFaceSizeAllowed){
        fprintf(stderr, "For these parameters the face size can be at most %d.\n", maximumFaceSizeAllowed);
        validConstraints = FALSE;
    } else if(maxFaceSize>maximumFaceSizeAllowed){
        maxFaceSize = maximumFaceSizeAllowed;
    }
    */
    boolean madeChanges = TRUE;
    while(madeChanges){
        madeChanges = FALSE;
        for(i=0; i<forbiddenFaceSizesCount; i++){
            if(minFaceSize==forbiddenFaceSizes[i]){
                minFaceSize++;
                madeChanges = TRUE;
            }
            if(maxFaceSize==forbiddenFaceSizes[i]){
                maxFaceSize--;
                madeChanges = TRUE;
            }
        }
    }
    madeChanges = TRUE;
    while(madeChanges){
        madeChanges = FALSE;
        for(i=0; i<forbiddenVertexDegreesCount; i++){
            if(minVertexDegree==forbiddenVertexDegrees[i]){
                minVertexDegree++;
                madeChanges = TRUE;
            }
            if(maxVertexDegree==forbiddenVertexDegrees[i]){
                maxVertexDegree--;
                madeChanges = TRUE;
            }
        }
    }
    
    //Adjust maximum number of face and vertex orbits based upon size of symbol
    if(maxFaceOrbitCount > maxVertexCount){
        maxFaceOrbitCount = maxVertexCount;
    }
    if(maxVertexOrbitCount > maxVertexCount){
        maxVertexOrbitCount = maxVertexCount;
    }
    
    //Adjust maximum face size and vertex degree based upon size of symbol
    int theoreticalMaxFaceSize = (maxVertexCount - minFaceOrbitCount + 1)*12;
    if(orientable) theoreticalMaxFaceSize/=2;
    if(maxFaceSize > theoreticalMaxFaceSize){
        maxFaceSize = theoreticalMaxFaceSize;
    }
    int theoreticalMaxVertexDegree = (maxVertexCount - minVertexOrbitCount + 1)*12;
    if(orientable) theoreticalMaxVertexDegree/=2;
    if(maxVertexDegree > theoreticalMaxVertexDegree){
        maxVertexDegree = theoreticalMaxVertexDegree;
    }
    return validConstraints;
}

void calculateSymbolSize(){
    int minVerticesNeeded, maxVerticesNeeded, i;
    
    //calculate minVertexCount
    minVerticesNeeded = (minFaceOrbitCount-requestedFaceSizesCount)*(minFaceSize/6);
    if(minFaceSize/6) minVerticesNeeded += (minFaceOrbitCount-requestedFaceSizesCount);
    for(i=0; i<requestedFaceSizesCount; i++){
        minVerticesNeeded += requestedFaceSizes[i]/6;
        if(requestedFaceSizes[i]%6) minVerticesNeeded++;
    }
    minVertexCount = MAX(minVertexCount, minVerticesNeeded);
    
    minVerticesNeeded = (minVertexOrbitCount-requestedVertexDegreesCount)*(minVertexDegree/6);
    if(minVertexDegree/6) minVerticesNeeded += (minVertexOrbitCount-requestedVertexDegreesCount);
    for(i=0; i<requestedVertexDegreesCount; i++){
        minVerticesNeeded += requestedVertexDegrees[i]/6;
        if(requestedVertexDegrees[i]%6) minVerticesNeeded++;
    }
    minVertexCount = MAX(minVertexCount, minVerticesNeeded);
    
    //calculate maxVertexCount
    maxVerticesNeeded = (maxFaceOrbitCount- requestedFaceSizesCount)*maxFaceSize*2;
    for(i=0; i<requestedFaceSizesCount; i++){
        maxVerticesNeeded += requestedFaceSizes[i]*2;
    }
    maxVertexCount = MIN(maxVertexCount, maxVerticesNeeded);
    
    maxVerticesNeeded = (maxVertexOrbitCount - requestedVertexDegreesCount)*maxVertexDegree*2;
    for(i=0; i<requestedVertexDegreesCount; i++){
        maxVerticesNeeded += requestedVertexDegrees[i]*2;
    }
    maxVertexCount = MIN(maxVertexCount, maxVerticesNeeded);
    
    minVertexCount = MAX(minVertexCount, minFaceOrbitCount);
    minVertexCount = MAX(minVertexCount, minVertexOrbitCount);
    minVertexCount = MAX(minVertexCount, minEdgeOrbitCount);
    
    maxVertexCount = MIN(maxVertexCount, maxEdgeOrbitCount*4);
    
    maxFaceOrbitCount = MIN(maxFaceOrbitCount, maxVertexCount);
    maxVertexOrbitCount = MIN(maxVertexOrbitCount, maxVertexCount);
    maxEdgeOrbitCount = MIN(maxEdgeOrbitCount, maxVertexCount);
}

//creates a table that keeps track for each face size whether it is forbidden or not
void createForbiddenTable(){
    int i;
    for(i=0; i<6*MAXN; i++){
        forbiddenFaceSizesTable[i]=FALSE;
        forbiddenVertexDegreesTable[i]=FALSE;
    }
    
    for(i=0;i<forbiddenFaceSizesCount; i++){
        forbiddenFaceSizesTable[forbiddenFaceSizes[i]]=TRUE;
    }
    
    for(i=0; i<forbiddenVertexDegreesCount; i++){
        forbiddenVertexDegreesTable[forbiddenVertexDegrees[i]]=TRUE;
    }
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
    static struct option long_options[] = {
        {"maxfacecount", required_argument, NULL, 0},
        {"minfacecount", required_argument, NULL, 0},
        {"maxvertexcount", required_argument, NULL, 0},
        {"minvertexcount", required_argument, NULL, 0},
        {"minfacesize", required_argument, NULL, 0},
        {"maxfacesize", required_argument, NULL, 0},
        {"minvertexdegree", required_argument, NULL, 0},
        {"maxvertexdegree", required_argument, NULL, 0},
        {"verbose", no_argument, NULL, 0},
        {"restrictionsonly", no_argument, NULL, 0},
        {"maxedgecount", required_argument, NULL, 0},
        {"minedgecount", required_argument, NULL, 0},
        {"statistics", no_argument, NULL, 0},
        {"filter", no_argument, NULL, 0},
        {"maxfaceR", required_argument, NULL, 0},
        {"minfaceR", required_argument, NULL, 0},
        {"maxvertexR", required_argument, NULL, 0},
        {"minvertexR", required_argument, NULL, 0},
        {"intermediate", no_argument, NULL, 0},
        {"help", no_argument, NULL, 'h'},
        {"lists", no_argument, NULL, 'L'},
        {"marked", no_argument, NULL, 't'},
        {"coloured", no_argument, NULL, 'c'},
        {"symbols", no_argument, NULL, 's'},
        {"listfile", required_argument, NULL, 'l'},
        {"output", required_argument, NULL, 'o'},
        {"modulo", required_argument, NULL, 'm'},
        {"requiredface", required_argument, NULL, 'R'},
        {"forbiddenface", required_argument, NULL, 'F'},
        {"requiredvertex", required_argument, NULL, 'r'},
        {"forbiddenvertex", required_argument, NULL, 'f'},
        {"minvertices", required_argument, NULL, 'n'},
        {"maxvertices", required_argument, NULL, 'N'},
        {"bipartite", no_argument, NULL, 'b'},
        {"orientable", no_argument, NULL, 'O'}
    };
    int option_index = 0;

    boolean failAfterArgumentParsing = FALSE;
    while ((c = getopt_long(argc, argv, "hl:Ltcso:m:R:F:r:f:n:N:bO", long_options, &option_index)) != -1) {
        switch (c) {
            case 0:
                //handle long option with no alternative
                switch(option_index) {
                    case 0:
                        maxFaceOrbitCount = atoi(optarg);
                        if(!checkIntegerValue(long_options[option_index].name, maxFaceOrbitCount, 1, MAXN)){
                            failAfterArgumentParsing = TRUE;
                        }
                        break;
                    case 1:
                        minFaceOrbitCount = atoi(optarg);
                        if(!checkIntegerValue(long_options[option_index].name, minFaceOrbitCount, 1, MAXN)){
                            failAfterArgumentParsing = TRUE;
                        }
                        break;
                    case 2:
                        maxVertexOrbitCount = atoi(optarg);
                        if(!checkIntegerValue(long_options[option_index].name, maxVertexOrbitCount, 1, MAXN)){
                            failAfterArgumentParsing = TRUE;
                        }
                        break;
                    case 3:
                        minVertexOrbitCount = atoi(optarg);
                        if(!checkIntegerValue(long_options[option_index].name, minVertexOrbitCount, 1, MAXN)){
                            failAfterArgumentParsing = TRUE;
                        }
                        break;
                    case 4:
                        minFaceSize = atoi(optarg);
                        if(!checkIntegerValue(long_options[option_index].name, minFaceSize, 3, 6*MAXN)){
                            failAfterArgumentParsing = TRUE;
                        }
                        break;
                    case 5:
                        maxFaceSize = atoi(optarg);
                        if(!checkIntegerValue(long_options[option_index].name, maxFaceSize, 3, 6*MAXN)){
                            failAfterArgumentParsing = TRUE;
                        }
                        break;
                    case 6:
                        minVertexDegree = atoi(optarg);
                        if(!checkIntegerValue(long_options[option_index].name, minVertexDegree, 3, 6*MAXN)){
                            failAfterArgumentParsing = TRUE;
                        }
                        break;
                    case 7:
                        maxVertexDegree = atoi(optarg);
                        if(!checkIntegerValue(long_options[option_index].name, maxVertexDegree, 3, 6*MAXN)){
                            failAfterArgumentParsing = TRUE;
                        }
                        break;
                    case 8:
                        verbose = TRUE;
                        break;
                    case 9:
                        restrictionsOnly = TRUE;
                        break;
                    case 10:
                        maxEdgeOrbitCount = atoi(optarg);
                        if(!checkIntegerValue(long_options[option_index].name, maxEdgeOrbitCount, 1, MAXN)){
                            failAfterArgumentParsing = TRUE;
                        }
                        break;
                    case 11:
                        minEdgeOrbitCount = atoi(optarg);
                        if(!checkIntegerValue(long_options[option_index].name, minEdgeOrbitCount, 1, MAXN)){
                            failAfterArgumentParsing = TRUE;
                        }
                        break;
                    case 12:
                        giveStatistics = TRUE;
                        break;
                    case 13:
                        filterDelaneyDressGraphs = TRUE;
                        break;
                    case 14:
                        maxFaceR = atoi(optarg);
                        if(!checkIntegerValue(long_options[option_index].name, maxFaceR, 1, 2*MAXN)){
                            failAfterArgumentParsing = TRUE;
                        }
                        break;
                    case 15:
                        minFaceR = atoi(optarg);
                        if(!checkIntegerValue(long_options[option_index].name, minFaceR, 1, 2*MAXN)){
                            failAfterArgumentParsing = TRUE;
                        }
                        break;
                    case 16:
                        maxVertexR = atoi(optarg);
                        if(!checkIntegerValue(long_options[option_index].name, maxVertexR, 1, 2*MAXN)){
                            failAfterArgumentParsing = TRUE;
                        }
                        break;
                    case 17:
                        minVertexR = atoi(optarg);
                        if(!checkIntegerValue(long_options[option_index].name, minVertexR, 1, 2*MAXN)){
                            failAfterArgumentParsing = TRUE;
                        }
                        break;
                    case 18:
                        intermediateStructures = TRUE;
                        break;
                    default:
                        fprintf(stderr, "Illegal option index %d.\n", option_index);
                        usage(name);
                        return EXIT_FAILURE;
                }
                break;
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
            case 'b':
                bipartite = TRUE;
                break;
            case 'O':
                bipartite = TRUE;
                orientable = TRUE;
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
            case 'R':
                requestedFaceSizes[requestedFaceSizesCount++] = atoi(optarg);
                if(!checkIntegerValue("requiredface", requestedFaceSizes[requestedFaceSizesCount-1], 3, 6*MAXN)){
                    failAfterArgumentParsing = TRUE;
                }
                break;
            case 'F':
                forbiddenFaceSizes[forbiddenFaceSizesCount++] = atoi(optarg);
                if(!checkIntegerValue("forbiddenface", forbiddenFaceSizes[forbiddenFaceSizesCount-1], 3, 6*MAXN)){
                    failAfterArgumentParsing = TRUE;
                }
                break;
            case 'r':
                requestedVertexDegrees[requestedVertexDegreesCount++] = atoi(optarg);
                if(!checkIntegerValue("requiredvertex", requestedVertexDegrees[requestedVertexDegreesCount-1], 3, 6*MAXN)){
                    failAfterArgumentParsing = TRUE;
                }
                break;
            case 'f':
                forbiddenVertexDegrees[forbiddenVertexDegreesCount++] = atoi(optarg);
                if(!checkIntegerValue("forbiddenvertex", forbiddenVertexDegrees[forbiddenVertexDegreesCount-1], 3, 6*MAXN)){
                    failAfterArgumentParsing = TRUE;
                }
                break;
            case 'n':
                minVertexCount = atoi(optarg);
                if(!checkIntegerValue("minvertices", minVertexCount, 1, MAXN)){
                    failAfterArgumentParsing = TRUE;
                }
                break;
            case 'N':
                maxVertexCount = atoi(optarg);
                if(!checkIntegerValue("maxvertices", maxVertexCount, 1, MAXN)){
                    failAfterArgumentParsing = TRUE;
                }
                break;
            case '?':
                usage(name);
                return EXIT_FAILURE;
            default:
                fprintf(stderr, "Illegal option %c.\n", c);
                usage(name);
                return EXIT_FAILURE;
        }
    }
    
    if(failAfterArgumentParsing){
        return EXIT_FAILURE;
    }

    // check the non-option arguments
    if (argc - optind != 1 && listFilename==NULL && !symbols) {
        usage(name);
        return EXIT_FAILURE;
    }

    /*=========== initialization ===========*/
    struct tms TMS;
    unsigned int oldtime = 0;

#ifdef _PROFILING
    {
        int i;
        for(i=0; i<MAXN/2; i++){
            graphsWithTrivialSymmetry[i]=0;
            graphsWithNonTrivialSymmetry[i]=0;
            graphsFromClosedGraphsWithTrivialSymmetry[i]=0;
            graphsWithTrivialSymmetryForRemainingConnections[i]=0;
            closedGraphsWithTrivialSymmetry[i]=0;
            closedGraphsWithNonTrivialSymmetry[i]=0;
            closedGraphsWithTrivialSymmetryForRemainingConnections[i]=0;
        }
    }
#endif
    
    //adjust some constant in case certain restrictions were chosen
    if(edgeColourable){
        Q1TypeComponentsSmallestCase[3] = 2;
    }
    
    if(symbols){
        int i;
        
        //validate restrictions for Delaney-Dress symbols and for tilings
        adjustSymbolConstraints();
        while(constraintsChanged()){
            if(!validateSymbolConstraints()){
                return EXIT_FAILURE;
            }
            //recalculate size limits for Delaney-Dress graph
            calculateSymbolSize();
        }
        
        createForbiddenTable();
        
        //start generation
        fprintf(stderr, "Generating Delaney-Dress symbols with %d to %d vertices.\n", minVertexCount, maxVertexCount);
        if(listFilename!=NULL){
            fprintf(stderr, "Using the component lists in %s as input.\n", listFilename);
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "Improved parameter bounds\n-------------------------\n");
        fprintf(stderr, "Number of face orbits lies in [%d,%d].\n", minFaceOrbitCount, maxFaceOrbitCount);
        fprintf(stderr, "Number of vertex orbits lies in [%d,%d].\n", minVertexOrbitCount, maxVertexOrbitCount);
        fprintf(stderr, "Number of edge orbits lies in [%d,%d].\n", minEdgeOrbitCount, maxEdgeOrbitCount);
        fprintf(stderr, "Face sizes lie in [%d,%d].\n", minFaceSize, maxFaceSize);
        fprintf(stderr, "Vertex degrees lie in [%d,%d].\n", minVertexDegree, maxVertexDegree);
        fprintf(stderr, "Required face sizes are [");
        for(i=0; i<requestedFaceSizesCount-1; i++){
            fprintf(stderr, "%d,", requestedFaceSizes[i]);
        }
        if(requestedFaceSizesCount>0){
            fprintf(stderr, "%d", requestedFaceSizes[requestedFaceSizesCount-1]);
        }
        fprintf(stderr, "].\n");
        fprintf(stderr, "Forbidden face sizes are [");
        for(i=0; i<forbiddenFaceSizesCount-1; i++){
            fprintf(stderr, "%d,", forbiddenFaceSizes[i]);
        }
        if(forbiddenFaceSizesCount>0){
            fprintf(stderr, "%d", forbiddenFaceSizes[forbiddenFaceSizesCount-1]);
        }
        fprintf(stderr, "].\n");
        fprintf(stderr, "Required vertex degrees are [");
        for(i=0; i<requestedVertexDegreesCount-1; i++){
            fprintf(stderr, "%d,", requestedVertexDegrees[i]);
        }
        if(requestedVertexDegreesCount>0){
            fprintf(stderr, "%d", requestedVertexDegrees[requestedVertexDegreesCount-1]);
        }
        fprintf(stderr, "].\n");
        fprintf(stderr, "Forbidden vertex degrees are [");
        for(i=0; i<forbiddenVertexDegreesCount-1; i++){
            fprintf(stderr, "%d,", forbiddenVertexDegrees[i]);
        }
        if(forbiddenVertexDegreesCount>0){
            fprintf(stderr, "%d", forbiddenVertexDegrees[forbiddenVertexDegreesCount-1]);
        }
        fprintf(stderr, "].\n");
        fprintf(stderr, "\n");
        
        if(!restrictionsOnly){
            if(listFilename!=NULL){
                startFromListFile(listFilename, minVertexCount, maxVertexCount);
            } else {
                startMultipleGenerations(minVertexCount, maxVertexCount);
            }
        }
        
        
    } else if(listFilename!=NULL){
        if(onlyLists){
            fprintf(stderr, "Generating component lists for Delaney-Dress graphs based on component lists in %s.\n",
                listFilename);
        } else if(colouredEdges){
            fprintf(stderr, "Generating edge-coloured Delaney-Dress graphs based on component lists in %s.\n",
                listFilename);
        } else {
            fprintf(stderr, "Generating Delaney-Dress graphs%s  based on component lists in %s.\n",
                markedTwoFactors ? (char *)" with marked 2-factors" : (char *)"",
                listFilename);
        }
        startFromListFile(listFilename, 0, MAXN);
    } else {

        //parse the order
        int vertexCount = strtol(argv[optind], NULL, 10);
        DEBUGDUMP(vertexCount, "%d")
                
        if(onlyLists){
            fprintf(stderr, "Generating component lists for Delaney-Dress graphs with %d %s.\n",
                vertexCount,
                vertexCount==1 ? (char *)"vertex" : (char *)"vertices");
        } else if(colouredEdges){
            fprintf(stderr, "Generating edge-coloured Delaney-Dress graphs with %d %s.\n",
                vertexCount,
                vertexCount==1 ? (char *)"vertex" : (char *)"vertices");
        } else {
            fprintf(stderr, "Generating Delaney-Dress graphs with%s %d %s.\n",
                markedTwoFactors ? (char *)" marked 2-factors and" : (char *)"",
                vertexCount,
                vertexCount==1 ? (char *)"vertex" : (char *)"vertices");
        }

        startGeneration(vertexCount);
    }



    times(&TMS);
    unsigned int savetime = oldtime + (unsigned int) TMS.tms_utime;
    fprintf(stderr, "CPU time: %.1f seconds.\n", (double) savetime / time_factor);
    
    if(giveStatistics && symbols && symbolsCount){
        fprintf(stderr, "\nDelaney-Dress symbol statistics:\n");
        
        fprintf(stderr, "Delaney-Dress graph order lies in [%d,%d].\n", minOrderStatistic, maxOrderStatistic);
        fprintf(stderr, "Number of face orbits lies in [%d,%d].\n", minFaceOrbitCountStatistic, maxFaceOrbitCountStatistic);
        fprintf(stderr, "Number of vertex orbits lies in [%d,%d].\n", minVertexOrbitCountStatistic, maxVertexOrbitCountStatistic);
        fprintf(stderr, "Number of edge orbits lies in [%d,%d].\n", minEdgeOrbitCountStatistic, maxEdgeOrbitCountStatistic);
        fprintf(stderr, "Face sizes lie in [%d,%d].\n", minFaceSizeStatistic, maxFaceSizeStatistic);
        fprintf(stderr, "Vertex degrees lie in [%d,%d].\n", minVertexDegreeStatistic, maxVertexDegreeStatistic);
    }
    
#ifdef _PROFILING
    
    if(!onlyLists){
        fprintf(stderr, "Extra profiling info:\n");
        fprintf(stderr, "Connections rejected\n");
        fprintf(stderr, "     based on colour  : %7llu\n", rejectedByColour);
        fprintf(stderr, "     by nauty         : %7llu\n", rejectedByNauty);
        fprintf(stderr, "Connections accepted\n");
        fprintf(stderr, "     because only one : %7llu (%llu)\n", acceptedBecauseOnlyOne, skippedNautyBecauseOnlyOne);
        fprintf(stderr, "     based on colour  : %7llu (%llu)\n", acceptedBecauseOnlyOneMinimalColour, skippedNautyBecauseOnlyOneMinimalColour);
        fprintf(stderr, "     by nauty         : %7llu\n\n", connectionsAccepted);

        {
            int i = 0;
            fprintf(stderr, "A : graphs with trivial symmetry.\n");
            fprintf(stderr, "B : graphs with non-trivial symmetry.\n");
            fprintf(stderr, "C : graphs derived from a closed graph with trivial symmetry.\n");
            fprintf(stderr, "C : graphs symmetry group that acts trivial on the remaining connections.\n");
            fprintf(stderr, "+-------------------------------------------------------------+\n");
            fprintf(stderr, "|                     |    A    |    B    |    C    |    D    |\n");
            fprintf(stderr, "|---------------------+---------+---------+---------+---------|\n");
            while(i<MAXN/2 && (graphsWithTrivialSymmetry[i] || graphsWithNonTrivialSymmetry[i])){
                fprintf(stderr, "|After %2d connections | %7llu | %7llu | %7llu | %7llu |\n",
                        i, graphsWithTrivialSymmetry[i], graphsWithNonTrivialSymmetry[i],
                        graphsFromClosedGraphsWithTrivialSymmetry[i], graphsWithTrivialSymmetryForRemainingConnections[i]);
                i++;
            }
            fprintf(stderr, "+-------------------------------------------------------------+\n\n");
        }
        {
            int i = 0;
            fprintf(stderr, "A : closed graphs with trivial symmetry.\n");
            fprintf(stderr, "B : closed graphs with non-trivial symmetry.\n");
            fprintf(stderr, "C : closed graphs symmetry group that acts trivial on the remaining connections.\n");
            fprintf(stderr, "+-------------------------------------------------------------+\n");
            fprintf(stderr, "|                     |    A    |    B    |    C    |   C-A   |\n");
            fprintf(stderr, "|---------------------+---------+---------+---------+---------|\n");
            while(i<MAXN/2 && (closedGraphsWithTrivialSymmetry[i] || closedGraphsWithNonTrivialSymmetry[i])){
                fprintf(stderr, "|After %2d connections | %7llu | %7llu | %7llu | %7llu |\n",
                        i, closedGraphsWithTrivialSymmetry[i], closedGraphsWithNonTrivialSymmetry[i],
                        closedGraphsWithTrivialSymmetryForRemainingConnections[i],
                        closedGraphsWithTrivialSymmetryForRemainingConnections[i] - closedGraphsWithTrivialSymmetry[i]);
                i++;
            }
            fprintf(stderr, "+-------------------------------------------------------------+\n");
        }
    }
    
    if(symbols){
        fprintf(stderr, "Number of accepted lists: %llu\n", numberOfListsAcceptedForSymbols);
        fprintf(stderr, "Number of rejected lists\n");
        fprintf(stderr, "    because pearl chain too long            : %llu\n", rejectedListsBecausePearlChainTooLong);
        fprintf(stderr, "    because locked pearl chain too long     : %llu\n", rejectedListsBecauseLockedPearlChainTooLong);
        fprintf(stderr, "    because too many closing semi-edges (1) : %llu\n", rejectedListsBecauseTooManySemiEdgesForVertexOrFaceOrbitCount);
        fprintf(stderr, "    because too many closing semi-edges (2) : %llu\n", rejectedListsBecauseTooManySemiEdgesForCombinedOrbitCount);
        fprintf(stderr, "    because too few edges with colour 1     : %llu\n", rejectedListsBecauseTooFewColour1Edges);
        fprintf(stderr, "Number of marked graphs: %llu\n", graphsCount);
        fprintf(stderr, "Number of edge-coloured graphs: %llu\n", edgeColouredGraphsCount);
        fprintf(stderr, "Number of accepted edge-coloured graphs: %llu\n", acceptedColouredGraphs);
        fprintf(stderr, "Number of rejected edge-coloured graphs\n");
        fprintf(stderr, "    because wrong number face orbits        : %llu\n", rejectedColouredGraphBecauseWrongNumberFaceOrbits);
        fprintf(stderr, "    because wrong number vertex orbits      : %llu\n", rejectedColouredGraphBecauseWrongNumberVertexOrbits);
        fprintf(stderr, "    because too many orbits                 : %llu\n", rejectedColouredGraphBecauseTooManyOrbits);
        fprintf(stderr, "    because too big face orbit              : %llu\n", rejectedColouredGraphBecauseFaceOrbitTooBig);
        fprintf(stderr, "    because too big vertex orbit            : %llu\n", rejectedColouredGraphBecauseVertexOrbitTooBig);
        fprintf(stderr, "    because incompatible parameters         : %llu\n", rejectedColouredGraphBecauseIncompatibleParameters);
        fprintf(stderr, "Number of possible symbols: %llu\n", possibleAssignments);
        fprintf(stderr, "Number of valid symbols: %llu\n", validAssignments);
        fprintf(stderr, "Number of Delaney-Dress graphs which were not used  : %llu\n", colouredDelaneyDressGraphsWithoutSymbol);
        fprintf(stderr, "Number of Delaney-Dress graphs which were used once : %llu\n", colouredDelaneyDressGraphsWithOneSymbol);
        fprintf(stderr, "Number of Delaney-Dress graphs which were used more : %llu\n", colouredDelaneyDressGraphsWithMultipleSymbols);
    }
#endif 
    
    return EXIT_SUCCESS;
}
