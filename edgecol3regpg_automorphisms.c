// cc -O4 -Wall -o edgecol3regpg_automorphisms edgecol3regpg_automorphisms.c

//reads a list of edge-coloured cubic pregraphs and determines the number of automorphisms

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <limits.h>
#include <getopt.h>

#define SEMIEDGE INT_MAX
#define MAXORDER 200

typedef int boolean;

#define TRUE 1
#define FALSE 0

short endian = LITTLE_ENDIAN; // defines which endian should be used while exporting pregraph code

typedef struct _graph {
    int order;
    int degrees[MAXORDER];
    int adjacency[MAXORDER][3];
    int colours[MAXORDER][3];
} GRAPH;

int automorphismCounts[MAXORDER];

boolean filterOnAutomorphismCount;
int filterOnAutomorphismCountValue = -1;

/****************************Read pregraphcolor_code************************/

/* Reads pregraphcolor_code and returns EOF if the end of the file is reached
 * and 1 in case of success.
 */

int read_pgc_code(FILE *f, GRAPH *g) {
    int order, a, b, i, j;
    int adjacency[4];
    int adjacencyIndex=-1;
    int colours[4];
    int coloursIndex=-1;

    if ((order = getc(f)) == EOF) return EOF;
    if (order == 0) {
        fprintf(stderr, "Version for graphs of this order is not yet implemented.\n");
        exit(0);
    }

    if (order == '>') {
        // could be a header or a 62
        a = getc(f);
        b = getc(f);
        if ((a == '>') && (b == 'p')) {
            //now we are sure that this is a header
            unsigned char temp;
            while ((temp = getc(f)) != '<');
            //almost at the end of the header
            temp = getc(f);
            if (temp != '<') {
                fprintf(stderr, "Problems with header -- single '<'\n");
                exit(1);
            }
            if ((order = getc(f)) == EOF) return EOF; //empty file

            g->order = order;
            for (i = 0; i < order; i++) g->degrees[i] = 0;

        } else {
            //otherwise no header

            g->order = order;
            for (i = 0; i < order; i++) g->degrees[i] = 0;

            //a and b can not be 0 in this case
            adjacency[0] = a;
            adjacency[1] = b;
            adjacencyIndex = 1;
        }
    } else {
        g->order = order;
        for (i = 0; i < order; i++) g->degrees[i] = 0;
    }

    for (i = 0; i < order; i++) {
        //read adjacency information
        while ((adjacency[++adjacencyIndex] = getc(f)) != 0);
        //read colour information
        while ((colours[++coloursIndex] = getc(f)) != 0);

        if (adjacencyIndex != coloursIndex) {
            fprintf(stderr, "Problem with file: adjacency list and colour list have different length.\n");
            exit(1);
        }

        for (j = 0; j < adjacencyIndex; j++) {
            if (adjacency[j] == order + 1) {
                int currentDegree = g->degrees[i];
                g->adjacency[i][currentDegree] = SEMIEDGE;
                g->colours[i][currentDegree] = colours[j] - 1;
                g->degrees[i]++;
            } else {
                int v = adjacency[j] - 1;
                g->adjacency[i][g->degrees[i]] = v;
                g->colours[i][g->degrees[i]] = colours[j] - 1;
                g->adjacency[v][g->degrees[v]] = i;
                g->colours[v][g->degrees[v]] = colours[j] - 1;
                g->degrees[i]++;
                g->degrees[v]++;
            }
        }
        adjacencyIndex = coloursIndex = -1;
    }

    return 1;
}

/************************************Write graph*****************************/

char write_2byte_number(FILE *f, unsigned short n, short writeEndian) {
    if (writeEndian == BIG_ENDIAN) {
        fprintf(f, "%c%c", n / 256, n % 256);
    } else {
        fprintf(f, "%c%c", n % 256, n / 256);
    }
    return (ferror(f) ? 2 : 1);
}

char write_pgc_code(FILE *f, GRAPH *g){
    static boolean first = TRUE;
    
    unsigned short i, j;
    unsigned short semiEdge = g->order + 1;
    if (first) { //if first graph
        fprintf(f, ">>pregraphcolor_code %s<<", (endian == LITTLE_ENDIAN ? "le" : "be"));
        first = FALSE;
    }
    if (g->order + 1 <= UCHAR_MAX) {
        fprintf(f, "%c", (unsigned char) g->order);
    } else {
        fprintf(f, "%c", 0);
        /* big graph */
        if (write_2byte_number(f, (unsigned short) g->order, endian) == 2) {
            return (2);
        }
    }

    int colours[3];
    int adjacencyListSize;

    for (i = 0; i < g->order; i++) {
        adjacencyListSize = 0;
        for (j = 0; j < 3; j++) {
            int neighbour = g->adjacency[i][j];
            if(neighbour==SEMIEDGE){
                if (g->order + 1 <= UCHAR_MAX) {
                    fprintf(f, "%c", (unsigned char)semiEdge);
                } else {
                    if (write_2byte_number(f, semiEdge, endian) == 2) {
                        return (2);
                    }
                }
                colours[adjacencyListSize] = g->colours[i][j]+1;
                adjacencyListSize++;
            } else {
                if(neighbour>i){
                    //only include adjacency information for vertices with a larger index
                    if (g->order + 1 <= UCHAR_MAX) {
                        fprintf(f, "%c", (unsigned char) (neighbour + 1));
                    } else {
                        if (write_2byte_number(f, neighbour + 1, endian) == 2) {
                            return (2);
                        }
                    }
                    colours[adjacencyListSize] = g->colours[i][j]+1;
                    adjacencyListSize++;
                }
            }
        }
        //closing 0
        if (g->order + 1 <= UCHAR_MAX) {
            fprintf(f, "%c", 0);
        } else {
            if (write_2byte_number(f, 0, endian) == 2) {
                return (2);
            }
        }
        //colour list
        for(j=0; j<adjacencyListSize; j++){
            if (g->order + 1 <= UCHAR_MAX) {
                fprintf(f, "%c", colours[j]);
            } else {
                if (write_2byte_number(f, colours[j], endian) == 2) {
                    return (2);
                }
            }
        }
        //closing 0
        if (g->order + 1 <= UCHAR_MAX) {
            fprintf(f, "%c", 0);
        } else {
            if (write_2byte_number(f, 0, endian) == 2) {
                return (2);
            }
        }
    }
    return (ferror(f) ? 2 : 1);
}

void printGraph(FILE *f, GRAPH *g){
    int i, j;
    for(i=0; i<g->order; i++){
        fprintf(f, "%3d) ", i);
        for(j = 0; j<g->degrees[i]; j++){
            if(g->adjacency[i][j]==SEMIEDGE){
                fprintf(f, "  S (%d)", g->colours[i][j]);
            } else {
                fprintf(f, "%3d (%d)", g->adjacency[i][j], g->colours[i][j]);
            }
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n");
}

//*******************************************************************************


void handle_new_graph(GRAPH *g){
    int i, j, k;
    int automorphismCount = 1; //we'll always have the identity
    
    int certificate[3*MAXORDER];
    int alternateCertificate[3*MAXORDER];
    int relabelling[MAXORDER];
    int inverseRelabelling[MAXORDER];
    boolean labelled[MAXORDER];
    
    for(i = 0; i < g->order; i++){
        for(j = 0; j < 3; j++){
            certificate[3*i + j] = g->adjacency[i][j];
        }
    }
    
    for(i = 1; i < g->order; i++){
        for(j = 0; j < MAXORDER; j++){
            labelled[j] = FALSE;
        }
        //relabel with i -> 0
        relabelling[i] = 0;
        labelled[i] = TRUE;
        int nextLabel = 1;
        int queue[MAXORDER];
        int tail, head;
        queue[0] = i;
        head = 1;
        tail = 0;
        
        while(tail!=head){
            int currentVertex = queue[tail++]; //pop vertex from queue
            for(j = 0; j < 3; j++){
                int nextVertex = g->adjacency[currentVertex][j];
                if(nextVertex < g->order && !labelled[nextVertex]){
                    labelled[nextVertex] = TRUE;
                    relabelling[nextVertex] = nextLabel++;
                    queue[head++] = nextVertex;
                }
            }
        }
        
        for(j = 0; j < g->order; j++){
            inverseRelabelling[relabelling[j]] = j;
        }
        //check alternate certificate
        for(j = 0; j < g->order; j++){
            for(k = 0; k < 3; k++){
                if(g->adjacency[inverseRelabelling[j]][k] < g->order){
                    alternateCertificate[3*j+k] = relabelling[g->adjacency[inverseRelabelling[j]][k]];
                } else {
                    //semi-edge
                    alternateCertificate[3*j+k] = g->adjacency[inverseRelabelling[j]][k];
                }
            }
        }
        
        if(memcmp(certificate, alternateCertificate, sizeof(int)*3*g->order)==0){
            //found automorphism
            automorphismCount++;
        }
    }
    automorphismCounts[automorphismCount]++;
    //fprintf(stderr, "Size of automorphism group: %d\n", automorphismCount);
    if(filterOnAutomorphismCount && automorphismCount == filterOnAutomorphismCountValue){
        write_pgc_code(stdout, g);
    }
}

void normaliseGraph(GRAPH *g){
    int i, j;
    int adj[3];
    for(i = 0; i < g->order; i++){
        for(j = 0; j < 3; j++){
            adj[g->colours[i][j]] = g->adjacency[i][j];
        }
        for(j = 0; j < 3; j++){
            g->adjacency[i][j] = adj[j];
            g->colours[i][j] = j;
        }
    }
}

//*******************************************************************************

void help(char *name){
    fprintf(stderr, "The program %s reads edge-coloured cubic pregraphs from stdin\n", name);
    fprintf(stderr, "and determines the number of automorphisms in each one.\n");
    fprintf(stderr, "Usage: %s [options]\n\n", name);
    fprintf(stderr, "Valid options:\n");
    fprintf(stderr, "  -h    : Print this help and return.\n");
}

void usage(char *name){
    fprintf(stderr, "Usage: %s [options]\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

int main(int argc, char *argv[]) {
    
    int i;

    GRAPH g;
    int countTotal = 0;

    int c;
    char *name = argv[0];

    while ((c = getopt(argc, argv, "hf:")) != -1) {
        switch (c) {
            case 'f':
                filterOnAutomorphismCount = TRUE;
                filterOnAutomorphismCountValue = atoi(optarg);
                break;
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            default:
                usage(name);
                return EXIT_FAILURE;
        }
    }
    
    for(i = 0; i < MAXORDER; i++){
        automorphismCounts[i] = 0;
    }

    while (read_pgc_code(stdin, &g) != EOF) {
        countTotal++;
        normaliseGraph(&g);
        handle_new_graph(&g);

    }

    fprintf(stderr, "Read %d graphs\n", countTotal);
    
    //print frequency table
    for(i = 0; i < MAXORDER; i++){
        if(automorphismCounts[i])
            fprintf(stderr, "%2d: %d\n", i, automorphismCounts[i]);
    }
    
    return (0);

}
