
// cc -O4 -Wall -o validateddgraph validateddgraph.c

//reads a list of edge-coloured cubic pregraphs and determines whether it are valid Delaney-Dress graph

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <limits.h>
#include <unistd.h>
#include <sys/times.h>

#define SEMIEDGE INT_MAX
#define MAXORDER 200
#define TRUE 1
#define FALSE 0
#define time_factor sysconf(_SC_CLK_TCK)

short endian = LITTLE_ENDIAN; // defines which endian should be used while exporting pregraph code

typedef int boolean;

typedef struct _graph {
    int order;
    int degrees[MAXORDER];
    int adjacency[MAXORDER][3];
    int colours[MAXORDER][3];
} GRAPH;

boolean check_vertex(GRAPH *g, int v){
    //fprintf(stderr, "Vertex: %d\n", v);
    int p0=0, p2=0;

    while(g->colours[v][p0]!=0) p0++;
    while(g->colours[v][p2]!=2) p2++;

    //fprintf(stderr, "Colour 0 at position %d\n", p0);
    //fprintf(stderr, "Colour 2 at position %d\n", p2);

    if(g->adjacency[v][p0]==g->adjacency[v][p2]){
        //fprintf(stderr, "Found q2 or q4");
        //q2 and q4
        return TRUE;
    }

    if(g->adjacency[v][p0]==SEMIEDGE){
        //fprintf(stderr, "Semi-edge has colour 0.\n");
        //should be q3
        int n2 = g->adjacency[v][p2];
        int np0 = 0;
        while(g->colours[n2][np0]!=0) np0++;
        return g->adjacency[n2][np0]==SEMIEDGE;
    } else if(g->adjacency[v][p2]==SEMIEDGE){
        //fprintf(stderr, "Semi-edge has colour 2.\n");
        //should be q3
        int n0 = g->adjacency[v][p0];
        int np2 = 0;
        while(g->colours[n0][np2]!=2) np2++;
        //fprintf(stderr, "  Neigbour at colour 0 is %d\n", n0);
        //fprintf(stderr, "  Position of colour 2 is %d\n", np2);
        return g->adjacency[n0][np2]==SEMIEDGE;
    } else {
        //fprintf(stderr, "Found q1.\n");
        //should be q1
        int n2 = g->adjacency[v][p2];
        int np0 = 0;
        while(g->colours[n2][np0]!=0) np0++;

        int n0 = g->adjacency[v][p0];
        int np2 = 0;
        while(g->colours[n0][np2]!=2) np2++;

        return g->adjacency[n0][np2] == g->adjacency[n2][np0];
    }
}

boolean validate_ddgraph(GRAPH *g){
    int i = 0;

    while(i<g->order && check_vertex(g, i)) i++;

    //fprintf(stderr, "\n\n");

    return i==g->order;
}

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

char write_pgc_code(FILE *f, GRAPH *g, boolean first){
    unsigned short i, j;
    unsigned short semiEdge = g->order + 1;
    if (first) { //if first graph
        fprintf(f, ">>pregraphcolor_code %s<<", (endian == LITTLE_ENDIAN ? "le" : "be"));
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

void help(char *name){
    fprintf(stderr, "The program %s reads edge-coloured cubic pregraphs from stdin\n", name);
    fprintf(stderr, "and checks whether the colouring corresponds to a valid Delaney-Dress graph.\n");
    fprintf(stderr, "Usage: %s [options]\n\n", name);
    fprintf(stderr, "Valid options:\n");
    fprintf(stderr, "  -h    : Print this help and return.\n");
    fprintf(stderr, "  -w    : Write the valid graphs in pregraphcolor_code to stdout.\n");
    fprintf(stderr, "  -p    : Print the valid graphs in human readable form to stderr.\n");
    fprintf(stderr, "  -n    : Causes -p and -w to output the non-valid graphs.\n");
}

void usage(char *name){
    fprintf(stderr, "Usage: %s [options]\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

int main(int argc, char *argv[]) {

    GRAPH g;
    int countTotal = 0;
    int validTotal = 0;
    boolean writethem = FALSE, printThem = FALSE;
    boolean first = TRUE;
    boolean negate = FALSE;

    int c;
    char *name = argv[0];

    while ((c = getopt(argc, argv, "hpwn")) != -1) {
        switch (c) {
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            case 'p':
                printThem = TRUE;
                break;
                writethem = TRUE;
                break;
            case 'w':
                writethem = TRUE;
                break;
            case 'n':
                negate = TRUE;
                break;
            default:
                usage(name);
                return EXIT_FAILURE;
        }
    }

    while (read_pgc_code(stdin, &g) != EOF) {
        countTotal++;
        
        if(negate){
            if(validate_ddgraph(&g)){
                validTotal++;
            } else {
                if(writethem){
                    write_pgc_code(stdout, &g, first);
                    first=FALSE;
                }
                if(printThem) {
                    printGraph(stderr, &g);
                }
            }
        } else {
            if(validate_ddgraph(&g)){
                validTotal++;
                if(writethem){
                    write_pgc_code(stdout, &g, first);
                    first=FALSE;
                }
                if(printThem) {
                    printGraph(stderr, &g);
                }
            }
        }

    }

    struct tms TMS;
    unsigned int oldtime = 0;

    fprintf(stderr, "Read %d graphs, of which %d had a valid colouring.\n", countTotal, validTotal);

    times(&TMS);
    unsigned int savetime = oldtime + (unsigned int) TMS.tms_utime;
    fprintf(stderr, "CPU time: %.1f seconds.\n", (double) savetime / time_factor);
    
    return (0);

}
