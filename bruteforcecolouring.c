// cc -O4 -Wall -o bfcolouring bruteforcecolouring.c

//reads a list of cubic pregraphs and does a brute-force colouring outputting all colourings
//without trying to filter out isomorphic colourings

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <unistd.h>
#include <sys/times.h>

#define SEMIEDGE INT_MAX
#define UNSET (INT_MAX - 1)
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

int countColoured = 0;

boolean writethem = FALSE;
boolean printThem = FALSE;

void colour_graph_impl(GRAPH *g, int vertex, int position);
char write_pgc_code(FILE *f, GRAPH *g, boolean first);
void printGraph(FILE *f, GRAPH *g);

void remove_colour(GRAPH *g, int colour){
    int i, j;
    for(i=0; i<g->order; i++){
        for(j=0; j<3; j++){
            if(g->colours[i][j]==colour){
                g->colours[i][j]=UNSET;
            }
        }
    }
}

/**************************Colour a new graph*************************************/

void do_next_edge(GRAPH *g, int vertex, int position){
    if(position==2 && vertex+1<g->order){
        colour_graph_impl(g, vertex + 1, 0);
    } else if(position<2){
        colour_graph_impl(g, vertex, position+1);
    } else {
        countColoured++;
        if(writethem){
            write_pgc_code(stdout, g, countColoured==1);
        }
        if(printThem){
            printGraph(stderr, g);
        }
    }
    return;
}

boolean try_colour(GRAPH *g, int vertex, int position, int colour){
    int i;
    int uncolouredEdgesCount = 0;
    int sumColours = 0;
    int target = g->adjacency[vertex][position];
    for(i = 0; i<3; i++){
        if(g->colours[target][i]==UNSET){
            uncolouredEdgesCount++;
        } else {
            sumColours+=g->colours[target][i];
        }
    }

    if((uncolouredEdgesCount==2 && sumColours==colour) ||
            (uncolouredEdgesCount==1 && sumColours + colour != 3)){
        return FALSE;
    }

    i=0;
    while(g->adjacency[target][i]!=vertex || g->colours[target][i]!=UNSET) i++;

    g->colours[vertex][position] = colour;
    g->colours[target][i] = colour;

    return TRUE;
}

void unset_colour(GRAPH *g, int vertex, int position){
    int i = 2;
    int target = g->adjacency[vertex][position];

    while(g->adjacency[target][i]!=vertex || g->colours[target][i]==UNSET) i--;

    g->colours[vertex][position] = UNSET;
    g->colours[target][i] = UNSET;
}

void colour_graph_impl(GRAPH *g, int vertex, int position){
    //skip to next edge if this edge is already coloured
    if(g->colours[vertex][position]!=UNSET){
        do_next_edge(g, vertex, position);
        return;
    }

    if(g->adjacency[vertex][position]==SEMIEDGE){
        int i;
        int uncolouredEdgesCount = 0;
        int sumColours = 0;
        for(i = 0; i<3; i++){
            if(g->colours[vertex][i]==UNSET){
                uncolouredEdgesCount++;
            } else {
                sumColours+=g->colours[vertex][i];
            }
        }
        if(uncolouredEdgesCount==3){
            g->colours[vertex][position] = 0;
            do_next_edge(g, vertex, position);
            g->colours[vertex][position] = UNSET;

            g->colours[vertex][position] = 1;
            do_next_edge(g, vertex, position);
            g->colours[vertex][position] = UNSET;

            g->colours[vertex][position] = 2;
            do_next_edge(g, vertex, position);
            g->colours[vertex][position] = UNSET;
        } else if(uncolouredEdgesCount==2){
            if(sumColours!=0) {
                g->colours[vertex][position] = 0;
                do_next_edge(g, vertex, position);
                g->colours[vertex][position] = UNSET;
            }
            if(sumColours!=1){
                g->colours[vertex][position] = 1;
                do_next_edge(g, vertex, position);
                g->colours[vertex][position] = UNSET;
            }
            if(sumColours!=2){
                g->colours[vertex][position] = 2;
                do_next_edge(g, vertex, position);
                g->colours[vertex][position] = UNSET;
            }
        } else /*uncolouredEdgesCount==1*/ {
            if(sumColours==1){
                g->colours[vertex][position] = 2;
                do_next_edge(g, vertex, position);
                g->colours[vertex][position] = UNSET;
            } else if(sumColours==2){
                g->colours[vertex][position] = 1;
                do_next_edge(g, vertex, position);
                g->colours[vertex][position] = UNSET;
            } else /*sumColours==3*/ {
                g->colours[vertex][position] = 0;
                do_next_edge(g, vertex, position);
                g->colours[vertex][position] = UNSET;
            }
        }
    } else {
        int i;
        int uncolouredEdgesCount = 0;
        int sumColours = 0;
        for(i = 0; i<3; i++){
            if(g->colours[vertex][i]==UNSET){
                uncolouredEdgesCount++;
            } else {
                sumColours+=g->colours[vertex][i];
            }
        }
        if(uncolouredEdgesCount==3){
            for(i=0; i<3; i++){
                if(try_colour(g, vertex, position, i)){
                    do_next_edge(g, vertex, position);
                    unset_colour(g, vertex, position);
                }
            }
        } else if(uncolouredEdgesCount==2){
            for(i=0; i<3; i++){
                if(sumColours!=i && try_colour(g, vertex, position, i)){
                    do_next_edge(g, vertex, position);
                    unset_colour(g, vertex, position);
                }
            }
        } else /*uncolouredEdgesCount==1*/ {
            if(sumColours==1){
                if(try_colour(g, vertex, position, 2)){
                    do_next_edge(g, vertex, position);
                    unset_colour(g, vertex, position);
                }
            } else if(sumColours==2){
                if(try_colour(g, vertex, position, 1)){
                    do_next_edge(g, vertex, position);
                    unset_colour(g, vertex, position);
                }
            } else /*sumColours==3*/ {
                if(try_colour(g, vertex, position, 0)){
                    do_next_edge(g, vertex, position);
                    unset_colour(g, vertex, position);
                }
            }
        }
    }
}

void colour_graph(GRAPH *g){
    colour_graph_impl(g, 0, 0);
}

/****************************Read pregraphcolor_code************************/

/* Reads pregraph_code and returns EOF if the end of the file is reached
 * and 1 in case of success.
 */

int read_pg_code(FILE *f, GRAPH *g) {
    int order, a, b, i, j;
    int adjacency[4];
    int adjacencyIndex=-1;

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

        for (j = 0; j < adjacencyIndex; j++) {
            if (adjacency[j] == order + 1) {
                int currentDegree = g->degrees[i];
                g->adjacency[i][currentDegree] = SEMIEDGE;
                g->degrees[i]++;
            } else {
                int v = adjacency[j] - 1;
                g->adjacency[i][g->degrees[i]] = v;
                g->adjacency[v][g->degrees[v]] = i;
                g->degrees[i]++;
                g->degrees[v]++;
            }
        }
        adjacencyIndex = -1;

        for(j=0; j<3; j++) g->colours[i][j] = UNSET;
    }

    return 1;
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

void printGraph(FILE *f, GRAPH *g){
    int i, j;
    for(i=0; i<g->order; i++){
        fprintf(f, "%3d) ", i);
        for(j = 0; j<g->degrees[i]; j++){
            if(g->adjacency[i][j]==SEMIEDGE){
                fprintf(f, "  S (%2d)", g->colours[i][j]==UNSET ? -1 : g->colours[i][j]);
            } else {
                fprintf(f, "%3d (%2d)", g->adjacency[i][j], g->colours[i][j]==UNSET ? -1 : g->colours[i][j]);
            }
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n");
}

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

void help(char *name){
    fprintf(stderr, "The program %s reads cubic pregraphs from stdin\n", name);
    fprintf(stderr, "and determines all possible colourings.\n");
    fprintf(stderr, "Usage: %s [options]\n\n", name);
    fprintf(stderr, "Valid options:\n");
    fprintf(stderr, "  -h     : Print this help and return.\n");
    fprintf(stderr, "  -w     : Write the graphs in pregraphcolor_code to stdout.\n");
    fprintf(stderr, "  -p     : Print the graphs in human readable tables to stderr.\n");
    fprintf(stderr, "  -e     : Read a partially coloured pregraph and complete that colouring.\n");
    fprintf(stderr, "  -c col : Specify the colour that needs to be ignored in case option e\n");
    fprintf(stderr, "           is used. This defaults to colour 3.\n");
}

void usage(char *name){
    fprintf(stderr, "Usage: %s [options]\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

int main(int argc, char *argv[]) {

    GRAPH g;
    int countRead = 0;
    boolean extend = FALSE;
    int dummyColour = 3;

    int c;
    char *name = argv[0];

    while ((c = getopt(argc, argv, "hwpec:")) != -1) {
        switch (c) {
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            case 'w':
                writethem = TRUE;
                break;
            case 'p':
                printThem = TRUE;
                break;
            case 'e':
                extend = TRUE;
                break;
            case 'c':
                dummyColour = strtol(optarg, NULL, 10);
                break;
            default:
                usage(name);
                return EXIT_FAILURE;
        }
    }

    struct tms TMS;
    unsigned int oldtime = 0;

    if(extend){
        while (read_pgc_code(stdin, &g) != EOF) {
            remove_colour(&g, dummyColour);
            countRead++;
            colour_graph(&g);
        }
    } else {
        while (read_pg_code(stdin, &g) != EOF) {
            countRead++;
            colour_graph(&g);
        }
    }

    fprintf(stderr, "Read %d graphs, found %d coloured graphs.\n", countRead, countColoured);

    times(&TMS);
    unsigned int savetime = oldtime + (unsigned int) TMS.tms_utime;
    fprintf(stderr, "CPU time: %.1f seconds.\n", (double) savetime / time_factor);
    
    return (0);

}
