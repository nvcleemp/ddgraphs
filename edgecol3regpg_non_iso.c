// cc -O4 -Wall -o edgecol3regpg_non_iso edgecol3regpg_non_iso.c nautil.c naugraph.c nauty.c

//reads a list of edge-coloured cubic pregraphs and determines whether it contains isomorphic copies

#include<stdio.h>
#include<stdlib.h>
#include<memory.h>
#include <limits.h>
#include "nauty.h"

#define SEMIEDGE INT_MAX
#define MAXORDER 200

short endian = LITTLE_ENDIAN; // defines which endian should be used while exporting pregraph code

typedef struct _graph {
    int order;
    int degrees[MAXORDER];
    int adjacency[MAXORDER][3];
    int colours[MAXORDER][3];
} GRAPH;

typedef graph *NAUTYGRAPH;

typedef struct le {
    NAUTYGRAPH nautyg;
    int order;
    struct le *smaller;
    struct le *larger;
} LISTENTRY;

LISTENTRY *list;

/**********Translate****************************************************/

void translate(GRAPH *g, NAUTYGRAPH nautyg, int m, int newOrder) {

    int i, j, k, dummy;

    dummy = g->order;

    for (i = 0; i < newOrder * m; i++) nautyg[i] = 0;

    for (i = 0; i < g->order; i++) {
        for (j = 0; j < 3; j++) {
            int neighbour = g->adjacency[i][j];
            if (neighbour == SEMIEDGE) neighbour = dummy++;
            if (neighbour > i) {
                if (g->colours[i][j] == 0) {
                    ADDELEMENT(nautyg + (m * i), neighbour);
                    ADDELEMENT(nautyg + (m * neighbour), i);
                } else {
                    int next = dummy++;
                    ADDELEMENT(nautyg + (m * i), next);
                    ADDELEMENT(nautyg + (m * next), i);
                    for (k = 0; k < g->colours[i][j] - 1; k++) {
                        ADDELEMENT(nautyg + (m * dummy), next);
                        ADDELEMENT(nautyg + (m * next), dummy);
                        next = dummy++;
                    }
                    ADDELEMENT(nautyg + (m * neighbour), next);
                    ADDELEMENT(nautyg + (m * next), neighbour);
                }
            }
        }
    }
}

int getNewOrder(GRAPH *g) {
    int newOrder, i, j;

    newOrder = g->order;

    for (i = 0; i < g->order; i++) {
        for (j = 0; j < 3; j++) {
            newOrder += g->colours[i][j];
            if (g->adjacency[i][j] == SEMIEDGE) {
                newOrder++;
            }
        }
    }
    return newOrder;
}

/**************************Construct list********************************/
void construct_list(LISTENTRY *el, NAUTYGRAPH canong, int m, int order) {

    el->nautyg = canong;
    el->order = order;
    el->smaller = NULL;
    el->larger = NULL;
}

/**************************Add to list********************************/

void add_to_list(LISTENTRY *el, NAUTYGRAPH canong, int m, int order, int *test) {
    int compare;

    compare = order - el->order;
    if (compare == 0) {
        compare = memcmp(canong, el->nautyg, m * order * sizeof (unsigned long));
    }

    if (compare == 0) {
        *test = 0;
    } else if (compare < 0) {
        if (el->smaller == NULL) {
            *test = 1;
            el->smaller = (LISTENTRY *) malloc(sizeof (LISTENTRY));
            if (el->smaller == NULL) {
                fprintf(stderr, "Can not get more memory !\n");
                exit(99);
            }
            construct_list(el->smaller, canong, m, order);
        } else {
            add_to_list(el->smaller, canong, m, order, test);
        }
    } else /* compare > 0 */ {
        if (el->larger == NULL) {
            *test = 1;
            el->larger = (LISTENTRY *) malloc(sizeof (LISTENTRY));
            if (el->larger == NULL) {
                fprintf(stderr, "Can not get more memory !\n");
                exit(99);
            }
            construct_list(el->larger, canong, m, order);
        } else {
            add_to_list(el->larger, canong, m, order, test);
        }
    }

}

/**************************Handle a new graph*************************************/

int handle_new_graph(GRAPH *g, LISTENTRY **list) {
    NAUTYGRAPH nautyg, canong;
    nvector lab[MAXORDER], ptn[MAXORDER], orbits[MAXORDER];
    static DEFAULTOPTIONS(options);
    statsblk(stats);
    setword workspace[100 * MAXORDER];
    int m, test, order;

    order = getNewOrder(g);
    options.getcanon = TRUE;
    options.writeautoms = FALSE;
    options.writemarkers = FALSE;

    if ((order % 32) == 0) m = order / 32;
    else m = order / 32 + 1;

    nautyg = (NAUTYGRAPH) calloc(m*order, sizeof (unsigned long));
    canong = (NAUTYGRAPH) calloc(m*order, sizeof (unsigned long));
    if ((nautyg == NULL) || (canong == NULL)) {
        fprintf(stderr, "Can not get more memory (1)!\n");
        exit(99);
    }

    translate(g, nautyg, m, order);
    nauty(nautyg, lab, ptn, NULL, orbits, &options, &stats, workspace, 100, m, order, canong);

    if ((*list) == NULL) {
        test = 1;
        *list = (LISTENTRY *) malloc(sizeof (LISTENTRY));
        if (*list == NULL) {
            fprintf(stderr, "Can not get more memory (0)!\n");
            exit(99);
        }
        construct_list(*list, canong, m, order);
    } else {
        add_to_list(*list, canong, m, order, &test);
    }

    free(nautyg);
    if (test == 0) free(canong); /* sonst wurde er in die Liste geschrieben */

    return (test);
}

/****************************Determine list size****************************/

int listSize(LISTENTRY *listLocal){
    if(listLocal==NULL){
        return 0;
    } else {
        return 1 + listSize(listLocal->larger) + listSize(listLocal->smaller);
    }
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
    fprintf(stderr, "and determines the number of non-isomorphic graphs it has read.\n");
    fprintf(stderr, "Usage: %s [options]\n\n", name);
    fprintf(stderr, "Valid options:\n");
    fprintf(stderr, "  -h    : Print this help and return.\n");
    fprintf(stderr, "  -w    : Write one graph for each isomorphism class.\n");
    fprintf(stderr, "          The graphs are written in the order they are read.\n");
    fprintf(stderr, "  -i    : Write info about which graphs are new and which aren't.\n");
    fprintf(stderr, "  -n    : Write info about which graphs are new.\n");
    fprintf(stderr, "  -o    : Write info about which graphs aren't new.\n");
}

void usage(char *name){
    fprintf(stderr, "Usage: %s [options]\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

int main(int argc, char *argv[]) {

    GRAPH g;
    int countTotal = 0;
    boolean writethem = FALSE, printThem = FALSE, printNew = FALSE, printOld = FALSE;
    boolean first = TRUE;

    list = NULL;

    int c;
    char *name = argv[0];

    while ((c = getopt(argc, argv, "hiwnop")) != -1) {
        switch (c) {
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            case 'i':
                printNew = TRUE;
                printOld = TRUE;
                break;
            case 'n':
                printNew = TRUE;
                break;
            case 'o':
                printOld = TRUE;
                break;
            case 'p':
                printThem = TRUE;
                break;
            case 'w':
                writethem = TRUE;
                break;
            default:
                usage(name);
                return EXIT_FAILURE;
        }
    }

    while (read_pgc_code(stdin, &g) != EOF) {
        countTotal++;
        int isNew = handle_new_graph(&g, &list);

        if(isNew && writethem) {
            write_pgc_code(stdout, &g, first);
            first = FALSE;
        }

        if (isNew) {
            if(printNew) fprintf(stderr, "Graph %d is new.\n", countTotal);
        } else {
            if(printOld) fprintf(stderr, "Graph %d is not new.\n", countTotal);
        }

        if(isNew && printThem) {
            printGraph(stderr, &g);
        }

    }

    int countNonIso = listSize(list);

    fprintf(stderr, "Read %d graphs, of which %d pairwise not isomorph.\n", countTotal, countNonIso);

    return (0);

}
