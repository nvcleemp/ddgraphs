// cc -O4 -Wall -o ddsymbol_non_iso ddsymbol_non_iso.c

//reads a list of Delaney-Dress symbols and determines whether it contains isomorphic copies

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <limits.h>
#include <getopt.h>
#include <sys/times.h>

#define MEMORYCHECK(ptr) if((ptr)==NULL){\
                fprintf(stderr, "Sorry, but I really need more memory to continue.\n");\
                exit(1);\
            }

#define MAXORDER 100
#define DIMENSION 2
#define TRUE 1
#define FALSE 0

typedef int boolean;

typedef struct _symbol {
    int order;
    int adjacency[MAXORDER][DIMENSION+1];
    int m[MAXORDER][DIMENSION];
} SYMBOL;

typedef struct le {
    int *certificate;
    int certificateLength;
    
    int timesSeen;
    
    struct le *smaller;
    struct le *larger;
} LISTENTRY;

LISTENTRY *list;

/*************************** I/O ************************************/

void printCertificate(SYMBOL *symbol, int *certificate){
    int i;
    for(i=0; i<symbol->order*(2*DIMENSION+1); i++){
        fprintf(stderr, "%d ", certificate[i]);
    }
    fprintf(stderr, "\n");
}

void printSymbol(SYMBOL *symbol){
    int i, j;
    
    for(i=0; i<symbol->order; i++){
        fprintf(stderr, "%d) ", i+1);
        for(j=0; j<=DIMENSION; j++){
            fprintf(stderr, "%d ", symbol->adjacency[i][j]+1);
        }
        fprintf(stderr, "| ");
        for(j=0; j<DIMENSION; j++){
            fprintf(stderr, "%d ", symbol->m[i][j]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

/* Slightly less save version of getdelim
 */
ssize_t mygetdelim (char **lineptr, size_t *n, int delimiter, FILE *fp) {
  ssize_t result = EOF;
  size_t cur_len = 0;

  if (lineptr == NULL || n == NULL || fp == NULL){
      return EOF;
  }

  if (*lineptr == NULL || *n == 0){
      *n = 120;
      *lineptr = (char *) malloc (*n);
      if (*lineptr == NULL){
          return EOF;
      }
  }

  for (;;){
      int i;

      i = getc (fp);
      if (i == EOF){
          break;
      }

      /* Make enough space for len+1 (for final NUL) bytes.  */
      if (cur_len + 1 >= *n) {
          size_t needed = 2 * *n + 1;   /* Be generous. */
          char *new_lineptr;

          new_lineptr = (char *) realloc (*lineptr, needed);
          if (new_lineptr == NULL){
              return -1;
          }

          *lineptr = new_lineptr;
          *n = needed;
      }

      (*lineptr)[cur_len] = i;
      cur_len++;

      if (i == delimiter) break;
  }
  (*lineptr)[cur_len] = '\0';
  result = cur_len ? cur_len : result;

  return result;
}

void setOrbitNumber(SYMBOL *ddsymbol, int startVertex, int mValue, int mNumber){
    int vertex = startVertex;
    int s1 = mNumber;
    int s2 = mNumber + 1;
    while(ddsymbol->m[vertex][mNumber]<0){
        ddsymbol->m[vertex][mNumber]=mValue;
        vertex = ddsymbol->adjacency[vertex][s1];
        ddsymbol->m[vertex][mNumber]=mValue;
        vertex = ddsymbol->adjacency[vertex][s2];
    }
    vertex = ddsymbol->adjacency[startVertex][s2];
    while(ddsymbol->m[vertex][mNumber]<0){
        ddsymbol->m[vertex][mNumber]=mValue;
        vertex = ddsymbol->adjacency[vertex][s1];
        ddsymbol->m[vertex][mNumber]=mValue;
        vertex = ddsymbol->adjacency[vertex][s2];
    }
}

/* Reads a Delaney-Dress symbol and returns EOF if the end of the file is reached
 * and 1 in case of success.
 */
int readSymbol(FILE *f, SYMBOL *ddsymbol, char **string) {
    int i,j,nextFreeVertex;
    
    //Read symbol from file
    size_t nbytes = 200;
    char *symbolString = (char *)malloc(nbytes);
    
    if(mygetdelim (&symbolString, &nbytes, '\n', stdin)==EOF){
        return EOF;
    }
    
    *string = symbolString;
    
    //Split into parts
    //char *part1 = strchr(symbolString, '<');
    char *part2 = strchr(symbolString, ':') + 1;
    char *part3 = strchr(part2, ':') + 1;
    char *part4 = strchr(part3, ':') + 1;
    
    //Parse order and dimension and validate
    int order, dimension = 2;
    if(sscanf(part2, "%d %d:", &order, &dimension)==EOF){
        fprintf(stderr, "Reached end of file\n");
        exit(1);
    }
    
    if(dimension!=DIMENSION){
        fprintf(stderr, "This program needs to be recompiled to handle this dimension.\n");
        exit(1);
    }
    
    if(order>MAXORDER){
        fprintf(stderr, "This program needs to be recompiled to handle this order.\n");
        exit(1);
    }
    ddsymbol->order=order;
    
    //Clean symbol
    for(i=0; i<=DIMENSION; i++){
        for(j=0; j<order; j++){
            ddsymbol->adjacency[j][i]=-1;
        }
    }
    for(i=0; i<DIMENSION; i++){
        for(j=0; j<order; j++){
            ddsymbol->m[j][i]=-1;
        }
    }
    
    //Parse adjacency information
    for(i=0; i<DIMENSION; i++){
        nextFreeVertex = 0;
        char *nextS = strchr(part3, ',');
        while(nextFreeVertex<order){
            int neighbour;
            if(part3==NULL || part3>=nextS){
                fprintf(stderr, "Error while reading symbol (Adjacency information for s%d)\n", i);
                exit(1);
            }
            if(sscanf(part3, "%d", &neighbour)==EOF){
                fprintf(stderr, "Reached end of file\n");
                exit(1);
            }
            ddsymbol->adjacency[nextFreeVertex][i]=neighbour-1;
            ddsymbol->adjacency[neighbour-1][i]=nextFreeVertex;
            while(nextFreeVertex<order && ddsymbol->adjacency[nextFreeVertex][i]>=0) nextFreeVertex++;
            part3 = strchr(part3, ' ');
            if(part3!=NULL) part3++;
        }
        part3 = nextS + 1;
    }
    
    nextFreeVertex = 0;
    while(nextFreeVertex<order){
        int neighbour;
        if(part3==NULL || part3>=part4){
            fprintf(stderr, "Error while reading symbol (Adjacency information for s%d)\n", DIMENSION);
            exit(1);
        }
        if(sscanf(part3, "%d", &neighbour)==EOF){
            fprintf(stderr, "Reached end of file\n");
            exit(1);
        }
        ddsymbol->adjacency[nextFreeVertex][DIMENSION]=neighbour-1;
        ddsymbol->adjacency[neighbour-1][DIMENSION]=nextFreeVertex;
        while(nextFreeVertex<order && ddsymbol->adjacency[nextFreeVertex][DIMENSION]>=0) nextFreeVertex++;
        part3 = strchr(part3, ' ');
        if(part3!=NULL) part3++;
    }
    
    //Parse orbit numbers
    for(i=0; i<DIMENSION-1; i++){
        nextFreeVertex = 0;
        char *nextM = strchr(part4, ',');
        while(nextFreeVertex<order){
            int m;
            if(part4==NULL || part4>=nextM){
                fprintf(stderr, "Error while reading symbol (Orbit numbers for m%d%d)\n", i, i+1);
                exit(1);
            }
            if(sscanf(part4, "%d", &m)==EOF){
                fprintf(stderr, "Reached end of file\n");
                exit(1);
            }
            setOrbitNumber(ddsymbol, nextFreeVertex, m, i);
            while(nextFreeVertex<order && ddsymbol->m[nextFreeVertex][i]>=0) nextFreeVertex++;
            part4 = strchr(part4, ' ');
            if(part4!=NULL) part4++;
        }
        part4 = nextM + 1;
    }
    
    nextFreeVertex = 0;
    char *end = strchr(part4, '>');
    while(nextFreeVertex<order){
        int m;
        if(part4==NULL || part4>=end){
            fprintf(stderr, "Error while reading symbol (Orbit number for m%d%d)\n", DIMENSION-1 , DIMENSION);
            exit(1);
        }
        if(sscanf(part4, "%d", &m)==EOF){
            fprintf(stderr, "Reached end of file\n");
            exit(1);
        }
        setOrbitNumber(ddsymbol, nextFreeVertex, m, DIMENSION-1);
        while(nextFreeVertex<order && ddsymbol->m[nextFreeVertex][DIMENSION-1]>=0) nextFreeVertex++;
        part4 = strchr(part4, ' ');
        if(part4!=NULL) part4++;
    }
    
    return 1;
}

/********************** Canonical form ****************************/

void getCertificate(SYMBOL *ddsymbol, int *certificate){
    int i,j;
    int old2new[MAXORDER];
    int new2old[MAXORDER];
    int enqueued[MAXORDER];
    int tempCertificate[MAXORDER*(2*DIMENSION+1)];
    for(i=0;i<MAXORDER;i++) enqueued[i]=FALSE;
            
    int queue[MAXORDER];
    int head = 0, tail = 0;
    int currentLabel = 0;

    //first we calculate the certificate if we label vertex 0 with 0 and look
    //for the canonical labeling
    queue[tail++] = 0;
    enqueued[0]=TRUE;
    while(head<tail){
        int v = queue[head++];
        new2old[currentLabel]=v;
        old2new[v]=currentLabel++;
        for(i=0; i<(DIMENSION+1); i++){
            if(!enqueued[ddsymbol->adjacency[v][i]]){
                queue[tail++] = ddsymbol->adjacency[v][i];
                enqueued[ddsymbol->adjacency[v][i]] = TRUE;
            }
        }
    }
    
    for(i=0; i<ddsymbol->order; i++){
        int old = new2old[i];
        for(j=0; j<(DIMENSION+1); j++){
            certificate[i*(DIMENSION+1) + j] = old2new[ddsymbol->adjacency[old][j]];
        }
        for(j=0; j<DIMENSION; j++){
            certificate[(ddsymbol->order)*(DIMENSION+1) + i*DIMENSION + j] = ddsymbol->m[old][j];
        }
    }
    
    int start;
    for(start=1; start<ddsymbol->order; start++){
        //we try to give each vertex the label 0
        head = 0, tail = 0;
        queue[tail++] = start;
        currentLabel = 0;
        for(i=0;i<MAXORDER;i++) enqueued[i]=FALSE;
        enqueued[start]=TRUE;
        
        while(head<tail){
            int v = queue[head++];
            new2old[currentLabel]=v;
            old2new[v]=currentLabel++;
            for(i=0; i<3; i++){
                if(!enqueued[ddsymbol->adjacency[v][i]]){
                    queue[tail++] = ddsymbol->adjacency[v][i];
                    enqueued[ddsymbol->adjacency[v][i]] = TRUE;
                }
            }
        }

        for(i=0; i<ddsymbol->order; i++){
            int old = new2old[i];
            for(j=0; j<(DIMENSION+1); j++){
                tempCertificate[i*(DIMENSION+1) + j] = old2new[ddsymbol->adjacency[old][j]];
            }
            for(j=0; j<DIMENSION; j++){
                tempCertificate[(ddsymbol->order)*(DIMENSION+1) + i*DIMENSION + j] = ddsymbol->m[old][j];
            }
        }
        
        int i=0;
        while(i<(2*DIMENSION+1)*ddsymbol->order && certificate[i]==tempCertificate[i]) i++;
        if(i<(2*DIMENSION+1)*(ddsymbol->order)){
            if(certificate[i]>tempCertificate[i]){
                //copy tempCertificate to certificate
                //TODO: use memcpy
                for(j=i; j<(2*DIMENSION+1)*ddsymbol->order; j++){
                    certificate[j] = tempCertificate[j];
                }
            }
        }
    }
}

/**************************** List operations *******************************/

/* Construct a list entry for the given symbol
 */
void constructListEntry(LISTENTRY *el, int *certificate, int certificateLength){
    el->certificate = certificate;
    el->certificateLength = certificateLength;
    el->timesSeen = 1;
    el->smaller = NULL;
    el->larger = NULL;
}

/* Recursively searches the list for the appropriate place to add this symbol.
 * In case the symbol is not yet in the list it is added and TRUE is returned.
 * If the symbol is in the list, nothing is done and FALSE is returned.
 */
boolean addToList(LISTENTRY *el, int *certificate, int certificateLength){
    int compare = certificateLength - el->certificateLength;
    
    if(compare == 0){
        //compare = memcmp(certificate, el->certificate, certificateLength);
        int i=0;
        while(i<certificateLength && compare==0){
            compare = certificate[i] - el->certificate[i];
            i++;
        }
    }
    
    if(compare == 0){
        el->timesSeen++;
        return FALSE;
    } else if(compare < 0){
        if(el->smaller == NULL){
            el->smaller = (LISTENTRY *)malloc(sizeof(LISTENTRY));
            MEMORYCHECK(el->smaller)
            constructListEntry(el->smaller, certificate, certificateLength);
            return TRUE;
        } else {
            return addToList(el->smaller, certificate, certificateLength);
        }
    } else /*(compare > 0)*/ {
        if(el->larger == NULL){
            el->larger = (LISTENTRY *)malloc(sizeof(LISTENTRY));
            MEMORYCHECK(el->larger)
            constructListEntry(el->larger, certificate, certificateLength);
            return TRUE;
        } else {
            return addToList(el->larger, certificate, certificateLength);
        }
    }
}

int listSize(LISTENTRY *listLocal){
    if(listLocal == NULL){
        return 0;
    } else {
        return 1 + listSize(listLocal->larger) + listSize(listLocal->smaller);
    }
}

/* Handle a symbol once it has been read, i.e., determine its certificate
 * and add it to the list when it is not present at the moment.
 * This method returns TRUE if the symbol was added to the list and FALSE
 * if it was already in the list.
 */
boolean handleSymbol(SYMBOL *symbol, LISTENTRY **list){
    int *certificate = malloc((2*DIMENSION+1)*MAXORDER*sizeof(int));
    MEMORYCHECK(certificate)
    
    getCertificate(symbol, certificate);
    
    if((*list)==NULL){
        *list = (LISTENTRY *)malloc(sizeof(LISTENTRY));
        MEMORYCHECK((*list))
        constructListEntry(*list, certificate, symbol->order*(2*DIMENSION+1));
        return TRUE;
    } else {
        return addToList(*list, certificate, symbol->order*(2*DIMENSION+1));
    }
}
/*
print a usage message. name is the name of the current program.
 */
void usage(char *name) {
    fprintf(stderr, "Usage: %s [options] n\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

/*
print a help message. name is the name of the current program.
 */
void help(char *name) {
    fprintf(stderr, "The program %s reads a list of Delaney-Dress symbols and checks it for isomorphic copies.\n", name);
    fprintf(stderr, "Usage: %s [options] n \n\n", name);
    fprintf(stderr, "Valid options:\n");
    fprintf(stderr, "  -h    : Print this help and return.\n");
    fprintf(stderr, "  -N    : Print symbols when they are read for the first time.\n");
    fprintf(stderr, "  -O    : Print symbols when they are read for a second, third, ... time.\n");
}

int main(int argc, char *argv[]){
    
    int c;
    char *name = argv[0];
    boolean printOld = FALSE;
    boolean printNew = FALSE;
    

    while ((c = getopt(argc, argv, "hON")) != -1) {
        switch (c) {
            case 'O': //(defaults to FALSE)
                printOld = TRUE;
                break;
            case 'N': //(defaults to FALSE)
                printNew = TRUE;
                break;
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            default:
                fprintf(stderr, "Illegal option %c.\n", c);
                usage(name);
                return EXIT_FAILURE;
        }
    }
    
    SYMBOL s;
    char *symbolString;
    int structuresRead = 0;
    int uniques = 0;
    
    while(readSymbol(stdin, &s, &symbolString)!=EOF){
        structuresRead++;
        boolean isNew = handleSymbol(&s, &list);
        if(isNew){
            if(printNew){
                fprintf(stderr, "Symbol %d is new: %s", structuresRead, symbolString);
            }
            uniques++;
        } else if(printOld){
            fprintf(stderr, "Symbol %d is old: %s", structuresRead, symbolString);
        }
    }
    fprintf(stderr, "Read %d symbol(s).\nFound %d unique symbol(s).\n", structuresRead, uniques);
    
    return EXIT_SUCCESS;
}