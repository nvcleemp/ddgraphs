// cc -O4 -Wall -o ddsymbol_summary ddsymbol_summary.c

//reads a list of Delaney-Dress symbols and gives a summary of the list

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <limits.h>
#include <getopt.h>
#include <sys/times.h>

#define MAXORDER 100
#define DIMENSION 2
#define TRUE 1
#define FALSE 0

typedef int boolean;

typedef struct _symbol {
    int order;
    int adjacency[MAXORDER][DIMENSION+1];
    int m[MAXORDER][DIMENSION];
    int orbitCount[DIMENSION];
    int orbitM[MAXORDER][DIMENSION];
} SYMBOL;

int minimumOrder = INT_MAX;
int maximumOrder = 0;
int minimumNumberFaceOrbits = INT_MAX;
int maximumNumberFaceOrbits = 0;
int minimumNumberVertexOrbits = INT_MAX;
int maximumNumberVertexOrbits = 0;
int minimumFaceSize = INT_MAX;
int maximumFaceSize = 0;
int minimumVertexDegree = INT_MAX;
int maximumVertexDegree = 0;
int requiredFaces[6*MAXORDER];
int requiredVertices[6*MAXORDER];
boolean forbiddenFaces[6*MAXORDER];
boolean forbiddenVertices[6*MAXORDER];
boolean semiEdgeSigma[DIMENSION+1];

/*************************** I/O ************************************/

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
    ddsymbol->orbitM[ddsymbol->orbitCount[mNumber]][mNumber] = mValue;
    ddsymbol->orbitCount[mNumber]++;
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
        ddsymbol->orbitCount[i] = 0;
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

/* Handle a symbol once it has been read.
 */
void handleSymbol(SYMBOL *symbol){
    if(maximumOrder < symbol->order){
        maximumOrder = symbol->order;
    }
    if(minimumOrder > symbol->order){
        minimumOrder = symbol->order;
    }
    if(maximumNumberFaceOrbits < symbol->orbitCount[0]){
        maximumNumberFaceOrbits = symbol->orbitCount[0];
    }
    if(minimumNumberFaceOrbits > symbol->orbitCount[0]){
        minimumNumberFaceOrbits = symbol->orbitCount[0];
    }
    if(maximumNumberVertexOrbits < symbol->orbitCount[1]){
        maximumNumberVertexOrbits = symbol->orbitCount[1];
    }
    if(minimumNumberVertexOrbits > symbol->orbitCount[1]){
        minimumNumberVertexOrbits = symbol->orbitCount[1];
    }
    
    int localFaceSizes[6*MAXORDER];
    int localVertexDegrees[6*MAXORDER];
    int i, j;
    
    for(i = 0; i < 6*MAXORDER; i++){
        localFaceSizes[i] = 0;
        localVertexDegrees[i] = 0;
    }
    
    for(i = 0; i < symbol->orbitCount[0]; i++){
        int currentFaceSize = symbol->orbitM[i][0];
        localFaceSizes[currentFaceSize]++;
        if(currentFaceSize > maximumFaceSize){
            maximumFaceSize = currentFaceSize;
        }
        if(currentFaceSize < minimumFaceSize){
            minimumFaceSize = currentFaceSize;
        }
    }
    
    for(i = 0; i < 6*MAXORDER; i++){
        if(localFaceSizes[i] && forbiddenFaces[i]){
            forbiddenFaces[i] = FALSE;
        }
        if(localFaceSizes[i] < requiredFaces[i]){
            requiredFaces[i] = localFaceSizes[i];
        }
    }
    
    for(i = 0; i < symbol->orbitCount[1]; i++){
        int currentVertexDegree = symbol->orbitM[i][1];
        localVertexDegrees[currentVertexDegree]++;
        if(currentVertexDegree > maximumVertexDegree){
            maximumVertexDegree = currentVertexDegree;
        }
        if(currentVertexDegree < minimumVertexDegree){
            minimumVertexDegree = currentVertexDegree;
        }
    }
    
    for(i = 0; i < 6*MAXORDER; i++){
        if(localVertexDegrees[i] && forbiddenVertices[i]){
            forbiddenVertices[i] = FALSE;
        }
        if(localVertexDegrees[i] < requiredVertices[i]){
            requiredVertices[i] = localVertexDegrees[i];
        }
    }
    
    for(i = 0; i < DIMENSION + 1; i++){
        for(j = 0; j < symbol->order; j++){
            if(symbol->adjacency[j][i] == j){
                semiEdgeSigma[i] = TRUE;
            }
        }
    }
}
/*
print a usage message. name is the name of the current program.
 */
void usage(char *name) {
    fprintf(stderr, "Usage: %s [options]\n", name);
    fprintf(stderr, "For more information type: %s -h\n\n", name);
}

/*
print a help message. name is the name of the current program.
 */
void help(char *name) {
    fprintf(stderr, "The program %s reads a list of Delaney-Dress symbols and prints a summary.\n", name);
    fprintf(stderr, "Usage: %s [options]\n\n", name);
    fprintf(stderr, "Valid options:\n");
    fprintf(stderr, "  -h    : Print this help and return.\n");
}

int main(int argc, char *argv[]){
    
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
    
    SYMBOL s;
    int i, j;
    char *symbolString;
    int structuresRead = 0;
    
    for(i = 0; i < 6*MAXORDER; i++){
        requiredFaces[i] = INT_MAX;
        requiredVertices[i] = INT_MAX;
        forbiddenFaces[i] = TRUE;
        forbiddenVertices[i] = TRUE;
    }
    for(i = 0; i < DIMENSION+1; i++){
        semiEdgeSigma[i] = FALSE;
    }
    
    while(readSymbol(stdin, &s, &symbolString)!=EOF){
        structuresRead++;
        handleSymbol(&s);
    }
    fprintf(stderr, "Read %d symbol(s).\n\n", structuresRead);
    
    fprintf(stderr, "Order lies in [%d,%d].\n", minimumOrder, maximumOrder);
    fprintf(stderr, "Number of face orbits lies in [%d,%d].\n", minimumNumberFaceOrbits, maximumNumberFaceOrbits);
    fprintf(stderr, "Number of vertex orbits lies in [%d,%d].\n", minimumNumberVertexOrbits, maximumNumberVertexOrbits);
    fprintf(stderr, "Face sizes lie in [%d,%d].\n", minimumFaceSize, maximumFaceSize);
    fprintf(stderr, "Vertex degrees lie in [%d,%d].\n", minimumVertexDegree, maximumVertexDegree);
    fprintf(stderr, "Required faces: ");
    for(i = minimumFaceSize; i <= maximumFaceSize; i++){
        for(j = 0; j < requiredFaces[i]; j++){
            fprintf(stderr, "%d ", i);
        }
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "Required vertices: ");
    for(i = minimumVertexDegree; i <= maximumVertexDegree; i++){
        for(j = 0; j < requiredVertices[i]; j++){
            fprintf(stderr, "%d ", i);
        }
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "Forbidden faces: ");
    for(i = minimumFaceSize; i <= maximumFaceSize; i++){
        if(forbiddenFaces[i]){
            fprintf(stderr, "%d ", i);
        }
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "Forbidden vertices: ");
    for(i = minimumVertexDegree; i <= maximumVertexDegree; i++){
        if(forbiddenVertices[i]){
            fprintf(stderr, "%d ", i);
        }
    }
    fprintf(stderr, "\n");
    
    for(i = 0; i <= DIMENSION; i++){
        if(semiEdgeSigma[i]){
            fprintf(stderr, "Symbol contains s%d semi edge.\n", i);
        } else {
            fprintf(stderr, "Symbol doesn't contain s%d semi edge.\n", i);
        }
    }
    
    return EXIT_SUCCESS;
}
