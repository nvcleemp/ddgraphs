// cc -O4 -Wall -o ddsymbol_filter ddsymbol_filter.c

//reads a list of Delaney-Dress symbols and filters it based upon given restrictions

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

int minimumOrder = 0;
int maximumOrder = INT_MAX;
int minimumNumberFaceOrbits = 0;
int maximumNumberFaceOrbits = INT_MAX;
int minimumNumberVertexOrbits = 0;
int maximumNumberVertexOrbits = INT_MAX;
int minimumFaceSize = 0;
int maximumFaceSize = INT_MAX;
int minimumVertexDegree = 0;
int maximumVertexDegree = INT_MAX;

boolean allownonplain = TRUE;

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

int getCurvatureSign(SYMBOL *symbol){
    int i, j;
    int numerator = 2;
    int curvature = 0;
    
    for(j = 0; j < DIMENSION; j++){
        for(i = 0; i < symbol->orbitCount[j]; i++){
            numerator = lcm(numerator, symbol->orbitM[i][j]);
        }
    }
    
    for(i = 0; i < symbol->order; i++){
        curvature += numerator/symbol->m[i][0] + numerator/symbol->m[i][1] - numerator/2;
    }
    
    return curvature;
}

/* Handle a symbol once it has been read.
 */
boolean acceptSymbol(SYMBOL *symbol){
    if(maximumOrder < symbol->order){
        return FALSE;
    }
    if(minimumOrder > symbol->order){
        return FALSE;
    }
    if(maximumNumberFaceOrbits < symbol->orbitCount[0]){
        return FALSE;
    }
    if(minimumNumberFaceOrbits > symbol->orbitCount[0]){
        return FALSE;
    }
    if(maximumNumberVertexOrbits < symbol->orbitCount[1]){
        return FALSE;
    }
    if(minimumNumberVertexOrbits > symbol->orbitCount[1]){
        return FALSE;
    }
    
    int i;
    for(i = 0; i < symbol->orbitCount[0]; i++){
        if(symbol->orbitM[i][0] > maximumFaceSize){
            return FALSE;
        }
        if(symbol->orbitM[i][0] < minimumFaceSize){
            return FALSE;
        }
    }
    for(i = 0; i < symbol->orbitCount[1]; i++){
        if(symbol->orbitM[i][1] > maximumVertexDegree){
            return FALSE;
        }
        if(symbol->orbitM[i][1] < minimumVertexDegree){
            return FALSE;
        }
    }
    
    if(!allownonplain){
        if(getCurvatureSign(symbol)){
            return FALSE;
        }
    }
    
    return TRUE;
}
/*
print a usage message. name is the name of the current program.
 */
void usage(char *name) {
    fprintf(stderr, "Usage: %s [options]\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

/*
print a help message. name is the name of the current program.
 */
void help(char *name) {
    fprintf(stderr, "The program %s reads a list of Delaney-Dress symbols and filters it based upon given restrictions.\n", name);
    fprintf(stderr, "Usage: %s [options]\n\n", name);
    fprintf(stderr, "Valid options:\n");
    fprintf(stderr, "  -h    : Print this help and return.\n");
    fprintf(stderr, "  -n n  : Only allow symbols of order at least n.\n");
    fprintf(stderr, "  -N n  : Only allow symbols of order at most n.\n");
    fprintf(stderr, "  -s n  : Only allow symbols with faces of size at least n.\n");
    fprintf(stderr, "  -S n  : Only allow symbols with faces of size at most n.\n");
    fprintf(stderr, "  -d n  : Only allow symbols with vertices of degree at least n.\n");
    fprintf(stderr, "  -D n  : Only allow symbols with vertices of degree at most n.\n");
    fprintf(stderr, "  -f n  : Only allow symbols with at least n face orbits.\n");
    fprintf(stderr, "  -F n  : Only allow symbols with at most n face orbits.\n");
    fprintf(stderr, "  -v n  : Only allow symbols with at least n vertex orbits.\n");
    fprintf(stderr, "  -V n  : Only allow symbols with at most n vertex orbits.\n");
    fprintf(stderr, "  -p    : Only allow symbols with curvature 0.\n");
}

int main(int argc, char *argv[]){
    
    int c;
    char *name = argv[0];
    

    while ((c = getopt(argc, argv, "n:N:f:F:v:V:ph")) != -1) {
        switch (c) {
            case 'n':
                minimumOrder = atoi(optarg);
                break;
            case 'N':
                maximumOrder = atoi(optarg);
                break;
            case 's':
                minimumFaceSize = atoi(optarg);
                break;
            case 'S':
                maximumFaceSize = atoi(optarg);
                break;
            case 'd':
                minimumVertexDegree = atoi(optarg);
                break;
            case 'D':
                maximumVertexDegree = atoi(optarg);
                break;
            case 'f':
                minimumNumberFaceOrbits = atoi(optarg);
                break;
            case 'F':
                maximumNumberFaceOrbits = atoi(optarg);
                break;
            case 'v':
                minimumNumberVertexOrbits = atoi(optarg);
                break;
            case 'V':
                maximumNumberVertexOrbits = atoi(optarg);
                break;
            case 'p':
                allownonplain = FALSE;
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
    int structuresWritten = 0;
    
    while(readSymbol(stdin, &s, &symbolString)!=EOF){
        structuresRead++;
        if(acceptSymbol(&s)){
            fprintf(stdout, "%s", symbolString);
            structuresWritten++;
        }
    }
    fprintf(stderr, "Read %d symbol(s). Written %d symbol(s).\n\n", structuresRead, structuresWritten);
    
    return EXIT_SUCCESS;
}

