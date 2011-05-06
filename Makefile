#
# Makefile for pregraphs
#

SHELL = /bin/sh

# Compiling executing program with DWORDSIZE=32 is slightly faster, 
# but limits the order of the graphs to 32.
CC32 = gcc
CC64 = gcc
CFLAGS = -Wall -DWORDSIZE=64 -DMAXN=64 -rdynamic
COMPLETE = ddgraphs ddgraphs-64 ddgraphs-profile ddgraphs-debug
SOURCES = ddgraphs.c ddgraphs.h util.h Makefile COPYRIGHT.txt LICENSE.txt
DDGRAPHS_SOURCES = ddgraphs.c nauty.c nautil.c nausparse.c

all : 32bit

complete: $(COMPLETE)

32bit: ddgraphs

64bit : ddgraphs-64

profile : ddgraphs-profile

debug : ddgraphs-debug

ddgraphs: $(DDGRAPHS_SOURCES)
	${CC32} $(CFLAGS) $(DDGRAPHS_SOURCES) -o ddgraphs

ddgraphs-64: $(DDGRAPHS_SOURCES)
	${CC64} $(CFLAGS) $(DDGRAPHS_SOURCES) -o ddgraphs-64

ddgraphs-profile: $(DDGRAPHS_SOURCES)
	${CC32} $(CFLAGS) -pg -g $(DDGRAPHS_SOURCES) -o ddgraphs-profile

ddgraphs-debug: $(DDGRAPHS_SOURCES)
	${CC32} $(CFLAGS) -rdynamic -g $(DDGRAPHS_SOURCES) -o ddgraphs-debug

sources: ddgraphs-sources.zip ddgraphs-sources.tar.gz

ddgraphs-sources.zip: $(SOURCES)
	zip ddgraphs-sources $(SOURCES)

ddgraphs-sources.tar.gz: $(SOURCES)
	tar czf ddgraphs-sources.tar.gz $(SOURCES)

clean:
	rm -f $(COMPLETE)
