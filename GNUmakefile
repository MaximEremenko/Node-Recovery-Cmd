#
# GNUmakefile for Galois field library
#
# The default flags do *not* have the SSE instructions enabled.
# Please cd to flag_tester and run which_compile_flags.sh to see which SSE instructions
# your machine and compiler support, and which flags you should include below.

CFLAGS = -O3 -msse4 -DINTEL_SSE4 -maes -mpclmul -DINTEL_PCLMUL
LDFLAGS = -O3 -msse4 -maes -mpclmul
#CFLAGS = -O3 
#LDFLAGS = -O3 

SRCS = gf.c gf_method.c gf_wgen.c gf_w4.c gf_w8.c gf_w16.c gf_w32.c \
       gf_w64.c gf_w128.c gf_rand.c gf_general.c node_recover_cmd.cpp

HDRS = gf_complete.h gf_int.h

EXECUTABLES = node_recover_cmd

RM = /bin/rm -f

LIBOBJS = gf.o gf_method.o gf_wgen.o gf_w4.o gf_w8.o gf_w16.o gf_w32.o \
          gf_w64.o gf_w128.o gf_rand.o gf_general.o

OBJS = $(addsuffix .o, $(basename $(SRCS)))

DEFAULT = $(EXECUTABLES) gf_complete.a

default: $(DEFAULT)

all: $(OBJS)

gf_complete.a: $(LIBOBJS)
	ar ru libgf_complete.a $(LIBOBJS)
	ranlib libgf_complete.a

node_recover_cmd: node_recover_cmd.o gf_complete.a

clean:
	$(RM) $(OBJS)

spotless: clean
	$(RM) *~ $(EXECUTABLES)
	$(RM) gf_complete.a

node_recover_cmd.o: gf_complete.h
gf_wgen.o: gf_int.h gf_complete.h
gf_w4.o: gf_int.h gf_complete.h
gf_w8.o: gf_int.h gf_complete.h
gf_w16.o: gf_int.h gf_complete.h
gf_w32.o: gf_int.h gf_complete.h
gf_w64.o: gf_int.h gf_complete.h
gf_general.o: gf_complete.h gf_int.h gf_general.h gf_rand.h
gf.o: gf_complete.h gf_int.h
gf_method.o: gf_complete.h