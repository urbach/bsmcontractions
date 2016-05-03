# ##########
srcdir=.

# Makefile


# ##########


CXX=g++
CXXFLAGS=-g -Wall -O3 -DF_
# CXXFLAGS=-O3

CCDEP = gcc 
DEPFLAGS = -MM

INCLUDE = -I/home/urbach/daten/workdir/head/lime-1.3.2/include/ -I..

#LIBS = -llime  -llapack -lblas -lm -lg2c
LIBS = -llime -llapack -lblas -lm -lgfortran
LDFLAGS = -L/home/urbach/daten/workdir/head/lime-1.3.2/lib/
LINK = $(CXX) -o $@ ${LDFLAGS}
COMPILE = ${CXX} $(INCLUDE) -o $@ ${CXXFLAGS}

MODULES = DML_crc32 dml fields io io_utils propagator_io contract_twopoint

PROGRAM = bsm_conn

all: dep bsm_conn


# ##########

$(addsuffix .d,$(MODULES)): %.d: ${srcdir}/%.cc Makefile
	  $(CCDEP) ${DEPFLAGS} ${INCLUDE} $< > $@

$(addsuffix .d,$(PROGRAM)): %.d: %.cc Makefile
	 @ $(CCDEP) ${DEPFLAGS} ${INCLUDE} $< > $@

dep: $(addsuffix .d,$(MODULES) ${PROGRAM})

$(addsuffix .o,${MODULES}): %.o: ${srcdir}/%.cc %.d Makefile
	${COMPILE} ${INCLUDE} ${OPTARGS} -c $< 

$(addsuffix .o,${PROGRAM}): %.o: %.cc %.d Makefile
	${COMPILE} ${INCLUDE} ${OPTARGS} -c $< 

${PROGRAM}: %: %.o $(addsuffix .o,${MODULES}) Makefile
	${LINK}  $(addsuffix .o,${MODULES}) $@.o $(LIBS)

# ##########


clean:
	rm -f *~ *.o *.d bsm_conn

.PHONY: clean

# ##########
