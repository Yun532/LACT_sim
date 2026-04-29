# Makefile for hessio hdata2root converter.
# Author: Konrad Bernloehr

# This Makefile requires GNU make.

PWD := .

SYSTEM := $(shell uname)
ifeq ($(SYSTEM),Darwin)
   SOLIBEXT = dylib
else
   SOLIBEXT = so
endif
lib_expand = $(patsubst %,lib%.$(SOLIBEXT),$(1))

OBJS3 = $(ROOT_OBJ)/hdata2root.o $(ROOT_OBJ)/toroot.o
override INCLUDES := -I../include
override CFLAGS := $(CFLAGS) $(INCLUDES) $(shell root-config --cflags)
override CXXFLAGS := $(CXXFLAGS) $(INCLUDES) $(shell root-config --cflags)
# ROOT itself badly fails to obey const-ness. Remove -Wcast-qual from CXXFLAGS.
override CXXFLAGS := $(CXXFLAGS:-Wcast-qual=)
override LDFLAGS := $(LDEXTRA) -Wl,-rpath -Wl,$(LIBDIR)

ROOT_LIBDIR := $(shell root-config --libdir 2>/dev/null)
LIBS += -lm -lc

HESSIOLIB  := ../lib/$(call lib_expand,hessio)

ifeq ($(ROOT_LIBDIR),)
nothing:
	@echo '*** Cannot build hdata2root. No ROOT found.'
endif
ifeq ($(ROOT_OBJ),)
   ROOT_OBJ = ../out
endif
ifeq ($(BINOUT),)
   BINOUT = ../bin
endif

all: $(BINOUT)/hdata2root

# Default rules which work with gcc/g77 and DEC cc/f77

$(ROOT_OBJ)/%.o: %.cc
	$(CXX) ${CXXFLAGS} -c -o $@ $<

$(ROOT_OBJ)/toroot.o: toroot.cc toroot.hh
	$(CXX) ${CXXFLAGS} -c -o $@ $<

$(BINOUT)/hdata2root: $(OBJS3) $(HESSIOLIB)
	$(CXX) $(LDFLAGS) $^ $(shell root-config --libs) $(LIBS) -o $@
