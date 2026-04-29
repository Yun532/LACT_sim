# Makefile for hessio hdata2hbook converter.
# Author: Konrad Bernloehr

# This Makefile requires GNU make.

PWD := .
# include $(HESSROOT)/makefiles/Makefile.Common
HBOOK_OBJ=.
BINOUT=../bin

CERNLIB_COMPILER := $(shell ./cernlib_compiler)
# $(warning CERNLIB_COMPILER = [$(CERNLIB_COMPILER)])
ifeq ($(CERNLIB_COMPILER),gfortran)
 ifneq ($(shell which gfortran 2>/dev/null),)
   F77 = gfortran -std=legacy
 endif
endif
ifeq ($(CERNLIB_COMPILER),g77)
   override F77 = g77
endif
# $(warning F77 = [$(F77)])

OBJS = $(HBOOK_OBJ)/cvt2.o $(HBOOK_OBJ)/tohbook.o $(HBOOK_OBJ)/hbookc.o
INCLUDES := -I../include
CFLAGS += $(INCLUDES)

ifeq ($(MACHINE),x86_64)
override FFLAGS += -fno-f2c
endif

CERNLIB_SAFE := -safe # Assume we have shared CERNlib libraries
CERN_PACKLIB := $(shell cernlib $(CERNLIB_SAFE) packlib 2>/dev/null)
ifeq ($(CERN_PACKLIB),)
ifneq ($(wildcard /cern/pro/lib),)
CERN_PACKLIB := -L/cern/pro/lib -lpacklib
endif
endif
ifeq ($(CERNLIB_COMPILER),gfortran)
   LIBS += -lgfortran -lnsl -lm -lc
else
   LIBS += -lnsl -lm -lc
endif
ifneq ($(findstring -D_REENTRANT,$(CFLAGS)),)
   LIBS += -lpthread
endif

REAL_F77 := $(shell which $(F77) 2>/dev/null)

ifeq ($(SYSTEM),Darwin)
   SOLIBEXT = dylib
else
   SOLIBEXT = so
endif

# If we want to link against hessio, check if our code would
# match the version with or without multithreading support.
ifneq ($(wildcard ../lib/libhessio*.$(ALIBEXT)),)
 ifeq ($(findstring -D_REENTRANT,$(CXXFLAGS)),)
    HESSIOLIB  := ../lib/libhessio.$(ALIBEXT)
 else
    HESSIOLIB  := ../lib/libhessio_r.$(ALIBEXT)
 endif
else
 ifeq ($(findstring -D_REENTRANT,$(CFLAGS)),)
    HESSIOLIB  := ../lib/libhessio.$(SOLIBEXT)
 else
    HESSIOLIB  := ../lib/libhessio_r.$(SOLIBEXT)
 endif
endif

ifeq ($(REAL_F77),)
nothing:
	@echo '*** Cannot build hdata2hbook. No FORTRAN compiler found.'
else
ifeq ($(CERN_PACKLIB),)
nothing:
	@echo '*** Cannot build hdata2hbook. No CERN libraries found.'
endif
endif

all: Makefile.depend $(BINOUT)/hdata2hbook

# Default rules which work with gcc/g77 and DEC cc/f77

$(HBOOK_OBJ)/%.o: %.f
	$(F77) ${FFLAGS} -c -o $@ $<
$(HBOOK_OBJ)/%.o: %.F
	$(F77) ${FFLAGS} -c -o $@ $<
$(HBOOK_OBJ)/%.o: %.c
	$(CC) ${CFLAGS} -c -o $@ $<

$(BINOUT)/cvt2: $(OBJS) $(HESSIOLIB)
	$(F77) $(LDFLAGS) $^ $(CERN_PACKLIB) $(LIBS) -o $@

$(BINOUT)/hdata2hbook: $(BINOUT)/cvt2
	cp -p $(BINOUT)/cvt2 $(BINOUT)/hdata2hbook

clean:
	rm -f $(HBOOK_OBJ)/*.o $(BINOUT)/cvt2 $(BINOUT)/hdata2hbook

distclean: clean
	rm -f Makefile.depend

depend:
	-rm -f Makefile.depend
	$(MAKE) Makefile.depend

.PHONY: clean distclean depend

Makefile.depend: Makefile
	-rm -f Makefile.depend
	@echo "# DO NOT DELETE THIS LINE -- make  depend  depends  on  it." >>Makefile.depend
ifeq ($(USE_MAKEDEPEND),Y)
	makedepend -Y -fMakefile.depend -- $(CFLAGS) -- \
           cvt2.c tohbook.c hbookc.F > /dev/null 2>&1
else
	for file in cvt2.c tohbook.c hbookc.F; do \
   $(CC) -MM $(CFLAGS) -I../include $$file \
   -MT `echo $${file} | sed 's/\(.*\)\.\(c\|cc\|cxx\|C\|F\)$$/\1.o/'` \
   >> Makefile.depend; done; exit 0
endif

ifeq ($(wildcard Makefile.depend),Makefile.depend)
include Makefile.depend
endif

