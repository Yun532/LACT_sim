# Makefile for the hessio library and tools as used with the
# exported 'hessioxxx' package (differs from the in-CVS-tree
# Makefile for hessio).

PROGRAMS := read_hess_nr listio testio
# Optional programs (C API).
ifneq ($(shell which root 2>/dev/null),)
   PROGRAMS += hdata2root
endif
ifneq ($(wildcard src/read_hess.c),)
   PROGRAMS += read_hess gen_lookup
endif
ifneq ($(wildcard src/read_iact.c),)
   PROGRAMS += read_iact
endif
ifneq ($(wildcard src/select_iact.c),)
   PROGRAMS += select_iact
endif
ifneq ($(wildcard src/extract_hess.c),)
   PROGRAMS += extract_hess
endif
ifneq ($(wildcard src/extract_simtel.c),)
   PROGRAMS += extract_simtel
endif
ifneq ($(wildcard src/split_hessio.c),)
   PROGRAMS += split_hessio
endif
ifneq ($(wildcard src/gen_trgmask.c),)
   PROGRAMS += gen_trgmask
endif
ifneq ($(wildcard src/check_trgmask.c),)
   PROGRAMS += check_trgmask
endif
ifneq ($(wildcard src/merge_simtel.c),)
   PROGRAMS += merge_simtel
endif
ifneq ($(wildcard src/fcat.c),)
   PROGRAMS += fcat
endif
ifneq ($(wildcard src/list_histograms.c),)
   PROGRAMS += list_histograms
endif
ifneq ($(wildcard src/add_histograms.c),)
   PROGRAMS += add_histograms
endif

SYSTEM := $(shell uname)
MACHINE:= $(shell uname -m)
WITH_GCC := $(shell ($(CC) --version 2>/dev/null | fgrep GCC) && echo 1 || echo 0)
ifeq ($(WITH_GCC),1)
   GCC_VERSION := $(shell $(CC) -dumpversion 2>/dev/null)
endif

# Optional programs (C++ API).
ifeq ($(NOCXX),)
ifneq ($(wildcard src/TestIO.cc),)
 ifeq ($(SYSTEM),Darwin)
# On Mac OS X with HFS+ filesystem we must avoid a clash between testio and TestIO.
   PROGRAMS += TestIO++
 else
   PROGRAMS += TestIO
 endif
endif
ifneq ($(wildcard src/statio.cc),)
   PROGRAMS += statio
endif
ifneq ($(wildcard src/filterio.cc),)
   PROGRAMS += filterio
endif
endif
LIBDIR := $(shell /bin/pwd)/lib

ifneq ($(RPATHDIR),)
   WITH_RPATH=1
else
   RPATHDIR=$(LIBDIR)
endif

ifeq ($(SYSTEM),Darwin)
$(warning Support for Mac OS X is experimental.)
$(info You need to export DYLD_LIBRARY_PATH=$(LIBDIR))
   SOLIBEXT = dylib
   SOLIBFLAGS = -dynamiclib -single_module \
    -flat_namespace -Wl,-install_name -Wl,$(LIBDIR)/lib$(1).$(SOLIBEXT) \
    -multiply_defined suppress
   # SOEXEFLAGS = -Wl,-rpath,$(LIBDIR) ## Apparently not working with Mac OS X
else
   SOCOMPILE = -fPIC -DPIC
   SOLIBEXT = so
   SOLIBFLAGS = -shared -Wl,-soname -Wl,lib$(1).$(SOLIBEXT) 
   ifneq ($(WITH_RPATH),)
      SOEXEFLAGS =  -Wl,-rpath -Wl,$(RPATHDIR)
   endif
endif
lib_expand = $(patsubst %,lib%.$(SOLIBEXT),$(1))

# Further libraries needed.
HESSIO_LIB=-L$(LIBDIR) -lhessio
HESSIOPP_LIB=-L$(LIBDIR) -lhessio++

# Listing of all the files used for exports
LIBHESSIO_SOURCES_NOPATH := \
    atmprof.c \
    atmprof.h \
    current.c \
    current.h \
    dhsort.c \
    dhsort.h \
    eventio.c \
    eventio_registry.c \
    eventio_registry.h \
    fileopen.c \
    fileopen.h \
    histogram.c \
    histogram.h \
    hconfig.c \
    hconfig.h \
    initial.h \
    moments.c \
    io_basic.h \
    io_histogram.c \
    io_histogram.h \
    io_history.c \
    io_history.h \
    io_simtel.c \
    mc_atmprof.c \
    mc_atmprof.h \
    mc_tel.h \
    straux.c \
    straux.h \
    warning.c \
    warning.h \
    io_hess.c \
    io_hess.h
# Listing of all the files used for exports
LIBHESSIOPP_SOURCES_NOPATH := \
    EventIO.cc \
    EventIO.hh

LIBHESSIO_SOURCES := $(patsubst %.c,src/%.c,$(patsubst %.h,include/%.h,$(LIBHESSIO_SOURCES_NOPATH)))
LIBHESSIO_LO := $(patsubst %.c,out/%.lo,$(filter %.c,$(LIBHESSIO_SOURCES_NOPATH)))

LIBHESSIOPP_SOURCES := $(patsubst %.cc,src/%.cc,$(patsubst %.hh,include/%.hh,$(LIBHESSIOPP_SOURCES_NOPATH)))
ifeq ($(SYSTEM),Darwin)
# On Mac OS X with HFS+ filesystem we must avoid a clash between eventio.lo and EventIO.lo.
LIBHESSIOPP_LO := $(subst EventIO.lo,EventIO++.lo,$(patsubst %.cc,out/%.lo,$(filter %.cc,$(LIBHESSIOPP_SOURCES_NOPATH))))
else
LIBHESSIOPP_LO := $(patsubst %.cc,out/%.lo,$(filter %.cc,$(LIBHESSIOPP_SOURCES_NOPATH)))
endif

DEFINES = -DCTA $(EXTRA_DEFINES)
CDEBUGFLAGS = -g
ifeq ($(WITH_GCC),1)
 ifneq ($(shell echo $(GCC_VERSION) | egrep '4\.[2-9]\.'),)
  OPT_mach := -march=native -mtune=native
 endif
endif
OPTIMIZE = -O2 $(OPT_mach) # Enable optimization when compiling for production.
OPTIONS = -Wall -Wmissing-prototypes -Wstrict-prototypes -Wshadow \
   -Wunused -Wuninitialized -Wnested-externs -Wcast-qual \
   -Wcast-align $(EXTRA_OPTIONS)
NOTCXXOPTIONS = -Wnested-externs -Wmissing-prototypes -Wstrict-prototypes

# Add option for transition support to 64 bit file access
# (files > 2 GB on 32 bit system), as supported by the eventio code. 
# Opening large files may need fopen64() instead of fopen().
# For automatic 64 bit support use '-D_FILE_OFFSET_BITS=64'.
# Then no change to fopen() calls is required.
LFS_OPTION = -D_LARGEFILE64_SOURCE
# LFS_OPTION = -D_FILE_OFFSET_BITS=64
OPTIONS += $(LFS_OPTION)

CFLAGS = -Iinclude $(DEFINES) $(CDEBUGFLAGS) $(OPTIONS) $(OPTIMIZE)
CXXFLAGS = -Iinclude $(DEFINES) $(CDEBUGFLAGS) \
   $(filter-out $(NOTCXXOPTIONS),$(OPTIONS)) $(OPTIMIZE)

ifeq ($(findstring -O,$(CXXFLAGS)),)
 ifeq ($(findstring -O1,$(CXXFLAGS)),)
  ifeq ($(findstring -O2,$(CXXFLAGS)),)
   ifeq ($(findstring -O3,$(CXXFLAGS)),)
     CFLAGS := $(CFLAGS:-Wuninitialized=)
     CXXFLAGS := $(CXXFLAGS:-Wuninitialized=)
   endif
  endif
 endif
endif

LDFLAGS = $(LDEBUGFLAGS) $(LDEXTRA)

.PHONY: all static semistatic install

all: $(patsubst %,bin/%,$(PROGRAMS))

static:
	$(MAKE) LINK_STATIC=1 WITH_STATIC_LIBS=1 -f Makefile.static

semistatic:
	$(MAKE) WITH_STATIC_LIBS=1 -f Makefile.static

install: all
	@echo "No install step being used here."

out/%.o: src/%.c
	$(CC) ${CFLAGS} -c -o $@ $<
out/%.o: src/%.cc
	$(CXX) ${CXXFLAGS} -c -o $@ $<
out/%++.o: src/%.cc
	$(CXX) ${CXXFLAGS} -c -o $@ $<
out/%.lo:   src/%.c
	$(CC) ${CFLAGS} $(SOCOMPILE) -c -o $@ $<
out/%.lo:   src/%.cc
	$(CXX) ${CXXFLAGS} $(SOCOMPILE) -c -o $@ $<
out/%++.lo:  src/%.cc
	$(CXX) ${CXXFLAGS} $(SOCOMPILE) -c -o $@ $<

lib/$(call lib_expand,hessio): $(LIBHESSIO_LO)
	$(CC) $(call SOLIBFLAGS,hessio) $(LIBHESSIO_LO) -o $@
lib/$(call lib_expand,hessio++): $(LIBHESSIOPP_LO) $(LIBHESSIO_LO)
	$(CXX) $(call SOLIBFLAGS,hessio++) $(LIBHESSIOPP_LO) $(LIBHESSIO_LO) -o $@

bin/read_iact:  out/read_iact.o out/fileopen.o  \
   out/io_simtel.o out/mc_atmprof.o \
   out/io_history.o out/current.o \
   out/eventio.o out/straux.o out/warning.o
	$(CC) $(LDFLAGS) $(SOEXEFLAGS) $^ -lm -o $@

bin/select_iact:  out/select_iact.o out/fileopen.o  \
   out/io_simtel.o out/mc_atmprof.o \
   out/eventio.o out/straux.o out/warning.o
	$(CC) $(LDFLAGS) $(SOEXEFLAGS) $^ -lm -o $@

bin/read_hess_nr: out/read_hess_nr.o out/rec_tools_nr.o \
           out/camera_image.o \
           lib/$(call lib_expand,hessio)
	$(CC) $(LDFLAGS) $(SOEXEFLAGS) \
           $(filter %.o,$^) \
           $(HESSIO_LIB) -lm  \
           -o $@
	-(cd `dirname $@` && rm -f read_simtel_nr && ln -s read_hess_nr read_simtel_nr)
	-(cd `dirname $@` && rm -f read_cta_nr && ln -s read_hess_nr read_cta_nr)

bin/read_hess: out/read_hess.o out/rec_tools.o out/reconstruct.o \
           out/camera_image.o out/user_analysis.o out/basic_ntuple.o \
           out/io_trgmask.o \
           lib/$(call lib_expand,hessio)
	$(CC) $(LDFLAGS) $(SOEXEFLAGS) \
           $(filter %.o,$^) \
           $(HESSIO_LIB) -lm  \
           -o $@
	-(cd `dirname $@` && rm -f read_simtel && ln -sf read_hess read_simtel)
	-(cd `dirname $@` && rm -f read_cta && ln -sf read_hess read_cta)

bin/merge_simtel: out/merge_simtel.o  \
           out/basic_ntuple.o out/io_trgmask.o \
           lib/$(call lib_expand,hessio)
	$(CC) $(LDFLAGS) $(SOEXEFLAGS) $(filter %.o,$^) \
           $(HESSIO_LIB) -lm  \
           -o $@

bin/split_hessio: out/split_hessio.o  \
           lib/$(call lib_expand,hessio)
	$(CC) $(LDFLAGS) $(SOEXEFLAGS) $(filter %.o,$^) \
           $(HESSIO_LIB) -lm  \
           -o $@

bin/best_of: src/best_of.cc lib/$(call lib_expand,hessio)
	if [ ! -L include/hessio ]; then (cd include && ln -s . hessio); fi
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(SOEXEFLAGS) \
           $(filter %.cc,$^) \
           $(HESSIO_LIB) -lm  \
           -o $@

# Use different name mainly for systems without proper handling of mixed-case filenames.
bin/TestIO++: src/TestIO.cc lib/$(call lib_expand,hessio) lib/$(call lib_expand,hessio++)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(SOEXEFLAGS) \
           $(filter %.cc,$^) \
           $(HESSIOPP_LIB) $(HESSIO_LIB) -lm \
           -o $@
ifneq ($(SYSTEM),Darwin)
# Do not build under original file name for Mac OS X as it overwrites testio.
bin/TestIO: src/TestIO.cc lib/$(call lib_expand,hessio) lib/$(call lib_expand,hessio++)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(SOEXEFLAGS) \
           $(filter %.cc,$^) \
           $(HESSIOPP_LIB) $(HESSIO_LIB) -lm \
           -o $@
endif

bin/statio: src/statio.cc lib/$(call lib_expand,hessio) lib/$(call lib_expand,hessio++)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(SOEXEFLAGS) \
           $(filter %.cc,$^) \
           $(HESSIOPP_LIB) $(HESSIO_LIB) -lm \
           -o $@

bin/filterio: src/filterio.cc lib/$(call lib_expand,hessio) lib/$(call lib_expand,hessio++)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(SOEXEFLAGS) \
           $(filter %.cc,$^) \
           $(HESSIOPP_LIB) $(HESSIO_LIB) -lm \
           -o $@

bin/gen_trgmask: out/gen_trgmask.o out/io_trgmask.o lib/$(call lib_expand,hessio)
	$(CC) $(LDFLAGS) $(SOEXEFLAGS) \
           $(filter %.o,$^) \
           $(HESSIO_LIB) -lm  \
           -o $@

bin/check_trgmask: out/check_trgmask.o out/io_trgmask.o lib/$(call lib_expand,hessio)
	$(CC) $(LDFLAGS) $(SOEXEFLAGS) \
           $(filter %.o,$^) \
           $(HESSIO_LIB) -lm  \
           -o $@

ifneq ($(EXTRAPROG),)
bin/$(EXTRAPROG): out/$(EXTRAPROG).o lib/$(call lib_expand,hessio)
	$(CC) $(LDFLAGS) $(SOEXEFLAGS) \
           $(filter %.o,$^) \
           $(HESSIO_LIB) -lm  \
           -o $@
endif

src/show_prepro_c.c: src/show_prepro.cc
	(cd src && ln -s show_prepro.cc show_prepro_c.c)

bin/show_types: out/show_types.o
	$(CC) $(LDFLAGS) $(SOEXEFLAGS) $^ -lm -o $@

bin/show_prepro: out/show_prepro.o
	$(CC) $(LDFLAGS) $(SOEXEFLAGS) $^ -lm -o $@

bin/show_prepro_c: out/show_prepro_c.o
	$(CC) $(LDFLAGS) $(SOEXEFLAGS) $^ -lm -o $@

bin/hdata2hbook:
	(cd hbook && make)

bin/extract_hess: out/extract_hess.o lib/$(call lib_expand,hessio)
	$(CC) $(LDFLAGS) $(SOEXEFLAGS) \
           $(filter %.o,$^) \
           $(HESSIO_LIB) -lm  \
           -o $@
	-(cd `dirname $@` && ln -s extract_hess extract_calibevent)

bin/extract_simtel: out/extract_simtel.o lib/$(call lib_expand,hessio) \
    out/io_trgmask.o
	$(CC) $(LDFLAGS) $(SOEXEFLAGS) \
           $(filter %.o,$^) \
           $(HESSIO_LIB) -lm  \
           -o $@

bin/%: out/%.o lib/$(call lib_expand,hessio)
	$(CC) $(LDFLAGS) $(SOEXEFLAGS) \
           $(filter %.o,$^) \
           $(HESSIO_LIB) -lm  \
           -o $@

bin/hdata2root: root/hdata2root.cc root/toroot.cc root/toroot.hh lib/$(call lib_expand,hessio)
	-@$(MAKE) -C root CXX="$(CXX)" CXXFLAGS="$(patsubst -I%,-I../%,$(CXXFLAGS))" \
LDEXTRA="$(LDFLAGS)" LIBDIR=$(LIBDIR) BINOUT=../bin ROOT_OBJ=../out

# Alternate build through cmake: very limited configuration possibilities
# but may (or may not) work on architectures not directly supported.

.PHONY: cmake
cmake:
	if [ ! -d build-cmake ]; then mkdir build-cmake; fi; \
	cd build-cmake || exit 1; \
	../cmake-run.sh $(EXTRA_DEFINES)

# Update the version whenever one of the sources changed:

include/eventio_version.h: $(LIBHESSIO_SOURCES)
	egrep 'RELEASE|ifdef|undef|endif' $@ | sed 's/ EVENTIO_RELEASE/ EVENTIO_BASE_RELEASE/' > $@.tmp
	echo "#define EVENTIO_VERSION \"`date '+%F %T %Z'` (`whoami`@`hostname`)\"" > $@
	cat $@.tmp >>$@


.PHONY: clean distclean depend
clean:
	-rm -f bin/* out/* lib/* test.ps test.root *.bck */*.bck hbook/*.o root/*.o

distclean: clean
	-rm -f Makefile.depend

depend:
	-rm -f Makefile.depend
	$(MAKE) Makefile.depend
#	(cd $(HBOOK_SRC) && $(MAKE) depend)
#	(cd $(ROOT_SRC) && $(MAKE) depend)

Makefile.depend: Makefile $(ALL_SOURCES)
	-rm -f Makefile.depend
	@echo "# DO NOT DELETE THIS LINE -- make  depend  depends  on  it." >>Makefile.depend
	for prg in $(PROGRAMS); do g="$$(/bin/ls src/$$prg{.c,.cc} 2>/dev/null)"; \
   if [ ! -z "$$g" ]; then $(CXX) -MM -MT $$prg $(CXXFLAGS) -Iinclude $$g \
   >> Makefile.depend; fi done; exit 0
	for file in src/*.c; do \
   $(CC) -MM -MT out/$$(basename $$file .c)".o" $(CXXFLAGS) -Iinclude $$file \
   >> Makefile.depend; done; exit 0
	for file in src/*.cc; do \
   $(CXX) -MM -MT out/$$(basename $$file .cc)".o" $(CXXFLAGS) -Iinclude $$file \
   >> Makefile.depend; done; exit 0

ifeq ($(wildcard Makefile.depend),Makefile.depend)
include Makefile.depend
endif
