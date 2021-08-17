# Change this to your local version of ITensor
LIBRARY_DIR=../../itensor

ifdef app
APP=$(app)
else
APP=cDMRG
endif

CCFILES=$(APP).cc cio.cc readvec.cc

#################################################################
#################################################################
#################################################################
#################################################################

include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

TENSOR_HEADERS=$(LIBRARY_DIR)/itensor/all.h
CCFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(CPPFLAGS) $(OPTIMIZATIONS)
CCGFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(DEBUGFLAGS)
LIBFLAGS=-L$(ITENSOR_LIBDIR) $(ITENSOR_LIBFLAGS)
LIBGFLAGS=-L$(ITENSOR_LIBDIR) $(ITENSOR_LIBGFLAGS)

# Mappings

OBJECTS=$(patsubst %.cc,%.o, $(CCFILES))
GOBJECTS=$(patsubst %,.debug_objs/%, $(OBJECTS))

# Rules

%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

# Targets

.PHONY: build
build: $(APP) 

.PHONY: debug
debug: $(APP)-g

$(APP): $(OBJECTS) $(ITENSOR_LIBS)
	$(CCCOM) $(CCFLAGS) $(OBJECTS) -o $(APP) $(LIBFLAGS)

$(APP)-g: mkdebugdir $(GOBJECTS) $(ITENSOR_GLIBS)
	$(CCCOM) $(CCGFLAGS) $(GOBJECTS) -o $(APP)-g $(LIBGFLAGS)

cDMRG.o: cSite.h cio.h readvec.h cDMRGeps.h localmpoproj.h

cio.o: cio.h

readvec.o: readvec.h

.PHONY: clean
clean:
	rm -fr .debug_objs *.o $(APP) $(APP)-g

.PHONY: mkdebugdir
mkdebugdir:
	mkdir -p .debug_objs
