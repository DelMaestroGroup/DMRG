LIBRARY_DIR=$(HOME)/local/src/ITensor

APP=dmrgSacc

# 4. Add any headers your program depends on here. The make program
#    will auto-detect if these headers have changed and recompile your app.
HEADERS=

# 5. For any additional .cc (source) files making up your project,
#    add their full filenames here.
CCFILES=$(APP).cc

#################################################################
#################################################################
#################################################################


include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

CCFLAGS+= -I$(HOME)/local/include/
CCFLAGS+= -I/usr/local/include/

CCGFLAGS+= -I$(HOME)/local/include/
CCGFLAGS+= -I/usr/local/include/

TENSOR_HEADERS=$(LIBRARY_DIR)/itensor/core.h

#Mappings --------------
OBJECTS=$(patsubst %.cc,%.o, $(CCFILES))
GOBJECTS=$(patsubst %,.debug_objs/%, $(OBJECTS))

#Rules ------------------

%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

build: $(APP)
debug: $(APP)-g

$(APP): $(OBJECTS) $(ITENSOR_LIBS)
	$(CCCOM) $(CCFLAGS) $(OBJECTS) -o $(APP) $(LIBFLAGS)

$(APP)-g: mkdebugdir $(GOBJECTS) $(ITENSOR_GLIBS)
	$(CCCOM) $(CCGFLAGS) $(GOBJECTS) -o $(APP)-g $(LIBGFLAGS)

clean:
	rm -fr .debug_objs *.o $(APP)m $(APP)-g

mkdebugdir:
	mkdir -p .debug_objs

