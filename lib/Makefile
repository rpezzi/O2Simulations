include $(ROOTSYS)/etc/Makefile.arch

SOURCES=src/MFTAnaSim.cxx
HEADERS=include/MFTAnaSim.h

all: libMFTAnaSim.so

DictMFTAnaSim.cxx: $(HEADERS) src/MFTAnaSimLinkDef.h
	rootcling -f $@ -c $(CXXFLAGS) -I${ROOT_INCLUDE_PATH} -p $^

libMFTAnaSim.so: DictMFTAnaSim.cxx $(SOURCES)
	g++ -shared -o$@ `root-config --ldflags` $(CXXFLAGS) -I$(ROOTSYS)/include -I${O2_ROOT}/include -I${O2_ROOT}/include/GPU $^

clean:
	rm -f libMFTAnaSim.so DictMFTAnaSim.cxx DictMFTAnaSim_rdict.pcm
