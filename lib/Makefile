include $(ROOTSYS)/etc/Makefile.arch

SOURCES=src/MFTAnaHits.cxx
HEADERS=include/MFTAnaHits.h

all: libMFTAnaHits.so

DictMFTAnaHits.cxx: $(HEADERS) src/MFTAnaHitsLinkDef.h
	rootcling -f $@ -c $(CXXFLAGS) -I${ROOT_INCLUDE_PATH} -p $^

libMFTAnaHits.so: DictMFTAnaHits.cxx $(SOURCES)
	g++ -shared -o$@ `root-config --ldflags` $(CXXFLAGS) -I$(ROOTSYS)/include -I${O2_ROOT}/include -I${O2_ROOT}/include/GPU $^

clean:
	rm -f libMFTAnaHits.so DictMFTAnaHits.cxx DictMFTAnaHits_rdict.pcm
