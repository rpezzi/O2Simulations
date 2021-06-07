SOURCES=src/MFTAnaSim.cxx src/MFTAnaSimTrack.cxx
HEADERS=include/MFTAnaSim.h include/MFTAnaSimTrack.h 

#CXXFLAGS=-O2 -Wall -fPIC -pthread -std=c++17 -m64 -I/home/vulpescu/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-08-alice1-6/include
CXXFLAGS=-Wall -fPIC -pthread -std=c++17 -m64 -I/home/vulpescu/alice/sw/ubuntu1804_x86-64/ROOT/v6-20-08-alice1-6/include

all: libMFTAnaSim.so

DictMFTAnaSim.cxx: $(HEADERS) src/MFTAnaSimLinkDef.h
	rootcling -f $@ $(CXXFLAGS) -I${ROOT_INCLUDE_PATH} $^

libMFTAnaSim.so: DictMFTAnaSim.cxx $(SOURCES)
	g++ -shared -o $@ $(CXXFLAGS) -I$(ROOTSYS)/include -I${O2_ROOT}/include -I${O2_ROOT}/include/GPU -I/home/vulpescu/alice/sw/ubuntu1804_x86-64/ms_gsl/2.1.0-1/include -I/home/vulpescu/alice/sw/ubuntu1804_x86-64/FairLogger/v1.9.1-1/include -I/home/vulpescu/alice/sw/ubuntu1804_x86-64/fmt/7.1.0-3/include -I/home/vulpescu/alice//sw/ubuntu1804_x86-64/FairRoot/v18.4.1-7/include $^

clean:
	echo $(CXXFLAGS)
	rm -f libMFTAnaSim.so DictMFTAnaSim.cxx DictMFTAnaSim_rdict.pcm
