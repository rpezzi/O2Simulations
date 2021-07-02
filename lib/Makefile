SOURCES=src/MFTAnaSim.cxx src/MFTAnaSimTrack.cxx src/MFTAnaSimSATrack.cxx

HEADERS=include/MFTAnaSim.h include/MFTAnaSimHit.h include/MFTAnaSimCluster.h include/MFTAnaSimTrack.h include/MFTAnaSimMCTrack.h include/MFTAnaSimSATrack.h 

SW_ROOT=/home/vulpescu/alice/sw/ubuntu1804_x86-64

#CXXFLAGS=-O2 -Wall -fPIC -pthread -std=c++17 -m64 -I${SW_ROOT}/ROOT/v6-20-08-alice1-6/include
CXXFLAGS=-Wall -fPIC -std=c++17 -I${SW_ROOT}/ROOT/v6-20-08-alice1-6/include

CXXFLAGS+= -g
CXXFLAGS+= -fopenmp

all: libMFTAnaSim.so

DictMFTAnaSim.cxx: $(HEADERS) src/MFTAnaSimLinkDef.h
	rootcling -f $@ $(CXXFLAGS) -I${ROOT_INCLUDE_PATH} $^

libMFTAnaSim.so: DictMFTAnaSim.cxx $(SOURCES)
	g++ -shared -o $@ $(CXXFLAGS) -I$(ROOTSYS)/include -I${O2_ROOT}/include -I${O2_ROOT}/include/GPU -I${SW_ROOT}/ms_gsl/2.1.0-1/include -I${SW_ROOT}/FairLogger/v1.9.1-1/include -I${SW_ROOT}/fmt/7.1.0-3/include -I${SW_ROOT}/FairRoot/v18.4.1-7/include $^

clean:
	echo $(CXXFLAGS)
	rm -f libMFTAnaSim.so DictMFTAnaSim.cxx DictMFTAnaSim_rdict.pcm
