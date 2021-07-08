SOURCES=src/MFTAnaSim.cxx src/MFTAnaSimTrack.cxx src/MFTAnaSimSATrack.cxx

HEADERS=include/MFTAnaSim.h include/MFTAnaSimHit.h include/MFTAnaSimCluster.h include/MFTAnaSimTrack.h include/MFTAnaSimMCTrack.h include/MFTAnaSimSATrack.h

SW_ROOT=${ALIBUILD_WORK_DIR}/`aliBuild architecture`

#CXXFLAGS=-O2 -Wall -fPIC -pthread -std=c++17 -m64 -I${SW_ROOT}/ROOT/v6-20-08-alice1-6/include
CXXFLAGS=-Wall -fPIC -std=c++17 -I${SW_ROOT}/ROOT/latest/include

CXXFLAGS+= -g
CXXFLAGS+= -fopenmp

all: libMFTAnaSim.so

DictMFTAnaSim.cxx: $(HEADERS) src/MFTAnaSimLinkDef.h
	rootcling -f $@ $(CXXFLAGS) -I${ROOT_INCLUDE_PATH} $^

libMFTAnaSim.so: DictMFTAnaSim.cxx $(SOURCES)
	g++ -shared -o $@ $(CXXFLAGS) -I${ROOTSYS}/include -I${O2_ROOT}/include -I${O2_ROOT}/include/GPU -I${SW_ROOT}/ms_gsl/latest/include -I${SW_ROOT}/FairLogger/latest/include -I${SW_ROOT}/FairLogger/latest/include/fairlogger -I${SW_ROOT}/fmt/latest/include -I${SW_ROOT}/FairRoot/latest/include -I${SW_ROOT}/boost/latest/include $^

install:
	mkdir -p ${ALIBUILD_WORK_DIR}/MFTAna
	cp libMFTAnaSim.* DictMFTAnaSim.cxx DictMFTAnaSim_rdict.pcm ${ALIBUILD_WORK_DIR}/MFTAna

clean:
	echo $(CXXFLAGS)
	rm -f libMFTAnaSim.* DictMFTAnaSim.cxx DictMFTAnaSim_rdict.pcm
