SOURCES=src/MFTAnaSim.cxx src/MFTAnaSimTrack.cxx
HEADERS=include/MFTAnaSim.h include/MFTAnaSimTrack.h 

SW_ROOT=${ALIBUILD_WORK_DIR}/ubuntu2004_x86-64

#CXXFLAGS=-O2 -Wall -fPIC -pthread -std=c++17 -m64 -I${SW_ROOT}/ROOT/latest/include/  -Dcxx17=On
CXXFLAGS=-Wall -fPIC -std=c++17 -I${SW_ROOT}/ROOT/latest/include/ 


#all: libMFTAnaSim.so

DictMFTAnaSim.cxx: $(HEADERS) src/MFTAnaSimLinkDef.h
	rootcling -f $@ $(CXXFLAGS) -I${ROOT_INCLUDE_PATH} $^

#libMFTAnaSim.so: DictMFTAnaSim.cxx $(SOURCES)
#	clang++ -shared -o $@ $(CXXFLAGS) -I$(ROOTSYS)/include -I${O2_ROOT}/include -I${O2_ROOT}/include/GPU -I${SW_ROOT}/ms_gsl/latest/include -I${SW_ROOT}/FairLogger/latest/include -I${SW_ROOT}/fmt/latest/include -I${SW_ROOT}/FairRoot/latest/include $^

clean:
	echo $(CXXFLAGS)
	rm -f libMFTAnaSim.so DictMFTAnaSim.cxx DictMFTAnaSim_rdict.pcm
