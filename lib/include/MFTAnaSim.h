#ifndef MFT_ANA_SIM
#define MFT_ANA_SIM

#include <TFile.h>
#include <TTree.h>

#include "../include/MFTAnaSimTrack.h"
#include "../include/MFTAnaSimSATrack.h"
#include "../include/MFTAnaSimCluster.h"
#include "../include/MFTAnaSimHit.h"

#include "ITSMFTSimulation/Hit.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "MFTBase/Constants.h"
#include "MFTBase/GeometryTGeo.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"

namespace o2::mftana
{

struct MCPart {
  //MCPart() = default;
  //~MCPart() = default;
  //MCPart& operator=(const MCPart&) = default;
  int mEvent = 0;
  int mPDGCode = -1;
  std::string mPDGName = "";
  int mCount = 0;
};

class MFTAnaSim
{
 public:
  MFTAnaSim();
  ~MFTAnaSim() = default;

  int mNEvents = 0;

  TTree* mHitTree = nullptr;
  TTree* mKineTree = nullptr;
  TTree* mClusTree = nullptr;
  TTree* mTrackTree = nullptr;
  
  o2::mft::GeometryTGeo* mGeoManager;
  o2::itsmft::ChipMappingMFT mChipMapper;

  bool initialize(int maxMCTracks);
  void initEvent(int event, int nMCTracks, int particleSource = 0);
  bool doParticles();
  bool doHits();
  bool doMCTracks();
  bool doSATracks();
  void finishEvent();
  void finish();
  void findMCTrackHits(int trkID, int& firstIndex, int& lastIndex);
  void findMCTrackClusters(int trkID, int& firstIndex, int& lastIndex);
  void setVerboseLevel(int vl) { mVerboseLevel = vl; }
  const std::vector<MCPart>& getParticles() { return mParticles; }

 private:
  bool mPrimary = kTRUE, mSecondary = kFALSE, mAll = kFALSE;
  enum ParticleSource {kPrimary, kSecondary, kAll};
  void filterPDGCode(int& pdgCode);
  void countParticle(int pdgCode);
  
  std::vector<o2::MCTrack> mMCTrkVec, *mMCTrkVecP = &mMCTrkVec;
  std::vector<o2::itsmft::Hit> mHitVec, *mHitVecP = &mHitVec;
  std::vector<o2::itsmft::CompClusterExt> mClusVec, *mClusVecP = &mClusVec;
  std::vector<unsigned char> mClusPatterns, *mClusPatternsP = &mClusPatterns;
  std::vector<o2::mft::TrackMFT> mTrackVec, *mTrackVecP = &mTrackVec;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* mClusLabels = nullptr;
  std::vector<o2::MCCompLabel> mTrackLabels, *mTrackLabelsInp = &mTrackLabels;
  std::vector<int> mTrackExtClsVec, *mTrackExtClsVecP = &mTrackExtClsVec;

  std::vector<std::array<bool, o2::mft::constants::DisksNumber>> mMCTrackHasHitsInDisk;
  std::vector<std::array<bool, o2::mft::constants::LayersNumber>> mMCTrackHasHitsInLayer;
  std::vector<MCPart> mParticles;
  std::vector<MFTAnaSimTrack> mAnaSimTracks;
  std::vector<MFTAnaSimCluster> mAnaSimClusters;
  std::vector<MFTAnaSimHit> mAnaSimHits;
  std::vector<std::pair<int, int>> mTrackIDtoMC;

  o2::itsmft::TopologyDictionary mTopoDict;
 
  TFile* mOutFile = nullptr;
  
  TTree* mOutTree1 = nullptr;   // MFTAnsSimTrack
  TTree* mOutTree2 = nullptr;   // MFTAnaSimHit
  TTree* mOutTree3 = nullptr;   // MFTAnaSimCluster 

  int mCurrEvent = 0;
  int mNMCTracks = 0;
  int mNSATracks = 0;
  int mMaxMCTracks = 0;
  int mNHitsInEvent = 0;
  int mNClusters = 0;
  int mVerboseLevel = 0;
};
  
};

#endif // MFT_ANA_SIM

