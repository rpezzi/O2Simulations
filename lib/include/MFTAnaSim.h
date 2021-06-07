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
  Int_t mEvent = 0;
  Int_t mPDGCode = -1;
  std::string mPDGName = "";
  Int_t mCount = 0;
};

class MFTAnaSim
{
 public:
  MFTAnaSim();
  ~MFTAnaSim() = default;

  Int_t mNEvents = 0;

  TTree* mHitTree = nullptr;
  TTree* mKineTree = nullptr;
  TTree* mClusTree = nullptr;
  TTree* mTrackTree = nullptr;
  
  o2::mft::GeometryTGeo* mGeoManager;
  o2::itsmft::ChipMappingMFT mChipMapper;

  Bool_t initialize(Int_t maxMCTracks);
  void initEvent(Int_t event, Int_t nMCTracks, Int_t particleSource = 0);
  Bool_t doParticles();
  Bool_t doHits();
  Bool_t doMCTracks();
  Bool_t doSATracks();
  void finishEvent();
  void finish();
  void findMCTrackHits(Int_t trkID, Int_t& firstIndex, Int_t& lastIndex);
  void findMCTrackClusters(Int_t trkID, Int_t& firstIndex, Int_t& lastIndex);
  void setVerboseLevel(Int_t vl) { mVerboseLevel = vl; }
  const std::vector<MCPart>& getParticles() { return mParticles; }

 private:
  Bool_t mPrimary = kTRUE, mSecondary = kFALSE, mAll = kFALSE;
  enum ParticleSource {kPrimary, kSecondary, kAll};
  void filterPDGCode(Int_t& pdgCode);
  void countParticle(Int_t pdgCode);
  
  std::vector<o2::MCTrack> mMCTrkVec, *mMCTrkVecP = &mMCTrkVec;
  std::vector<o2::itsmft::Hit> mHitVec, *mHitVecP = &mHitVec;
  std::vector<o2::itsmft::CompClusterExt> mClusVec, *mClusVecP = &mClusVec;
  std::vector<o2::mft::TrackMFT> mTrackVec, *mTrackVecP = &mTrackVec;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* mClusLabels = nullptr;
  std::vector<o2::MCCompLabel> mTrackLabels, *mTrackLabelsInp = &mTrackLabels;
  std::vector<Int_t> mTrackExtClsVec, *mTrackExtClsVecP = &mTrackExtClsVec;

  std::vector<std::array<bool, o2::mft::constants::DisksNumber>> mMCTrackHasHitsInDisk;
  std::vector<std::array<bool, o2::mft::constants::LayersNumber>> mMCTrackHasHitsInLayer;
  std::vector<MCPart> mParticles;
  std::vector<MFTAnaSimTrack> mAnaSimTracks;
  std::vector<MFTAnaSimCluster> mAnaSimClusters;
  std::vector<MFTAnaSimHit> mAnaSimHits;
  std::vector<std::pair<Int_t, Int_t>> mTrackIDtoMC;

  o2::itsmft::TopologyDictionary mTopoDict;
 
  TFile* mOutFile = nullptr;
  
  TTree* mOutTree1 = nullptr;   // MFTAnsSimTrack
  TTree* mOutTree2 = nullptr;   // MFTAnaSimHit
  TTree* mOutTree3 = nullptr;   // MFTAnaSimCluster 

  Int_t mCurrEvent = 0;
  Int_t mNMCTracks = 0;
  Int_t mNSATracks = 0;
  Int_t mMaxMCTracks = 0;
  Int_t mNHitsInEvent = 0;
  Int_t mNClusters = 0;
  Int_t mVerboseLevel = 0;
};
  
};

#endif // MFT_ANA_SIM

