#ifndef ALICEO2_MFTANAHITS_H
#define ALICEO2_MFTANAHITS_H

#include <TFile.h>
#include <TTree.h>

#include "MFTBase/Constants.h"
#include "ITSMFTSimulation/Hit.h"
#include "SimulationDataFormat/MCTrack.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"

namespace o2::mftana
{

class MFTAnaHits
{
 public:
  MFTAnaHits();
  ~MFTAnaHits() = default;

  Int_t mNrEvents = 0;
  TTree* mHitTree = nullptr;
  TTree* mKineTree = nullptr;
  o2::itsmft::ChipMappingMFT mChipMapper;

  void initialize(Int_t maxMCTracks);
  void initEvent(Int_t event, Int_t nMCTracks);
  Bool_t doHits();
  Bool_t doMCTracks();

 private:
  std::vector<o2::itsmft::Hit> mHitVec, *mHitVecP = &mHitVec;
  std::vector<o2::MCTrack> mMCTrkVec, *mMCTrkVecP = &mMCTrkVec;
  std::vector<std::array<bool, o2::mft::constants::DisksNumber>> mMCTrackHasHitsInDisk;
  Int_t mCurrEvent = 0;
  Int_t mNrMCTracks = 0;   
};
  
};

#endif // ALICEO2_MFTANAHITS_H

