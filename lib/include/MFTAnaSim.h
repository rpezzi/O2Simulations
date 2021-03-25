#ifndef ALICEO2_MFTANASIM_H
#define ALICEO2_MFTANASIM_H

#include <TFile.h>
#include <TTree.h>

#include "MFTBase/Constants.h"
#include "ITSMFTSimulation/Hit.h"
#include "SimulationDataFormat/MCTrack.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"

namespace o2::mftana
{

struct MCPart {
  //MCPart() = default;
  //MCPart& operator=(const MCPart&) = default;
  Int_t mPDGCode = -1;
  std::string mPDGName = "";
  Int_t mCount = 0;
};

class MFTAnaSim
{
 public:
  MFTAnaSim();
  ~MFTAnaSim() = default;

  Int_t mNrEvents = 0;
  TTree* mHitTree = nullptr;
  TTree* mKineTree = nullptr;
  o2::itsmft::ChipMappingMFT mChipMapper;

  void initialize(Int_t maxMCTracks);
  void initEvent(Int_t event, Int_t nMCTracks, Int_t particleSource = 0);
  Bool_t doParticles();
  Bool_t doHits();
  Bool_t doMCTracks();

  const std::vector<MCPart>& getParticles() { return mParticles; }

 private:
  enum ParticleSource {kPrimary, kSecondary, kAll};
  void countParticle(Int_t pdgCode);
  Bool_t mPrimary = kTRUE, mSecondary = kFALSE, mAll = kFALSE;
  std::vector<o2::itsmft::Hit> mHitVec, *mHitVecP = &mHitVec;
  std::vector<o2::MCTrack> mMCTrkVec, *mMCTrkVecP = &mMCTrkVec;
  std::vector<std::array<bool, o2::mft::constants::DisksNumber>> mMCTrackHasHitsInDisk;
  std::vector<MCPart> mParticles;
  Int_t mCurrEvent = 0;
  Int_t mNrMCTracks = 0;
};
  
};

#endif // ALICEO2_MFTANASIM_H

