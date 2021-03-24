#include "../include/MFTAnaHits.h"

namespace o2::mftana
{

using o2::itsmft::Hit;

//_____________________________________________________________________________
MFTAnaHits::MFTAnaHits()
{
}

//_____________________________________________________________________________
void MFTAnaHits::initialize(Int_t maxMCTracks)
{
  mKineTree->SetBranchAddress("MCTrack",&mMCTrkVecP);
  mNrEvents = mKineTree->GetEntries();
  printf("Number of generated events: %d \n", mNrEvents);
  mHitTree->SetBranchAddress("MFTHit",&mHitVecP);
  mMCTrackHasHitsInDisk = std::vector<std::array<bool, o2::mft::constants::DisksNumber>>(maxMCTracks, {0, 0, 0, 0, 0});
}

//_____________________________________________________________________________
void MFTAnaHits::initEvent(Int_t event, Int_t nMCTracks)
{
  mCurrEvent = event;
  mNrMCTracks = nMCTracks;
  for (auto i = 0; i < mMCTrackHasHitsInDisk.size(); i++) {
    for (auto d = 0; d < o2::mft::constants::DisksNumber; d++) {
      mMCTrackHasHitsInDisk.at(i)[d] = false;
    }
  }
}

//_____________________________________________________________________________
Bool_t MFTAnaHits::analyzeHits()
{
  mHitTree->GetEntry(mCurrEvent);
  Int_t nHits = mHitVec.size();
  printf("In event %d found %d hits.\n", mCurrEvent, nHits);

  // identify trackable tracks
  for (Int_t n_hit = 0 ; n_hit < nHits; n_hit++) {
    Hit* hitp = &(mHitVec).at(n_hit);
    Int_t trkID = hitp->GetTrackID(); // ID of the tracks having given the hit
    Float_t z = hitp->GetZ(); // Z position of the hit => Identify MFT disk
    mMCTrackHasHitsInDisk.at(trkID)[mChipMapper.chip2Layer(hitp->GetDetectorID()) / 2] = true;
  }
  
  return kTRUE;
}
  
} // end namespace o2::mftana
