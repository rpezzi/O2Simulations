#include "../include/MFTAnaSim.h"

namespace o2::mftana
{

using o2::itsmft::Hit;
using o2::MCTrack;

//_____________________________________________________________________________
MFTAnaSim::MFTAnaSim()
{
}

//_____________________________________________________________________________
void MFTAnaSim::initialize(Int_t maxMCTracks)
{
  mKineTree->SetBranchAddress("MCTrack",&mMCTrkVecP);
  mNrEvents = mKineTree->GetEntries();
  printf("Number of generated events: %d \n", mNrEvents);
  mHitTree->SetBranchAddress("MFTHit",&mHitVecP);
  mMCTrackHasHitsInDisk = std::vector<std::array<bool, o2::mft::constants::DisksNumber>>(maxMCTracks, {0, 0, 0, 0, 0});
}

//_____________________________________________________________________________
void MFTAnaSim::initEvent(Int_t event, Int_t nMCTracks)
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
Bool_t MFTAnaSim::doHits()
{
  mHitTree->GetEntry(mCurrEvent);
  Int_t nHits = mHitVec.size();
  printf("In event %d found %d hits.\n", mCurrEvent, nHits);

  // identify trackable tracks
  for (Int_t n_hit = 0 ; n_hit < nHits; n_hit++) {
    Hit* hitp = &(mHitVec).at(n_hit);
    Int_t trkID = hitp->GetTrackID(); // ID of the tracks having given the hit
    Float_t z = hitp->GetZ(); // Z position of the hit => Identify MFT disk
    mMCTrackHasHitsInDisk[trkID][mChipMapper.chip2Layer(hitp->GetDetectorID()) / 2] = true;
  }
  
  return kTRUE;
}
  
//_____________________________________________________________________________
Bool_t MFTAnaSim::doMCTracks()
{
  for (Int_t trkID = 0 ; trkID < mNrMCTracks; trkID++) {
    MCTrack* mcTrack =  &(mMCTrkVec)[trkID];
    if (!mcTrack->isPrimary()) {
      continue;
    }
    Int_t nMFTDisksHasHits = 0;
    for(auto disk : {0, 1, 2, 3, 4}) {
      nMFTDisksHasHits += (Int_t)(mMCTrackHasHitsInDisk[trkID][disk]);
    }
  }
  
  return kTRUE;
}

} // end namespace o2::mftana
