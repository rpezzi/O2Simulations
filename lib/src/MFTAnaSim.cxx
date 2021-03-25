#include "../include/MFTAnaSim.h"

#include "TDatabasePDG.h"

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

  // maximum MC tracks per event, from a scan of all events in the file
  mMCTrackHasHitsInDisk = std::vector<std::array<bool, o2::mft::constants::DisksNumber>>(maxMCTracks, {0, 0, 0, 0, 0});
}

//_____________________________________________________________________________
void MFTAnaSim::initEvent(Int_t event, Int_t nMCTracks, Int_t particleSource)
{
  mCurrEvent = event;
  mNrMCTracks = nMCTracks;
  for (auto i = 0; i < mMCTrackHasHitsInDisk.size(); i++) {
    for (auto d = 0; d < o2::mft::constants::DisksNumber; d++) {
      mMCTrackHasHitsInDisk.at(i)[d] = false;
    }
  }
  mParticles.clear();
  switch (particleSource) {
  case kPrimary:
    mPrimary = kTRUE;
    mSecondary = kFALSE;
    mAll = kFALSE;
    break;
  case kSecondary:
    mPrimary = kFALSE;
    mSecondary = kTRUE;
    mAll = kFALSE;
    break;
  case kAll:
    mPrimary = kFALSE;
    mSecondary = kFALSE;
    mAll = kTRUE;
    break;
  default:
    mPrimary = kFALSE;
    mSecondary = kFALSE;
    mAll = kFALSE;
    break;
  };
	
}

//_____________________________________________________________________________
Bool_t MFTAnaSim::doParticles()
{
  MCPart particle;
  Int_t pdgCode;
  for (Int_t trkID = 0 ; trkID < mNrMCTracks; trkID++) {
    MCTrack* mcTrack =  &(mMCTrkVec)[trkID];
    if (!mAll) {
      if (mPrimary && !mcTrack->isPrimary()) {
	continue;
      }
      if (mSecondary && mcTrack->isPrimary()) {
	continue;
      }
      if (!mPrimary && !mSecondary) {
	continue;
      }
    }
    pdgCode = mcTrack->GetPdgCode();
    // skip pomeron and reggeon
    if (pdgCode == 110 || pdgCode == 990) {
      continue;
    }
    // skip ...
    if (pdgCode > 990000000) {
      continue;
    }
    // diffractive states
    if (pdgCode > 9900000) {
      pdgCode -= 9900000;
    }
    //printf("MCTrack %d \n", trkID);
    //printf("pdg code: %d \n", pdgCode);
    //printf("name: %s \n", TDatabasePDG::Instance()->GetParticle(pdgCode)->GetName());
    
    countParticle(pdgCode);
    /*
    if (mCurrEvent == 0) {
      printf("MCTrack %4d   PDG %4d   name %s   isSec %d   E %7.3f \n",
	     trkID,
	     pdgCode,
	     TDatabasePDG::Instance()->GetParticle(pdgCode)->GetName(),
	     mcTrack->isSecondary(),
	     mcTrack->GetEnergy());
      printf("Particle name: %s \n", particle.mPDGName.c_str());
    }
    */
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

//_____________________________________________________________________________
void MFTAnaSim::countParticle(Int_t pdgCode)
{
  MCPart newPart;  
  Bool_t counted = kFALSE;
  for (auto i = 0; i < mParticles.size(); i++) {
    auto& part = mParticles.at(i);
    if (part.mPDGCode == pdgCode) {
      part.mCount++;
      counted = kTRUE;
      break;
    }
  }
  if (!counted) {
    newPart.mPDGCode = pdgCode;
    newPart.mPDGName = std::string(TDatabasePDG::Instance()->GetParticle(pdgCode)->GetName());
    newPart.mCount = 1;
    mParticles.push_back(newPart);
  }
}
  
} // end namespace o2::mftana
