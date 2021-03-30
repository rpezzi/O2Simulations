#include "../include/MFTAnaSim.h"

#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "MathUtils/Cartesian.h"

namespace o2::mftana
{

using o2::itsmft::Hit;
using o2::MCTrack;

//_____________________________________________________________________________
MFTAnaSim::MFTAnaSim()
{
  mOutFile = new TFile("MFTAnaSimTracks.root", "recreate");
  
  mOutTree1 = new TTree("MFTAnaSimTrack", "MFTAnaSimTrack");
  mOutTree1->Branch("MFTAnaSimTrack", &mAnaSimTracks);
  
  mOutTree2 = new TTree("MFTAnaSimHit", "MFTAnaSimHit");
  mOutTree2->Branch("MFTAnaSimHit", &mAnaSimHits);

  mOutTree3 = new TTree("MFTAnaSimCluster", "MFTAnaSimCluster");
  mOutTree3->Branch("MFTAnaSimCluster", &mAnaSimClusters);
}

//_____________________________________________________________________________
Bool_t MFTAnaSim::initialize(Int_t maxMCTracks)
{
  mKineTree->SetBranchAddress("MCTrack", &mMCTrkVecP);
  mNrEvents = mKineTree->GetEntries();
  if (mVerboseLevel > 0) {
    printf("Number of generated events: %d \n", mNrEvents);
  }
  mHitTree->SetBranchAddress("MFTHit", &mHitVecP);
  mClusTree->SetBranchAddress("MFTClusterComp", &mClusVecP);
  if (mClusTree->GetBranch("MFTClusterMCTruth")) {
    mClusTree->SetBranchAddress("MFTClusterMCTruth", &mClusLabels);
  } else {
    printf("No Monte-Carlo information in this file\n");
    return kFALSE;
  }
  mClusTree->GetEntry(0);
  mNClusters = mClusVec.size();

  // maximum MC tracks per event, from a scan of all events in the file
  mMCTrackHasHitsInDisk = std::vector<std::array<bool, o2::mft::constants::DisksNumber>>(maxMCTracks, {0, 0, 0, 0, 0});
  mMCTrackHasHitsInLayer = std::vector<std::array<bool, o2::mft::constants::LayersNumber>>(maxMCTracks, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
  mMaxMCTracks = maxMCTracks;

  std::string dictfile = "MFTdictionary.bin";
  std::ifstream file(dictfile.c_str());
  if (file.good()) {
    printf("Running with dictionary: %s \n", dictfile.c_str());
    mTopoDict.readBinaryFile(dictfile);
  } else {
    printf("Can not run without dictionary !\n");
    return kFALSE;
  }
  
  return kTRUE;
}

//_____________________________________________________________________________
void MFTAnaSim::finish()
{
  mOutFile->WriteTObject(mOutTree1);
  mOutFile->WriteTObject(mOutTree2);
  mOutFile->WriteTObject(mOutTree3);
  mOutFile->Close();
}
  
//_____________________________________________________________________________
void MFTAnaSim::finishEvent()
{
  mOutTree1->Fill();
  mOutTree2->Fill();
  mOutTree3->Fill();
}

//_____________________________________________________________________________
void MFTAnaSim::initEvent(Int_t event, Int_t nMCTracks, Int_t particleSource)
{
  mCurrEvent = event;
  mNrMCTracks = nMCTracks;
  mNHitsInEvent = 0;
  for (UInt_t i = 0; i < mMCTrackHasHitsInDisk.size(); i++) {
    for (auto di = 0; di < o2::mft::constants::DisksNumber; di++) {
      mMCTrackHasHitsInDisk.at(i)[di] = false;
    }
  }
  for (UInt_t i = 0; i < mMCTrackHasHitsInLayer.size(); i++) {
    for (auto la = 0; la < o2::mft::constants::LayersNumber; la++) {
      mMCTrackHasHitsInLayer.at(i)[la] = false;
    }
  }

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
  
  mParticles.clear();
  mAnaSimTracks.clear();	
  mAnaSimClusters.clear();	
  mAnaSimHits.clear();	
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
    filterPDGCode(pdgCode);
    if (pdgCode < 0) {
      continue;
    }
    if (TDatabasePDG::Instance()->GetParticle(pdgCode) == nullptr) {
      continue;
    }
    //printf("MCTrack %d \n", trkID);
    //printf("pdg code: %d \n", pdgCode);
    //printf("name: %s \n", TDatabasePDG::Instance()->GetParticle(pdgCode)->GetName());
    
    countParticle(pdgCode);
    
    if (mCurrEvent == 0) {
      printf("MCTrack %4d   PDG %4d   name %s   isSec %d   E %7.3f \n",
	     trkID,
	     pdgCode,
	     TDatabasePDG::Instance()->GetParticle(pdgCode)->GetName(),
	     mcTrack->isSecondary(),
	     mcTrack->GetEnergy());
      printf("Particle name: %s \n", particle.mPDGName.c_str());
    }
    
  }
}

//_____________________________________________________________________________
void MFTAnaSim::filterPDGCode(Int_t& pdgCode)
{
  // skip pomeron and reggeon
  if (pdgCode == 110 || pdgCode == 990) {
    pdgCode = -1;
  }
  // diffractive states
  if (pdgCode > 9900000) {
    pdgCode -= 9900000;
  }
}
  
//_____________________________________________________________________________
Bool_t MFTAnaSim::doHits()
{
  mHitTree->GetEntry(mCurrEvent);
  Int_t nHits = mHitVec.size();
  mNHitsInEvent = nHits;
  if (mVerboseLevel > 0) {
    printf("In event %d (%d maxMCTracks) found %d hits.\n", mCurrEvent, mMaxMCTracks, nHits);
  }
  // identify trackable tracks
  for (Int_t n_hit = 0 ; n_hit < nHits; n_hit++) {
    Hit* hitp = &(mHitVec).at(n_hit);
    Int_t trkID = hitp->GetTrackID();
    Int_t chipID = hitp->GetDetectorID();
    mMCTrackHasHitsInDisk[trkID][mChipMapper.chip2Layer(chipID) / 2] = true;
    mMCTrackHasHitsInLayer[trkID][mChipMapper.chip2Layer(chipID)] = true;
  }
  
  return kTRUE;
}
  
//_____________________________________________________________________________
Bool_t MFTAnaSim::doMCTracks()
{
  if (mNHitsInEvent == 0) {
    mHitTree->GetEntry(mCurrEvent);
    mNHitsInEvent = mHitVec.size();
  }
  
  MFTAnaSimTrack asTrack;
  Int_t pdgCode, nMFTHasLayers, nMFTHasDisks;
  Int_t firstHit = -1, lastHit = -1;
  Int_t firstCluster = -1, lastCluster = -1;
  for (Int_t trkID = 0 ; trkID < mNrMCTracks; trkID++) {
    MCTrack* mcTrack =  &(mMCTrkVec)[trkID];

    pdgCode = mcTrack->GetPdgCode();
    filterPDGCode(pdgCode);
    if (pdgCode < 0) continue;
    if (TDatabasePDG::Instance()->GetParticle(pdgCode) == nullptr) {
      continue;
    }

    nMFTHasDisks = 0;
    for(auto disk : {0, 1, 2, 3, 4}) {
      nMFTHasDisks += (Int_t)(mMCTrackHasHitsInDisk[trkID][disk]);
    }
    if (nMFTHasDisks > 0) {
      nMFTHasLayers = 0;
      for(auto layer : {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}) {
	if (mMCTrackHasHitsInLayer[trkID][layer]) {
	  asTrack.setLayer(nMFTHasLayers, layer);
	  nMFTHasLayers++;
	}
      }
      asTrack.setNDisks(nMFTHasDisks);
      asTrack.setNLayers(nMFTHasLayers);
      asTrack.setVertexPxPyPz(mcTrack->Px(), mcTrack->Py(), mcTrack->Pz());
      asTrack.setVertexXYZ(mcTrack->GetStartVertexCoordinatesX(), mcTrack->GetStartVertexCoordinatesY(), mcTrack->GetStartVertexCoordinatesZ());
      asTrack.setPDGCode(pdgCode);
      asTrack.setIsPrimary(mcTrack->isPrimary());
      asTrack.setMotherTrackId(mcTrack->getMotherTrackId());
      asTrack.setSecondMotherTrackId(mcTrack->getSecondMotherTrackId());
      asTrack.setFirstDaughterTrackId(mcTrack->getFirstDaughterTrackId());
      asTrack.setLastDaughterTrackId(mcTrack->getLastDaughterTrackId());

      // associate hits
      firstHit = lastHit = -1;
      findMCTrackHits(trkID, firstHit, lastHit);
      asTrack.setNHits(lastHit - firstHit + 1);
      asTrack.setFirstHitIndex(firstHit);
      asTrack.setLastHitIndex(lastHit);
      
      // associate clusters
      firstCluster = lastCluster = -1;
      findMCTrackClusters(trkID, firstCluster, lastCluster);
      asTrack.setNClusters(lastCluster - firstCluster + 1);
      asTrack.setFirstClusterIndex(firstCluster);
      asTrack.setLastClusterIndex(lastCluster);
      
      mAnaSimTracks.push_back(asTrack);
    }
  }
  
  return kTRUE;
}

//_____________________________________________________________________________
void MFTAnaSim::findMCTrackHits(Int_t trkID, Int_t& firstIndex, Int_t& lastIndex)
{
  // write the hits associated to the MC track
  MFTAnaSimHit asHit;
  for (Int_t n_hit = 0 ; n_hit < mNHitsInEvent; n_hit++) {
    Hit* hitp = &(mHitVec).at(n_hit);
    if (hitp->GetTrackID() != trkID) {
      continue;
    }
    if (firstIndex < 0) {
      firstIndex = mAnaSimHits.size();
    }
    asHit.setXYZ(hitp->GetStartX(), hitp->GetStartY(), hitp->GetStartZ());
    mAnaSimHits.push_back(asHit);
    lastIndex = mAnaSimHits.size() - 1;
  }
}

//_____________________________________________________________________________
void MFTAnaSim::findMCTrackClusters(Int_t trkID, Int_t& firstIndex, Int_t& lastIndex)
{
  // write the clusters associated to the MC track
  Int_t clusSrcID, clusTrkID, clusEvnID, npix;
  Bool_t fake;
  o2::math_utils::Point3D<Float_t> locC;
  MFTAnaSimCluster asCluster;
  for (Int_t n_cls = 0 ; n_cls < mNClusters; n_cls++) {
    auto cluster = mClusVec[n_cls];
    auto& label = (mClusLabels->getLabels(n_cls))[0];
    if (label.isNoise()) {
      continue;
    }
    label.get(clusTrkID, clusEvnID, clusSrcID, fake);
    if (clusTrkID != trkID) {
	continue;
    }
    auto chipID = cluster.getChipID(); 
    auto pattID = cluster.getPatternID();
    npix = 0;
    if (pattID == o2::itsmft::CompCluster::InvalidPatternID || mTopoDict.isGroup(pattID)) {
      // temporary fix ...
      locC = mTopoDict.getClusterCoordinates(cluster);
      
      //o2::itsmft::ClusterPattern patt(pattIt);
      //locC = mTopoDict.getClusterCoordinates(cluster, patt);
    } else {
      locC = mTopoDict.getClusterCoordinates(cluster);
      npix = mTopoDict.getNpixels(pattID);
    }
    // Transformation to the local --> global
    auto gloC = mGeoManager->getMatrixL2G(chipID) * locC;
    if (firstIndex < 0) {
      firstIndex = mAnaSimClusters.size();
    }
    asCluster.setSensorID(chipID);
    asCluster.setX(gloC.X());
    asCluster.setY(gloC.Y());
    asCluster.setZ(gloC.Z());
    mAnaSimClusters.push_back(asCluster);
    lastIndex = mAnaSimClusters.size() - 1;
  }
}

//_____________________________________________________________________________
void MFTAnaSim::countParticle(Int_t pdgCode)
{
  MCPart newPart;  
  Bool_t counted = kFALSE;
  for (UInt_t i = 0; i < mParticles.size(); i++) {
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
  
}; // end namespace o2::mftana
