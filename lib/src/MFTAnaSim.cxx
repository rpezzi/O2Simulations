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
bool MFTAnaSim::initialize(int maxMCTracks)
{
  // kinematics
  mKineTree->SetBranchAddress("MCTrack", &mMCTrkVecP);
  mNEvents = mKineTree->GetEntries();
  if (mVerboseLevel > 0) {
    printf("Number of generated events: %d \n", mNEvents);
  }
  // loading is per MC event

  // hits
  mHitTree->SetBranchAddress("MFTHit", &mHitVecP);
  // loading is per MC event

  // clusters
  mClusTree->SetBranchAddress("MFTClusterComp", &mClusVecP);
  if (mClusTree->GetBranch("MFTClusterMCTruth")) {
    mClusTree->SetBranchAddress("MFTClusterMCTruth", &mClusLabels);
  } else {
    printf("No Monte-Carlo information in this file\n");
    return kFALSE;
  }
  auto pattBranch = mClusTree->GetBranch("MFTClusterPatt");
  if (pattBranch) {
    pattBranch->SetAddress(&mClusPatternsP);
  }
  // loading is in a single vector (for all MC events)
  mClusTree->GetEntry(0);
  mNClusters = mClusVec.size();

  // reconstructed tracks, MFT standalone (SA)
  mTrackTree->SetBranchAddress("MFTTrack", &mTrackVecP);
  if (mTrackTree->GetBranch("MFTTrackMCTruth")) {
    mTrackTree->SetBranchAddress("MFTTrackMCTruth", &mTrackLabelsInp);
  } else {
    printf("No Monte-Carlo information in this file\n");
    return kFALSE;
  }
  mTrackTree->SetBranchAddress("MFTTrackClusIdx", &mTrackExtClsVecP);
  // loading is in a single vector (for all MC events)
  mTrackTree->GetEntry(0);
  mNSATracks = mTrackVec.size();
  if (mVerboseLevel >= 0) {
    printf("Number of SA reconstructed tracks: %d \n", mNSATracks);
  }
  
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
  
  mAnaSimClusters.clear();	

  return kTRUE;
}

//_____________________________________________________________________________
void MFTAnaSim::finish()
{
  mOutTree3->Fill();  // MFTAnaSimCluster

  mOutFile->WriteTObject(mOutTree1);
  mOutFile->WriteTObject(mOutTree2);
  mOutFile->WriteTObject(mOutTree3);
  
  mOutFile->Close();
}
  
//_____________________________________________________________________________
void MFTAnaSim::finishEvent()
{
  mOutTree1->Fill();  // MFTAnaSimTrack
  mOutTree2->Fill();  // MFTAnaSimHit
}

//_____________________________________________________________________________
void MFTAnaSim::initEvent(int event, int nMCTracks, int particleSource)
{
  mCurrEvent = event;
  mNMCTracks = nMCTracks;
  mNHitsInEvent = 0;
  for (int i = 0; i < mMCTrackHasHitsInDisk.size(); i++) {
    for (auto di = 0; di < o2::mft::constants::DisksNumber; di++) {
      mMCTrackHasHitsInDisk.at(i)[di] = false;
    }
  }
  for (int i = 0; i < mMCTrackHasHitsInLayer.size(); i++) {
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
  mAnaSimHits.clear();	
}

//_____________________________________________________________________________
bool MFTAnaSim::doParticles()
{
  mKineTree->GetEvent(mCurrEvent);

  int pdgCode, trkID = -1;
  while (++trkID < mNMCTracks) {
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
    }    
  }
}

//_____________________________________________________________________________
void MFTAnaSim::filterPDGCode(int& pdgCode)
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
bool MFTAnaSim::doHits()
{
  mHitTree->GetEntry(mCurrEvent);
  int nHits = mHitVec.size();
  mNHitsInEvent = nHits;
  if (mVerboseLevel > 0) {
    printf("In event %d (%d maxMCTracks) found %d hits.\n", mCurrEvent, mMaxMCTracks, nHits);
  }
  // identify trackable tracks
  for (int n_hit = 0; n_hit < nHits; n_hit++) {
    Hit* hitp = &(mHitVec).at(n_hit);
    int trkID = hitp->GetTrackID();
    int chipID = hitp->GetDetectorID();
    mMCTrackHasHitsInDisk[trkID][mChipMapper.chip2Layer(chipID) / 2] = true;
    mMCTrackHasHitsInLayer[trkID][mChipMapper.chip2Layer(chipID)] = true;
  }
  
  return kTRUE;
}
  
//_____________________________________________________________________________
bool MFTAnaSim::doMCTracks()
{
  if (mNHitsInEvent == 0) {
    mHitTree->GetEntry(mCurrEvent);
    mNHitsInEvent = mHitVec.size();
  }
  
  MFTAnaSimTrack asTrack;
  asTrack.setEvent(mCurrEvent);
  int pdgCode, nMFTHasLayers, nMFTHasDisks;
  int firstHit = -1, lastHit = -1;
  int firstCluster = -1, lastCluster = -1;
  for (int trkID = 0; trkID < mNMCTracks; trkID++) {
    MCTrack* mcTrack =  &(mMCTrkVec)[trkID];

    pdgCode = mcTrack->GetPdgCode();
    filterPDGCode(pdgCode);
    if (pdgCode < 0) continue;
    if (TDatabasePDG::Instance()->GetParticle(pdgCode) == nullptr) {
      continue;
    }

    nMFTHasDisks = 0;
    for(auto disk : {0, 1, 2, 3, 4}) {
      nMFTHasDisks += (int)(mMCTrackHasHitsInDisk[trkID][disk]);
    }
    if (nMFTHasDisks == 0) {
      continue;
    }
    
    nMFTHasLayers = 0;
    for(auto layer : {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}) {
      if (mMCTrackHasHitsInLayer[trkID][layer]) {
	asTrack.setLayer(nMFTHasLayers, layer);
	nMFTHasLayers++;
      }
    }
    asTrack.setMCTrackID(trkID);
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
  
  return kTRUE;
}

//_____________________________________________________________________________
void MFTAnaSim::findMCTrackHits(int trkID, int& firstIndex, int& lastIndex)
{
  // write the hits associated to the MC track
  MFTAnaSimHit asHit;
  for (int n_hit = 0; n_hit < mNHitsInEvent; n_hit++) {
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
void MFTAnaSim::findMCTrackClusters(int trkID, int& firstIndex, int& lastIndex)
{
  // write the clusters associated to the MC track
  int clusSrcID, clusTrkID, clusEvnID, nPixels;
  bool fake;
  o2::math_utils::Point3D<float> locC;
  MFTAnaSimCluster asCluster;
  auto pattIt = mClusPatternsP->cbegin();
  for (int n_cls = 0; n_cls < mNClusters; n_cls++) {
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
    nPixels = 0;
    if (pattID == o2::itsmft::CompCluster::InvalidPatternID || mTopoDict.isGroup(pattID)) {
      o2::itsmft::ClusterPattern patt(pattIt);
      locC = mTopoDict.getClusterCoordinates(cluster, patt, false);
    } else {
      locC = mTopoDict.getClusterCoordinates(cluster);
      nPixels = mTopoDict.getNpixels(pattID);
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
void MFTAnaSim::countParticle(int pdgCode)
{
  MCPart newPart;  
  bool counted = kFALSE;
  for (int i = 0; i < mParticles.size(); i++) {
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
  
//_____________________________________________________________________________
bool MFTAnaSim::doSATracks()
{
  MFTAnaSimSATrack asSATrack;
  
  int srcID, trkID, evnID, nPixels;
  bool fake;
  o2::math_utils::Point3D<float> locC;
  int iTrack = 0;
  auto pattIt = mClusPatternsP->cbegin();
  printf("Found %d SA tracks \n", mTrackVec.size());
  for (auto &track : mTrackVec) {
    auto trkX = track.getX();
    auto trkY = track.getY();
    auto trkZ = track.getZ();
    auto trkLabel = mTrackLabels.at(iTrack);
    //trkLabel.print();
    auto eventID = trkLabel.getEventID();
    auto outParam = track.getOutParam();
    auto trkOutX = outParam.getX();
    auto trkOutY = outParam.getY();
    auto trkOutZ = outParam.getZ();
    printf("Track %3d   isCA %1d   x,y,z-in  %7.3f  %7.3f  %7.3f  x,y,z-out  %7.3f  %7.3f  %7.3f   ev %2d   labels ", iTrack, track.isCA(), trkX, trkY, trkZ, trkOutX, trkOutY, trkOutZ, eventID);
    auto trkID = trkLabel.getTrackID();
    printf("%4d   \n", trkID);
    auto nClus = track.getNumberOfPoints();
    auto offset = track.getExternalClusterIndexOffset();
    printf("Number of clusters %2d \n", nClus);
    for (int i_cls = 0; i_cls < nClus; i_cls++) {
      auto clusEntry = mTrackExtClsVec[offset + i_cls];
      auto cluster = mClusVec[clusEntry];
      auto& clusLabel = (mClusLabels->getLabels(clusEntry))[0];
      if (clusLabel.isNoise()) {
	continue;
      }
      clusLabel.get(trkID, evnID, srcID, fake);
      auto chipID = cluster.getChipID(); 
      auto pattID = cluster.getPatternID();
      //printf("Cluster %5d  chip %03d \n", i_cls, chipID);
      nPixels = 0;
      if (pattID != itsmft::CompCluster::InvalidPatternID) {
	if (mTopoDict.isGroup(pattID)) {
	  o2::itsmft::ClusterPattern patt(pattIt);
	  locC = mTopoDict.getClusterCoordinates(cluster, patt, false);
	} else {
	  locC = mTopoDict.getClusterCoordinates(cluster);
	  nPixels = mTopoDict.getNpixels(pattID);
	}
      } else {
	o2::itsmft::ClusterPattern patt(pattIt);
	locC = mTopoDict.getClusterCoordinates(cluster, patt);
      }
      // Transformation to the local --> global
      auto gloC = mGeoManager->getMatrixL2G(chipID) * locC;
      printf("Cluster %5d   chip ID %03d   evn %2d   mctrk %4d   x,y,z  %7.3f  %7.3f  %7.3f \n", i_cls, cluster.getChipID(), evnID, trkID, gloC.X(), gloC.Y(), gloC.Z());
    }
    iTrack++;
  }
}

}; // end namespace o2::mftana
