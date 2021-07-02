#include "../include/MFTAnaSim.h"

#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "MathUtils/Cartesian.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace o2::mftana
{

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

  mOutTree4 = new TTree("MFTAnaSimSATrack", "MFTAnaSimSATrack");
  mOutTree4->Branch("MFTAnaSimSATrack", &mAnaSimSATracks);
}

//_____________________________________________________________________________
bool MFTAnaSim::initialize()
{
  // kinematics
  mKineTree->SetBranchAddress("MCTrack", &mMCTrkVecP);
  mNEvents = mKineTree->GetEntries();
  if (mVerboseLevel >= 0) {
    printf("Number of generated events: %d \n", mNEvents);
  }
  
  // hits
  mHitTree->SetBranchAddress("MFTHit", &mHitVecP);

  // calculate the total number of MC tracks and hits and the maximum number of MC tracks per event
  int nMCTracks, maxMCTracks = 0, totMCTracks = 0, totHits = 0;
  for (int event = 0; event < mNEvents ; event++) {
    mKineTree->GetEntry(event);
    nMCTracks = mMCTrkVec.size();
    totMCTracks += nMCTracks;
    maxMCTracks = std::max(maxMCTracks, nMCTracks);
    mHitTree->GetEntry(event);
    totHits += mHitVec.size();
  }
  /*
  if (mVerboseLevel >= 0) {
    printf("Initialize reserve %d total MC tracks. \n", totMCTracks);
  }
  mAnaSimTracks.reserve(totMCTracks);
  */
  if (mVerboseLevel >= 0) {
    printf("Initialize reserve %d for MC tracks with hits. \n", totHits);
  }
  mAnaSimTracks.reserve(totHits);
  
  if (mVerboseLevel >= 0) {
    printf("Initialize reserve %d total hits. \n", totHits);
  }
  mAnaSimHits.reserve(totHits);
  
  if (mVerboseLevel >= 0) {
    printf("Maximum number of MC tracks per event %d. \n", maxMCTracks);
  }
  mMCTrackHasHitsInDisk = std::vector<std::array<bool, o2::mft::constants::DisksNumber>>(maxMCTracks, {0, 0, 0, 0, 0});
  mMCTrackHasHitsInLayer = std::vector<std::array<bool, o2::mft::constants::LayersNumber>>(maxMCTracks, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
  mMaxMCTracks = maxMCTracks;

  // pattern dictionary for the clusters
  std::string dictfile = "MFTdictionary.bin";
  std::ifstream file(dictfile.c_str());
  if (file.good()) {
    printf("Running with dictionary: %s \n", dictfile.c_str());
    mTopoDict.readBinaryFile(dictfile);
  } else {
    printf("Can not run without dictionary !\n");
    return false;
  }
  
  // clusters
  mClusTree->SetBranchAddress("MFTClusterComp", &mClusVecP);
  if (mClusTree->GetBranch("MFTClusterMCTruth")) {
    mClusTree->SetBranchAddress("MFTClusterMCTruth", &mClusLabels);
  } else {
    printf("No Monte-Carlo information in this file\n");
    return false;
  }
  // cluster patterns
  auto pattBranch = mClusTree->GetBranch("MFTClusterPatt");
  if (pattBranch) {
    pattBranch->SetAddress(&mClusPatternsP);
  } else {
    printf("No patterns!\n");
    return false;
  }
  // loading is in a single vector (for all MC events)
  mClusTree->GetEntry(0);
  mNClusters = mClusVec.size();
  if (mVerboseLevel >= 0) {
    printf("Number of clusters %d \n", mNClusters);
  }
  mAnaSimClusters.reserve(mNClusters);
  mAnaSimClustersSort.reserve(mNClusters);
  extractClusters();
  
  // reconstructed tracks, MFT standalone (SA)
  mTrackTree->SetBranchAddress("MFTTrack", &mTrackVecP);
  if (mTrackTree->GetBranch("MFTTrackMCTruth")) {
    mTrackTree->SetBranchAddress("MFTTrackMCTruth", &mTrackLabelsInp);
  } else {
    printf("No Monte-Carlo information in this file\n");
    return false;
  }
  mTrackTree->SetBranchAddress("MFTTrackClusIdx", &mTrackExtClsVecP);
  // loading is in a single vector (for all MC events)
  mTrackTree->GetEntry(0);
  mNSATracks = mTrackVec.size();
  if (mVerboseLevel >= 0) {
    printf("Number of SA reconstructed tracks: %d \n", mNSATracks);
  }
  mAnaSimSATracks.reserve(mNSATracks);

  mEventClusterRange.reserve(mNEvents);
  
#ifdef _OPENMP
  omp_set_num_threads(mNThreads);
#endif
  
  return true;
}

//_____________________________________________________________________________
void MFTAnaSim::finish()
{
  mOutTree1->Fill();  // MFTAnaSimTrack
  mOutTree2->Fill();  // MFTAnaSimHit
  mOutTree3->Fill();  // MFTAnaSimCluster
  mOutTree4->Fill();  // MFTAnaSimSATrack

  mOutFile->WriteTObject(mOutTree1);
  mOutFile->WriteTObject(mOutTree2);
  mOutFile->WriteTObject(mOutTree3);
  mOutFile->WriteTObject(mOutTree4);
  
  mOutFile->Close();
}
  
//_____________________________________________________________________________
void MFTAnaSim::finishEvent()
{
}

//_____________________________________________________________________________
void MFTAnaSim::initEvent(int event, int particleSource)
{
  mCurrEvent = event;
  mNHitsInEvent = 0;
  for (long unsigned int i = 0; i < mMCTrackHasHitsInDisk.size(); i++) {
    for (auto di = 0; di < o2::mft::constants::DisksNumber; di++) {
      mMCTrackHasHitsInDisk.at(i)[di] = false;
    }
  }
  for (long unsigned int i = 0; i < mMCTrackHasHitsInLayer.size(); i++) {
    for (auto la = 0; la < o2::mft::constants::LayersNumber; la++) {
      mMCTrackHasHitsInLayer.at(i)[la] = false;
    }
  }

  switch (particleSource) {
  case kPrimary:
    mPrimary = true;
    mSecondary = false;
    mAll = false;
    break;
  case kSecondary:
    mPrimary = false;
    mSecondary = true;
    mAll = false;
    break;
  case kAll:
    mPrimary = false;
    mSecondary = false;
    mAll = true;
    break;
  default:
    mPrimary = false;
    mSecondary = false;
    mAll = false;
    break;
  };
  
}

//_____________________________________________________________________________
bool MFTAnaSim::doParticles()
{
  int pdgCode;
  for (int trkID = 0; trkID < mNMCTracks; trkID++) {
    if (!trackHasHits(trkID)) {
      continue;
    }
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
    if (!filterPDGCode(pdgCode)) {
      continue;
    }
    if (TDatabasePDG::Instance()->GetParticle(pdgCode) == nullptr) {
      printf("WARN: no particle with this PDG code in the database %d \n", pdgCode);
      continue;
    }
    //printf("MCTrack %d \n", trkID);
    //printf("pdg code: %d \n", pdgCode);
    //printf("name: %s \n", TDatabasePDG::Instance()->GetParticle(pdgCode)->GetName());
    
    countParticle(pdgCode, mCurrEvent);
    /*
    if (mCurrEvent == 0) {
      printf("MCTrack %4d   PDG %4d   name %10s   isSec %d   E %7.3f \n",
	     trkID,
	     pdgCode,
	     TDatabasePDG::Instance()->GetParticle(pdgCode)->GetName(),
	     mcTrack->isSecondary(),
	     mcTrack->GetEnergy());
    } 
    */   
  }
}

//_____________________________________________________________________________
bool MFTAnaSim::doHits()
{
  if (mNHitsInEvent == 0) {
    mHitTree->GetEntry(mCurrEvent);
    mNHitsInEvent = mHitVec.size();
  }
  
  if (mVerboseLevel > 0) {
    printf("In event %d (%d maxMCTracks) found %d hits.\n", mCurrEvent, mMaxMCTracks, mNHitsInEvent);
  }
  
  // identify trackable tracks
  for (int i_hit = 0; i_hit < mNHitsInEvent; i_hit++) {
    Hit* hitp = &(mHitVec).at(i_hit);
    int trkID = hitp->GetTrackID();
    int chipID = hitp->GetDetectorID();
    mMCTrackHasHitsInDisk[trkID][mChipMapper.chip2Layer(chipID) / 2] = true;
    mMCTrackHasHitsInLayer[trkID][mChipMapper.chip2Layer(chipID)] = true;
  }

  // calculate how many MC tracks have hits
  mKineTree->GetEvent(mCurrEvent);
  mNMCTracks = mMCTrkVec.size();
  if (mVerboseLevel > 0) {
    printf("Number of MC tracks %d \n", mNMCTracks);
  }
  int nMCTracksWHits = 0;
  for (int trkID = 0; trkID < mNMCTracks; trkID++) {
    if (!trackHasHits(trkID)) {
      continue;
    }
    nMCTracksWHits++;
  }
  mNMCTracksWHits += nMCTracksWHits;
  if (mVerboseLevel > 0) {
    printf("In event %d only %d MC tracks have hits (total %d).\n", mCurrEvent, nMCTracksWHits, mNMCTracksWHits);
  }
  
  return true;
}
  
//_____________________________________________________________________________
bool MFTAnaSim::doMCTracks()
{
  mEvnClsIndexMin = mEventClusterRange.at(mCurrEvent).first;
  mEvnClsIndexMax = mEventClusterRange.at(mCurrEvent).second;
  MFTAnaSimTrack asTrack;
  asTrack.setEvent(mCurrEvent);
  int pdgCode, nMFTHasLayers, nMFTHasDisks;
  int firstHit = -1, lastHit = -1;
  int firstCluster = -1, lastCluster = -1;
  int nMCTracksWClusters = 0;
  int nextEvnIndex = 0;
  for (int trkID = 0; trkID < mNMCTracks; trkID++) {
    if (!trackHasHits(trkID)) {
      continue;
    }
    MCTrack* mcTrack =  &(mMCTrkVec)[trkID];
    pdgCode = mcTrack->GetPdgCode();
    if (!filterPDGCode(pdgCode)) {
      continue;
    }
    nMFTHasDisks = 0;
    for(auto disk : {0, 1, 2, 3, 4}) {
      nMFTHasDisks += (int)(mMCTrackHasHitsInDisk[trkID][disk]);
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
    findMCTrackClusters(asTrack);
    if (asTrack.getNClusters() > 0) {
      nMCTracksWClusters++;
    }

    mAnaSimTracks.push_back(asTrack);
  }
  if (mVerboseLevel > 0) {
    printf("MFTAnaSim::doMCTracks found %d tracks with clusters in event %d \n", nMCTracksWClusters, mCurrEvent);
  }
  
  return true;
}

//_____________________________________________________________________________
void MFTAnaSim::findMCTrackHits(int trkID, int& firstIndex, int& lastIndex)
{
  // write the hits associated to the MC track
  MFTAnaSimHit asHit;
  for (int i_hit = 0; i_hit < mNHitsInEvent; i_hit++) {
    Hit* hitp = &(mHitVec).at(i_hit);
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
  if (mVerboseLevel > 0) {
    printf("findMCTrackHits in range %d %d \n", firstIndex, lastIndex);
  }
}

//_____________________________________________________________________________
void MFTAnaSim::findMCTrackClusters(MFTAnaSimTrack& asTrack)
{
  if (mVerboseLevel > 1) {
    printf("findMCTrackClusters in range   %d   %d \n", mEvnClsIndexMin, mEvnClsIndexMax);
  }
  int i_cls, nClusters = 0;
  std::vector<std::pair<int, float>> clusIndexZ;
  for (i_cls = mEvnClsIndexMin; i_cls <= mEvnClsIndexMax; i_cls++) {
    auto& asCluster = mAnaSimClustersSort.at(i_cls);
    if (asCluster.getMCTrackID() != asTrack.getMCTrackID()) {
      continue;
    }
    assert(nClusters < (o2::mftana::MCSplitCluster * o2::mft::constants::LayersNumber));
    clusIndexZ.push_back(std::pair<int, float>(i_cls, asCluster.getZ()));
    nClusters++;
  }
  asTrack.setNClusters(nClusters);
  // sort in decreasing "z" position of the clusters
  sort(clusIndexZ.begin(), clusIndexZ.end(), [](std::pair<int, float>& clsiz1, std::pair<int, float>& clsiz2) { return (clsiz1.second > clsiz2.second); });
  i_cls = 0;
  for (auto& clsiz : clusIndexZ) {
    auto& asCluster = mAnaSimClustersSort.at(clsiz.first);
    if (mVerboseLevel > 1) {
      asCluster.print();
    }
    asTrack.setIntClusIndex(i_cls, asCluster.getID());
    i_cls++;
  }

}

//_____________________________________________________________________________
void MFTAnaSim::countParticle(int pdgCode, int currEvent)
{
  MCPart newPart;  
  bool counted = false;
  for (long unsigned int i = 0; i < mParticles.size(); i++) {
    auto& part = mParticles.at(i);
    if (part.mEvent != currEvent) {
      continue;
    }
    if (part.mPDGCode == pdgCode) {
      part.mCount++;
      counted = true;
      break;
    }
  }
  if (!counted) {
    newPart.mEvent = currEvent;
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
  
  int nPixels, nDisks, nLayers;
  bool hasHitsInDisk[o2::mft::constants::DisksNumber], hasHitsInLayer[o2::mft::constants::LayersNumber];
  o2::math_utils::Point3D<float> locC;
  auto pattIt = mClusPatternsP->cbegin();
  
  int iTrack = 0;
  for (auto& track : mTrackVec) {
    asSATrack.copy(track);
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
    auto trkID = trkLabel.getTrackID();
    auto nPoints = track.getNumberOfPoints();
    asSATrack.setNPoints(nPoints);
    if (mVerboseLevel > 0) {
      printf("Track %3d   isCA %1d   x,y,z-in  %7.3f  %7.3f  %7.3f  x,y,z-out  %7.3f  %7.3f  %7.3f   ev %2d   label %4d   points %d \n", iTrack, track.isCA(), trkX, trkY, trkZ, trkOutX, trkOutY, trkOutZ, eventID, trkID, nPoints);
    }
    for (int j = 0; j < o2::mft::constants::DisksNumber; j++) {
      hasHitsInDisk[j] = false;
    }
    for (int j = 0; j < o2::mft::constants::LayersNumber; j++) {
      hasHitsInLayer[j] = false;
    }
    
    auto offset = track.getExternalClusterIndexOffset();
    for (int i_cls = 0; i_cls < nPoints; i_cls++) {
      auto clusEntry = mTrackExtClsVec[offset + i_cls];
      auto cluster = mAnaSimClusters.at(clusEntry);
      asSATrack.setLayer(i_cls, cluster.getLayer());
      asSATrack.setEvent(i_cls, cluster.getEvent());
      asSATrack.setMCTrackID(i_cls, cluster.getMCTrackID());
      asSATrack.setIntClusIndex(i_cls, clusEntry);
      hasHitsInLayer[cluster.getLayer()] = true;
      hasHitsInDisk[cluster.getLayer() / 2] = true;
      if (mVerboseLevel > 0) {
	printf("Cluster %5d (%5d)  chip ID %03d   evn %2d   mctrk %4d   x,y,z  %7.3f  %7.3f  %7.3f  layer %d \n", i_cls, clusEntry, cluster.getSensorID(), cluster.getEvent(), cluster.getMCTrackID(), cluster.getX(), cluster.getY(), cluster.getZ(), cluster.getLayer());
      }
    }
    nDisks = 0;
    for(auto disk : {0, 1, 2, 3, 4}) {
      nDisks += (int)(hasHitsInDisk[disk]);
    }
    nLayers = 0;
    for(auto layer : {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}) {
      nLayers += (int)(hasHitsInLayer[layer]);
    }
    asSATrack.setNDisks(nDisks);
    asSATrack.setNLayers(nLayers);

    mAnaSimSATracks.push_back(asSATrack);
    
    iTrack++;
  }
  
}

//_____________________________________________________________________________
void MFTAnaSim::extractClusters()
{
  auto pattIt = mClusPatternsP->cbegin();
  o2::itsmft::ClusterPattern patt(pattIt);
  int i_cls, clusSrcID, clusTrkID, clusEvnID, nPixels;
  bool fake;
  o2::math_utils::Point3D<float> locC;
  MFTAnaSimCluster asCluster;
  for (i_cls = 0; i_cls < mNClusters; i_cls++) {
    asCluster.setID(i_cls);
    auto cluster = mClusVec[i_cls];
    auto chipID = cluster.getChipID(); 
    asCluster.setSensorID(chipID);
    auto layer = mChipMapper.chip2Layer(chipID);
    asCluster.setLayer(layer);
    auto& label = (mClusLabels->getLabels(i_cls))[0];
    clusTrkID = clusEvnID = clusSrcID = -1;
    if (label.isNoise()) {
      asCluster.setIsNoise(true);
    } else {
      asCluster.setIsNoise(false);
      label.get(clusTrkID, clusEvnID, clusSrcID, fake);
    }
    
    asCluster.setEvent(clusEvnID);
    asCluster.setMCTrackID(clusTrkID);
    auto pattID = cluster.getPatternID();
    nPixels = 0;   
    if (pattID != o2::itsmft::CompCluster::InvalidPatternID) {
      if (mTopoDict.isGroup(pattID)) {
	locC = mTopoDict.getClusterCoordinates(cluster, patt);
      } else {
	locC = mTopoDict.getClusterCoordinates(cluster);
	nPixels = mTopoDict.getNpixels(pattID);
      }
    } else {
      locC = mTopoDict.getClusterCoordinates(cluster, patt, false);
    }
    
    // Transformation to the local --> global
    auto gloC = mGeoManager->getMatrixL2G(chipID) * locC;
    asCluster.setX(gloC.X());
    asCluster.setY(gloC.Y());
    asCluster.setZ(gloC.Z());
    asCluster.setNPixels(nPixels);
    if (mVerboseLevel > 3) {
      printf("Extract cluster  %5d  %3d  %3d (%3d)  %1d  %7.3f  %7.3f  %7.3f event %d \n", i_cls, chipID, pattID, o2::itsmft::CompCluster::InvalidPatternID, mTopoDict.isGroup(pattID), gloC.X(), gloC.Y(), gloC.Z(), clusEvnID);
    }
    mAnaSimClusters.push_back(asCluster);
    // for attaching to MC tracks
    if (clusTrkID >= 0) {
      mAnaSimClustersSort.push_back(asCluster);
    }
  }

  // sort the cluster in event ID incfreasing order (starting with noise)
  sort(mAnaSimClustersSort.begin(), mAnaSimClustersSort.end(), [](MFTAnaSimCluster& cls1, MFTAnaSimCluster& cls2) {
      return (cls1.getEvent() < cls2.getEvent()); });
  // and calculate the index ranges for the MC tracks
  int devn, minClsIndex, maxClsIndex, prevEvnID = -1;
  int nClustersSort = mAnaSimClustersSort.size();
  for (i_cls = 0; i_cls < nClustersSort; i_cls++) {
    auto& asCluster = mAnaSimClustersSort.at(i_cls);
    clusEvnID = asCluster.getEvent();
    if (clusEvnID > prevEvnID) {
      if (prevEvnID < 0) {
	minClsIndex = i_cls;
	prevEvnID = clusEvnID;
      } else {
	// events without clusters
	devn = clusEvnID - prevEvnID - 1;
	for (int ievn = 0; ievn < devn; ievn++) {
	  if (mVerboseLevel > 1) {
	    printf("extractClusters evn %d cls range   %d   %d \n", (prevEvnID + ievn), 0, -1);
	  }
	  mEventClusterRange.push_back(std::pair(0, -1));
	  prevEvnID++;
	}
	maxClsIndex = i_cls - 1;
	if (mVerboseLevel > 1) {
	  printf("extractClusters evn %d cls range   %d   %d \n", prevEvnID, minClsIndex, maxClsIndex);
	}
	mEventClusterRange.push_back(std::pair(minClsIndex, maxClsIndex));
	minClsIndex = i_cls;
	prevEvnID = clusEvnID;
      }
    }
  }
  maxClsIndex = i_cls - 1;
  if (mVerboseLevel > 1) {
    printf("extractClusters evn %d cls range   %d   %d \n", prevEvnID, minClsIndex, maxClsIndex);
  }
  mEventClusterRange.push_back(std::pair(minClsIndex, maxClsIndex));
}
  
//_____________________________________________________________________________
void MFTAnaSim::linkTracks()
{
#ifdef _OPENMP
  int thread_nr, n_threads, thread_tracks, thread_start;
  int nMCTracks = mAnaSimTracks.size();
#pragma omp parallel default(shared) private(thread_nr, n_threads, thread_tracks, thread_start)
  {
    thread_nr = omp_get_thread_num();
    n_threads =  omp_get_num_threads();
    thread_tracks = nMCTracks / n_threads;
    thread_start = thread_nr * thread_tracks;
    if (thread_nr == (n_threads - 1)) {
      thread_tracks = nMCTracks - thread_start;
    }
    for (int iMCTrack = thread_start; iMCTrack < (thread_start + thread_tracks); iMCTrack++) {
      auto& mcTrack = mAnaSimTracks.at(iMCTrack);
      for (auto iSATrack = 0; iSATrack < mAnaSimSATracks.size(); iSATrack++) {
	auto& saTrack = mAnaSimSATracks.at(iSATrack);
	for (int ip = 0; ip < saTrack.getNPoints(); ip++) {
	  if (saTrack.getEvent(ip) == mcTrack.getEvent() && saTrack.getMCTrackID(ip) == mcTrack.getMCTrackID()) {
	    mcTrack.addIntSATrackIndex(iSATrack);
	    saTrack.addIntMCTrackIndex(mcTrack.getEvent(), iMCTrack);
	  }
	}
      }
    }
  }
#else
  for (auto iMCTrack = 0; iMCTrack < mAnaSimTracks.size(); iMCTrack++) {
    auto& mcTrack = mAnaSimTracks.at(iMCTrack);
    for (auto iSATrack = 0; iSATrack < mAnaSimSATracks.size(); iSATrack++) {
      auto& saTrack = mAnaSimSATracks.at(iSATrack);
      for (int ip = 0; ip < saTrack.getNPoints(); ip++) {
	if (saTrack.getEvent(ip) == mcTrack.getEvent() && saTrack.getMCTrackID(ip) == mcTrack.getMCTrackID()) {
	  mcTrack.addIntSATrackIndex(iSATrack);
	  saTrack.addIntMCTrackIndex(mcTrack.getEvent(), iMCTrack);
	}
      }
    }
  }
#endif // _OPENMP
}

}; // end namespace o2::mftana
