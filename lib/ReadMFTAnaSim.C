#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "include/MFTAnaSim.h"

#include "SimulationDataFormat/MCEventHeader.h"

#endif

using MFTAnaSimSATrack = o2::mftana::MFTAnaSimSATrack;
using MFTAnaSimTrack = o2::mftana::MFTAnaSimTrack;
using MFTAnaSimCluster = o2::mftana::MFTAnaSimCluster;
using MFTAnaSimHit = o2::mftana::MFTAnaSimHit;

void ReadMFTAnaSim()
{
  TFile inFile("MFTAnaSimTracks_36.root");
  
  TTree *tree1 = (TTree*)inFile.Get("MFTAnaSimTrack");
  std::vector<MFTAnaSimTrack> anaSimTracks, *anaSimTracksP = &anaSimTracks;
  tree1->SetBranchAddress("MFTAnaSimTrack", &anaSimTracksP);
  
  TTree *tree2 = (TTree*)inFile.Get("MFTAnaSimHit");
  std::vector<MFTAnaSimHit> anaSimHits, *anaSimHitsP = &anaSimHits;
  tree2->SetBranchAddress("MFTAnaSimHit", &anaSimHitsP);
  
  TTree *tree3 = (TTree*)inFile.Get("MFTAnaSimCluster");
  std::vector<MFTAnaSimCluster> anaSimClusters, *anaSimClustersP = &anaSimClusters;
  tree3->SetBranchAddress("MFTAnaSimCluster", &anaSimClustersP);
  
  TTree *tree4 = (TTree*)inFile.Get("MFTAnaSimSATrack");
  std::vector<MFTAnaSimSATrack> anaSimSATracks, *anaSimSATracksP = &anaSimSATracks;
  tree4->SetBranchAddress("MFTAnaSimSATrack", &anaSimSATracksP);

  tree1->GetEntry(0);
  tree2->GetEntry(0);
  tree3->GetEntry(0);
  tree4->GetEntry(0);
  
  int nPoints, nClusters, iTrack = 0;
  for (auto& track : anaSimTracks) {
    
    printf("Event %d Track %d NDisks %d NLayers %d MC %d ", track.getEvent(), iTrack, track.getNDisks(), track.getNLayers(), track.getMCTrackID());
    printf("layers: ");
    for (int index = 0; index < track.getNLayers(); index++) {
      printf("%d ", track.getLayer(index));
    }
    printf("\n");
    
    // hits
    printf("Number of hits %d , hit range: %d %d \n", track.getNHits(), track.getFirstHitIndex(), track.getLastHitIndex());
    if (track.getFirstHitIndex() >= 0) {
      for (int iHit = track.getFirstHitIndex(); iHit <= track.getLastHitIndex(); iHit++) {
	auto& hit = anaSimHits.at(iHit);
	printf("%5d:   %7.3f   %7.3f   %7.3f \n", iHit, hit.getX(), hit.getY(), hit.getZ());
      }
    }
    
    // clusters
    nClusters = track.getNClusters();
    printf("Number of clusters %d \n", nClusters);
    if (nClusters > 0) {
      for (int i_cls = 0; i_cls < nClusters; i_cls++) {
	auto& cluster = anaSimClusters.at(track.getIntClusIndex(i_cls));
	printf("%5d:   %7.3f   %7.3f   %7.3f \n", i_cls, cluster.getX(), cluster.getY(), cluster.getZ());
      }
    } else {
      printf("... no clusters ...\n");
    }
    
    iTrack++;
  }
  
  int clusEntry, iSATrack = 0;
  for (auto& saTrack : anaSimSATracks) {
    printf("SATrack %d points %d disks %d layers %d \n", iSATrack, saTrack.getNPoints(), saTrack.getNDisks(), saTrack.getNLayers());
    nPoints = saTrack.getNPoints();
    printf("Number of points %d \n", nPoints);
    if (nPoints > 0) {
      for (int i_p = 0; i_p < nPoints; i_p++) {
	clusEntry = saTrack.getIntClusIndex(i_p);
	auto& cluster = anaSimClusters.at(clusEntry);
	printf("%5d  %5d:   %7.3f   %7.3f   %7.3f   %d \n", i_p, clusEntry, cluster.getX(), cluster.getY(), cluster.getZ(), cluster.getLayer());
      }
    } else {
      printf("... no clusters ...\n");
    }
    
    iSATrack++;
  }
    
}
