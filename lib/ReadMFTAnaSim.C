#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "include/MFTAnaSim.h"

#include "SimulationDataFormat/MCEventHeader.h"

#endif

using MFTAnaSimTrack = o2::mftana::MFTAnaSimTrack;
using MFTAnaSimCluster = o2::mftana::MFTAnaSimCluster;
using MFTAnaSimHit = o2::mftana::MFTAnaSimHit;

void ReadMFTAnaSim()
{
  TFile inFile("MFTAnaSimTracks.root");
  TTree *tree1 = (TTree*)inFile.Get("MFTAnaSimTrack");
  std::vector<MFTAnaSimTrack> anaSimTracks, *anaSimTracksP = &anaSimTracks;
  tree1->SetBranchAddress("MFTAnaSimTrack", &anaSimTracksP);
  TTree *tree2 = (TTree*)inFile.Get("MFTAnaSimHit");
  std::vector<MFTAnaSimHit> anaSimHits, *anaSimHitsP = &anaSimHits;
  tree2->SetBranchAddress("MFTAnaSimHit", &anaSimHitsP);
  TTree *tree3 = (TTree*)inFile.Get("MFTAnaSimCluster");
  std::vector<MFTAnaSimCluster> anaSimClusters, *anaSimClustersP = &anaSimClusters;
  tree3->SetBranchAddress("MFTAnaSimCluster", &anaSimClustersP);

  Int_t nEntries1 = tree1->GetEntries();
  Int_t nEntries2 = tree2->GetEntries();
  Int_t nEntries3 = tree3->GetEntries();
  for (Int_t entry = 0; entry < nEntries1; entry++) {
    tree1->GetEntry(entry);
    tree2->GetEntry(entry);
    tree3->GetEntry(entry);
    Int_t iTrack = 0;
    for (auto& track : anaSimTracks) {
      
      printf("Event %d Track %d NDisks %d NLayers %d ", entry, iTrack, track.getNDisks(), track.getNLayers());
      printf("layers: ");
      for (Int_t index = 0; index < track.getNLayers(); index++) {
	printf("%d ", track.getLayer(index));
      }
      printf("\n");
      
      // hits
      printf("NHits %d Hit range: %d %d \n", track.getNHits(), track.getFirstHitIndex(), track.getLastHitIndex());
      if (track.getFirstHitIndex() >= 0) {
	for (Int_t iHit = track.getFirstHitIndex(); iHit <= track.getLastHitIndex(); iHit++) {
	  auto& hit = anaSimHits.at(iHit);
	  printf("%5d:   %7.3f   %7.3f   %7.3f \n", iHit, hit.getX(), hit.getY(), hit.getZ());
	}
      }
      
      // clusters
      printf("NClusters %d Cluster range: %d %d \n", track.getNClusters(), track.getFirstClusterIndex(), track.getLastClusterIndex());
      if (track.getFirstClusterIndex() >= 0) {
	for (Int_t iCls = track.getFirstClusterIndex(); iCls <= track.getLastClusterIndex(); iCls++) {
	  auto& cluster = anaSimClusters.at(iCls);
	  printf("%5d:   %7.3f   %7.3f   %7.3f \n", iCls, cluster.getX(), cluster.getY(), cluster.getZ());
	}
      } else {
	printf("... no clusters ...\n");
      }
      
      iTrack++;
    }
  }
  
}
