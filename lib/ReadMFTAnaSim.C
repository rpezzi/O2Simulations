#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "include/MFTAnaSim.h"
#include "ReadMFTAnaSim.h"
#include "SimulationDataFormat/MCEventHeader.h"

#endif

using MFTAnaSimSATrack = o2::mftana::MFTAnaSimSATrack;
using MFTAnaSimTrack = o2::mftana::MFTAnaSimTrack;
using MFTAnaSimCluster = o2::mftana::MFTAnaSimCluster;
using MFTAnaSimHit = o2::mftana::MFTAnaSimHit;

//_____________________________________________________________________________
void ReadMFTAnaSim(std::string run = "36")
{
  int verbose = 0;
  
  std::string outfilename = "ReadMFTAnaSim_" + run + ".root";
  TFile outFile(outfilename.c_str(),"RECREATE");
  createHistograms();
  
  std::string inpfilename = "MFTAnaSimTracks_" + run + ".root";
  TFile inFile(inpfilename.c_str());
  
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

  // MC tracks
  int nPoints, nClusters, iTrack = 0;
  for (auto& track : anaSimTracks) {

    printf("===============================================================\n");
    printf("Event %d Track %d NDisks %d NLayers %d MC %d ", track.getEvent(), iTrack, track.getNDisks(), track.getNLayers(), track.getMCTrackID());
    printf("layers: ");
    for (int index = 0; index < track.getNLayers(); index++) {
      printf("%d ", track.getLayer(index));
    }
    printf("\n");
    
    // hits
    printf("---------------------------------------------------------------\n");
    printf("Number of hits %d , hit range: %d %d \n", track.getNHits(), track.getFirstHitIndex(), track.getLastHitIndex());
    if (track.getFirstHitIndex() >= 0) {
      for (int iHit = track.getFirstHitIndex(); iHit <= track.getLastHitIndex(); iHit++) {
	auto& hit = anaSimHits.at(iHit);
	printf("%5d:   %7.3f   %7.3f   %7.3f \n", iHit, hit.getX(), hit.getY(), hit.getZ());
      }
    }
    
    TH1Histos[MCNrOfHits]->Fill(track.getNHits());
    TH1Histos[MCHasRecPt]->Fill(track.getPt());
    
    // clusters
    nClusters = track.getNClusters();
    printf("---------------------------------------------------------------\n");
    printf("Number of clusters %d \n", nClusters);
    if (nClusters > 0) {
      for (int i_cls = 0; i_cls < nClusters; i_cls++) {
	auto& cluster = anaSimClusters.at(track.getIntClusIndex(i_cls));
	printf("%5d:   %7.3f   %7.3f   %7.3f \n", i_cls, cluster.getX(), cluster.getY(), cluster.getZ());
      }
    } else {
      printf("... no clusters ...\n");
    }
    
    printf("---------------------------------------------------------------\n");
    for (int iSATrack = 0; iSATrack < track.getNSATracks(); iSATrack++) {
      printf("SATrack %d (%d) \n", track.getIntSATrackIndex(iSATrack), track.getIntSATrackMult(iSATrack));  
    }
    
    iTrack++;
  }

  // SA tracks
  int clusEntry, iSATrack = 0;
  for (auto& saTrack : anaSimSATracks) {
    auto trkInX = saTrack.getX();
    auto trkInY = saTrack.getY();
    auto trkInZ = saTrack.getZ();
    auto outParam = saTrack.getOutParam();
    auto trkOutX = outParam.getX();
    auto trkOutY = outParam.getY();
    auto trkOutZ = outParam.getZ();
    printf("---------------------------------------------------------------\n");
    printf("SATrack %d points %d disks %d layers %d  x,y,z-in  %7.3f  %7.3f  %7.3f  x,y,z-out  %7.3f  %7.3f  %7.3f  pt  %f\n", iSATrack, saTrack.getNPoints(), saTrack.getNDisks(), saTrack.getNLayers(), trkInX, trkInY, trkInZ, trkOutX, trkOutY, trkOutZ, outParam.getPt());
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
    
    printf("---------------------------------------------------------------\n");
    for (int iMCTrack = 0; iMCTrack < saTrack.getNMCTracks(); iMCTrack++) {
      printf("MCTrack %d (%d) evn %d \n", saTrack.getIntMCTrackIndex(iMCTrack), saTrack.getIntMCTrackMult(iMCTrack), saTrack.getIntMCTrackEvent(iMCTrack));  
    }
    
    iSATrack++;
  }
  inFile.Close();
   
  outFile.cd();
  for (auto& h : TH1Histos) {
    h->Write();
  }
  outFile.Close();
  
}
