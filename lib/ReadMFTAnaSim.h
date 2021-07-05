#include <vector>
#include <map>

#include "include/MFTAnaSimHit.h"
#include "include/MFTAnaSimCluster.h"
#include "include/MFTAnaSimTrack.h"
#include "include/MFTAnaSimMCTrack.h"
#include "include/MFTAnaSimSATrack.h"

#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>

enum TH1HistoCodes {
  MCNrOfHits,
  MCall_Pt,
  MCinSA_Pt,
  NHistograms
};

std::map<int, const char *> TH1Names {
  {MCNrOfHits, "MC number of hits"},
  {MCall_Pt, "Pt pf MC tracks with hits"},
  {MCinSA_Pt, "Pt of MC tracks with clusters in SA tracks"}
};

std::map<int, const char *> TH1Titles {
  {MCNrOfHits, "MC number of hits"},
  {MCall_Pt, "Pt pf MC tracks with hits"},
  {MCinSA_Pt, "Pt of MC tracks with clusters in SA tracks"}
};

std::map<int, std::array<double, 3>> TH1Binning {
  {MCNrOfHits, {50, 0., 50.}},
  {MCall_Pt, {100, 0., 5.}},
  {MCinSA_Pt, {100, 0., 5.}}
};

std::map<int, const char *> TH1XaxisTitles {
  {MCNrOfHits, "Nr of hits"},
  {MCall_Pt, "pt [GeV/c]"},
  {MCinSA_Pt, "pt [GeV/c]"}
};

std::vector<TH1F*> TH1Histos(NHistograms);

std::vector<o2::mftana::MFTAnaSimTrack> anaSimTracks, *anaSimTracksP = &anaSimTracks;
std::vector<o2::mftana::MFTAnaSimHit> anaSimHits, *anaSimHitsP = &anaSimHits;
std::vector<o2::mftana::MFTAnaSimCluster> anaSimClusters, *anaSimClustersP = &anaSimClusters;
std::vector<o2::mftana::MFTAnaSimSATrack> anaSimSATracks, *anaSimSATracksP = &anaSimSATracks;

//_____________________________________________________________________________
void createHistograms() {  
  auto nHisto = 0;
  for (auto& h : TH1Histos) {
    h = new TH1F(TH1Names[nHisto], TH1Titles[nHisto], (int)TH1Binning[nHisto][0], TH1Binning[nHisto][1], TH1Binning[nHisto][2]);
    h->GetXaxis()->SetTitle(TH1XaxisTitles[nHisto]);
    ++nHisto;
  }
}

//_____________________________________________________________________________
void printSATrack(o2::mftana::MFTAnaSimSATrack* saTrack, int i) {
  
  int nPoints, clusEntry;
  auto trkInX = saTrack->getX();
  auto trkInY = saTrack->getY();
  auto trkInZ = saTrack->getZ();
  auto outParam = saTrack->getOutParam();
  auto trkOutX = outParam.getX();
  auto trkOutY = outParam.getY();
  auto trkOutZ = outParam.getZ();
  printf("---------------------------------------------------------------\n");
  printf("SATrack %d points %d disks %d layers %d  x,y,z-in  %7.3f  %7.3f  %7.3f  x,y,z-out  %7.3f  %7.3f  %7.3f  pt  %f\n", i, saTrack->getNPoints(), saTrack->getNDisks(), saTrack->getNLayers(), trkInX, trkInY, trkInZ, trkOutX, trkOutY, trkOutZ, outParam.getPt());
  nPoints = saTrack->getNPoints();
  printf("Number of points %d \n", nPoints);
  if (nPoints > 0) {
    for (int i_p = 0; i_p < nPoints; i_p++) {
      clusEntry = saTrack->getIntClusIndex(i_p);
      auto& cluster = anaSimClusters.at(clusEntry);
      printf("%5d  %5d:   %7.3f   %7.3f   %7.3f   %d \n", i_p, clusEntry, cluster.getX(), cluster.getY(), cluster.getZ(), cluster.getLayer());
    }
  } else {
    printf("... no clusters ...\n");
  }
  
  printf("---------------------------------------------------------------\n");
  for (int iMCTrack = 0; iMCTrack < saTrack->getNMCTracks(); iMCTrack++) {
    printf("MCTrack %d (%d) evn %d \n", saTrack->getIntMCTrackIndex(iMCTrack), saTrack->getIntMCTrackMult(iMCTrack), saTrack->getIntMCTrackEvent(iMCTrack));  
  }
  
}

//_____________________________________________________________________________
void printMCTrack(o2::mftana::MFTAnaSimTrack* mcTrack, int i) {

  int nClusters;
  printf("===============================================================\n");
  printf("Event %d Track %d NDisks %d NLayers %d MC %d ", mcTrack->getEvent(), i, mcTrack->getNDisks(), mcTrack->getNLayers(), mcTrack->getMCTrackID());
  printf("layers: ");
  for (int index = 0; index < mcTrack->getNLayers(); index++) {
    printf("%d ", mcTrack->getLayer(index));
  }
  printf("\n");
  
  // hits
  printf("---------------------------------------------------------------\n");
  printf("Number of hits %d , hit range: %d %d \n", mcTrack->getNHits(), mcTrack->getFirstHitIndex(), mcTrack->getLastHitIndex());
  if (mcTrack->getFirstHitIndex() >= 0) {
    for (int iHit = mcTrack->getFirstHitIndex(); iHit <= mcTrack->getLastHitIndex(); iHit++) {
      auto& hit = anaSimHits.at(iHit);
      printf("%5d:   %7.3f   %7.3f   %7.3f \n", iHit, hit.getX(), hit.getY(), hit.getZ());
    }
  }
  
  // clusters
  nClusters = mcTrack->getNClusters();
  printf("---------------------------------------------------------------\n");
  printf("Number of clusters %d \n", nClusters);
  if (nClusters > 0) {
    for (int i_cls = 0; i_cls < nClusters; i_cls++) {
      auto& cluster = anaSimClusters.at(mcTrack->getIntClusIndex(i_cls));
      printf("%5d:   %7.3f   %7.3f   %7.3f \n", i_cls, cluster.getX(), cluster.getY(), cluster.getZ());
    }
  } else {
    printf("... no clusters ...\n");
  }
  
  printf("---------------------------------------------------------------\n");
  for (int iSATrack = 0; iSATrack < mcTrack->getNSATracks(); iSATrack++) {
    printf("SATrack %d (%d) \n", mcTrack->getIntSATrackIndex(iSATrack), mcTrack->getIntSATrackMult(iSATrack));  
  }

}
