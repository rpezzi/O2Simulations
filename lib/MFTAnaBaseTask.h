#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <vector>
#include <map>

#include "include/MFTAnaSimHit.h"
#include "include/MFTAnaSimCluster.h"
#include "include/MFTAnaSimTrack.h"
#include "include/MFTAnaSimMCTrack.h"
#include "include/MFTAnaSimSATrack.h"

#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#endif

class MFTAnaBaseTask
{
public:
  void init(std::string run = "")
  {

    // file with histograms
    std::string outfilename = "ReadMFTAnaSim" + run + ".root";
    outFile = new TFile(outfilename.c_str(), "RECREATE");
    createHistograms();

    // input file (from MFTAna.C)
    std::string inpfilename = "MFTAnaSimTracks" + run + ".root";
    inFile = new TFile(inpfilename.c_str());

    tree1 = (TTree *)inFile->Get("MFTAnaSimTrack");
    tree1->SetBranchAddress("MFTAnaSimTrack", &anaSimTracksP);

    tree2 = (TTree *)inFile->Get("MFTAnaSimHit");
    tree2->SetBranchAddress("MFTAnaSimHit", &anaSimHitsP);

    tree3 = (TTree *)inFile->Get("MFTAnaSimCluster");
    tree3->SetBranchAddress("MFTAnaSimCluster", &anaSimClustersP);

    tree4 = (TTree *)inFile->Get("MFTAnaSimSATrack");
    tree4->SetBranchAddress("MFTAnaSimSATrack", &anaSimSATracksP);

    tree1->GetEntry(0);
    tree2->GetEntry(0);
    tree3->GetEntry(0);
    tree4->GetEntry(0);
    outFile->cd();
  }

  void finalize()
  {
    inFile->Close();

    // save histograms
    outFile->cd();

    for (auto &h : TH1HistoCodes)
    {
      TH1Histos[h]->Write();
    }

    for (auto &h : TH2HistoCodes)
    {
      TH2Histos[h]->Write();
    }

    for (auto &h : TH3HistoCodes)
    {
      TH3Histos[h]->Write();
    }
    outFile->Close();
  }
  void run(std::string run = "")
  {
    std::cout << "MFTAna task: init." << std::endl;
    init(run);
    std::cout << "MFTAna task: process." << std::endl;
    process();
    std::cout << "MFTAna task: finalize." << std::endl;
    finalize();
    std::cout << "MFTAna task: done." << std::endl;
  }

  virtual void process();

  std::vector<std::string> TH1HistoCodes{};
  void setTH1HistoCodes(std::vector<std::string> v) { TH1HistoCodes = v; };

  std::map<std::string, const char *> TH1Names{};
  void setTH1Names(std::map<std::string, const char *> v) { TH1Names = v; };

  std::map<std::string, const char *> TH1Titles{};
  void setTH1Titles(std::map<std::string, const char *> v) { TH1Titles = v; };

  std::map<std::string, std::array<double, 3>> TH1Binning{};
  void setTH1Binning(std::map<std::string, std::array<double, 3>> v) { TH1Binning = v; };

  std::map<std::string, const char *> TH1XaxisTitles{};
  void setTH1XaxisTitles(std::map<std::string, const char *> v) { TH1XaxisTitles = v; };

  std::vector<std::string> TH2HistoCodes{};
  void setTH2HistoCodes(std::vector<std::string> v) { TH2HistoCodes = v; };

  std::map<std::string, const char *> TH2Names{};
  void setTH2Names(std::map<std::string, const char *> v) { TH2Names = v; };

  std::map<std::string, const char *> TH2Titles{};
  void setTH2Titles(std::map<std::string, const char *> v) { TH2Titles = v; };

  std::map<std::string, std::array<double, 6>> TH2Binning{};
  void setTH2Binning(std::map<std::string, std::array<double, 6>> v) { TH2Binning = v; };

  std::map<std::string, const char *> TH2XaxisTitles{};
  void setTH2XaxisTitles(std::map<std::string, const char *> v) { TH2XaxisTitles = v; };

  std::map<std::string, const char *> TH2YaxisTitles{};
  void setTH2YaxisTitles(std::map<std::string, const char *> v) { TH2YaxisTitles = v; };

  std::vector<std::string> TH3HistoCodes{};
  void setTH3HistoCodes(std::vector<std::string> v) { TH3HistoCodes = v; };

  std::map<std::string, const char *> TH3Names{};
  void setTH3Names(std::map<std::string, const char *> v) { TH3Names = v; };

  std::map<std::string, const char *> TH3Titles{};
  void setTH3Titles(std::map<std::string, const char *> v) { TH3Titles = v; };

  std::map<std::string, std::array<double, 9>> TH3Binning{};
  void setTH3Binning(std::map<std::string, std::array<double, 9>> v) { TH3Binning = v; };

  std::map<std::string, const char *> TH3XaxisTitles{};
  void setTH3XaxisTitles(std::map<std::string, const char *> v) { TH3XaxisTitles = v; };

  std::map<std::string, const char *> TH3YaxisTitles{};
  void setTH3YaxisTitles(std::map<std::string, const char *> v) { TH3YaxisTitles = v; };

  std::map<std::string, const char *> TH3ZaxisTitles{};
  void setTH3ZaxisTitles(std::map<std::string, const char *> v) { TH3ZaxisTitles = v; };

  std::map<std::string, TH1F *> TH1Histos;
  std::map<std::string, TH2F *> TH2Histos;
  std::map<std::string, TH3F *> TH3Histos;

  std::vector<o2::mftana::MFTAnaSimTrack> anaSimTracks, *anaSimTracksP = &anaSimTracks;
  std::vector<o2::mftana::MFTAnaSimHit> anaSimHits, *anaSimHitsP = &anaSimHits;
  std::vector<o2::mftana::MFTAnaSimCluster> anaSimClusters, *anaSimClustersP = &anaSimClusters;
  std::vector<o2::mftana::MFTAnaSimSATrack> anaSimSATracks, *anaSimSATracksP = &anaSimSATracks;
  TTree *tree1, *tree2, *tree3, *tree4;
  TFile *outFile, *inFile;

  //_____________________________________________________________________________
  void createHistograms()
  {

    // 1D histograms
    for (auto &h : TH1HistoCodes)
    {
      std::cout << "Creating TH1: " << h << std::endl;
      TH1Histos[h] = new TH1F(TH1Names[h], TH1Titles[h], (int)TH1Binning[h][0], TH1Binning[h][1], TH1Binning[h][2]);
      TH1Histos[h]->GetXaxis()->SetTitle(TH1XaxisTitles[h]);
    }

    // 2D histograms
    for (auto &h : TH2HistoCodes)
    {
      std::cout << "Creating TH2: " << h << std::endl;
      TH2Histos[h] = new TH2F(TH2Names[h], TH2Titles[h],
                              (int)TH2Binning[h][0], TH2Binning[h][1], TH2Binning[h][2],
                              (int)TH2Binning[h][3], TH2Binning[h][4], TH2Binning[h][5]);
      TH2Histos[h]->GetXaxis()->SetTitle(TH2XaxisTitles[h]);
      TH2Histos[h]->GetYaxis()->SetTitle(TH2YaxisTitles[h]);
    }

    // 3D histograms
    for (auto &h : TH3HistoCodes)
    {
      std::cout << "Creating TH3: " << h << std::endl;
      TH3Histos[h] = new TH3F(TH3Names[h], TH3Titles[h],
                              (int)TH3Binning[h][0], TH3Binning[h][1], TH3Binning[h][2],
                              (int)TH3Binning[h][3], TH3Binning[h][4], TH3Binning[h][5],
                              (int)TH3Binning[h][6], TH3Binning[h][7], TH3Binning[h][8]);
      TH3Histos[h]->GetXaxis()->SetTitle(TH3XaxisTitles[h]);
      TH3Histos[h]->GetYaxis()->SetTitle(TH3YaxisTitles[h]);
      TH3Histos[h]->GetZaxis()->SetTitle(TH3ZaxisTitles[h]);
    }
  }

  //_____________________________________________________________________________
  void printSATrack(o2::mftana::MFTAnaSimSATrack *saTrack, int i)
  {

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
    if (nPoints > 0)
    {
      for (int i_p = 0; i_p < nPoints; i_p++)
      {
        clusEntry = saTrack->getIntClusIndex(i_p);
        auto &cluster = anaSimClusters.at(clusEntry);
        printf("%5d  %5d:   %7.3f   %7.3f   %7.3f   %d \n", i_p, clusEntry, cluster.getX(), cluster.getY(), cluster.getZ(), cluster.getLayer());
      }
    }
    else
    {
      printf("... no clusters ...\n");
    }

    printf("---------------------------------------------------------------\n");
    for (int iMCTrack = 0; iMCTrack < saTrack->getNMCTracks(); iMCTrack++)
    {
      printf("MCTrack %d (%d) evn %d \n", saTrack->getIntMCTrackIndex(iMCTrack), saTrack->getIntMCTrackMult(iMCTrack), saTrack->getIntMCTrackEvent(iMCTrack));
    }
  }

  //_____________________________________________________________________________
  void printMCTrack(o2::mftana::MFTAnaSimTrack *mcTrack, int i)
  {

    int nClusters;
    printf("===============================================================\n");
    printf("Event %d Track %d NDisks %d NLayers %d MC %d ", mcTrack->getEvent(), i, mcTrack->getNDisks(), mcTrack->getNLayers(), mcTrack->getMCTrackID());
    printf("layers: ");
    for (int index = 0; index < mcTrack->getNLayers(); index++)
    {
      printf("%d ", mcTrack->getLayer(index));
    }
    printf("\n");

    // hits
    printf("---------------------------------------------------------------\n");
    printf("Number of hits %d , hit range: %d %d \n", mcTrack->getNHits(), mcTrack->getFirstHitIndex(), mcTrack->getLastHitIndex());
    if (mcTrack->getFirstHitIndex() >= 0)
    {
      for (int iHit = mcTrack->getFirstHitIndex(); iHit <= mcTrack->getLastHitIndex(); iHit++)
      {
        auto &hit = anaSimHits.at(iHit);
        printf("%5d:   %7.3f   %7.3f   %7.3f \n", iHit, hit.getX(), hit.getY(), hit.getZ());
      }
    }

    // clusters
    nClusters = mcTrack->getNClusters();
    printf("---------------------------------------------------------------\n");
    printf("Number of clusters %d \n", nClusters);
    if (nClusters > 0)
    {
      for (int i_cls = 0; i_cls < nClusters; i_cls++)
      {
        auto &cluster = anaSimClusters.at(mcTrack->getIntClusIndex(i_cls));
        printf("%5d:   %7.3f   %7.3f   %7.3f \n", i_cls, cluster.getX(), cluster.getY(), cluster.getZ());
      }
    }
    else
    {
      printf("... no clusters ...\n");
    }

    printf("---------------------------------------------------------------\n");
    for (int iSATrack = 0; iSATrack < mcTrack->getNSATracks(); iSATrack++)
    {
      printf("SATrack %d (%d) \n", mcTrack->getIntSATrackIndex(iSATrack), mcTrack->getIntSATrackMult(iSATrack));
    }
  }
};

//_____________________________________________________________________________
inline void MFTAnaBaseTask::process()
{

std::cout << "Error! This virtual method should be overriden.." << std::endl;
return;
}
