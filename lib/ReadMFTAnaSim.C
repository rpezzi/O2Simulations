#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "ReadMFTAnaSim.h"
#include "SimulationDataFormat/MCEventHeader.h"

#endif

using MFTAnaSimSATrack = o2::mftana::MFTAnaSimSATrack;
using MFTAnaSimTrack = o2::mftana::MFTAnaSimTrack;
using MFTAnaSimCluster = o2::mftana::MFTAnaSimCluster;
using MFTAnaSimHit = o2::mftana::MFTAnaSimHit;

//_____________________________________________________________________________
void ReadMFTAnaSim(std::string run = "")
{
  int verbose = 0;

  // file with histograms
  std::string outfilename = "ReadMFTAnaSim" + run + ".root";
  TFile outFile(outfilename.c_str(),"RECREATE");
  createHistograms();

  // input file (from MFTAna.C)
  std::string inpfilename = "MFTAnaSimTracks" + run + ".root";
  TFile inFile(inpfilename.c_str());
  
  TTree *tree1 = (TTree*)inFile.Get("MFTAnaSimTrack");
  tree1->SetBranchAddress("MFTAnaSimTrack", &anaSimTracksP);
  
  TTree *tree2 = (TTree*)inFile.Get("MFTAnaSimHit");
  tree2->SetBranchAddress("MFTAnaSimHit", &anaSimHitsP);
  
  TTree *tree3 = (TTree*)inFile.Get("MFTAnaSimCluster");
  tree3->SetBranchAddress("MFTAnaSimCluster", &anaSimClustersP);
  
  TTree *tree4 = (TTree*)inFile.Get("MFTAnaSimSATrack");
  tree4->SetBranchAddress("MFTAnaSimSATrack", &anaSimSATracksP);

  tree1->GetEntry(0);
  tree2->GetEntry(0);
  tree3->GetEntry(0);
  tree4->GetEntry(0);

  // MC tracks
  int iMCTrack = 0;
  for (auto& mcTrack : anaSimTracks) {
    //printMCTrack(&mcTrack, iMCTrack++);
    
    TH1Histos[MCNrOfHits]->Fill(mcTrack.getNHits());
    TH1Histos[MCall_Pt]->Fill(mcTrack.getPt());
    if (mcTrack.getNSATracks() > 0) {
      TH1Histos[MCinSA_Pt]->Fill(mcTrack.getPt());
    }
    
  }

  // SA tracks
  int iSATrack = 0;
  for (auto& saTrack : anaSimSATracks) {
    //printSATrack(&saTrack, iSATrack++);
  }
  
  inFile.Close();

  // save histograms
  outFile.cd();
  for (auto& h : TH1Histos) {
    h->Write();
  }
  outFile.Close();
  
}
