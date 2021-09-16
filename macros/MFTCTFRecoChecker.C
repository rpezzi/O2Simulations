#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#endif


//_________________________________________________________________________________________________
void MFTCTFRecoChecker(const Char_t* trkFile = "mfttracks.root", const Char_t* clsFile = "mftclusters.root")
{

  // MFT Tracks
  TFile* trkFileIn = new TFile(trkFile);
  TTree* mftTrackTree = (TTree*)trkFileIn->Get("o2sim");

  // MFT Clusters
  TFile* clsFileIn = new TFile(clsFile);
  TTree* mftClsTree = (TTree*)clsFileIn->Get("o2sim");

  std::vector<o2::itsmft::CompClusterExt> compClusMFTVec, *compClusMFTVecP = &compClusMFTVec;
  mftClsTree->SetBranchAddress("MFTClusterComp", &compClusMFTVecP);


  auto c = new TCanvas();
  mftClsTree->Draw("MFTClustersROF.mROFEntry.mEntries>>MFTClusters","","");
  c->Print("MFTClustersROF.png");

  c = new TCanvas();
  mftClsTree->Draw("MFTClusterComp.mChipID>>ChipIDs","","");
  c->Print("MFTClustersChipIDs.png");

  c = new TCanvas();
  mftTrackTree->Draw("MFTTracksROF.mROFEntry.mEntries>>MFTTracks","","");
  c->Print("MFTTracksROFs.png");

}
