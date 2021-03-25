#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "include/MFTAnaSim.h"

#include "SimulationDataFormat/MCEventHeader.h"

#endif

void MFTAna(const Char_t *hitsFileName = "o2sim_HitsMFT.root",
	    const Char_t *kineFileName = "o2sim_Kine.root")
{

  auto anaSim = o2::mftana::MFTAnaSim();

  // kinematics, MC tracks
  TFile kineFile(kineFileName);
  TTree *kineTree = (TTree*)kineFile.Get("o2sim");
  o2::dataformats::MCEventHeader* eventHeader = nullptr;
  kineTree->SetBranchAddress("MCEventHeader.", &eventHeader);
  Int_t nrEvents = kineTree->GetEntries();
  anaSim.mKineTree = kineTree;
  
  // hits
  TFile hitsFile(hitsFileName);
  TTree* hitTree = (TTree*)hitsFile.Get("o2sim");
  anaSim.mHitTree = hitTree;

  Int_t nHits, nMCTracks, maxMCTracks = 0;
  for (Int_t event = 0; event < nrEvents ; event++) {
    kineTree->GetEntry(event);
    nMCTracks = eventHeader->getMCEventStats().getNKeptTracks();
    maxMCTracks = std::max(maxMCTracks, nMCTracks);
  }

  anaSim.initialize(maxMCTracks);
  
  for (Int_t event = 0; event < nrEvents ; event++) {
    kineTree->GetEntry(event);
    nMCTracks = eventHeader->getMCEventStats().getNKeptTracks();
    printf("Analyze event %d with %d MC tracks.\n", event, nMCTracks);
    anaSim.initEvent(event, nMCTracks);
    anaSim.doHits();
    anaSim.doMCTracks();
  }
}
