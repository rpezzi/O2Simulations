#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "include/MFTAnaSim.h"

#include "SimulationDataFormat/MCEventHeader.h"

#endif

void MFTAna(const Char_t *hitsFileName = "o2sim_HitsMFT.root",
	    const Char_t *kineFileName = "o2sim_Kine.root")
{

  enum ParticleSource {kPrimary, kSecondary, kAll};
  
  auto anaSim = o2::mftana::MFTAnaSim();

  // kinematics (MC tracks, particles)
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

  // maximum MC tracks per event, from a scan of all events in the file
  Int_t nMCTracks, maxMCTracks = 0;
  for (Int_t event = 0; event < nrEvents ; event++) {
    kineTree->GetEntry(event);
    nMCTracks = eventHeader->getMCEventStats().getNKeptTracks();
    maxMCTracks = std::max(maxMCTracks, nMCTracks);
  }
  anaSim.initialize(maxMCTracks);
  
  Int_t nHits;
  for (Int_t event = 0; event < nrEvents ; event++) {
    kineTree->GetEntry(event);
    nMCTracks = eventHeader->getMCEventStats().getNKeptTracks();
    printf("Analyze event %d with %d MC tracks.\n", event, nMCTracks);
    
    anaSim.initEvent(event, nMCTracks, kSecondary);
    
    anaSim.doParticles();
    if (kTRUE || event == 0) {
      auto particles = anaSim.getParticles();
      printf("and %zu particles.\n", particles.size());
      for (auto& part : particles) {
	printf("%6d   %16s   %5d \n", part.mPDGCode, part.mPDGName.c_str(), part.mCount);
      }
    }
    
    anaSim.doHits();

    anaSim.doMCTracks();
  }
}
