#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "include/MFTAnaSim.h"

#include "SimulationDataFormat/MCEventHeader.h"
#include "MFTBase/GeometryTGeo.h"

#endif

void MFTAna(const char *hitsFileName = "../tracking/31/o2sim_HitsMFT.root",
	    const char *kineFileName = "../tracking/31/o2sim_Kine.root",
	    const char *clusFileName = "../tracking/31/mftclusters.root",
	    const char *tracksFileName = "../tracking/31/mfttracks.root",
	    const char *geomFileName = "../tracking/31/o2sim_geometry.root")
{
  enum ParticleSource {kPrimary, kSecondary, kAll};
  
  // create the main analysis class
  auto anaSim = o2::mftana::MFTAnaSim();
  
  // geometry manager
  o2::base::GeometryManager::loadGeometry(geomFileName);
  auto gman = o2::mft::GeometryTGeo::Instance();
  gman->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
  anaSim.mGeoManager = gman;
  
  // kinematics (MC tracks, particles)
  TFile kineFile(kineFileName);
  TTree *kineTree = (TTree*)kineFile.Get("o2sim");
  o2::dataformats::MCEventHeader* eventHeader = nullptr;
  kineTree->SetBranchAddress("MCEventHeader.", &eventHeader);
  int nrEvents = kineTree->GetEntries();
  anaSim.mKineTree = kineTree;
  
  // maximum MC tracks per event, from a scan of all events in the file
  int nMCTracks, maxMCTracks = 0;
  for (int event = 0; event < nrEvents ; event++) {
    kineTree->GetEntry(event);
    nMCTracks = eventHeader->getMCEventStats().getNKeptTracks();
    maxMCTracks = std::max(maxMCTracks, nMCTracks);
  }
  
  // hits
  TFile hitsFile(hitsFileName);
  TTree* hitTree = (TTree*)hitsFile.Get("o2sim");
  anaSim.mHitTree = hitTree;
  
  // clusters
  TFile clusFile(clusFileName);
  TTree *clusTree = (TTree*)clusFile.Get("o2sim");
  anaSim.mClusTree = clusTree;
  
  // reconstructed tracks, MFT standalone (SA)
  TFile tracksFile(tracksFileName);
  TTree *trackTree = (TTree*)tracksFile.Get("o2sim");
  anaSim.mTrackTree = trackTree;
  
  // read the input trees and dimension internal vector containers
  if (!anaSim.initialize(maxMCTracks)) {
    printf("MFTAnaSim::initialize returns false!\n");
    return;
  }
  
  for (int event = 0; event < nrEvents ; event++) {
    kineTree->GetEntry(event);
    nMCTracks = eventHeader->getMCEventStats().getNKeptTracks();
    printf("Analyze event %d with %d MC tracks.\n", event, nMCTracks);
    
    anaSim.initEvent(event, nMCTracks, kSecondary);
    
    anaSim.doParticles();
    
    if (kTRUE || event == 0) {
      auto particles = anaSim.getParticles();
      printf("and %zu particle species.\n", particles.size());
      for (auto& part : particles) {
	printf("%6d   %16s   %5d \n", part.mPDGCode, part.mPDGName.c_str(), part.mCount);
      }
    }
    
    anaSim.doHits();
    anaSim.doMCTracks();
    anaSim.finishEvent();
  }
  
  anaSim.doSATracks();
  
  anaSim.finish();
}
