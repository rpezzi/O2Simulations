#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "include/MFTAnaSim.h"
#include "include/MFTAnaSimCluster.h"
#include "include/MFTAnaSimHit.h"
#include "include/MFTAnaSimTrack.h"

#include "SimulationDataFormat/MCEventHeader.h"
#include "MFTBase/GeometryTGeo.h"

#include "src/MFTAnaSimLinkDef.h"

#endif

void MFTAna(const char *sim_path = "./") {
  enum ParticleSource {kPrimary, kSecondary, kAll};

  std::string hitsFileName = sim_path + std::string("/o2sim_HitsMFT.root");
  std::string kineFileName = sim_path + std::string("/o2sim_Kine.root");
  std::string clusFileName = sim_path + std::string("/mftclusters.root");
  std::string tracksFileName = sim_path + std::string("/mfttracks.root");
  std::string geomFileName = sim_path + std::string("/o2sim");

  // create the main analysis class
  auto anaSim = o2::mftana::MFTAnaSim();

  // geometry manager
  o2::base::GeometryManager::loadGeometry(geomFileName.c_str());
  auto gman = o2::mft::GeometryTGeo::Instance();
  gman->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
  anaSim.mGeoManager = gman;
  
  // kinematics (MC tracks, particles)
  TFile kineFile(kineFileName.c_str());
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
  TFile hitsFile(hitsFileName.c_str());
  TTree* hitTree = (TTree*)hitsFile.Get("o2sim");
  anaSim.mHitTree = hitTree;
  
  // clusters
  TFile clusFile(clusFileName.c_str());
  TTree *clusTree = (TTree*)clusFile.Get("o2sim");
  anaSim.mClusTree = clusTree;
  
  // reconstructed tracks, MFT standalone (SA)
  TFile tracksFile(tracksFileName.c_str());
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
