#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "include/MFTAnaSim.h"

#include "SimulationDataFormat/MCEventHeader.h"
#include "MFTBase/GeometryTGeo.h"

#endif

void MFTAna(const char *hitsFileName = "../tracking/36/o2sim_HitsMFT.root",
	    const char *kineFileName = "../tracking/36/o2sim_Kine.root",
	    const char *clusFileName = "../tracking/36/mftclusters.root",
	    const char *tracksFileName = "../tracking/36/mfttracks.root",
	    const char *geomFileName = "../tracking/36/o2sim_geometry.root")
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
  
  anaSim.setVerboseLevel(1);

  // read the input trees and dimension internal vector containers
  if (!anaSim.initialize()) {
    printf("MFTAnaSim::initialize returns false!\n");
    return;
  }

  for (int event = 0; event < nrEvents ; event++) {
    //printf("Analyze event %d.\n", event);
    
    anaSim.initEvent(event, kAll);
    
    anaSim.doHits();
    anaSim.doParticles();
    anaSim.doMCTracks();
    
    anaSim.finishEvent();
  }
  
  if (false) {
    auto particles = anaSim.getParticles();
    printf("Particles %zu\n", particles.size());
    for (auto& part : particles) {
      printf("%6d   %16s   %5d   %3d \n", part.mPDGCode, part.mPDGName.c_str(), part.mCount, part.mEvent);
    }
  }
    
  anaSim.doSATracks();
  
  anaSim.finish();
}
