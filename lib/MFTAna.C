#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "include/MFTAnaSim.h"

#include "SimulationDataFormat/MCEventHeader.h"
#include "MFTBase/GeometryTGeo.h"

#include <sys/time.h>

#endif

//_____________________________________________________________________________
void MFTAna(std::string prefix = "", 
	    std::string hitsFileName = "o2sim_HitsMFT.root",
	    std::string kineFileName = "o2sim_Kine.root",
	    std::string clusFileName = "mftclusters.root",
	    std::string tracksFileName = "mfttracks.root",
	    std::string geomFileName = "o2sim_geometry.root",
	    std::string dictFileName = "MFTdictionary.bin")
{
  hitsFileName = prefix + hitsFileName;
  kineFileName = prefix + kineFileName;
  clusFileName = prefix + clusFileName;
  tracksFileName = prefix + tracksFileName;
  geomFileName = prefix + geomFileName;
  dictFileName = prefix + dictFileName;
  printf("Using: \n");
  printf("%s \n", hitsFileName.data());
  printf("%s \n", kineFileName.data());
  printf("%s \n", clusFileName.data());
  printf("%s \n", tracksFileName.data());
  printf("%s \n", geomFileName.data());
  printf("%s \n", dictFileName.data());

  struct timeval tvStart, tvEnd;

  enum ParticleSource {kPrimary, kSecondary, kAll};
  
  // create the main analysis class
  auto anaSim = o2::mftana::MFTAnaSim();

  anaSim.mNThreads = 4;

  anaSim.setTopoDictFileName(dictFileName);
  
  // geometry manager
  o2::base::GeometryManager::loadGeometry(geomFileName.data());
  auto gman = o2::mft::GeometryTGeo::Instance();
  gman->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
  anaSim.mGeoManager = gman;
  
  // kinematics (MC tracks, particles)
  TFile kineFile(kineFileName.data());
  TTree *kineTree = (TTree*)kineFile.Get("o2sim");
  o2::dataformats::MCEventHeader* eventHeader = nullptr;
  kineTree->SetBranchAddress("MCEventHeader.", &eventHeader);
  int nrEvents = kineTree->GetEntries();
  anaSim.mKineTree = kineTree;
  
  // hits
  TFile hitsFile(hitsFileName.data());
  TTree* hitTree = (TTree*)hitsFile.Get("o2sim");
  anaSim.mHitTree = hitTree;
  
  // clusters
  TFile clusFile(clusFileName.data());
  TTree *clusTree = (TTree*)clusFile.Get("o2sim");
  anaSim.mClusTree = clusTree;
  
  // reconstructed tracks, MFT standalone (SA)
  TFile tracksFile(tracksFileName.data());
  TTree *trackTree = (TTree*)tracksFile.Get("o2sim");
  anaSim.mTrackTree = trackTree;
  
  anaSim.setVerboseLevel(0);

  // execution timing
  gettimeofday(&tvStart, NULL);

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

  // print particles summary
  if (false) {
    auto particles = anaSim.getParticles();
    printf("Particles %zu\n", particles.size());
    for (auto& part : particles) {
      printf("%6d   %16s   %5d   %3d \n", part.mPDGCode, part.mPDGName.c_str(), part.mCount, part.mEvent);
    }
  }
    
  anaSim.doSATracks();

  anaSim.linkTracks();
 
  anaSim.finish();

  // execution timing
  gettimeofday(&tvEnd, NULL);
  std::cout << "Time " << (tvEnd.tv_usec - tvStart.tv_usec) / 1000000.0 + (tvEnd.tv_sec - tvStart.tv_sec) << " [seconds]" << std::endl;

  kineFile.Close();
  hitsFile.Close();
  clusFile.Close();
  tracksFile.Close();
}
