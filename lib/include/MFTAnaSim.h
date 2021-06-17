/// @brief Main class for the analysis of Muon Forward Tracker simulations
/// @details Class to create at the beginning of a steering macro:
/// @details auto anaSim = o2::mftana::MFTAnaSim();

#ifndef MFT_ANA_SIM
#define MFT_ANA_SIM

#include <TFile.h>
#include <TTree.h>

#include "../include/MFTAnaSimTrack.h"
#include "../include/MFTAnaSimSATrack.h"
#include "../include/MFTAnaSimCluster.h"
#include "../include/MFTAnaSimHit.h"

#include "ITSMFTSimulation/Hit.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "MFTBase/Constants.h"
#include "MFTBase/GeometryTGeo.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"

namespace o2::mftana
{

struct MCPart {
  //MCPart() = default;
  //~MCPart() = default;
  //MCPart& operator=(const MCPart&) = default;
  int mEvent = 0;
  int mPDGCode = -1;
  std::string mPDGName = "";
  int mCount = 0;
};

class MFTAnaSim
{
 public:
  MFTAnaSim();
  ~MFTAnaSim() = default;

  int mNEvents = 0;             ///< Number of events in the input files
  TTree* mHitTree = nullptr;    ///< Tree with hits, set from the steering macro
  TTree* mKineTree = nullptr;   ///< Tree with particle kinematics, set from the steering macro
  TTree* mClusTree = nullptr;   ///< Tree with clusters, set from the steering macro
  TTree* mTrackTree = nullptr;  ///< Tree with standalone reconstructed tracks, set from the steering macro
  
  o2::mft::GeometryTGeo* mGeoManager;      ///< Manager of the MFT geometry, set from the steering macro
  o2::itsmft::ChipMappingMFT mChipMapper;  ///< Information on the MFT chip mapping

  bool initialize(int maxMCTracks);   ///< Global initialization for all events
  void initEvent(int event, int nMCTracks, int particleSource = 0); ///< Event initialization
  bool doParticles();   ///< Extract information of all particles generated in one event, may select primaries and/or secondaries
  bool doHits();   ///< Count the hits per layer and disk
  bool doMCTracks();   ///< Extract information of the Monte-Carlo tracks in a vector of MFTAnaSimTrack
  bool doSATracks();   ///< Extract information of the standalone reconstructed tracks in a vector of MFTAnaSimSATrack
  void finishEvent();   ///< Called after each event
  void finish();   ///< Called at the end
  void findMCTrackHits(int trkID, int& firstIndex, int& lastIndex);   ///< Return the index range in the vector of all extracted hits (MFTAnaSimHit) for a given MC track
  void findMCTrackClusters(int trkID, int& firstIndex, int& lastIndex);   ///< Return the index range in the vector of all extracted clusters (MFTAnaSimCluster) for a given MC track
  void setVerboseLevel(int vl) { mVerboseLevel = vl; }
  const std::vector<MCPart>& getParticles() { return mParticles; }   ///< Return the vector of (MCPart) extracted particles in the event

  struct ClusterCoord {
    float xGlo = 0.;
    float yGlo = 0.;
    float zGlo = 0.;
    int nPixels = 0;
  };
  
 private:
  bool mPrimary = true, mSecondary = false, mAll = false;
  enum ParticleSource {kPrimary, kSecondary, kAll};
  void filterPDGCode(int& pdgCode);   ///< Skip some PDG codes
  void countParticle(int pdgCode);   ///< Count the particle species and store in a vector of MCPart
  void extractClustersCoord();
  
  std::vector<o2::MCTrack> mMCTrkVec, *mMCTrkVecP = &mMCTrkVec;   ///< Vector of extracted MCTrack
  std::vector<o2::itsmft::Hit> mHitVec, *mHitVecP = &mHitVec;   ///< Vector of hits
  std::vector<o2::itsmft::CompClusterExt> mClusVec, *mClusVecP = &mClusVec;   ///< vector of clusters
  std::vector<unsigned char>* mClusPatternsP = nullptr;   ///< vector of cluster patterns
  std::vector<o2::mft::TrackMFT> mTrackVec, *mTrackVecP = &mTrackVec;   ///< Vector of standalone reconstructed o2::mft::TrackMFT
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* mClusLabels = nullptr;   ///< MC labels for clusters
  std::vector<o2::MCCompLabel> mTrackLabels, *mTrackLabelsInp = &mTrackLabels;   ///< MC labels for the standalone reconstructed tracks
  std::vector<int> mTrackExtClsVec, *mTrackExtClsVecP = &mTrackExtClsVec;   ///< External index for the clusters associated to the standalone reconstructed tracks

  std::vector<std::array<bool, o2::mft::constants::DisksNumber>> mMCTrackHasHitsInDisk;   ///< Counting MFT disks with hits
  std::vector<std::array<bool, o2::mft::constants::LayersNumber>> mMCTrackHasHitsInLayer;   ///< Counting MFT layers with hits
  std::vector<MCPart> mParticles;   ///< Vector counting the generated particles 
  std::vector<MFTAnaSimTrack> mAnaSimTracks;   ///< Vector with extracted MC tracks
  std::vector<MFTAnaSimCluster> mAnaSimClusters;   ///< Vector with clusters associated to MFTAnaSimTrack
  std::vector<MFTAnaSimHit> mAnaSimHits;   ///< Vector with hits associated to MFTAnaSimTrack

  o2::itsmft::TopologyDictionary mTopoDict;   ///< Dictionary for the cluster topologies

  std::vector<ClusterCoord> mClustersCoord;   ///< Simple structure to store cluster coordinates
 
  TFile* mOutFile = nullptr;    ///< Output ROOT file
  
  TTree* mOutTree1 = nullptr;   ///< Tree of MFTAnaSimTrack in the output file
  TTree* mOutTree2 = nullptr;   ///< Tree of MFTAnaSimHit in the output file
  TTree* mOutTree3 = nullptr;   ///< Tree of MFTAnaSimCluster in the output file 

  int mCurrEvent = 0;     ///< Index of the current event
  int mNMCTracks = 0;     ///< Number of MC tracks in event
  int mNSATracks = 0;     ///< Number of standalone reconstructed tracks for all events in the input file
  int mMaxMCTracks = 0;   ///< Maximum number of MC tracks per event
  int mNHitsInEvent = 0;  ///< Number of hits in event
  int mNClusters = 0;     ///< Number of clusters for all events in the input file
  int mVerboseLevel = 0;
};
  
};

#include "../src/MFTAnaSim.cxx"

#endif // MFT_ANA_SIM

