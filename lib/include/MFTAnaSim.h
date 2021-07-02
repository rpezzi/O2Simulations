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
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"

namespace o2::mftana
{

using o2::itsmft::Hit;
using o2::MCTrack;

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

  int mNThreads;   ///< Number of threads to use when compiled with OpenMP
  int mNEvents = 0;             ///< Number of events in the input files
  TTree* mHitTree = nullptr;    ///< Tree with hits, set from the steering macro
  TTree* mKineTree = nullptr;   ///< Tree with particle kinematics, set from the steering macro
  TTree* mClusTree = nullptr;   ///< Tree with clusters, set from the steering macro
  TTree* mTrackTree = nullptr;  ///< Tree with standalone reconstructed tracks, set from the steering macro
  
  o2::mft::GeometryTGeo* mGeoManager;      ///< Manager of the MFT geometry, set from the steering macro
  o2::itsmft::ChipMappingMFT mChipMapper;  ///< Information on the MFT chip mapping

  bool initialize();   ///< Global initialization for all events
  void initEvent(int event, int particleSource = 2); ///< Event initialization; the source selection is only for the vector of particles
  bool trackHasHits(int trkID);   ///< Check if a track (by its MC ID) has hits in the detector
  bool doParticles();   ///< Extract information of all particles generated in one event, which left hits in the MFT; may select primaries and/or secondaries
  bool doHits();   ///< Mark MC tracks with hits per layer and per disk
  bool doMCTracks();   ///< Extract information of the Monte-Carlo tracks in a vector of MFTAnaSimTrack
  bool doSATracks();   ///< Extract information of the standalone reconstructed tracks in a vector of MFTAnaSimSATrack
  void finishEvent();   ///< Called after each event
  void finish();   ///< Called at the end
  void findMCTrackHits(int trkID, int& firstIndex, int& lastIndex);   ///< Return the index range in the vector of all extracted hits (MFTAnaSimHit) for a given MC track
  void findMCTrackClusters(MFTAnaSimTrack& asTrack);   ///< Associate indexes in the vector of MFTAnsSimCluster to this MC track
  void setVerboseLevel(int vl) { mVerboseLevel = vl; }   ///< Set the verbosity level of the printed messages
  const std::vector<MCPart>& getParticles() { return mParticles; }   ///< Return the vector of (MCPart) extracted particles from all events
  const std::vector<MFTAnaSimTrack>& getSimTracks() { return mAnaSimTracks; } ///< Return the vector of (MFTAnaSimTrack) extracted MC tracks from all events
  void linkTracks();   ///< Link SA tracks (indexes) to MC tracks (indexes) and viceversa
  
 private:
  bool mPrimary = true, mSecondary = false, mAll = false;
  enum ParticleSource {kPrimary, kSecondary, kAll};
  bool filterPDGCode(int pdgCode);   ///< Skip some PDG codes
  void countParticle(int pdgCode, int currEvent);   ///< Count the particle species and store in a vector of MCPart
  void extractClusters();   ///< Extract the clusters in a vector of MFTAnsSimCluster
  
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
  std::vector<MFTAnaSimCluster> mAnaSimClustersSort;   ///< Vector with clusters associated to MFTAnaSimTrack, sorted in event ID order
  std::vector<MFTAnaSimHit> mAnaSimHits;   ///< Vector with hits associated to MFTAnaSimTrack
  std::vector<MFTAnaSimSATrack> mAnaSimSATracks;   ///< Vector with standalone reconstructed tracks

  o2::itsmft::TopologyDictionary mTopoDict;   ///< Dictionary for the cluster topologies

  TFile* mOutFile = nullptr;    ///< Output ROOT file
  
  TTree* mOutTree1 = nullptr;   ///< Tree of MFTAnaSimTrack in the output file
  TTree* mOutTree2 = nullptr;   ///< Tree of MFTAnaSimHit in the output file
  TTree* mOutTree3 = nullptr;   ///< Tree of MFTAnaSimCluster in the output file 
  TTree* mOutTree4 = nullptr;   ///< Tree of MFTAnaSimSATrack in the output file

  int mCurrEvent = 0;     ///< Index of the current event
  int mNMCTracks = 0;     ///< Number of MC tracks in event
  int mNMCTracksWHits = 0;   ///< Number of MC tracks in event with hits
  int mNSATracks = 0;     ///< Number of standalone reconstructed tracks from all events in the input file
  int mMaxMCTracks = 0;   ///< Maximum number of MC tracks per event
  int mNHitsInEvent = 0;  ///< Number of hits in event
  int mNClusters = 0;     ///< Number of clusters from all events in the input file
  std::vector<std::pair<int, int>> mEventClusterRange;   ///< Split the cluster vector indexes by event
  int mEvnClsIndexMin = 0;   ///< Split the cluster vector: start index
  int mEvnClsIndexMax = 0;   ///< Split the cluster vector: end index
  int mVerboseLevel = 0;   ///< Verbosity level of the printed messages
};

//_____________________________________________________________________________
inline bool MFTAnaSim::trackHasHits(int trkID)
{
  for (int i_hit = 0; i_hit < mNHitsInEvent; i_hit++) {
    Hit* hitp = &(mHitVec).at(i_hit);
    if (hitp->GetTrackID() == trkID) {
      return true;
    }
  }
  return false;
}
 
//_____________________________________________________________________________
inline bool MFTAnaSim::filterPDGCode(int pdgCode)
{
  // skip reggeon and pomeron
  if (abs(pdgCode) == 110 || abs(pdgCode) == 990) {
    return false;
  }
  // diffractive states
  if (pdgCode > 9900000) {
    return false;
  }
  // nuclear code
  if (abs(pdgCode) > 1000000000) {
    return false;
  }
  return true;
}
  
};

#endif // MFT_ANA_SIM

