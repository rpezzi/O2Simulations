#ifndef MFT_ANA_SA_TRACK
#define MFT_ANA_SA_TRACK

#include "MFTBase/Constants.h"
#include "DataFormatsMFT/TrackMFT.h"

namespace o2::mftana
{

constexpr int SASplitCluster = 4;   ///< Maximum split of a cluster with the same MC track ID

class MFTAnaSimSATrack : public o2::mft::TrackMFT
{
 public:
  MFTAnaSimSATrack();
  ~MFTAnaSimSATrack() = default;
  MFTAnaSimSATrack& operator=(const MFTAnaSimSATrack&) = default;

  void setNDisks(int nd) { mNDisks = nd; }
  int getNDisks() const { return mNDisks; }
  
  void setNLayers(int nl) { mNLayers = nl; }
  int getNLayers() const { return mNLayers; }
  
  void setNPoints(int np) { mNPoints = np; }
  int getNPoints() const { return mNPoints; }
  
  void setLayer(int ipoint, int layer) {
    assert(ipoint < o2::mft::constants::LayersNumber);
    mLayers[ipoint] = layer;
  }
  int getLayer(int ipoint) const { return mLayers[ipoint]; }

  void setIntClusIndex(int ipoint, int index) {
    assert(ipoint < o2::mft::constants::LayersNumber);
    mIntClusIndex[ipoint] = index;
  }
  int getIntClusIndex(int ipoint) const {
    assert(ipoint < o2::mft::constants::LayersNumber);
    return mIntClusIndex[ipoint];
  }
  
  void setEvent(int ipoint, int evnID) {
    assert(ipoint < o2::mft::constants::LayersNumber);
    mEventID[ipoint] = evnID;
  }
  int getEvent(int ipoint) const {
    assert(ipoint < o2::mft::constants::LayersNumber);
    return mEventID[ipoint];
  }

  void setMCTrackID(int ipoint, int trkID) {
    assert(ipoint < o2::mft::constants::LayersNumber);
    mMCTrackID[ipoint] = trkID;
  }
  int getMCTrackID(int ipoint) const {
    assert(ipoint < o2::mft::constants::LayersNumber);
    return mMCTrackID[ipoint];
  }

  void copy(o2::mft::TrackMFT& track);
  
  int getNMCTracks() const { return mNMCTracks; }
  
  void addIntMCTrackIndex(int evn, int id);
  
  int getIntMCTrackEvent(int i) const {
    assert(i < (SASplitCluster * o2::mft::constants::LayersNumber));
    return mIntMCTrackEvent[i];
  }
  
  int getIntMCTrackIndex(int i) const {
    assert(i < (SASplitCluster * o2::mft::constants::LayersNumber));
    return mIntMCTrackIndex[i];
  }
  
  int getIntMCTrackMult(int i) const {
    assert(i < (SASplitCluster * o2::mft::constants::LayersNumber));
    return mIntMCTrackMult[i];
  }
  
 private: 
  int mNDisks = 0;   ///< Number of MFT disks
  int mNLayers = 0;   ///< Number of MFT layers
  int mNPoints = 0;   ///< Number of points in the track
  int mLayers[o2::mft::constants::LayersNumber];   ///< ID of the layers
  int mIntClusIndex[o2::mft::constants::LayersNumber];   ///< Internal index for the attached clusters
  int mEventID[o2::mft::constants::LayersNumber];   ///< ID of the events to which the clusters belong
  int mMCTrackID[o2::mft::constants::LayersNumber];   ///< ID of the MC tracks which contribute to the points
  int mIntMCTrackEvent[SASplitCluster * o2::mft::constants::LayersNumber];   ///< List of the event IDs which contribute with MC tracks
  int mIntMCTrackIndex[SASplitCluster * o2::mft::constants::LayersNumber];   ///< List of MC track indexes which contribute with clusters to this SA track
  int mIntMCTrackMult[SASplitCluster * o2::mft::constants::LayersNumber];   ///< Multiplicity of MC track indexes which contribute with clusters to this SA track
  int mNMCTracks = 0;   ///< Number of MC tracks which contribute with clusters to this SA track
};

//_____________________________________________________________________________
inline void MFTAnaSimSATrack::addIntMCTrackIndex(int evn, int id)
{
  int i;
  for (i = 0; i < mNMCTracks; i++) {
    assert(i < (SASplitCluster * o2::mft::constants::LayersNumber));
    if (mIntMCTrackEvent[i] == evn && mIntMCTrackIndex[i] == id ) {
      mIntMCTrackMult[i]++;
      return;
    }
  }
  assert(i < (SASplitCluster * o2::mft::constants::LayersNumber));
  mIntMCTrackEvent[i] = evn;
  mIntMCTrackIndex[i] = id;
  mIntMCTrackMult[i]++;
  mNMCTracks++;
}

};

#endif // MFT_ANA_SA_TRACK
