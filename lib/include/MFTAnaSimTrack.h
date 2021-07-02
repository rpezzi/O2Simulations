#ifndef MFT_ANA_SIM_TRACK
#define MFT_ANA_SIM_TRACK

#include "MFTAnaSimMCTrack.h"

#include "MFTBase/Constants.h"

namespace o2::mftana
{

constexpr int MCSplitCluster = 4;   ///< Maximum split of a cluster with the same MC track ID

class MFTAnaSimTrack : public MFTAnaSimMCTrack
{
 public:
  MFTAnaSimTrack();
  ~MFTAnaSimTrack() = default;
  MFTAnaSimTrack& operator=(const MFTAnaSimTrack&) = default;

  void setMCTrackID(int id) { mMCTrackID = id; }
  int getMCTrackID() const { return mMCTrackID; }
  
  void setNDisks(int nd) { mNDisks = nd; }
  int getNDisks() const { return mNDisks; }
  
  void setNLayers(int nl) { mNLayers = nl; }
  int getNLayers() const { return mNLayers; }
  
  void setNHits(int nh) { mNHits = nh; }
  int getNHits() const { return mNHits; }
  
  void setFirstHitIndex(int ih) { mFirstHitIndex = ih; }
  int getFirstHitIndex() const { return mFirstHitIndex; }
  
  void setLastHitIndex(int ih) { mLastHitIndex = ih; }
  int getLastHitIndex() const { return mLastHitIndex; }

  void setNClusters(int nh) { mNClusters = nh; }
  int getNClusters() const { return mNClusters; }
  
  void setEvent(int ev) { mEvent = ev; }
  int getEvent() const { return mEvent; }

  void setLayer(int icluster, int layer) {
    assert(icluster < o2::mft::constants::LayersNumber);
    mLayers[icluster] = layer;
  }
  int getLayer(int icluster) const {
    assert(icluster < o2::mft::constants::LayersNumber);
    return mLayers[icluster];
  }
  
  void setIntClusIndex(int icluster, int index) {
    assert(icluster < (MCSplitCluster * o2::mft::constants::LayersNumber));
    mIntClusIndex[icluster] = index;
  }
  int getIntClusIndex(int icluster) const {
    assert(icluster < (MCSplitCluster * o2::mft::constants::LayersNumber));
    return mIntClusIndex[icluster];
  }

  int getNSATracks() const { return mNSATracks; }
  
  void addIntSATrackIndex(int index);
  
  int getIntSATrackIndex(int i) const {
    assert(i < (MCSplitCluster * o2::mft::constants::LayersNumber));
    return mIntSATrackIndex[i];
  }
  
  int getIntSATrackMult(int i) const {
    assert(i < (MCSplitCluster * o2::mft::constants::LayersNumber));
    return mIntSATrackMult[i];
  }
  
 private:
  int mEvent = 0;   ///< Event to which this MC track belongs
  int mNDisks = 0;   ///< Number of MFT disks
  int mNLayers = 0;   ///< Number of MFT layers
  int mNHits = 0;   ///< Number of hits
  int mFirstHitIndex = 0;   ///< Hit index range, start
  int mLastHitIndex = 0;   ///< Hit index range, end
  int mNClusters = 0;   ///< Number of clusters
  int mMCTrackID = 0;   ///< ID of this MC track
  int mLayers[o2::mft::constants::LayersNumber];   ///< ID of the layers
  int mIntClusIndex[MCSplitCluster * o2::mft::constants::LayersNumber];   ///< Internal index for the attached clusters
  int mIntSATrackIndex[MCSplitCluster * o2::mft::constants::LayersNumber];   ///< List of SA track indexes to which this MC track contributes with clusters
  int mIntSATrackMult[MCSplitCluster * o2::mft::constants::LayersNumber];   ///< Multiplicity of SA track indexes to which this MC track contributes with clusters
  int mNSATracks = 0;   ///< Number of SA tracks to which this MC track contributes with clusters
};

//_____________________________________________________________________________
inline void MFTAnaSimTrack::addIntSATrackIndex(int index)
{
  int i;
  for (i = 0; i < mNSATracks; i++) {
    assert(i < (MCSplitCluster * o2::mft::constants::LayersNumber));
    if (mIntSATrackIndex[i] == index) {
      mIntSATrackMult[i]++;
      return;
    }
  }
  assert(i < (MCSplitCluster * o2::mft::constants::LayersNumber));
  mIntSATrackIndex[i] = index;
  mIntSATrackMult[i]++;
  mNSATracks++;
}

};

#endif // MFT_ANA_SIM_TRACK
