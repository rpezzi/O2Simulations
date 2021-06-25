#ifndef MFT_ANA_SIM_TRACK
#define MFT_ANA_SIM_TRACK

#include "MFTAnaSimMCTrack.h"

#include "MFTBase/Constants.h"

namespace o2::mftana
{

constexpr int SplitCluster = 4;

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
  
  void setIntClusIndex(int icluster, int iclus) {
    assert(icluster < (SplitCluster * o2::mft::constants::LayersNumber));
    mIntClusIndex[icluster] = iclus;
  }
  int getIntClusIndex(int icluster) const {
    assert(icluster < (SplitCluster * o2::mft::constants::LayersNumber));
    return mIntClusIndex[icluster];
  }
  
 private:
  int mEvent = 0;
  int mNDisks = 0;
  int mNLayers = 0;
  int mNHits = 0;
  int mFirstHitIndex = 0;
  int mLastHitIndex = 0;
  int mNClusters = 0;
  int mFirstClusterIndex = 0;
  int mLastClusterIndex = 0;
  int mMCTrackID = 0;
  int mLayers[o2::mft::constants::LayersNumber];
  int mIntClusIndex[SplitCluster * o2::mft::constants::LayersNumber];
};

};

#endif // MFT_ANA_SIM_TRACK
