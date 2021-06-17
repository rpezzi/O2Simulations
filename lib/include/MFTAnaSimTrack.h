#ifndef MFT_ANA_SIM_TRACK
#define MFT_ANA_SIM_TRACK

#include "MFTAnaSimMCTrack.h"

#include "MFTBase/Constants.h"

namespace o2::mftana
{

class MFTAnaSimTrack : public MFTAnaSimMCTrack
{
 public:
   MFTAnaSimTrack() {
     for (int i = 0; i < o2::mft::constants::LayersNumber; i++) {
       mLayers[i] = 0;
     }
   }
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
  
  void setFirstClusterIndex(int ih) { mFirstClusterIndex = ih; }
  int getFirstClusterIndex() const { return mFirstClusterIndex; }
  
  void setLastClusterIndex(int ih) { mLastClusterIndex = ih; }
  int getLastClusterIndex() const { return mLastClusterIndex; }

  void setLayer(int index, int layer) { mLayers[index] = layer; }
  int getLayer(int index) const { return mLayers[index]; }
  
  void setEvent(int ev) { mEvent = ev; }
  int getEvent() const { return mEvent; }

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
};

};

//#include "../src/MFTAnaSimTrack.cxx"

#endif // MFT_ANA_SIM_TRACK
