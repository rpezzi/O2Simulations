#ifndef MFT_ANA_SIM_TRACK
#define MFT_ANA_SIM_TRACK

#include "MFTAnaSimMCTrack.h"

#include "MFTBase/Constants.h"

namespace o2::mftana
{

class MFTAnaSimTrack : public MFTAnaSimMCTrack
{
 public:
  MFTAnaSimTrack();
  ~MFTAnaSimTrack() = default;
  MFTAnaSimTrack& operator=(const MFTAnaSimTrack&) = default;

  void setMCTrackID(Int_t id) { mMCTrackID = id; }
  Int_t getMCTrackID() const { return mMCTrackID; }
  
  void setNDisks(Int_t nd) { mNDisks = nd; }
  Int_t getNDisks() const { return mNDisks; }
  
  void setNLayers(Int_t nl) { mNLayers = nl; }
  Int_t getNLayers() const { return mNLayers; }
  
  void setNHits(Int_t nh) { mNHits = nh; }
  Int_t getNHits() const { return mNHits; }
  
  void setFirstHitIndex(Int_t ih) { mFirstHitIndex = ih; }
  Int_t getFirstHitIndex() const { return mFirstHitIndex; }
  
  void setLastHitIndex(Int_t ih) { mLastHitIndex = ih; }
  Int_t getLastHitIndex() const { return mLastHitIndex; }

  void setNClusters(Int_t nh) { mNClusters = nh; }
  Int_t getNClusters() const { return mNClusters; }
  
  void setFirstClusterIndex(Int_t ih) { mFirstClusterIndex = ih; }
  Int_t getFirstClusterIndex() const { return mFirstClusterIndex; }
  
  void setLastClusterIndex(Int_t ih) { mLastClusterIndex = ih; }
  Int_t getLastClusterIndex() const { return mLastClusterIndex; }

  void setLayer(Int_t index, Int_t layer) { mLayers[index] = layer; }
  Int_t getLayer(Int_t index) const { return mLayers[index]; }
  
  void setEvent(Int_t ev) { mEvent = ev; }
  Int_t getEvent() const { return mEvent; }

 private:
  Int_t mEvent = 0;
  Int_t mNDisks = 0;
  Int_t mNLayers = 0;
  Int_t mNHits = 0;
  Int_t mFirstHitIndex = 0;
  Int_t mLastHitIndex = 0;
  Int_t mNClusters = 0;
  Int_t mFirstClusterIndex = 0;
  Int_t mLastClusterIndex = 0;
  Int_t mMCTrackID = 0;
  Int_t mLayers[o2::mft::constants::LayersNumber];
};

};

#endif // MFT_ANA_SIM_TRACK
