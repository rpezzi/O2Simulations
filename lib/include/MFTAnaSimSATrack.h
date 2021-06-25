#ifndef MFT_ANA_SA_TRACK
#define MFT_ANA_SA_TRACK

namespace o2::mftana
{

class MFTAnaSimSATrack 
{
 public:
  MFTAnaSimSATrack() = default;
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
  
  void setEventID(int ipoint, int evnID) {
    assert(ipoint < o2::mft::constants::LayersNumber);
    mEventID[ipoint] = evnID;
  }
  int getEventID(int ipoint) const {
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

 private: 
  int mNDisks = 0;
  int mNLayers = 0;
  int mNPoints = 0;
  int mLayers[o2::mft::constants::LayersNumber];
  int mIntClusIndex[o2::mft::constants::LayersNumber];
  int mEventID[o2::mft::constants::LayersNumber];
  int mMCTrackID[o2::mft::constants::LayersNumber];
};

};

#endif // MFT_ANA_SA_TRACK
