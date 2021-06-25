#ifndef MFT_ANA_SIM_CLUSTER
#define MFT_ANA_SIM_CLUSTER

#include "MFTTracking/Cluster.h"

namespace o2::mftana
{

class MFTAnaSimCluster : public o2::mft::Cluster
{
 public:
  MFTAnaSimCluster()= default;
  ~MFTAnaSimCluster() = default;
  MFTAnaSimCluster& operator=(const MFTAnaSimCluster&) = default;

  void setEvent(int currEvent) { mEvent = currEvent; }
  int getEvent() const { return mEvent; }

  void setIsNoise(bool val) { mIsNoise = val; }
  bool getIsNoise() const { return mIsNoise; }

  void setNPixels(int npix) { mNPixels = npix; }
  int getNPixels() const { return mNPixels; }

  void setMCTrackID(int trkID) { mMCTrackID = trkID; }
  int getMCTrackID() const { return mMCTrackID; }

  void setLayer(int layer) { mLayer = layer; }
  int getLayer() const { return mLayer; }

  void print() const;
  
 private:
  bool mIsNoise = true;
  int mEvent = -1;
  int mMCTrackID = -1;
  int mNPixels = 0;
  int mLayer = -1;
};

//_____________________________________________________________________________
inline void MFTAnaSimCluster::print() const
{
  printf("Cluster in event %d MCtrackID %d Noise %d x,y,z  %7.3f  %7.3f  %7.3f \n", mEvent, mMCTrackID, mIsNoise, getX(), getY(), getZ());
}

};

#endif // MFT_ANA_SIM_CLUSTER
