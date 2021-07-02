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

  void setID(int id) { mID = id; }
  int getID() const { return mID; }

  void print() const;
  
 private:
  bool mIsNoise = true;   ///< True if the cluster is from noise digits
  int mEvent = -1;   ///< Event ID to which belongs the MC track which contributed to this cluster
  int mMCTrackID = -1;   ///< ID of the MC track which contributed to this clusters
  int mNPixels = 0;   ///< Number of pixels in the cluster
  int mLayer = -1;   ///< Layer ID
  int mID = -1;   ///< Internal ID of the cluster
};

//_____________________________________________________________________________
inline void MFTAnaSimCluster::print() const
{
  printf("Cluster in event %d MCtrackID %d Noise %d x,y,z  %7.3f  %7.3f  %7.3f \n", mEvent, mMCTrackID, mIsNoise, getX(), getY(), getZ());
}

};

#endif // MFT_ANA_SIM_CLUSTER
