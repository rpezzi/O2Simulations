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

 private:
};

};

#endif // MFT_ANA_SIM_CLUSTER
