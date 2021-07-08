#include "../include/MFTAnaSimTrack.h"

namespace o2::mftana
{

//_____________________________________________________________________________
MFTAnaSimTrack::MFTAnaSimTrack()
{
  for (int i = 0; i < o2::mft::constants::LayersNumber; i++) {
    mLayers[i] = -1;
  }
  for (int i = 0; i < (MCSplitCluster * o2::mft::constants::LayersNumber); i++) {
    mIntSATrackIndex[i] = -1;
    mIntSATrackMult[i] = 0;
  }
  mNSATracks = 0;
}
  
}; // end namespace o2::mftana
