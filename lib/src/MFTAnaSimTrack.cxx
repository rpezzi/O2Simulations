#include "../include/MFTAnaSimTrack.h"

namespace o2::mftana
{

//_____________________________________________________________________________
MFTAnaSimTrack::MFTAnaSimTrack()
{
  for (int i = 0; i < o2::mft::constants::LayersNumber; i++) {
    mLayers[i] = 0;
  }
}
  
}; // end namespace o2::mftana
