#include "../include/MFTAnaSimSATrack.h"

namespace o2::mftana
{

//_____________________________________________________________________________
MFTAnaSimSATrack::MFTAnaSimSATrack()
{
  for (int i = 0; i < o2::mft::constants::LayersNumber; i++) {
    mLayers[i] = -1;
  }
  for (int i = 0; i < (SASplitCluster * o2::mft::constants::LayersNumber); i++) {
    mIntMCTrackIndex[i] = -1;
    mIntMCTrackMult[i] = 0;
  }
  mNMCTracks = 0;
}
  
//_____________________________________________________________________________
void MFTAnaSimSATrack::copy(o2::mft::TrackMFT& track)
{
  this->setOutParam(track.getOutParam());
  this->setX(track.getX());
  this->setY(track.getY());
  this->setZ(track.getZ());
  this->setCharge(track.getCharge());
  this->setTrackChi2(track.getTrackChi2());
}
  
}; // end namespace o2::mftana
