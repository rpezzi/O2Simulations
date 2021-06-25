#ifndef MFT_ANA_MC_TRACK
#define MFT_ANA_MC_TRACK

#include <TMath.h>

#include "TDatabasePDG.h"

namespace o2::mftana
{

class MFTAnaSimMCTrack 
{
 public:
  MFTAnaSimMCTrack() = default;
  ~MFTAnaSimMCTrack() = default;
  MFTAnaSimMCTrack& operator=(const MFTAnaSimMCTrack&) = default;
  
  double getP() const
  {
    double mx(mStartVertexMomentumX);
    double my(mStartVertexMomentumY);
    double mz(mStartVertexMomentumZ);
    return std::sqrt(mx * mx + my * my + mz * mz);
  }

  double getPt() const
  {
    double mx(mStartVertexMomentumX);
    double my(mStartVertexMomentumY);
    return std::sqrt(mx * mx + my * my);
  }

  double getPhi() const
  {
    double mx(mStartVertexMomentumX);
    double my(mStartVertexMomentumY);
    return (TMath::Pi() + TMath::ATan2(-mx, -my));
  }

  double getEta() const
  {
    double pmom = getP();
    double mz(mStartVertexMomentumZ);
    if (pmom != TMath::Abs(mz)) {
      return 0.5 * std::log((pmom + mz) / (pmom - mz));
    } else {
      return 1.e30;
    }
  }

  double getTheta() const
  {
    double mz(mStartVertexMomentumZ);
    return (mz == 0) ? TMath::PiOver2() : TMath::ACos(mz / getP());
  }

  void setIsPrimary(bool isPrimary) { mIsPrimary = isPrimary; }
  bool isPrimary() const { return mIsPrimary; }

  void setVertexPxPyPz(double px, double py, double pz)
  {
    mStartVertexMomentumX = px;
    mStartVertexMomentumY = py;
    mStartVertexMomentumZ = pz;
  }
  void setVertexXYZ(double x, double y, double z)
  {
    mStartVertexCoordinatesX = x;
    mStartVertexCoordinatesY = y;
    mStartVertexCoordinatesZ = z;
  }

  double getVertexPx() const { return mStartVertexMomentumX; }
  double getVertexPy() const { return mStartVertexMomentumY; }
  double getVertexPz() const { return mStartVertexMomentumZ; }

  double getVertexX() const { return mStartVertexCoordinatesX; }
  double getVertexY() const { return mStartVertexCoordinatesY; }
  double getVertexZ() const { return mStartVertexCoordinatesZ; }

  void setPDGCode(int pdgCode)
  {
    mPDGCode = pdgCode;
    mPDGName = std::string(TDatabasePDG::Instance()->GetParticle(pdgCode)->GetName());
  }
  
  void setMotherTrackId(int trkid) { mMotherTrackId = trkid; }
  int getMotherTrackId() const { return mMotherTrackId; }

  void setSecondMotherTrackId(int trkid) { mSecondMotherTrackId = trkid; }
  int getSecondMotherTrackId() const { return mSecondMotherTrackId; }

  void setFirstDaughterTrackId(int trkid) { mFirstDaughterTrackId = trkid; }
  int getFirstDaughterTrackId() const { return mFirstDaughterTrackId; }

  void setLastDaughterTrackId(int trkid) { mLastDaughterTrackId = trkid; }
  int getLastDaughterTrackId() const { return mLastDaughterTrackId; }

 private:
  int mPDGCode = -1;
  std::string mPDGName = "";
  bool mIsPrimary = true;
  double mStartVertexCoordinatesX = 0.;
  double mStartVertexCoordinatesY = 0.;
  double mStartVertexCoordinatesZ = 0.;
  double mStartVertexMomentumX = 0.;
  double mStartVertexMomentumY = 0.;
  double mStartVertexMomentumZ = 0.;
  int mMotherTrackId = -1;
  int mSecondMotherTrackId = -1;
  int mFirstDaughterTrackId = -1;
  int mLastDaughterTrackId = -1;
};

};

#endif // MFT_ANA_MC_TRACK
