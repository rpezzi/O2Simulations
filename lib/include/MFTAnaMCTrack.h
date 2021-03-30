#ifndef MFT_ANA_MC_TRACK
#define MFT_ANA_MC_TRACK

#include <TMath.h>

#include "TDatabasePDG.h"

namespace o2::mftana
{

class MFTAnaMCTrack 
{
 public:
  MFTAnaMCTrack() = default;
  ~MFTAnaMCTrack() = default;
  MFTAnaMCTrack& operator=(const MFTAnaMCTrack&) = default;
  
  Double_t getP() const
  {
    double mx(mStartVertexMomentumX);
    double my(mStartVertexMomentumY);
    double mz(mStartVertexMomentumZ);
    return std::sqrt(mx * mx + my * my + mz * mz);
  }

  Double_t getPt() const
  {
    double mx(mStartVertexMomentumX);
    double my(mStartVertexMomentumY);
    return std::sqrt(mx * mx + my * my);
  }

  Double_t getPhi() const
  {
    double mx(mStartVertexMomentumX);
    double my(mStartVertexMomentumY);
    return (TMath::Pi() + TMath::ATan2(-mx, -my));
  }

  Double_t getEta() const
  {
    double_t pmom = getP();
    double mz(mStartVertexMomentumZ);
    if (pmom != TMath::Abs(mz)) {
      return 0.5 * std::log((pmom + mz) / (pmom - mz));
    } else {
      return 1.e30;
    }
  }

  Double_t getTheta() const
  {
    double mz(mStartVertexMomentumZ);
    return (mz == 0) ? TMath::PiOver2() : TMath::ACos(mz / getP());
  }

  void setIsPrimary(Bool_t isPrimary) { mIsPrimary = isPrimary; }
  Bool_t isPrimary() const { return mIsPrimary; }

  void setVertexPxPyPz(Double_t px, Double_t py, Double_t pz)
  {
    mStartVertexMomentumX = px;
    mStartVertexMomentumY = py;
    mStartVertexMomentumZ = pz;
  }
  void setVertexXYZ(Double_t x, Double_t y, Double_t z)
  {
    mStartVertexCoordinatesX = x;
    mStartVertexCoordinatesY = y;
    mStartVertexCoordinatesZ = z;
  }

  Double_t getVertexPx() const { return mStartVertexMomentumX; }
  Double_t getVertexPy() const { return mStartVertexMomentumY; }
  Double_t getVertexPz() const { return mStartVertexMomentumZ; }

  Double_t getVertexX() const { return mStartVertexCoordinatesX; }
  Double_t getVertexY() const { return mStartVertexCoordinatesY; }
  Double_t getVertexZ() const { return mStartVertexCoordinatesZ; }

  void setPDGCode(Int_t pdgCode)
  {
    mPDGCode = pdgCode;
    mPDGName = std::string(TDatabasePDG::Instance()->GetParticle(pdgCode)->GetName());
  }
  
  void setMotherTrackId(Int_t trkid) { mMotherTrackId = trkid; }
  Int_t getMotherTrackId() const { return mMotherTrackId; }
  void setSecondMotherTrackId(Int_t trkid) { mSecondMotherTrackId = trkid; }
  Int_t getSecondMotherTrackId() const { return mSecondMotherTrackId; }
  void setFirstDaughterTrackId(Int_t trkid) { mFirstDaughterTrackId = trkid; }
  Int_t getFirstDaughterTrackId() const { return mFirstDaughterTrackId; }
  void setLastDaughterTrackId(Int_t trkid) { mLastDaughterTrackId = trkid; }
  Int_t getLastDaughterTrackId() const { return mLastDaughterTrackId; }

 private: 
  Int_t mPDGCode = -1;
  std::string mPDGName = "";
  Bool_t mIsPrimary = kTRUE;
  Double_t mStartVertexCoordinatesX = 0.;
  Double_t mStartVertexCoordinatesY = 0.;
  Double_t mStartVertexCoordinatesZ = 0.;
  Double_t mStartVertexMomentumX = 0.;
  Double_t mStartVertexMomentumY = 0.;
  Double_t mStartVertexMomentumZ = 0.;
  Int_t mMotherTrackId = -1;
  Int_t mSecondMotherTrackId = -1;
  Int_t mFirstDaughterTrackId = -1;
  Int_t mLastDaughterTrackId = -1;
};

};

#endif // MFT_ANA_MC_TRACK
