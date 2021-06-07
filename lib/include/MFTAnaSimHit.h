#ifndef MFT_ANA_SIM_HIT
#define MFT_ANA_SIM_HIT

namespace o2::mftana
{

class MFTAnaSimHit
{
 public:
  MFTAnaSimHit() = default;
  ~MFTAnaSimHit() = default;
  MFTAnaSimHit& operator=(const MFTAnaSimHit&) = default;

  void setXYZ(Double_t x, Double_t y, Double_t z)
  {
    mX = x;
    mY = y;
    mZ = z;
  }
  Double_t getX() const { return mX; }
  Double_t getY() const { return mY; }
  Double_t getZ() const { return mZ; }
  
  void setEvent(Int_t ev) { mEvent = ev; }
  Int_t getEvent() const { return mEvent; }

 private:
  Int_t mEvent = 0;
  Double_t mX = 0.;
  Double_t mY = 0.;
  Double_t mZ = 0.;
};

};

#endif // MFT_ANA_SIM_HIT
