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
  
 private:
  Double_t mX;
  Double_t mY;
  Double_t mZ;
};

};

#endif // MFT_ANA_SIM_HIT
