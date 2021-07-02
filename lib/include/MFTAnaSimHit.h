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

  void setXYZ(double x, double y, double z)
  {
    mX = x;
    mY = y;
    mZ = z;
  }
  double getX() const { return mX; }
  double getY() const { return mY; }
  double getZ() const { return mZ; }
  
  void setEvent(int ev) { mEvent = ev; }
  int getEvent() const { return mEvent; }

 private:
  int mEvent = 0;   ///< ID of the event in which this hit was produced
  double mX = 0.;   ///< Hit x global coordinate
  double mY = 0.;   ///< Hit y global coordinate
  double mZ = 0.;   ///< Hit z global coordinate
};

};

#endif // MFT_ANA_SIM_HIT
