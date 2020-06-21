#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "DataFormatsMFT/TrackMFT.h"
#include <TMath.h>


#endif





//_________________________________________________________________________________________________
void extrapMFTTrackHelixToZ(o2::mft::TrackMFT& track, double zEnd, double Field)
{
   using TrackMFT = o2::mft::TrackMFT;

  /// Track extrapolated to the plane at "zEnd" considering a helix

  if (track.getZ() == zEnd) {
    return; // nothing to be done if same z
  }

  // Compute track parameters
  double dZ = (zEnd - track.getZ());
  double x0 = track.getX();
  double y0 = track.getY();
  double px0 = track.getPx();
  double py0 = track.getPy();
  double invtanl0 = 1.0 / track.getTanl();;
  double invqpt0 = track.getInvQPt();
  auto q = track.getCharge();
  auto Hz =  std::copysign(1, Field); 
  double k = TMath::Abs(3e-4 * Field);
  auto invk = 1.0 / k;
  double theta = -invqpt0 * dZ * k * invtanl0;
  double costheta, sintheta;
  o2::utils::sincos(theta, sintheta, costheta);
  double deltax = Hz * py0 * invk * (1.0 - costheta) - px0 * q * invk * sintheta;
  double deltay = -Hz * px0 * invk * (1.0 - costheta) - py0 * q * invk * sintheta;

  double x = x0 + deltax;
  double y = y0 + deltay;
  double phi = track.getPhi() + Hz * theta;

  track.setX(x);
  track.setY(y);
  track.setZ(zEnd);
  track.setPhi(phi);
}


//_________________________________________________________________________________________________
void extrapMFTTrackQuadraticToZ(o2::mft::TrackMFT& track, double zEnd, double Field)
{
   using TrackMFT = o2::mft::TrackMFT;

  /// Track extrapolated to the plane at "zEnd" considering a helix

  if (track.getZ() == zEnd) {
    return; // nothing to be done if same z
  }

  // Compute track parameters
  double dZ = (zEnd - track.getZ());
  double x0 = track.getX();
  double y0 = track.getY();
  double phi0 = track.getPhi();
  double cosphi0, sinphi0;
  o2::utils::sincos(phi0, sinphi0, cosphi0);
  double invtanl0 = 1.0 / track.getTanl();;
  double invqpt0 = track.getInvQPt();
  auto q = track.getCharge();
  auto Hz =  std::copysign(1, Field); 

  double n = dZ * invtanl0;
  double k = TMath::Abs(o2::constants::math::B2C * Field);
  double theta = -invqpt0 * dZ * k * invtanl0;
  double deltax = n * cosphi0 - 0.5 * n * theta * Hz * sinphi0;
  double deltay = n * sinphi0 + 0.5 * n * theta * Hz * cosphi0;

  double x = x0 + deltax;
  double y = y0 + deltay;
  double phi = phi0 + Hz * theta;

  track.setX(x);
  track.setY(y);
  track.setZ(zEnd);
  track.setPhi(phi);
  
}


//_________________________________________________________________________________________________
void extrapMFTTrackLinearToZ(o2::mft::TrackMFT& track, double zEnd)
{
   using TrackMFT = o2::mft::TrackMFT;

  /// Track extrapolated to the plane at "zEnd" considering a helix

  if (track.getZ() == zEnd) {
    return; // nothing to be done if same z
  }

  // Compute track parameters
  double dZ = (zEnd - track.getZ());
  double x0 = track.getX();
  double y0 = track.getY();
  double phi0 = track.getPhi();
  double invtanl0 = 1.0 / track.getTanl();;
  double cosphi0, sinphi0;
  o2::utils::sincos(phi0, sinphi0, cosphi0);
  double n = dZ * invtanl0;
  double x = x0 + n * cosphi0;
  double y = y0 + n * sinphi0;

  track.setX(x);
  track.setY(y);
  track.setZ(zEnd);
}
