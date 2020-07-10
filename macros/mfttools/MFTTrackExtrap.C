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
  auto dZ = (zEnd - track.getZ());
  auto x0 = track.getX();
  auto y0 = track.getY();
  auto phi0 = track.getPhi(); 
  auto invtanl0 = 1.0 / track.getTanl();;
  auto invqpt0 = track.getInvQPt();
  auto qpt0 = 1.0 / invqpt0;
  double cosphi0, sinphi0;
  o2::utils::sincos(phi0, sinphi0, cosphi0);
  auto q = track.getCharge();
  auto Hz =  std::copysign(1, Field); 
  auto k = TMath::Abs(3e-4 * Field);
  auto invk = 1.0 / k;
  auto theta = -invqpt0 * dZ * k * invtanl0;
  double costheta, sintheta;
  o2::utils::sincos(theta, sintheta, costheta);




  auto Y = sinphi0 * qpt0 * invk;
  auto X = cosphi0 * qpt0 * invk;
  auto YC = Y * costheta;
  auto YS = Y * sintheta;
  auto XC = X * costheta;
  auto XS = X * sintheta;

  // Extrapolate track parameters to "zEnd"
  auto x = x0 + Hz * (Y - YC) - XS;
  auto y = y0 + Hz * (-X + XC) - YS;
  auto phi = phi0 + Hz * theta;

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
  auto dZ = (zEnd - track.getZ());
  auto x0 = track.getX();
  auto y0 = track.getY();
  auto phi0 = track.getPhi();
  double cosphi0, sinphi0;
  o2::utils::sincos(phi0, sinphi0, cosphi0);
  auto invtanl0 = 1.0 / track.getTanl();;
  auto invqpt0 = track.getInvQPt();
  auto q = track.getCharge();
  auto Hz =  std::copysign(1, Field); 

  auto n = dZ * invtanl0;
  auto k = TMath::Abs(o2::constants::math::B2C * Field);
  auto theta = -invqpt0 * dZ * k * invtanl0;
  auto deltax = n * cosphi0 - 0.5 * n * theta * Hz * sinphi0;
  auto deltay = n * sinphi0 + 0.5 * n * theta * Hz * cosphi0;

  auto x = x0 + deltax;
  auto y = y0 + deltay;
  auto phi = phi0 + Hz * theta;

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
  auto dZ = (zEnd - track.getZ());
  auto x0 = track.getX();
  auto y0 = track.getY();
  auto phi0 = track.getPhi();
  auto invtanl0 = 1.0 / track.getTanl();;
  double cosphi0, sinphi0;
  o2::utils::sincos(phi0, sinphi0, cosphi0);
  auto n = dZ * invtanl0;
  auto x = x0 + n * cosphi0;
  auto y = y0 + n * sinphi0;

  track.setX(x);
  track.setY(y);
  track.setZ(zEnd);
}
