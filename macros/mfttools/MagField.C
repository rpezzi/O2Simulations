#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <TGeoGlobalMagField.h>

#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"

#endif

double getZField(double x, double y, double z)
{
  const auto grp = o2::parameters::GRPObject::loadFrom("o2sim_grp.root");
  std::unique_ptr<o2::parameters::GRPObject> mGRP = nullptr;
  mGRP.reset(grp);
  o2::base::Propagator::initFieldFromGRP(grp);
  auto field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
  
  double position[3] = {x, y, z};
  auto Bz = field->getBz(position);
  printf("B field z = %f [kGauss]\n", Bz);
  return Bz;
}
