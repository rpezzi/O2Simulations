
{

  TGeoManager::Import("o2sim_geometry.root");

  TGeoVolume *t = gGeoManager->GetVolume("MFT");
//  TGeoVolume *t = gGeoManager->GetVolume("Support_H0_D3");
  gGeoManager->SetVisLevel(10);
  //t->Draw("");
  t->Raytrace();
  //t->RandomPoints();

}
