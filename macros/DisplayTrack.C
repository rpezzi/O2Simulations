/// \file DisplayTrack.C
/// \brief Simple macro to display ITSU tracks

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <string>

#include <TEveGeoNode.h>
#include <TEveGeoShape.h>
#include <TEveGeoShapeExtract.h>
#include <TEveManager.h>
#include <TEvePointSet.h>
#include <TFile.h>
#include <TGLViewer.h>
#include <TGeoManager.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>

#include "MFTBase/GeometryTGeo.h"
#include "DataFormatsITSMFT/Cluster.h"
#include "ITSMFTSimulation/Hit.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "MathUtils/Cartesian3D.h"
#include "MathUtils/Utils.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#endif

void DisplayMFTTrack(Int_t event = 0, Int_t track = 0, std::string tracfile = "mfttracks.root", std::string clusfile = "mftclusters.root", std::string hitfile = "o2sim.root", std::string inputGeom = "O2geometry.root")
{
  using namespace o2::base;
  using namespace o2::mft;

  using o2::itsmft::Cluster;
  using o2::itsmft::Hit;

  TFile* f = nullptr;

  if (gEve == nullptr) {
    TEveManager::Create();
  }

  // Load geometry
  if (gGeoManager == nullptr) {
    o2::base::GeometryManager::loadGeometry(inputGeom, "FAIRGeom");
  }

//  gGeoManager->GetVolume("obSuppCyl")->SetInvisible();
//  gGeoManager->GetVolume("ibSuppCyl")->SetInvisible();
//  gGeoManager->GetVolume("ITSUStave0_StaveStruct")->SetInvisible();
//  gGeoManager->GetVolume("ITSUStave1_StaveStruct")->SetInvisible();
//  gGeoManager->GetVolume("ITSUStave2_StaveStruct")->SetInvisible();

//  gGeoManager->GetVolume("ITSUHalfStave0")->SetTransparency(50);
//  gGeoManager->GetVolume("ITSUHalfStave1")->SetTransparency(50);
//  gGeoManager->GetVolume("ITSUHalfStave2")->SetTransparency(50);
//  gGeoManager->GetVolume("ITSUHalfStave3")->SetTransparency(50);
//  gGeoManager->GetVolume("ITSUHalfStave4")->SetTransparency(50);
//  gGeoManager->GetVolume("ITSUHalfStave5")->SetTransparency(50);
//  gGeoManager->GetVolume("ITSV_2")->SetTransparency(80);
//  gGeoManager->GetVolume("HalfConeVolume_1")->SetTransparency(80);
  TGeoNode* tnode = gGeoManager->GetTopVolume()->FindNode("MFT_0");
  gGeoManager->GetTopVolume()->SetTransparency(80);
  TEveGeoTopNode* evenode = new TEveGeoTopNode(gGeoManager, tnode);
  evenode->SetVisLevel(4);
  //gEve->AddGlobalElement(evenode);

  TGLViewer* view = gEve->GetDefaultGLViewer();
  Double_t center[3]{0, 0, -61.4};
  view->CurrentCamera().Configure(3., 1200., center, 10., 5 * 3.14 / 180);
  view->CurrentCamera().Reset();
  //view->ResetCurrentCamera();
  //ViewerRedraw();
  gEve->Redraw3D();

  /*
  // Simplified geometry
  f = TFile::Open("simple_geom_ITS.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f->Get("ITS");
  TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
  gEve->AddElement(gsre,0);
  f->Close();
  */

  // Hits
  f = TFile::Open(hitfile.data());
  TTree* tree = (TTree*)gDirectory->Get("o2sim");

  string s{"event"};
  s += std::to_string(event);
  s += "_hits";
  s += std::to_string(track);
  TEvePointSet* points = new TEvePointSet(s.data());
  points->SetMarkerColor(kBlue);

  std::vector<Hit>* hitArr = nullptr;
  tree->SetBranchAddress("MFTHit", &hitArr);

  tree->GetEvent(event);

  Int_t nc = hitArr->size(), n = 0;
  while (nc--) {
    Hit& c = (*hitArr)[nc];
    points->SetNextPoint(c.GetX(), c.GetY(), c.GetZ());
      n++;
  }
  cout << "Number of hits: " << n << endl;

  gEve->AddElement(points, 0);
  f->Close();

  // Clusters
  f = TFile::Open(clusfile.data());
  tree = (TTree*)gDirectory->Get("o2sim");

  s = "event";
  s += std::to_string(event);
  s += "_clusters";
  s += std::to_string(track);
  points = new TEvePointSet(s.data());
  points->SetMarkerColor(kMagenta);

  std::vector<Cluster>* clusArr = nullptr;
  tree->SetBranchAddress("MFTCluster", &clusArr);
  // Cluster MC labels
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* clsLabArr = nullptr;
  tree->SetBranchAddress("MFTClusterMCTruth", &clsLabArr);

  int tf = 0;
  int lastTF = tree->GetEntries();
  for (; tf < lastTF; ++tf) {
    tree->GetEvent(tf);
    int nc = clusArr->size();
    std::cout << "Clusters nc = " << nc << std::endl;

    for (int i = 0; i < nc; i++) { // Find the TF containing this MC event
      auto mclab = (clsLabArr->getLabels(i))[0];
      auto id = mclab.getEventID();
      if (id == event)
        goto found;
    }
  }
  std::cout << "Time Frame containing the MC event " << event << " was not found" << std::endl;

found:
  std::cout << "MC event " << event << " found in the Time Frame #" << tf << std::endl;
  o2::mft::GeometryTGeo* gman = GeometryTGeo::Instance();
  gman->fillMatrixCache(o2::utils::bit2Mask(o2::TransformType::T2G)); // request cached transforms

  nc = clusArr->size();
  n = 0;
  std::cout << "Clusters nc = " << nc << std::endl;
  while (nc--) {
    Cluster& c = (*clusArr)[nc];
    auto xyz = c.getXYZGlo(*gman);
    points->SetNextPoint(xyz.x(), xyz.y(), xyz.z());
    n++;
//    }
  }
  cout << "Number of clusters: " << n << endl;

  gEve->AddElement(points, 0);
  f->Close();

  // Track
  f = TFile::Open(tracfile.data());
  tree = (TTree*)gDirectory->Get("o2sim");
  std::cout << "Break 1" << std::endl;
  s = "event";
  s += std::to_string(event);
  s += "_track";
  s += std::to_string(track);
  points = new TEvePointSet(s.data());
  points->SetMarkerColor(kGreen);
  std::cout << "Break 2" << std::endl;

  std::vector<TrackLTF>* trkArr = nullptr;
  std::vector<int>* clIdx = nullptr;
  tree->SetBranchAddress("MFTTrackLTF", &trkArr);
  //tree->SetBranchAddress("MFTTrackClusIdx", &clIdx);
  // Track MC labels
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* trkLabArr = nullptr;
  tree->SetBranchAddress("MFTTrackMCTruth", &trkLabArr);

  tree->GetEvent(tf);

  Int_t nt = trkArr->size();
  std::cout << "trkArr->size() nt = " << nt << std::endl;

  TEveTrackList* tracks = new TEveTrackList("tracks");
  auto prop = tracks->GetPropagator();
  prop->SetMagField(0.5);
  prop->SetMaxR(50.);
  prop->SetMaxZ(80.);


  std::vector<float> tXCoordinates, tYCoordinates, tZCoordinates;
  while (nt--) { //loop over tracks
    n = 0;
    const TrackLTF& t = (*trkArr)[nt];
    Int_t nc = t.getNPoints();
    tXCoordinates = t.getXCoordinates();
    tYCoordinates = t.getYCoordinates();
    tZCoordinates = t.getZCoordinates();

    TEveRecTrackD tEve;
    tEve.fP = {tXCoordinates[nc],tYCoordinates[nc],tZCoordinates[nc]};
    //t.fSign = (rec.getSign() < 0) ? -1 : 1;
    TEveTrack* track = new TEveTrack(&tEve, prop);
    track->SetLineColor(kGreen);
    tracks->AddElement(track);

    while (n < nc) { //loop over clusters on each track
      points->SetNextPoint(tXCoordinates[n],tYCoordinates[n],tZCoordinates[n]);
      //std::cout << "Track " << nt << " XYZ = " << tXCoordinates[n] << " , " << tYCoordinates[n] << " , " << tZCoordinates[n] << std::endl;
      n++;
    }
  //  break;
  }
  cout << "Number of attached tracks: " << trkArr->size() << endl;
  tracks->MakeTracks();
  gEve->AddElement(tracks, 0);
  gEve->AddElement(points, 0);
  f->Close();
}

