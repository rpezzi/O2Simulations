#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"


TFile* eventFile = nullptr; //! the file containing the persistent events
int eventsAvailable = 0;

void KinematicsGraphs(bool filterSecondaries = kFALSE,
                      const Char_t *kinFile = "Kinematics.root",
                      const Double_t pMin = 0,
                      const Double_t pMax = 20,
                      const Double_t etaMin = -10.0,
                      const Double_t etaMax = 10.0
                     )
{

  std::unique_ptr<TH1F> MCTrackspT = std::make_unique<TH1F> ("MC Tracks pT", "MC Tracks pT", 100, pMin, pMax);
  MCTrackspT->GetXaxis()->SetTitle("Transverse p");
  std::unique_ptr<TH1F> MCTracksp = std::make_unique<TH1F> ("MC Tracks p", "MC Tracks p", 100, pMin, pMax);
  MCTracksp->GetXaxis()->SetTitle("Total p");
  std::unique_ptr<TH1F> MCTrackEta = std::make_unique<TH1F> ("MC Tracks eta", "MC Tracks Pseudorapidity", 100, etaMin, etaMax);
  MCTrackEta->GetXaxis()->SetTitle("Pseudorapidity");

  std::unique_ptr<TH2F> VertexXZ = std::make_unique<TH2F> ("VertexXZ", "Vertexes of Tracks with hits in the MFT", 1000, -100, 100, 1000, -100, 100);
  VertexXZ->GetXaxis()->SetTitle("Z");
  VertexXZ->GetYaxis()->SetTitle("X");

  std::unique_ptr<TH2F> VertexYZ = std::make_unique<TH2F> ("VertexYZ", "Vertexes of Tracks with hits in the MFT", 1000, -100, 100, 1000, -100, 100);
  VertexYZ->GetXaxis()->SetTitle("Z");
  VertexYZ->GetYaxis()->SetTitle("Y");

  std::unique_ptr<TH2F> VertexXY = std::make_unique<TH2F> ("VertexXY", "Vertexes of Tracks with hits in the MFT", 1000, -10, 10, 1000, -10, 10);
  VertexXY->GetXaxis()->SetTitle("Z");
  VertexXY->GetYaxis()->SetTitle("Y");

  eventFile = TFile::Open(kinFile);
if (eventFile == nullptr) {
  std::cout << "EventFile " << kinFile << " not found \n";
  return;
}
// the kinematics will be stored inside a Tree "TreeK" with branch "Particles"
// different events are stored inside TDirectories

// we need to probe for the number of events
TObject* object = nullptr;
do {
  std::stringstream eventstringstr;
  eventstringstr << "Event" << eventsAvailable;
  // std::cout << "probing for " << eventstring << "\n";
  object = eventFile->Get(eventstringstr.str().c_str());
  // std::cout << "got " << object << "\n";
  if (object != nullptr)
    eventsAvailable++;
} while (object != nullptr);
std::cout << "Found " << eventsAvailable << " events in this file \n";



for (auto eventCounter = 0; eventCounter < eventsAvailable; eventCounter++) {

    // get the tree and the branch
    std::stringstream treestringstr;
    treestringstr << "Event" << eventCounter << "/TreeK";
    TTree* tree = (TTree*)eventFile->Get(treestringstr.str().c_str());
    if (tree == nullptr) {
      return kFALSE;
    }

    auto branch = tree->GetBranch("Particles");
    TParticle* particle = nullptr;
    branch->SetAddress(&particle);
    std::cout  << "Reading " << branch->GetEntries() << " particles from event " << eventCounter << "\n";

    // read the whole kinematics initially
    std::vector<TParticle> particles;
    for (int i = 0; i < branch->GetEntries(); ++i) {
      branch->GetEntry(i);
      if( filterSecondaries && particle->GetStatusCode() < 1 ) continue;
      particles.push_back(*particle);
    }

    for(std::vector<TParticle>::iterator particle = particles.begin(); particle != particles.end(); ++particle) {
      auto pdgid = particle->GetPdgCode();
      auto p = particle->P();
      auto pt = particle->Pt();
      auto eta = particle->Eta();
      auto vx = particle->Vx();
      auto vy = particle->Vy();
      auto vz = particle->Vz();
      auto parent = -1;
      auto e = particle->Energy();
      auto tof = particle->T();
      auto weight = particle->GetWeight();
      MCTrackspT->Fill(pt);
      MCTracksp->Fill(p);
      MCTrackEta->Fill(eta);
      VertexXZ->Fill(vz,vx);
      VertexYZ->Fill(vz,vy);
      VertexXY->Fill(vx,vy);
      //std::cout  << "Putting primary " << pdgid << " " << particle.GetStatusCode() << " " << particle.GetUniqueID() << std::endl;
    }
    std::cout << "Finished event " << eventCounter << std::endl;

}



//  TFile *simFileIn = new TFile(SimFile);
  TFile outFile("Kinematics_graphs.root","RECREATE");

MCTrackspT->Write();
MCTracksp->Write();
MCTrackEta->Write();
VertexXZ->Write();
VertexYZ->Write();
VertexXY->Write();

outFile.Close();

}
