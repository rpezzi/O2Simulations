#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

void HitTrackCorrelation(const Char_t *nameFile = "o2sim.root") {

  using o2::itsmft::Hit;
  using o2::MCTrackT;

  TFile *fileIn = new TFile(nameFile);
  TTree *tr = (TTree*) fileIn -> Get("o2sim");

  vector<Hit>* hit = nullptr;
  tr -> SetBranchAddress("MFTHit",&hit);
  vector<MCTrackT<float>>* mcTr = nullptr;
  tr -> SetBranchAddress("MCTrack",&mcTr);

  Int_t nbE = tr -> GetEntries();

  for (Int_t i=0; i<nbE ; i++) {

    tr -> GetEntry(i);

    Int_t nbTr = mcTr->size();
    Int_t nbH = hit->size();

    for (Int_t j=0 ; j<nbH; j++) {

      Hit* c = &(*hit)[j];
      Int_t myID = c->GetTrackID(); // ID of the tracks having given the hit
      MCTrackT<float>* myTr =  &(*mcTr)[myID];

      // here you can analyse the information of you MC track

     }
  }

} 
