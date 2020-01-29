#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

//constexpr Double_t MFTLayerZ[] = {-45.3, -46.7, -48.6, -50.0, -52.4, -53.8, -67.7, -69.1, -76.1, -77.5};

constexpr Double_t pMax = 5;
constexpr Double_t pMin = 0.1;
constexpr Double_t etaMin = -6.0;
constexpr Double_t etaMax = 0.0;

void MFTHitDistributions(const Char_t *SimFile = "o2sim.root") {

  using o2::itsmft::Hit;
  using o2::MCTrackT;
  o2::itsmft::ChipMappingMFT mftChipMapper;
  std::vector<Int_t> trackables_per_event;
  Int_t nMFTTrackables=0;

  Char_t histname[64], histtitle[64];
  Int_t ih = 0;
  std::array<std::unique_ptr<TH1F>, 6> hitRadialDistrib; // 0 to 4 -> hits in corresponding disks. 5-> hits in all disks
  for (auto& h : hitRadialDistrib) {
  snprintf(histname, 64, "hitRadialDistrib_%d", ih);
  snprintf(histtitle, 64, "Hit radial Distribution Disk %d", ih);
  h = std::make_unique<TH1F> (histname, histtitle, 30, 0, 20);
  ++ih;
}

ih = 0;
std::array<std::unique_ptr<TH2F>, 6> hitDistrib; // 0 to 4 -> hits in corresponding disks. 5-> hits in all disks
for (auto& h : hitDistrib) {
snprintf(histname, 64, "hitDistrib_%d", ih);
snprintf(histtitle, 64, "Hit Distribution in Disk %d", ih);
h = std::make_unique<TH2F> (histname, histtitle, 80, -20, 20, 80, -20, 20);
++ih;
}

  std::array<Int_t, 5> numberOfHitperDisk = {0,0,0,0,0};
  Int_t numberOfMFTHits = 0;

  using trackHasHitsinMFTDisks = std::array<bool,5>; // Disks with hits from a MFT track

  std::unique_ptr<TH1F> MCTrackspT = std::make_unique<TH1F> ("MC Tracks pT", "MC Tracks pT", 100, 0, pMax);
  MCTrackspT->GetXaxis()->SetTitle("Transverse p");

  std::unique_ptr<TH2F> hitVertex = std::make_unique<TH2F> ("HitVertexes", "Vertexes of Tracks with hits in the MFT", 1000, -100, 100, 1000, -100, 100);
  hitVertex->GetXaxis()->SetTitle("Z");
  hitVertex->GetYaxis()->SetTitle("X");


  std::unique_ptr<TH1F> MCTracksp = std::make_unique<TH1F> ("MC Tracks p", "MC Tracks p", 100, 0, pMax);
  MCTracksp->GetXaxis()->SetTitle("Total p");

  std::unique_ptr<TH1F> MCTracksEta = std::make_unique<TH1F> ("MC Tracks eta", "MC Tracks Pseudorapidity", 100, etaMin, etaMax);
  MCTracksEta->GetXaxis()->SetTitle("Pseudorapidity");


  std::unique_ptr<TH1I> Trackablility = std::make_unique<TH1I> ("Trackablility", "In how many disks the tracks has hits", 6, 0, 6);
  Trackablility->GetXaxis()->SetTitle("Number of disks");


  std::unique_ptr<TH1I> MultiplicityDistrib = std::make_unique<TH1I> ("Multiplicity", "MFT Trackables Distribution", 10000, 0, 10000);
  MultiplicityDistrib->GetXaxis()->SetTitle("Trackable Multiplicity");




  //Histos for Trackables
  std::unique_ptr<TH1F> TrackablepT = std::make_unique<TH1F> ("Trackables Tracks pT", "Trackables Tracks pT", 100, 0, pMax);
  TrackablepT->GetXaxis()->SetTitle("Transverse p");
  std::unique_ptr<TH1F> Trackablep = std::make_unique<TH1F> ("Trackables Tracks p", "Trackables Tracks p", 100, 0, pMax);
  Trackablep->GetXaxis()->SetTitle("Total p");

  std::unique_ptr<TH1F> MFTTrackablesEta = std::make_unique<TH1F> ("Trackables Tracks eta", "Trackables Pseudorapidity", 100, etaMin, etaMax);
  MFTTrackablesEta->GetXaxis()->SetTitle("Pseudorapidity");





  TFile *simFileIn = new TFile(SimFile);
  //TFile *trkFileIn = new TFile(trkFile);
  TFile outFile("MFTHitDistributions.root","RECREATE");


  TTree *o2SimTree = (TTree*) simFileIn -> Get("o2sim");
  //TTree *mftTrackTree = (TTree*) trkFileIn -> Get("o2sim");

  Int_t numberOfEvents = o2SimTree -> GetEntries();
  std::cout << "numberOfEvents = " << numberOfEvents << std::endl;

  //Int_t nCleanTracksLTF = 0, nCleanTracksCA = 0, nInvalidTracksLTF =0, nInvalidTracksCA = 0;

  vector<Hit>* mfthit = nullptr;
  o2SimTree -> SetBranchAddress("MFTHit",&mfthit);
  vector<MCTrackT<float>>* mcTr = nullptr;
  o2SimTree -> SetBranchAddress("MCTrack",&mcTr);


  o2::dataformats::MCEventHeader* eventHeader = nullptr;
  o2SimTree -> SetBranchAddress("MCEventHeader.",&eventHeader);

//  mftTrackTree->GetEntry(0);
  o2SimTree -> GetEntry(0);
  //Int_t numberOfTracksPerEvent = mcTr->size(); // Number of tracks in first MCEvent (assumed to be the same for all events)
//  vector<eventFoundTracks> allFoundTracksLTF(numberOfEvents), allFoundTracksCA(numberOfEvents); // True for reconstructed tracks - one vector of bool per event


  // Part 1: MC hits and tracks
  //   - Identify trackable tracks (clusters in at least 4 disks)
  //   - Fill histograms



  // Part 1: MC hits and tracks
  // std::cout << "Starting Part 2: MC hits and tracks!" << std::endl;
  for (Int_t event=0; event<numberOfEvents ; event++) { // Loop over events in o2sim
    // std::cout << "Loop over events in o2sim. Event = " << event << std::endl;
    o2SimTree -> GetEntry(event);
    MCTrackT<float>* thisTrack;// =  &(*mcTr)[trID];
    Int_t nMFTHits = mfthit->size(); // Number of mft hits in this event
    Int_t trackables_in_this_event=0;

    // std::cout << "Event " << event << " has " << eventHeader->getMCEventStats().getNKeptTracks() << " tracks and " << nMFTHits << " hits\n";

    std::vector<trackHasHitsinMFTDisks> mcTrackHasHitsInMFTDisks(eventHeader->getMCEventStats().getNKeptTracks(),{0,0,0,0,0}); //

    // std::cout << "Loop over " << nMFTHits << " mfthits to identify trackable MFT tracks in event " <<  event << std::endl;
    for (Int_t n_hit=0 ; n_hit < nMFTHits; n_hit++) { // Loop over mfthits to identify trackable tracks
      Hit* hitp = &(*mfthit).at(n_hit);
      Int_t trID = hitp->GetTrackID(); // ID of the tracks having given the hit
      //std::cout << "n_hit = " << n_hit << " ** trID = " << trID << std::endl;

      Float_t x = hitp->GetX(); // X position of the hit
      Float_t y = hitp->GetY(); // Y position of the hit
      Float_t z = hitp->GetZ(); // Z position of the hit
      Float_t r = sqrt(x*x+y*y);
      thisTrack = &(*mcTr)[trID];

      auto vz = thisTrack->GetStartVertexCoordinatesZ();
      auto vx = thisTrack->GetStartVertexCoordinatesX();
      hitVertex->Fill(vz,vx);
      auto this_disk = mftChipMapper.chip2Layer(hitp->GetDetectorID())/2;
      mcTrackHasHitsInMFTDisks[trID][this_disk] = true;
      hitRadialDistrib[this_disk]->Fill(r);
      hitRadialDistrib[5]->Fill(r);
      numberOfHitperDisk[this_disk]++;
      numberOfMFTHits++;
      hitDistrib[this_disk]->Fill(x,y);
      hitDistrib[5]->Fill(x,y);
      }

    for (Int_t trID=0 ; trID < eventHeader->getMCEventStats().getNKeptTracks(); trID++) { // Loop on MC tracks
      //std::cout << "Loop on tracks to build histos. Track " << trID << " at event " << event << " -> " ;

      //fill MC histograms
      thisTrack =  &(*mcTr)[trID];
      auto z = thisTrack->GetStartVertexCoordinatesZ();
      auto p = thisTrack->GetP();
      auto eta = atanh (thisTrack->GetStartVertexMomentumZ()/p); // eta;
      MCTrackspT->Fill(thisTrack->GetPt());
      MCTracksp->Fill(p);
      MCTracksEta->Fill(eta);


      // Count disks "touched" by the track
      int nMFTDisksHasHits = 0;
      for(auto disk: {0,1,2,3,4}) nMFTDisksHasHits+= int(mcTrackHasHitsInMFTDisks[trID][disk]);
      Trackablility->Fill(nMFTDisksHasHits);
      //std::cout << "nMFTDisksHasHits = " << nMFTDisksHasHits; // << std::endl;

      if(nMFTDisksHasHits>=4) {   //Track is trackable if has left hits on at least 4 disks

        MFTTrackablesEta->Fill(eta);
        trackables_in_this_event++;
}

    //std::cout << " Finished Track " << trID << std::endl;

    } // end Loop on tracks
    trackables_per_event.push_back(trackables_in_this_event);
    nMFTTrackables+=trackables_in_this_event;
    MultiplicityDistrib->Fill(trackables_in_this_event);

    //std::cout << "Finished event " << event << std::endl;
    } // end loop over events


// Write histograms to file
//std::cout << "Writting histograms to file..." << std::endl;

MCTrackspT->Write();
MCTracksp->Write();
MCTracksEta->Write();


Trackablility->Write();
MultiplicityDistrib->Write();

hitVertex->Write();

MFTTrackablesEta->Write();

for (auto& h : hitRadialDistrib) {
h->Write();
}

for (auto& h : hitDistrib) {
h->Write();
}



outFile.Close();

//Int_t totalRecoMFTTracks = nCleanTracksLTF + nCleanTracksCA + nInvalidTracksLTF + nInvalidTracksCA;
std::cout << "Number of MFT trackables = " << nMFTTrackables << std::endl;

auto this_disk=0;
for (auto nHits: numberOfHitperDisk) {
  std::cout << nHits << " hits in MFT disk " << this_disk << std::endl;
  this_disk++;
}
std::cout << numberOfMFTHits << " hits all MFT disks" << std::endl;

new TBrowser;
}
