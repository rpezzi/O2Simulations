#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

//constexpr Double_t MFTLayerZ[] = {-45.3, -46.7, -48.6, -50.0, -52.4, -53.8, -67.7, -69.1, -76.1, -77.5};

constexpr Double_t pMax = 5;
constexpr Double_t pMin = 0.1;
constexpr Double_t etaMin = -3.9;
constexpr Double_t etaMax = -2.0;

void MultiplicityEstimator(const Char_t *SimFile = "o2sim.root", const Char_t *trkFile = "mfttracks.root") {

  using o2::itsmft::Hit;
  using o2::MCTrackT;
  using o2::mft::TrackCA;
  using o2::mft::TrackLTF;
  o2::itsmft::ChipMappingMFT mftChipMapper;
  using eventFoundTracks = std::vector<bool>;

  std::vector<Int_t> trackables_per_event;

  using trackHasHitsinMFTDisks = std::array<bool,5>; // Disks with hits from a MFT track

  std::unique_ptr<TH1F> MCTrackspT = std::make_unique<TH1F> ("MC Tracks pT", "MC Tracks pT", 100, pMin, pMax);
  MCTrackspT->GetXaxis()->SetTitle("Transverse p");
  std::unique_ptr<TH1F> MCTracksp = std::make_unique<TH1F> ("MC Tracks p", "MC Tracks p", 100, pMin, pMax);
  MCTracksp->GetXaxis()->SetTitle("Total p");
  std::unique_ptr<TH1F> MCTrackEta = std::make_unique<TH1F> ("MC Tracks eta", "MC Tracks Pseudorapidity", 100, etaMin, etaMax);
  MCTrackEta->GetXaxis()->SetTitle("Pseudorapidity");

  std::unique_ptr<TH1F> MFTTrackspT = std::make_unique<TH1F> ("MFT Tracks pT", "MFT Tracks pT", 100, pMin, pMax);
  MFTTrackspT->GetXaxis()->SetTitle("Transverse p");
  std::unique_ptr<TH1F> MFTTracksp = std::make_unique<TH1F> ("MFT Tracks p", "MFT Tracks p", 100, pMin, pMax);
  MFTTracksp->GetXaxis()->SetTitle("Total p");
  std::unique_ptr<TH1F> MFTTrackEta = std::make_unique<TH1F> ("MFT Tracks eta", "MFT Tracks Pseudorapidity", 100, etaMin, etaMax);
  MFTTrackEta->GetXaxis()->SetTitle("Pseudorapidity");

  std::unique_ptr<TH1F> LTFTrackspT = std::make_unique<TH1F> ("LTF Tracks pT", "LTF Tracks pT", 100, pMin, pMax);
  LTFTrackspT->GetXaxis()->SetTitle("Transverse p");
  std::unique_ptr<TH1F> LTFTracksp = std::make_unique<TH1F> ("LTF Tracks p", "LTF Tracks p", 100, pMin, pMax);
  LTFTracksp->GetXaxis()->SetTitle("Total p");
  std::unique_ptr<TH1F> LTFTrackEta = std::make_unique<TH1F> ("LTF Tracks eta", "LTF Tracks Pseudorapidity", 100, etaMin, etaMax);
  LTFTrackEta->GetXaxis()->SetTitle("Pseudorapidity");

  std::unique_ptr<TH1F> CATrackspT = std::make_unique<TH1F> ("CA Tracks pT", "CA Tracks pT", 100, pMin, pMax);
  CATrackspT->GetXaxis()->SetTitle("Transverse p");
  std::unique_ptr<TH1F> CATracksp = std::make_unique<TH1F> ("CA Tracks p", "CA Tracks p", 100, pMin, pMax);
  CATracksp->GetXaxis()->SetTitle("Total p");
  std::unique_ptr<TH1F> CATrackEta = std::make_unique<TH1F> ("CA Tracks eta", "LTF Tracks Pseudorapidity", 100, etaMin, etaMax);
  CATrackEta->GetXaxis()->SetTitle("Pseudorapidity");


  std::unique_ptr<TH1F> MCTracksEta5 = std::make_unique<TH1F> ("MC Tracks 5 eta", "-5 cm < zVertex < 5 cm", 100, etaMin, etaMax);
  MCTracksEta5->GetXaxis()->SetTitle("Pseudorapidity");
  std::unique_ptr<TH1F> MCTracksEta5_10pos = std::make_unique<TH1F> ("MC Tracks -5 -10 eta", "-10 cm < zVertex < -5 cm", 100, etaMin, etaMax);
  MCTracksEta5_10pos->GetXaxis()->SetTitle("Pseudorapidity");
  std::unique_ptr<TH1F> MCTracksEta5_10neg = std::make_unique<TH1F> ("MC Tracks 5 10 eta", "5 cm < zVertex < 10 cm", 100, etaMin, etaMax);
  MCTracksEta5_10neg->GetXaxis()->SetTitle("Pseudorapidity");

  std::unique_ptr<TH1F> MFTTracksEta5 = std::make_unique<TH1F> ("MFT Tracks 5 eta", "-5 cm < zVertex < 5 cm", 100, etaMin, etaMax);
  MFTTracksEta5->GetXaxis()->SetTitle("Pseudorapidity");
  std::unique_ptr<TH1F> MFTTracksEta5_10pos = std::make_unique<TH1F> ("MFT Tracks -5 -10 eta", "-10 cm < zVertex < -5 cm", 100, etaMin, etaMax);
  MFTTracksEta5_10pos->GetXaxis()->SetTitle("Pseudorapidity");
  std::unique_ptr<TH1F> MFTTracksEta5_10neg = std::make_unique<TH1F> ("MFT Tracks 5 10 eta", "5 cm < zVertex < 10 cm", 100, etaMin, etaMax);
  MFTTracksEta5_10neg->GetXaxis()->SetTitle("Pseudorapidity");


  std::unique_ptr<TH1F> MCTracksp5 = std::make_unique<TH1F> ("MC Tracks 5 p", "-5 cm < zVertex < 5 cm", 100, pMin, pMax);
  MCTracksp5->GetXaxis()->SetTitle("Total p");
  std::unique_ptr<TH1F> MCTracksp5_10pos = std::make_unique<TH1F> ("MC Tracks -5 -10 p", "-10 cm < zVertex < -5 cm", 100, pMin, pMax);
  MCTracksp5_10pos->GetXaxis()->SetTitle("Total p");
  std::unique_ptr<TH1F> MCTracksp5_10neg = std::make_unique<TH1F> ("MC Tracks 5 10 p", "5 cm < zVertex < 10 cm", 100, pMin, pMax);
  MCTracksp5_10neg->GetXaxis()->SetTitle("Total p");

  std::unique_ptr<TH1F> MFTTracksp5 = std::make_unique<TH1F> ("MFT Tracks 5 p", "-5 cm < zVertex < 5 cm", 100, pMin, pMax);
  MFTTracksp5->GetXaxis()->SetTitle("Total p");
  std::unique_ptr<TH1F> MFTTracksp5_10pos = std::make_unique<TH1F> ("MFT Tracks -5 -10 p", "-10 cm < zVertex < -5 cm", 100, pMin, pMax);
  MFTTracksp5_10pos->GetXaxis()->SetTitle("Total p");
  std::unique_ptr<TH1F> MFTTracksp5_10neg = std::make_unique<TH1F> ("MFT Tracks 5 10 p", "5 cm < zVertex < 10 cm", 100, pMin, pMax);
  MFTTracksp5_10neg->GetXaxis()->SetTitle("Total p");

  std::unique_ptr<TH1I> Trackablility = std::make_unique<TH1I> ("MFT Trackablility", "In how many disks the tracks has hits", 6, 0, 6);
  Trackablility->GetXaxis()->SetTitle("Number of disks");


  std::unique_ptr<TH1I> TrackablesDistrib = std::make_unique<TH1I> ("Trackables", "MFT Trackables Distribution", 1000, 0, 30000);
  TrackablesDistrib->GetXaxis()->SetTitle("Trackable Multiplicity");



  //Histos for Missed (missed tracks that could be tracked)
  std::unique_ptr<TH1F> MissedlepT = std::make_unique<TH1F> ("Missed Tracks pT", "Missed Tracks pT", 100, pMin, pMax);
  std::unique_ptr<TH1F> Missedp = std::make_unique<TH1F> ("Missed Tracks p", "Missed Tracks p", 100, pMin, pMax);
  std::unique_ptr<TH1F> MissedEta = std::make_unique<TH1F> ("Missed Tracks eta", "Missed Pseudorapidity", 100, etaMin, etaMax);

  //2D Histos
  std::unique_ptr<TH2F> MFTTrackedEtaZ = std::make_unique<TH2F> ("MFT_Tracked_eta_z", "Reconstructed Tracks", 31, -15, 16, 25, etaMin, etaMax);
  MFTTrackedEtaZ->GetXaxis()->SetTitle("Vertex PosZ [cm]");
  std::unique_ptr<TH2F> MFTTrackablesEtaZ = std::make_unique<TH2F> ("MFT_Trackables_eta_z", "MFT Trackables:", 31, -15, 16, 25, etaMin, etaMax);
  MFTTrackablesEtaZ->GetXaxis()->SetTitle("Vertex PosZ [cm]");
  std::unique_ptr<TH2F> MCTracksEtaZ = std::make_unique<TH2F> ("MCTracks_eta_z", "MC Tracks: Pseudorapidity vs zVertex", 31, -15, 16, 25, etaMin, etaMax);
  MCTracksEtaZ->GetXaxis()->SetTitle("Vertex PosZ [cm]");



  TFile *simFileIn = new TFile(SimFile);
  TFile *trkFileIn = new TFile(trkFile);
  TFile outFile("MultiplicityEstimator.root","RECREATE");


  TTree *o2SimTree = (TTree*) simFileIn -> Get("o2sim");
  TTree *mftTrackTree = (TTree*) trkFileIn -> Get("o2sim");

  Int_t numberOfEvents = o2SimTree -> GetEntries();
  std::cout << "numberOfEvents = " << numberOfEvents << std::endl;

  Int_t nCleanTracksLTF = 0, nCleanTracksCA = 0, nInvalidTracksLTF =0, nInvalidTracksCA = 0;

  vector<Hit>* mfthit = nullptr;
  o2SimTree -> SetBranchAddress("MFTHit",&mfthit);
  vector<MCTrackT<float>>* mcTr = nullptr;
  o2SimTree -> SetBranchAddress("MCTrack",&mcTr);
  std::vector<o2::mft::TrackCA> trackCAVec, *trackCAVecP = &trackCAVec;
  mftTrackTree->SetBranchAddress("MFTTrackCA", &trackCAVecP);
  std::vector<o2::mft::TrackLTF> trackLTFVec, *trackLTFVecP = &trackLTFVec;
  mftTrackTree->SetBranchAddress("MFTTrackLTF", &trackLTFVecP);

  o2::dataformats::MCEventHeader* eventHeader = nullptr;
  o2SimTree -> SetBranchAddress("MCEventHeader.",&eventHeader);

  mftTrackTree->GetEntry(0);
  o2SimTree -> GetEntry(0);
  //Int_t numberOfTracksPerEvent = mcTr->size(); // Number of tracks in first MCEvent (assumed to be the same for all events)
  vector<eventFoundTracks> allFoundTracksLTF(numberOfEvents*2), allFoundTracksCA(numberOfEvents*2); // True for reconstructed tracks - one vector of bool per event


  for (auto event = 0 ; event < numberOfEvents ; event++) { // Resize vector to accomodate found status of all tracks in all events
    o2SimTree -> GetEntry(event);
    //o2::dataformats::MCEventStats& eventstats = eventHeader->getMCEventStats();
    auto numberOfTracksThisEvent = eventHeader->getMCEventStats().getNKeptTracks();
    //std::cout << "Resizing allFoundTracks for event " << event <<  " with ntracks = " << numberOfTracksThisEvent << std::endl;
    allFoundTracksLTF[event].resize(numberOfTracksThisEvent,false);
    allFoundTracksCA[event].resize(numberOfTracksThisEvent,false);
  }


  // Part 1: Quality of reconstructed MFT tracks
  //   - Loop over reconstructed tracks to identify clean and mixed/noise tracks
  //   - Clean tracks have at least 80% of its clusters from the same track
  //   - If track is not clear it is a mixed/noise track
  // Part 2: MC hits and tracks
  //   - Identify trackable tracks (clusters in at least 4 disks)
  //   - Identify successfully reconstructed tracks
  // Part 3: Calculate Efficiencies


  // Part 1: Reconstructed MFT Tracks
  //std::cout << "Starting Part 1: Reconstructed MFT Tracks!" << std::endl;

  // TracksLTF - Identify reconstructed tracks
  for (const auto &trackLTF : trackLTFVec) {
  auto thisTrackMCCompLabels = trackLTF.getMCCompLabels();
  std::map<Int_t, Int_t> trkIDs;
  for (auto iLabel = 0; iLabel < trackLTF.getNPoints(); iLabel++) { // Count trackIDs of this track
    auto id = thisTrackMCCompLabels[iLabel].getTrackID();
    trkIDs[id]=trkIDs[id]+1;
    //std::cout << thisTrackMCCompLabels[iLabel].getTrackID() << " ";
  }
  //std::cout << std::endl;

  Int_t thisTrkID=-1, thisTEventIDLabel=-1;
  for (auto iLabel = 0; iLabel < trackLTF.getNPoints(); iLabel++) { // Decide if track was successfully tracked
    auto id = thisTrackMCCompLabels[iLabel].getTrackID();
    //std::cout << "ID ; count ; NPoints ; count/NPoints : "<< id << " " << trkIDs[id] << " " << trackLTF.getNPoints()  << " " << 1.0*trkIDs[id]/trackLTF.getNPoints() << std::endl;
    if(1.0*trkIDs[id]/trackLTF.getNPoints() >= 0.8) { // Must have at least 80% of its clusters from the same MC Track
      thisTrkID = id;
      thisTEventIDLabel = iLabel;
    }
  }
  auto eventID =  thisTrackMCCompLabels[thisTEventIDLabel].getEventID();
  //td::cout << "This TrackLTF ID = " << thisTrkID << " in eventID = " << eventID << std::endl;

  if( (thisTrkID > 0 & thisTrkID != 0x7FFFFFF) & (eventID <= numberOfEvents) ) { // If is good match and not noise ...
   allFoundTracksLTF[eventID][thisTrkID]=true;
   nCleanTracksLTF++;
   }
   else {
    //std::cout << "Noise or Mixed TrackLTF!" << std::endl;
    nInvalidTracksLTF++;
  }

  } // Loop on TracksLTF


  // TracksCA - Identify reconstructed tracks
  for (const auto &trackCA : trackCAVec) {
  auto thisTrackMCCompLabels = trackCA.getMCCompLabels();
  std::map<Int_t, Int_t> trkIDs;
  for (auto iLabel = 0; iLabel < trackCA.getNPoints(); iLabel++) { // Count trackIDs of this track
    auto id = thisTrackMCCompLabels[iLabel].getTrackID();
    trkIDs[id]=trkIDs[id]+1;
    //std::cout << thisTrackMCCompLabels[iLabel].getTrackID() << " ";
  }
  //std::cout << std::endl;

  Int_t thisTrkID=-1, thisTEventIDLabel=-1;
  for (auto iLabel = 0; iLabel < trackCA.getNPoints(); iLabel++) { //Decide if track was successfully tracked
    auto id = thisTrackMCCompLabels[iLabel].getTrackID();
    //std::cout << "ID ; count ; NPoints ; count/NPoints : "<< id << " " << trkIDs[id] << " " << trackCA.getNPoints()  << " " << 1.0*trkIDs[id]/trackCA.getNPoints() << std::endl;
    if(1.0*trkIDs[id]/trackCA.getNPoints() >= 0.8) {
      thisTrkID = id;
      thisTEventIDLabel = iLabel;
    }
  }

  auto eventID =  thisTrackMCCompLabels[thisTEventIDLabel].getEventID();
  //std::cout << "This TrackCA ID = " << thisTrkID << " in eventID = " << eventID << std::endl;
  if( (thisTrkID > 0 & thisTrkID != 0x7FFFFFF) & (eventID <= numberOfEvents) ) { // If is good match and not noise ...
    allFoundTracksCA[eventID][thisTrkID]=true;
    nCleanTracksCA++;
  }
  else {
    //std::cout << "Noise or Mixed TrackCA!" << std::endl;
    nInvalidTracksCA++;
  }

  } // Loop on TracksCA


  // Part 2: MC hits and tracks
   //std::cout << "Starting Part 2: MC hits and tracks!" << std::endl;
  for (Int_t event=0; event<numberOfEvents ; event++) { // Loop over events in o2sim
    // std::cout << "Loop over events in o2sim. Event = " << event << std::endl;
    o2SimTree -> GetEntry(event);
    Int_t nMFTHits = mfthit->size(); // Number of mft hits in this event
    Int_t trackables_in_this_event=0;

     //std::cout << "Event " << event << " has " << eventHeader->getMCEventStats().getNKeptTracks() << " tracks and " << nMFTHits << " hits\n";

    std::vector<trackHasHitsinMFTDisks> mcTrackHasHitsInMFTDisks(eventHeader->getMCEventStats().getNKeptTracks(),{0,0,0,0,0}); //

     std::cout << "Loop over " << nMFTHits << " mfthits to identify trackable MFT tracks in event " <<  event << std::endl;
    for (Int_t n_hit=0 ; n_hit < nMFTHits; n_hit++) { // Loop over mfthits to identify trackable tracks
      Hit* hitp = &(*mfthit).at(n_hit);
      Int_t trID = hitp->GetTrackID(); // ID of the tracks having given the hit
      //std::cout << "n_hit = " << n_hit << " ** trID = " << trID << std::endl;

      Float_t z = hitp->GetZ(); // Z position of the hit => Identify MFT disk
      mcTrackHasHitsInMFTDisks[trID][mftChipMapper.chip2Layer(hitp->GetDetectorID())/2] = true;
      }

    for (Int_t trID=0 ; trID < eventHeader->getMCEventStats().getNKeptTracks(); trID++) { // Loop on MC tracks
      //std::cout << "Loop on tracks to build histos. Track " << trID << " at event " << event << " -> " ;

      //fill MC histograms
      MCTrackT<float>* thisTrack =  &(*mcTr)[trID];
      auto z = thisTrack->GetStartVertexCoordinatesZ();
      auto p = thisTrack->GetP();
      auto eta = atanh (thisTrack->GetStartVertexMomentumZ()/p); // eta;
      MCTrackspT->Fill(thisTrack->GetPt());
      MCTracksp->Fill(p);
      MCTrackEta->Fill(eta);
      MCTracksEtaZ->Fill(z,eta);
      if( (z>-5) & (z<5) ) {
        MCTracksEta5->Fill(eta);
        MCTracksp5->Fill(p);
      }
      if( (z>5) & (z<10) ) {
        MCTracksEta5_10pos->Fill(eta);
        MCTracksp5_10pos->Fill(p);
      }
      if( (z>-10) & (z<-5) ) {
        MCTracksEta5_10neg->Fill(eta);
        MCTracksp5_10neg->Fill(p);
      }

      // Count disks "touched" by the track
      int nMFTDisksHasHits = 0;
      for(auto disk: {0,1,2,3,4}) nMFTDisksHasHits+= int(mcTrackHasHitsInMFTDisks[trID][disk]);
      Trackablility->Fill(nMFTDisksHasHits);
      //std::cout << "nMFTDisksHasHits = " << nMFTDisksHasHits; // << std::endl;

      if(nMFTDisksHasHits>=4) {   //Track is trackable if has left hits on at least 4 disks

        MFTTrackablesEtaZ->Fill(z,eta);
        trackables_in_this_event++;

        bool wasFound = allFoundTracksLTF[event][trID] | allFoundTracksCA[event][trID];

        if(wasFound) {
          MFTTrackspT->Fill(thisTrack->GetPt());
          MFTTracksp->Fill(p);
          MFTTrackEta->Fill(eta);


          if( (z>-5) & (z<5) ) {
            MFTTracksEta5->Fill(eta);
            MFTTracksp5->Fill(p);
          }
          if( (z>5) & (z<10) ) {
            MFTTracksEta5_10pos->Fill(eta);
            MFTTracksp5_10pos->Fill(p);
          }
          if( (z>-10) & (z<-5) ) {
            MFTTracksEta5_10neg->Fill(eta);
            MFTTracksp5_10neg->Fill(p);
          }


          if(allFoundTracksLTF[event][trID]) {
            LTFTrackspT->Fill(thisTrack->GetPt());
            LTFTracksp->Fill(p);
            LTFTrackEta->Fill(eta);
            MFTTrackedEtaZ->Fill(thisTrack->GetStartVertexCoordinatesZ(),eta);


          }
          if(allFoundTracksCA[event][trID]) {
            CATrackspT->Fill(thisTrack->GetPt());
            CATracksp->Fill(p);
            CATrackEta->Fill(eta);
          }
        }
      } else {  // Fill histograms for Missed Tracks
        MissedlepT->Fill(thisTrack->GetPt());
        Missedp->Fill(p);
        MissedEta->Fill(eta);
      }
    //std::cout << " Finished Track " << trID << std::endl;

    } // end Loop on tracks
    trackables_per_event.push_back(trackables_in_this_event);
    TrackablesDistrib->Fill(trackables_in_this_event);

    std::cout << "Finished event " << event << std::endl;
    } // end loop over events

// Part 3: Calculate Efficiencies
std::cout << "Building efficiencies histos..." << std::endl;
TH1F MFTEfficiencypT = (*MFTTrackspT)/ (*MCTrackspT);
TH1F MFTTEfficiencyp = (*MFTTracksp) / (*MCTracksp);
TH1F MFTEfficiencyEta = (*MFTTrackEta) / (*MCTrackEta);
TH2F MFTTrackerEfficiency = (*MFTTrackedEtaZ) / (*MFTTrackablesEtaZ);
TH2F MFTEfficiency2D = (*MFTTrackedEtaZ) / (*MCTracksEtaZ);
TH2F MFTAcceptance = (*MFTTrackablesEtaZ) / (*MCTracksEtaZ);

TH1F MFTEffsEta5 = (*MFTTracksEta5)/(*MCTracksEta5);
TH1F MFTEffEta5_10pos = (*MFTTracksEta5_10pos)/(*MCTracksEta5_10pos);
TH1F MFTEffEta5_10neg = (*MFTTracksEta5_10neg)/(*MCTracksEta5_10neg);


TH1F MFTEffsp5 = (*MFTTracksp5)/(*MCTracksp5);
TH1F MFTEffp5_10pos = (*MFTTracksp5_10pos)/(*MCTracksp5_10pos);
TH1F MFTEffp5_10neg = (*MFTTracksp5_10neg)/(*MCTracksp5_10neg);


MFTEfficiencypT.SetNameTitle("MFT Efficiency pT", "MFT Efficiency pT");
MFTTEfficiencyp.SetNameTitle("MFT Efficiency p", "MFT Efficiency p");
MFTEfficiencyEta.SetNameTitle("MFT Efficiency eta", "MFT Efficiency Pseudorapidity");
MFTTrackerEfficiency.SetNameTitle("MFT Tracker Efficiency", "MFT Tracker Efficiency");
MFTTrackerEfficiency.GetXaxis()->SetTitle("Vertex PosZ [cm]");

MFTEfficiency2D.SetNameTitle("MFT Efficiency", "MFT Efficiency");
MFTEfficiency2D.GetXaxis()->SetTitle("Vertex PosZ [cm]");

MFTAcceptance.SetNameTitle("MFT Acceptance", "MFT Acceptance");
MFTAcceptance.GetXaxis()->SetTitle("Vertex PosZ [cm]");

MFTEffsEta5.SetNameTitle("MFT Eta Efficiency5_5", "-5 cm < z < 5 cm");
MFTEffsEta5.GetXaxis()->SetTitle("Pseudorapidity");
MFTEffEta5_10pos.SetNameTitle("MFT Eta Efficiency_5_10pos", "5 cm < z < 10 cm");
MFTEffEta5_10pos.GetXaxis()->SetTitle("Pseudorapidity");
MFTEffEta5_10neg.SetNameTitle("MFT Eta Efficiency_10_5neg", "-10 cm < z < -5 cm");
MFTEffEta5_10neg.GetXaxis()->SetTitle("Pseudorapidity");

MFTEffsp5.SetNameTitle("MFT P Efficiency5_5", "-5 cm < z < 5 cm");
MFTEffsp5.GetXaxis()->SetTitle("P (GeV)");
MFTEffp5_10pos.SetNameTitle("MFT P Efficiency_5_10pos", "5 cm < z < 10 cm");
MFTEffp5_10pos.GetXaxis()->SetTitle("P (GeV)");
MFTEffp5_10neg.SetNameTitle("MFT P Efficiency_10_5neg", "-10 cm < z < -5 cm");
MFTEffp5_10neg.GetXaxis()->SetTitle("P (GeV)");

MFTEffsEta5.SetNameTitle("MFT Eta Efficiency5_5", "-5 cm < z < 5 cm");
MFTEffsEta5.GetXaxis()->SetTitle("Pseudorapidity");
MFTEffEta5_10pos.SetNameTitle("MFT Eta Efficiency_5_10pos", "5 cm < z < 10 cm");
MFTEffEta5_10pos.GetXaxis()->SetTitle("Pseudorapidity");
MFTEffEta5_10neg.SetNameTitle("MFT Eta Efficiency_10_5neg", "-10 cm < z < -5 cm");
MFTEffEta5_10neg.GetXaxis()->SetTitle("Pseudorapidity");


// Stacks
THStack mftEtaStack("PStack","MFT Tracks");
MFTEffsEta5.SetLineColor(kBlack);
mftEtaStack.Add(&MFTEffsEta5,"nostack");
MFTEffEta5_10pos.SetLineColor(kBlue);
mftEtaStack.Add(&MFTEffEta5_10pos,"nostack");
MFTEffEta5_10neg.SetLineColor(kRed);
mftEtaStack.Add(&MFTEffEta5_10neg,"nostack");



// Write histograms to file
std::cout << "Writting histograms to file..." << std::endl;

MCTrackspT->Write();
MCTracksp->Write();
MCTrackEta->Write();

MFTTrackspT->Write();
MFTTracksp->Write();
MFTTrackEta->Write();

LTFTrackspT->Write();
LTFTracksp->Write();
LTFTrackEta->Write();

CATrackspT->Write();
CATracksp->Write();
CATrackEta->Write();

MFTEfficiencypT.Write();
MFTTEfficiencyp.Write();
MFTEfficiencyEta.Write();

MCTracksEta5->Write();
MCTracksEta5_10pos->Write();
MCTracksEta5_10neg->Write();

MFTTracksEta5->Write();
MFTTracksEta5_10pos->Write();
MFTTracksEta5_10neg->Write();

MCTracksp5->Write();
MCTracksp5_10pos->Write();
MCTracksp5_10neg->Write();


MFTTracksp5->Write();
MFTTracksp5_10pos->Write();
MFTTracksp5_10neg->Write();

MFTEffsEta5.Write();
MFTEffEta5_10pos.Write();
MFTEffEta5_10neg.Write();

MFTEffsp5.Write();
MFTEffp5_10pos.Write();
MFTEffp5_10neg.Write();

mftEtaStack.Write();

MissedlepT->Write();
Missedp->Write();
MissedEta->Write();

Trackablility->Write();
TrackablesDistrib->Write();

MCTracksEtaZ->SetOption("CONT4");
MCTracksEtaZ->Write();

MFTTrackablesEtaZ->SetOption("CONT4");
MFTTrackablesEtaZ->Write();

MFTTrackedEtaZ->SetOption("CONT4");
MFTTrackedEtaZ->Write();


MFTTrackerEfficiency.SetOption("CONT4");
MFTTrackerEfficiency.Write();

MFTEfficiency2D.SetOption("CONT4");
MFTEfficiency2D.Write();

MFTAcceptance.SetOption("CONT4");
MFTAcceptance.Write();


outFile.Close();

Int_t totalRecoMFTTracks = nCleanTracksLTF + nCleanTracksCA + nInvalidTracksLTF + nInvalidTracksCA;
Int_t totalCleanMFTTracks = nCleanTracksLTF + nCleanTracksCA;
Int_t totalInvalidMFTTracks = nInvalidTracksLTF + nInvalidTracksCA;
std::cout << "---------------------------------------------------" << std::endl;
std::cout << "Number of reconstructed MFT Tracks = " << totalRecoMFTTracks << std::endl;
std::cout << "Number of clean MFT Tracks = " << totalCleanMFTTracks << " (" << 100.f*totalCleanMFTTracks/(totalRecoMFTTracks) << " %)" << std::endl;
std::cout << "Number of mixed MFT Tracks = " << totalInvalidMFTTracks << " (" << 100.f*totalInvalidMFTTracks/(totalRecoMFTTracks) << " %)" << std::endl;
std::cout << "---------------------------------------------------" << std::endl;
std::cout << "nCleanTracksLTF = " << nCleanTracksLTF << std::endl;
std::cout << "nCleanTracksCA = " << nCleanTracksCA << std::endl;
std::cout << "nInvalidTracksLTF = " << nInvalidTracksLTF  << " (" << 100.f*nInvalidTracksLTF/(nCleanTracksLTF+nInvalidTracksLTF) << " % of LTF tracks)" << std::endl;
std::cout << "nInvalidTracksCA = " << nInvalidTracksCA << " (" << 100.f*nInvalidTracksCA/(nCleanTracksCA+nInvalidTracksCA) << " % of CA tracks)" << std::endl;
std::cout << "---------------------------------------------------" << std::endl;

//new TBrowser;

}
