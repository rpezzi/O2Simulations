#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

constexpr Double_t MFTLayerZ[] = {-45.3, -46.7, -48.6, -50.0, -52.4, -53.8, -67.7, -69.1, -76.1, -77.5};

constexpr Double_t ITSLayerRMin[] = {2.24, 3.01, 3.78};
constexpr Double_t ITSLayerRMax[] = {2.67, 3.46, 4.21};


constexpr Double_t pMax = 4;
constexpr Double_t pMin = 0.1;
constexpr Double_t etaMin = -3.9;
constexpr Double_t etaMax = -2.0;

void CheckITSMFTCommonTracks(const Char_t *SimFile = "o2sim.root", const Char_t *trkFile = "mfttracks.root") {

  using o2::itsmft::Hit;
  using o2::MCTrackT;
  using o2::mft::TrackCA;
  using o2::mft::TrackLTF;
  using eventFoundTracks = std::vector<bool>;

  using trackHasHitsinMFTDisks = std::array<bool,5>; // Disks with hits from a MFT track
  using trackHasHitsinITSLayers = std::array<bool,3>; // Disks with hits from a MFT track


  std::unique_ptr<TH1F> MCTrackspT = std::make_unique<TH1F> ("MC Tracks pT", "MC Tracks pT", 100, 0, pMax);
  MCTrackspT->GetXaxis()->SetTitle("Transverse p");
  std::unique_ptr<TH1F> MCTracksp = std::make_unique<TH1F> ("MC Tracks p", "MC Tracks p", 100, 0, pMax);
  MCTracksp->GetXaxis()->SetTitle("Total p");
  std::unique_ptr<TH1F> MCTrackRap = std::make_unique<TH1F> ("MC Tracks eta", "MC Tracks Pseudorapidity", 100, etaMin, etaMax);
  MCTrackRap->GetXaxis()->SetTitle("Pseudorapidity");

  std::unique_ptr<TH1F> MFTTrackspT = std::make_unique<TH1F> ("MFT Tracks pT", "MFT Tracks pT", 100, 0, pMax);
  MFTTrackspT->GetXaxis()->SetTitle("Transverse p");
  std::unique_ptr<TH1F> MFTTracksp = std::make_unique<TH1F> ("MFT Tracks p", "MFT Tracks p", 100, 0, pMax);
  MFTTracksp->GetXaxis()->SetTitle("Total p");
  std::unique_ptr<TH1F> MFTTrackRap = std::make_unique<TH1F> ("MFT Tracks eta", "MFT Tracks Pseudorapidity", 100, etaMin, etaMax);
  MFTTrackRap->GetXaxis()->SetTitle("Pseudorapidity");


  std::unique_ptr<TH1I> MFTTrackablility = std::make_unique<TH1I> ("MFTTrackablility", "In how many MFT disks the tracks has hits", 6, 0, 6);
  MFTTrackablility->GetXaxis()->SetTitle("Number of MFT disks");
  std::unique_ptr<TH1I> ITSTrackablility = std::make_unique<TH1I> ("ITSTrackablility", "In how many ITS Layers the tracks has hits", 4, 0, 4);
  ITSTrackablility->GetXaxis()->SetTitle("Number of ITS layers");



  //Histos for Trackables
  std::unique_ptr<TH1F> TrackablepT = std::make_unique<TH1F> ("Trackables Tracks pT", "Trackables Tracks pT", 100, 0, pMax);
  TrackablepT->GetXaxis()->SetTitle("Transverse p");
  std::unique_ptr<TH1F> Trackablep = std::make_unique<TH1F> ("Trackables Tracks p", "Trackables Tracks p", 100, 0, pMax);
  Trackablep->GetXaxis()->SetTitle("Total p");


  //2D Histos
  std::unique_ptr<TH2F> MFTTrackedEtaZ = std::make_unique<TH2F> ("MFT_Tracked_eta_z", "Reconstructed Tracks: Pseudorapidity vs zVertex", 31, -15, 16, 25, etaMin, etaMax);
  MFTTrackedEtaZ->GetXaxis()->SetTitle("Vertex PosZ [cm]");
  std::unique_ptr<TH2F> MFTAccepEtaZ = std::make_unique<TH2F> ("MFT_Acceptance_eta_z", "MFT Acceptance (Trackables): Pseudorapidity vs zVertex", 31, -15, 16, 25, etaMin, etaMax);
  MFTAccepEtaZ->GetXaxis()->SetTitle("Vertex PosZ [cm]");
  std::unique_ptr<TH2F> MCTracksEtaZ = std::make_unique<TH2F> ("MCTracks_eta_z", "MC Tracks: Pseudorapidity vs zVertex", 31, -15, 16, 25, etaMin, etaMax);
  MCTracksEtaZ->GetXaxis()->SetTitle("Vertex PosZ [cm]");

  std::unique_ptr<TH2F> MFTandITSLInnerBarrel = std::make_unique<TH2F> ("MFT and ITSInner", "MFT Trackables with hits in all ITS Inner Layers", 25, -5, 20, 25, etaMin, etaMax);
  MFTandITSLInnerBarrel->GetXaxis()->SetTitle("Vertex PosZ [cm]");
  std::unique_ptr<TH2F> MFTandTwoITSLInnerLayers = std::make_unique<TH2F> ("MFT and 2 ITSInner", "MFT Trackables with 2 or more hits in ITS Inner Layers", 25, -5, 20, 25, etaMin, etaMax);
  MFTandTwoITSLInnerLayers->GetXaxis()->SetTitle("Vertex PosZ [cm]");
  std::unique_ptr<TH2F> MFTandOneITSLInnerLayer = std::make_unique<TH2F> ("MFT and 1 ITSInner", "MFT Trackables with 1 or more hits in ITS Inner Layers", 25, -5, 20, 25, etaMin, etaMax);
  MFTandOneITSLInnerLayer->GetXaxis()->SetTitle("Vertex PosZ [cm]");


  TFile *simFileIn = new TFile(SimFile);
  TFile *trkFileIn = new TFile(trkFile);
  TFile outFile("ITSMFTCommonTracks.root","RECREATE");


  TTree *o2SimTree = (TTree*) simFileIn -> Get("o2sim");
  TTree *mftTrackTree = (TTree*) trkFileIn -> Get("o2sim");

  Int_t numberOfEvents = o2SimTree -> GetEntries();
  std::cout << "numberOfEvents = " << numberOfEvents << std::endl;

  Int_t nCleanTracksLTF = 0, nCleanTracksCA = 0, nInvalidTracksLTF =0, nInvalidTracksCA = 0;

  vector<Hit>* mfthit = nullptr;
  vector<Hit>* itshit = nullptr;
  o2SimTree -> SetBranchAddress("MFTHit",&mfthit);
  o2SimTree -> SetBranchAddress("ITSHit",&itshit);
  vector<MCTrackT<float>>* mcTr = nullptr;
  o2SimTree -> SetBranchAddress("MCTrack",&mcTr);
  std::vector<o2::mft::TrackCA> trackCAVec, *trackCAVecP = &trackCAVec;
  mftTrackTree->SetBranchAddress("MFTTrackCA", &trackCAVecP);
  std::vector<o2::mft::TrackLTF> trackLTFVec, *trackLTFVecP = &trackLTFVec;
  mftTrackTree->SetBranchAddress("MFTTrackLTF", &trackLTFVecP);

  mftTrackTree->GetEntry(0);
  o2SimTree -> GetEntry(0);
  Int_t numberOfTracksPerEvent = mcTr->size(); // Number of tracks in first MCEvent (assumed to be the same for all events)
  vector<eventFoundTracks> allFoundTracksLTF(numberOfEvents+10), allFoundTracksCA(numberOfEvents+10); // True for reconstructed tracks - one vector of bool per event
  for (auto event = 0 ; event < numberOfEvents ; event++) { // Resize vector to accomodate found status of all tracks in all events
    //std::cout << "Resizing allFoundTracks for event " << event <<  " with ntracks = " << numberOfTracksPerEvent << std::endl;
    allFoundTracksLTF[event].resize(numberOfTracksPerEvent+10,false);
    allFoundTracksCA[event].resize(numberOfTracksPerEvent+10,false);
  }

  // Part 1: Reconstructed MFT Tracks
  //  Loop over reconstructed tracks to identify clean and mixed/noise tracks
  //   1.1 Clean tracks have at least 80% of its clusters from the same track
  //   1.2 If track is not clear it is a invalid (mixed/noise) track
  // Part 2: MC hits and tracks
  //   2.1 Identify trackable tracks (clusters in at least 4 disks)
  //   2.2 Identify successfully reconstructed tracks
  //   2.3 Fill Histograms


  // Part 1: Reconstructed MFT Tracks
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
  for (auto iLabel = 0; iLabel < trackLTF.getNPoints(); iLabel++) { //Decide if track was successfully tracked
    auto id = thisTrackMCCompLabels[iLabel].getTrackID();
    //std::cout << "ID ; count ; NPoints ; count/NPoints : "<< id << " " << trkIDs[id] << " " << trackLTF.getNPoints()  << " " << 1.0*trkIDs[id]/trackLTF.getNPoints() << std::endl;
    if(1.0*trkIDs[id]/trackLTF.getNPoints() >= 0.8) {
      thisTrkID = id;
      thisTEventIDLabel = iLabel;
    }
  }
  auto eventID =  thisTrackMCCompLabels[thisTEventIDLabel].getEventID();
  //std::cout << "This TrackLTF ID = " << thisTrkID << " in eventID = " << eventID << std::endl;

  if( (thisTrkID > 0 & thisTrkID != 0x7FFFFFF) & (eventID <= numberOfEvents) ) {
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
  if( (thisTrkID > 0 & thisTrkID != 0x7FFFFFF) & (eventID <= numberOfEvents) ) {
    allFoundTracksCA[eventID][thisTrkID]=true;
    nCleanTracksCA++;
  }
  else {
    //std::cout << "Noise or Mixed TrackCA!" << std::endl;
    nInvalidTracksCA++;
  }

  } // Loop on TracksCA


  // Part 2: MC hits and tracks
  for (Int_t event=0; event<numberOfEvents ; event++) { // Loop over events in o2sim
    //std::cout << "Loop over events in o2sim. Event = " << event << std::endl;
    o2SimTree -> GetEntry(event);
    Int_t nMFTHits = mfthit->size(); // Number of mft hits in this event
    Int_t nITSHits = itshit->size(); // Number of its hits in this event

    //std::cout << "Event " << event << " has " << numberOfTracksPerEvent << " tracks and " << nMFTHits << " hits\n";

    std::vector<trackHasHitsinMFTDisks> mcTrackHasHitsInMFTDisks(numberOfTracksPerEvent+10,{0,0,0,0,0}); //
    std::vector<trackHasHitsinITSLayers> mcTrackHasHitsInITSLayers(numberOfTracksPerEvent+10,{0,0,0}); //


    //std::cout << "Loop over " << nMFTHits << " mfthits to identify trackable MFT tracks in event " <<  event << std::endl;
    for (Int_t n_hit=0 ; n_hit < nMFTHits; n_hit++) { // Loop over mfthits to identify trackable tracks
      Hit* hitp = &(*mfthit).at(n_hit);
      Int_t trID = hitp->GetTrackID(); // ID of the tracks having given the hit
      //std::cout << "n_hit = " << n_hit << " ** trID = " << trID << std::endl;

      Float_t z = hitp->GetZ(); // Z position of the hit => Identify MFT disk.
      for(auto disk: {0,1,2,3,4}) if( z < MFTLayerZ[disk*2] + .3  & z > MFTLayerZ[disk*2+1] -.3 ) mcTrackHasHitsInMFTDisks[trID][disk] = true;
      }

      //std::cout << "Loop over " << nITSHits << " itshits to identify trackable ITS tracks in event " <<  event << std::endl;
      for (Int_t n_hit=0 ; n_hit < nITSHits; n_hit++) { // Loop over itshits to identify trackable tracks
        Hit* hitp = &(*itshit).at(n_hit);
        Int_t trID = hitp->GetTrackID(); // ID of the tracks having given the hit
        //std::cout << "n_hit = " << n_hit << " ** trID = " << trID << std::endl;

        Float_t x = hitp->GetX(); // X position of the hit => Identify ITS Layer.
        Float_t y = hitp->GetY(); // Y position of the hit => Identify ITS Layer.
        Float_t r = sqrt(x*x+y*y);
          //std::cout << std::endl << "R = " << r << " -> ";
        for(auto layer: {0,1,2}) if( r < ITSLayerRMax[layer] + .1  & r > ITSLayerRMin[layer] -.1 )
        {
          //std::cout << "ITS Layer: " << layer << std::endl;
          mcTrackHasHitsInITSLayers[trID][layer] = true;
        }
        }

    for (Int_t trID=0 ; trID < numberOfTracksPerEvent; trID++) { // Loop on tracks
      //std::cout << "Loop on tracks to build histos. Track " << trID << " at event " << event << " -> " ;

      //fill MC histograms
      MCTrackT<float>* thisTrack =  &(*mcTr)[trID];
      auto z = thisTrack->GetStartVertexCoordinatesZ();
      auto p = thisTrack->GetP();
      auto eta = atanh (thisTrack->GetStartVertexMomentumZ()/p); // eta;
      MCTrackspT->Fill(thisTrack->GetPt());
      MCTracksp->Fill(p);
      MCTrackRap->Fill(eta);
      MCTracksEtaZ->Fill(z,eta);



      // Count MFT disks "touched" by the track
      int nMFTDisksHasHits = 0;
      for(auto disk: {0,1,2,3,4}) nMFTDisksHasHits+= int(mcTrackHasHitsInMFTDisks[trID][disk]);
      MFTTrackablility->Fill(nMFTDisksHasHits);
      //std::cout << "nMFTDisksHasHits = " << nMFTDisksHasHits; // << std::endl;

      // Count ITS Layers "touched" by the track
      int nITSLayersHasHits = 0;
      for(auto layer: {0,1,2}) nITSLayersHasHits+= int(mcTrackHasHitsInITSLayers[trID][layer]);
      ITSTrackablility->Fill(nITSLayersHasHits);


      if(nMFTDisksHasHits>=4) {   //Track is trackable if has left hits on at least 4 disks

        MFTAccepEtaZ->Fill(z,eta);
        if(nITSLayersHasHits==3) MFTandITSLInnerBarrel->Fill(z,eta);
        if(nITSLayersHasHits>=2) MFTandTwoITSLInnerLayers->Fill(z,eta);
        if(nITSLayersHasHits>=1) MFTandOneITSLInnerLayer->Fill(z,eta);



        bool wasFound = allFoundTracksLTF[event][trID] | allFoundTracksCA[event][trID];

        if(wasFound) {
          MFTTrackspT->Fill(thisTrack->GetPt());
          MFTTracksp->Fill(p);
          MFTTrackRap->Fill(eta);

          if(allFoundTracksLTF[event][trID]) {
            MFTTrackedEtaZ->Fill(thisTrack->GetStartVertexCoordinatesZ(),eta);


          }

        }
      }
    //std::cout << " Finished Track " << trID << std::endl;

    } // end Loop on tracks
    //std::cout << "Finished event " << event << std::endl;
    } // end loop over events




// Write histograms to file
//std::cout << "Writting histograms to file..." << std::endl;

MCTrackspT->Write();
MCTracksp->Write();
MCTrackRap->Write();

MFTTrackspT->Write();
MFTTracksp->Write();
MFTTrackRap->Write();


MFTTrackablility->Write();
ITSTrackablility->Write();


MCTracksEtaZ->SetOption("CONT4");
MCTracksEtaZ->Write();

MFTAccepEtaZ->SetOption("CONT4");
MFTAccepEtaZ->Write();

MFTTrackedEtaZ->SetOption("CONT4");
MFTTrackedEtaZ->Write();



MFTandITSLInnerBarrel->SetOption("CONT4");
MFTandITSLInnerBarrel->Write();
MFTandTwoITSLInnerLayers->SetOption("CONT4");
MFTandTwoITSLInnerLayers->Write();
MFTandOneITSLInnerLayer->SetOption("CONT4");
MFTandOneITSLInnerLayer->Write();


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

}
