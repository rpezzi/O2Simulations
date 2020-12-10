#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "ITSMFTSimulation/Hit.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsMFT/TrackMFT.h"

#endif

bool DEBUG_VERBOSE = false;

void MFTTrackerChecker(Double_t pMax = 6.0,
		       Double_t pMin = 0.0,
		       Double_t etaMin = -3.8,
		       Double_t etaMax = -2.2,
		       const Char_t *kineFileName = "o2sim_Kine.root",
		       const Char_t *hitsFileName = "o2sim_HitsMFT.root",
		       const Char_t *clsFileName = "mftclusters.root",
		       const Char_t *trkFileName = "mfttracks.root")
{
  using o2::itsmft::Hit;
  using o2::itsmft::CompClusterExt;
  using o2::MCTrack;
  
  using eventFoundTracks = std::vector<bool>;
  vector<eventFoundTracks> allFoundTracksLTF, allFoundTracksCA; // True for reconstructed tracks - one vector of bool per event
  using trackHasHitsInMFTDisks = std::array<bool,5>; // Disks with hits from a MFT track

  o2::itsmft::ChipMappingMFT mftChipMapper;
 
  std::unique_ptr<TH1F> MCTrackspT = std::make_unique<TH1F> ("MC Tracks pT", "MC Tracks pT", 100, pMin, 0.2 * pMax);
  MCTrackspT->GetXaxis()->SetTitle("Transverse p");
  std::unique_ptr<TH1F> MCTracksp = std::make_unique<TH1F> ("MC Tracks p", "MC Tracks p", 100, pMin, pMax);
  MCTracksp->GetXaxis()->SetTitle("Total p");
  std::unique_ptr<TH1F> MCTrackEta = std::make_unique<TH1F> ("MC Tracks eta", "MC Tracks Pseudorapidity", 100, etaMin, etaMax);
  MCTrackEta->GetXaxis()->SetTitle("Pseudorapidity");

  std::unique_ptr<TH1F> MFTTrackspT = std::make_unique<TH1F> ("MFT Tracks pT", "MFT Tracks pT", 100, pMin, 0.2 * pMax);
  MFTTrackspT->GetXaxis()->SetTitle("Transverse p");
  std::unique_ptr<TH1F> MFTTracksp = std::make_unique<TH1F> ("MFT Tracks p", "MFT Tracks p", 100, pMin, pMax);
  MFTTracksp->GetXaxis()->SetTitle("Total p");
  std::unique_ptr<TH1F> MFTTrackEta = std::make_unique<TH1F> ("MFT Tracks eta", "MFT Tracks Pseudorapidity", 100, etaMin, etaMax);
  MFTTrackEta->GetXaxis()->SetTitle("Pseudorapidity");

  std::unique_ptr<TH1F> LTFTrackspT = std::make_unique<TH1F> ("LTF Tracks pT", "LTF Tracks pT", 100, pMin, 0.2 * pMax);
  LTFTrackspT->GetXaxis()->SetTitle("Transverse p");
  std::unique_ptr<TH1F> LTFTracksp = std::make_unique<TH1F> ("LTF Tracks p", "LTF Tracks p", 100, pMin, pMax);
  LTFTracksp->GetXaxis()->SetTitle("Total p");
  std::unique_ptr<TH1F> LTFTrackEta = std::make_unique<TH1F> ("LTF Tracks eta", "LTF Tracks Pseudorapidity", 100, etaMin, etaMax);
  LTFTrackEta->GetXaxis()->SetTitle("Pseudorapidity");

  std::unique_ptr<TH1F> CATrackspT = std::make_unique<TH1F> ("CA Tracks pT", "CA Tracks pT", 100, pMin, 0.2 * pMax);
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

  std::unique_ptr<TH1I> Trackablility = std::make_unique<TH1I> ("Trackablility", "In how many disks the tracks has hits", 6, 0, 6);
  Trackablility->GetXaxis()->SetTitle("Number of disks");

  //Histos for Missed (missed tracks that could be tracked)
  std::unique_ptr<TH1F> MissedlepT = std::make_unique<TH1F> ("Missed Tracks pT", "Missed Tracks pT", 100, pMin, 0.2 * pMax);
  std::unique_ptr<TH1F> Missedp = std::make_unique<TH1F> ("Missed Tracks p", "Missed Tracks p", 100, pMin, pMax);
  std::unique_ptr<TH1F> MissedEta = std::make_unique<TH1F> ("Missed Tracks eta", "Missed Pseudorapidity", 100, etaMin, etaMax);

  //Histos for Trackables
  std::unique_ptr<TH1F> TrackablepT = std::make_unique<TH1F> ("Trackables Tracks pT", "Trackables Tracks pT", 100, pMin, 0.2 * pMax);
  TrackablepT->GetXaxis()->SetTitle("Transverse p");
  std::unique_ptr<TH1F> Trackablep = std::make_unique<TH1F> ("Trackables Tracks p", "Trackables Tracks p", 100, pMin, pMax);
  Trackablep->GetXaxis()->SetTitle("Total p");
  std::unique_ptr<TH1F> TrackableEta = std::make_unique<TH1F> ("Trackables Tracks eta", "Trackables Pseudorapidity", 100, etaMin, etaMax);
  TrackableEta->GetXaxis()->SetTitle("Pseudorapidity");

  //2D Histos
  std::unique_ptr<TH2F> MFTTrackedEtaZ = std::make_unique<TH2F> ("MFT_Tracked_eta_z", "Reconstructed Tracks", 31, -15, 16, 25, etaMin, etaMax);
  MFTTrackedEtaZ->GetXaxis()->SetTitle("Vertex PosZ [cm]");
  std::unique_ptr<TH2F> MFTTrackablesEtaZ = std::make_unique<TH2F> ("MFT_Trackables_eta_z", "MFT Trackables:", 31, -15, 16, 25, etaMin, etaMax);
  MFTTrackablesEtaZ->GetXaxis()->SetTitle("Vertex PosZ [cm]");
  std::unique_ptr<TH2F> MCTracksEtaZ = std::make_unique<TH2F> ("MCTracks_eta_z", "MC Tracks: Pseudorapidity vs zVertex", 31, -15, 16, 25, etaMin, etaMax);
  MCTracksEtaZ->GetXaxis()->SetTitle("Vertex PosZ [cm]");
  
  // MC tracks
  TFile kineFile(kineFileName);
  TTree *kineTree = (TTree*)kineFile.Get("o2sim");
  std::vector<o2::MCTrack> mcTrkVec, *mcTrkVecP = &mcTrkVec;
  kineTree->SetBranchAddress("MCTrack",&mcTrkVecP);
  o2::dataformats::MCEventHeader* eventHeader = nullptr;
  kineTree->SetBranchAddress("MCEventHeader.", &eventHeader);

  // hits
  TFile hitsFile(hitsFileName);
  TTree* hitTree = (TTree*)hitsFile.Get("o2sim");
  std::vector<Hit> hitVec, *hitVecP = &hitVec;
  hitTree->SetBranchAddress("MFTHit", &hitVecP);

  Int_t nEvents = hitTree->GetEntries();
  printf("Number of events %d \n", nEvents);

  // clusters
  TFile clsFile(clsFileName);
  TTree *clsTree = (TTree*)clsFile.Get("o2sim");
  std::vector<CompClusterExt> clsVec, *clsVecP = &clsVec;
  clsTree->SetBranchAddress("MFTClusterComp", &clsVecP);
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* clsLabels = nullptr;
  if (clsTree->GetBranch("MFTClusterMCTruth")) {
    clsTree->SetBranchAddress("MFTClusterMCTruth", &clsLabels);
  } else {
    printf("No Monte-Carlo information in this file\n");
    return;
  }
  
  clsTree->GetEntry(0);
  Int_t nClusters = clsVec.size();
  printf("Number of clusters %d \n", nClusters);
  
  // tracks
  TFile trkFile(trkFileName);
  TTree *trkTree = (TTree*)trkFile.Get("o2sim");
  std::vector<o2::mft::TrackMFT> trackVec, *trackVecP = &trackVec;
  trkTree->SetBranchAddress("MFTTrack", &trackVecP);
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* trkLabels = nullptr;
  if (trkTree->GetBranch("MFTTrackMCTruth")) {
    trkTree->SetBranchAddress("MFTTrackMCTruth", &trkLabels);
  } else {
    printf("No Monte-Carlo information in this file\n");
    return;
  }
  std::vector<int> trackExtClsVec, *trackExtClsVecP = &trackExtClsVec;
  trkTree->SetBranchAddress("MFTTrackClusIdx", &trackExtClsVecP);

  trkTree->GetEntry(0);

  // MFTTrackerChecker.C
  // Resize vector to accomodate found status of all tracks in all events
  for (auto event = 0 ; event < nEvents ; event++) { 
    kineTree->GetEntry(event);
    hitTree->GetEntry(event);
    auto nTracks = eventHeader->getMCEventStats().getNKeptTracks();
    if(DEBUG_VERBOSE) std::cout << "Resizing allFoundTracks for event " << event <<  " with ntracks = " << nTracks << std::endl;
    eventFoundTracks tempFoundTracks(nTracks, false);
    allFoundTracksLTF.push_back(tempFoundTracks); // reserve size and initialize
    allFoundTracksCA.push_back(tempFoundTracks); // reserve size and initialize
  }

  // MFTTrackerChecker.C
  // Part 1: Quality of reconstructed MFT tracks
  //   - Loop over reconstructed tracks to identify clean and mixed/noise tracks
  //   - Clean tracks have at least 80% of its clusters from the same track
  //   - If track is not clear it is a mixed/noise track
  // Part 2: MC hits and tracks
  //   - Identify trackable tracks (clusters in at least 4 disks)
  //   - Identify successfully reconstructed tracks
  // Part 3: Calculate Efficiencies

  // Part 1: Quality of reconstructed MFT tracks

  Int_t nCleanTracksLTF = 0, nCleanTracksCA = 0, nInvalidTracksLTF = 0, nInvalidTracksCA = 0, nMFTTrackable = 0;

  int srcID, trkID, evnID;
  bool fake;
  Int_t iTrack = 0;
  for (auto &track : trackVec) {
    
    auto ncls = track.getNumberOfPoints();
    auto offset = track.getExternalClusterIndexOffset();
    std::map<Int_t, Int_t> trkIDs;
    for (int icls = 0; icls < ncls; ++icls) {
      auto clsEntry = trackExtClsVec[offset + icls];
      auto cluster = clsVec[clsEntry];
      auto& clsLabel = (clsLabels->getLabels(clsEntry))[0];
      clsLabel.get(trkID, evnID, srcID, fake);
      if (!clsLabel.isNoise()) {
	trkIDs[trkID] = trkIDs[trkID] + 1;
      }
    }
    
    Int_t thisEvnID = -1, thisSrcID = -1, thisTrkID = -1, thisEventIDLabel = -1, nPoints = track.getNumberOfPoints();
    for (int icls = 0; icls < ncls; ++icls) {
      auto clsEntry = trackExtClsVec[offset + icls];
      auto cluster = clsVec[clsEntry];
      auto& clsLabel = (clsLabels->getLabels(clsEntry))[0];
      clsLabel.get(trkID, evnID, srcID, fake);
      if (!clsLabel.isNoise()) {
	if (((Float_t)(trkIDs[trkID]) / (Float_t)(nPoints)) >= 0.8) { // Must have at least 80% of its clusters from the same MC Track
	  thisTrkID = trkID;
	  thisSrcID = srcID;
	  thisEvnID = evnID;
	  thisEventIDLabel = icls;
	}
      }
    }
    
    auto eventID = thisEvnID;
    
    if ((thisTrkID >= 0) & (thisTrkID != 0x7FFFFFF) & (eventID < nEvents)) { // If is good match and not noise ...
      if (!track.isCA()) {
	allFoundTracksLTF[eventID][thisTrkID] = true;
	nCleanTracksLTF++;
      } else {
	allFoundTracksCA[eventID][thisTrkID] = true;
	nCleanTracksCA++;
      }
    } else {
      if (!track.isCA()) {
	nInvalidTracksLTF++;
      } else {
	nInvalidTracksCA++;
      }
    }
    if(DEBUG_VERBOSE) std::cout << "This TrackLTF ID = " << thisTrkID << " from sourceID = " << thisSrcID << " in eventID = " << eventID << " nPoints = " << nPoints << " nCleanTracksLTF = " << nCleanTracksLTF << std::endl;
   
    ++iTrack;
  }

  // Part 2: MC hits and tracks

  for (Int_t event = 0; event < nEvents ; event++) {

    hitTree->GetEntry(event);
    kineTree->GetEntry(event);

    Int_t nMFTHits = hitVec.size();
    std::vector<trackHasHitsInMFTDisks> mcTrackHasHitsInMFTDisks(eventHeader->getMCEventStats().getNKeptTracks(), {0, 0, 0, 0, 0});

    if(DEBUG_VERBOSE) std::cout << "Loop over " << nMFTHits << " mfthits to identify trackable MFT tracks in event " <<  event << std::endl;

    for (Int_t n_hit = 0 ; n_hit < nMFTHits; n_hit++) { // Loop over mfthits to identify trackable tracks
      
      Hit* hitp = &(hitVec).at(n_hit);
      Int_t trkID = hitp->GetTrackID(); // ID of the tracks having given the hit
      Float_t z = hitp->GetZ(); // Z position of the hit => Identify MFT disk
      mcTrackHasHitsInMFTDisks.at(trkID)[mftChipMapper.chip2Layer(hitp->GetDetectorID()) / 2] = true;
    }
    
    for (Int_t trkID = 0 ; trkID < eventHeader->getMCEventStats().getNKeptTracks(); trkID++) { // Loop on MC tracks

      //fill MC histograms
      MCTrack* thisTrack =  &(mcTrkVec)[trkID];
      if (!thisTrack->isPrimary()) {
	continue;
      }
      
      auto z = thisTrack->GetStartVertexCoordinatesZ();
      auto pt = thisTrack->GetPt();
      auto p = thisTrack->GetP();
      auto eta = atanh (thisTrack->GetStartVertexMomentumZ()/p); // eta;
      MCTrackspT->Fill(pt);
      MCTracksp->Fill(p);
      MCTrackEta->Fill(eta);
      MCTracksEtaZ->Fill(z,eta);
      if( (z > -5) & (z < 5) ) {
        MCTracksEta5->Fill(eta);
        MCTracksp5->Fill(p);
      }
      if( (z > 5) & (z < 10) ) {
        MCTracksEta5_10pos->Fill(eta);
        MCTracksp5_10pos->Fill(p);
      }
      if( (z > -10) & (z < -5) ) {
        MCTracksEta5_10neg->Fill(eta);
        MCTracksp5_10neg->Fill(p);
      }
      // Count disks "touched" by the track
      int nMFTDisksHasHits = 0;
      for(auto disk: {0, 1, 2, 3, 4}) {
	nMFTDisksHasHits+= int(mcTrackHasHitsInMFTDisks[trkID][disk]);
      }
      Trackablility->Fill(nMFTDisksHasHits);
      //std::cout << "nMFTDisksHasHits = " << nMFTDisksHasHits; // << std::endl;
      if (nMFTDisksHasHits >= 4) {   //Track is trackable if has left hits on at least 4 disks
        nMFTTrackable++;
        MFTTrackablesEtaZ->Fill(z,eta);
        TrackablepT->Fill(pt);
        Trackablep->Fill(p);
        TrackableEta->Fill(eta);
        bool wasFound = allFoundTracksLTF[event][trkID] | allFoundTracksCA[event][trkID];
        if(wasFound) {	  
          MFTTrackspT->Fill(pt);
          MFTTracksp->Fill(p);
          MFTTrackEta->Fill(eta);
          if( (z > -5) & (z < 5) ) {
            MFTTracksEta5->Fill(eta);
            MFTTracksp5->Fill(p);
          }
          if( (z > 5) & (z < 10) ) {
            MFTTracksEta5_10pos->Fill(eta);
            MFTTracksp5_10pos->Fill(p);
          }
          if( (z > -10) & (z < -5) ) {
            MFTTracksEta5_10neg->Fill(eta);
            MFTTracksp5_10neg->Fill(p);
          }
          if(allFoundTracksLTF[event][trkID]) {
            LTFTrackspT->Fill(thisTrack->GetPt());
            LTFTracksp->Fill(thisTrack->GetP());
            LTFTrackEta->Fill(eta);
            MFTTrackedEtaZ->Fill(thisTrack->GetStartVertexCoordinatesZ(),eta);
          }
          if(allFoundTracksCA[event][trkID]) {
            CATrackspT->Fill(thisTrack->GetPt());
            CATracksp->Fill(thisTrack->GetP());
            CATrackEta->Fill(eta);
          }
        }
      } else {  // Fill histograms for Missed Tracks
        MissedlepT->Fill(thisTrack->GetPt());
        Missedp->Fill(thisTrack->GetP());
        MissedEta->Fill(eta);
      }
    } // end loop on tracks   
  } // end loop on events

  // Part 3: Calculate Efficiencies
  //std::cout << "Building efficiencies histos..." << std::endl;
  TH1F MFTEfficiencypT = (*MFTTrackspT)/ (*TrackablepT);
  TH1F MFTTEfficiencyp = (*MFTTracksp) / (*Trackablep);
  TH1F MFTEfficiencyEta = (*MFTTrackEta) / (*TrackableEta);
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
  //std::cout << "Writting histograms to file..." << std::endl;
  TFile outFile("MFTTrackerCheck.root","RECREATE");
  
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
  
  TrackablepT->Write();
  Trackablep->Write();
  TrackableEta->Write();
  
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
  std::cout << "---------------------------------------------------" << std::endl;
  
  std::cout << "Number of MFT trackables = " << nMFTTrackable << std::endl;
  std::cout << "Number of reconstructed MFT Tracks = " << totalRecoMFTTracks << std::endl;
  std::cout << "Number of clean MFT Tracks = " << nCleanTracksLTF + nCleanTracksCA << std::endl;
  std::cout << "Number of mixed MFT Tracks = " << nInvalidTracksLTF + nInvalidTracksCA << std::endl;
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "nCleanTracksLTF = " << nCleanTracksLTF << std::endl;
  std::cout << "nCleanTracksCA = " << nCleanTracksCA << std::endl;
  std::cout << "nInvalidTracksLTF = " << nInvalidTracksLTF  << " (" << 100.f*nInvalidTracksLTF/(nCleanTracksLTF+nInvalidTracksLTF) << " %)" << std::endl;
  std::cout << "nInvalidTracksCA = " << nInvalidTracksCA << " (" << 100.f*nInvalidTracksCA/(nCleanTracksCA+nInvalidTracksCA) << " %)" << std::endl;
  std::cout << "---------------------------------------------------" << std::endl;

}
