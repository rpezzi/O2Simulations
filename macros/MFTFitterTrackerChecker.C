#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <TMath.h>


//constexpr Double_t MFTLayerZ[] = {-45.3, -46.7, -48.6, -50.0, -52.4, -53.8, -67.7, -69.1, -76.1, -77.5};
using o2::itsmft::Hit;
using o2::MCTrackT;
using o2::mft::TrackMFT;
o2::itsmft::ChipMappingMFT mftChipMapper;
using eventFoundTracks = std::vector<bool>;
vector<eventFoundTracks> allFoundTracksMFT; // True for reconstructed tracks - one vector of bool per event

using trackHasHitsinMFTDisks = std::array<bool,5>; // Disks with hits from a MFT track

bool DEBUG_VERBOSE = false;

bool DISABLE_PART2 = true;


//_________________________________________________________________________________________________
void extrapMFTTrackHelixToZ(o2::mft::TrackMFT& track, double zEnd, double Field)
{
   using TrackMFT = o2::mft::TrackMFT;

  /// Track extrapolated to the plane at "zEnd" considering a helix
  /// On return, results from the extrapolation are returned as a new TrackMFT.

  if (track.getZ() == zEnd) {
    return; // nothing to be done if same z
  }

  // Compute track parameters
  double dZ = (zEnd - track.getZ()); // Propagate in meters
  double x0 = track.getX();
  double y0 = track.getY();
  double phi0 = track.getPhi();
  double cosphi0 = TMath::Cos(phi0);
  double sinphi0 = TMath::Sin(phi0);
  double tanl0 = track.getTanl();
  double invqpt0 = track.getInvQPt();

  double k = - Field * 0.299792458e-3;
  double deltax = (dZ * cosphi0 / tanl0 - dZ * dZ * k * invqpt0 * sinphi0 / (2. * tanl0 * tanl0));
  double deltay = (dZ * sinphi0 / tanl0 + dZ * dZ * k * invqpt0 * cosphi0 / (2. * tanl0 * tanl0));

  double x = x0 + deltax;
  double y = y0 + deltay;
  double deltaphi = +dZ * k * invqpt0 / tanl0;

  float phi = phi0 + deltaphi;
  double tanl = tanl0;
  double invqpt = invqpt0;
  track.setX(x);
  track.setY(y);
  track.setZ(zEnd);
  track.setPhi(phi);
  track.setTanl(tanl);
  track.setInvQPt(invqpt);
}

//_________________________________________________________________________________________________
void MFTFitterTrackerChecker( Double_t pMax = 40.0,
                              Double_t pMaxFit = 1000.0,
                              Double_t pMin = 0.0,
                              Double_t etaMin = -.2,//-4.0,
                              Double_t etaMax = +.2,// +4.0,
                              Double_t phiMin = -.2, //-3.15,
                              Double_t phiMax = .2, //+3.15,
                              const Char_t *o2sim_KineFile = "o2sim_Kine.root",
                              const Char_t *HitsMFTFile = "o2sim_HitsMFT.root",
                              const Char_t *trkFile = "mfttracks.root") {

  std::unique_ptr<TH1F> MCTrackspT = std::make_unique<TH1F> ("MC Tracks pT", "MC Tracks pT", 100, pMin, pMax);
  MCTrackspT->GetXaxis()->SetTitle("Transverse p");
  std::unique_ptr<TH1F> MCTracksp = std::make_unique<TH1F> ("MC Tracks p", "MC Tracks p", 100, pMin, pMax);
  MCTracksp->GetXaxis()->SetTitle("Total p");
  std::unique_ptr<TH1F> MCTrackEta = std::make_unique<TH1F> ("MC Tracks eta", "MC Tracks Pseudorapidity", 100, etaMin, etaMax);
  MCTrackEta->GetXaxis()->SetTitle("\\eta ");

  std::unique_ptr<TH1F> MFTTracksMCpT = std::make_unique<TH1F> ("MFT Tracks MC pT", "MFT Tracks MC pT", 100, pMin, pMax);
  MFTTracksMCpT->GetXaxis()->SetTitle("Transverse p");
  std::unique_ptr<TH1F> MFTTracksMCp = std::make_unique<TH1F> ("MFT Tracks MC p", "MFT Tracks MC p", 100, pMin, pMax);
  MFTTracksMCp->GetXaxis()->SetTitle("Total p");
  std::unique_ptr<TH1F> MFTTrackMCEta = std::make_unique<TH1F> ("MFT Tracks MC eta", "MFT Tracks MC Pseudorapidity", 100, etaMin, etaMax);
  MFTTrackMCEta->GetXaxis()->SetTitle("\\eta ");


  std::unique_ptr<TH1F> MCTracksEta5 = std::make_unique<TH1F> ("MC Tracks 5 MC eta", "-5 cm < zVertex < 5 cm", 100, etaMin, etaMax);
  MCTracksEta5->GetXaxis()->SetTitle("\\eta ");
  std::unique_ptr<TH1F> MCTracksEta5_10pos = std::make_unique<TH1F> ("MC Tracks -5 -10 MC eta", "-10 cm < zVertex < -5 cm", 100, etaMin, etaMax);
  MCTracksEta5_10pos->GetXaxis()->SetTitle("\\eta ");
  std::unique_ptr<TH1F> MCTracksEta5_10neg = std::make_unique<TH1F> ("MC Tracks 5 10 MC eta", "5 cm < zVertex < 10 cm", 100, etaMin, etaMax);
  MCTracksEta5_10neg->GetXaxis()->SetTitle("\\eta ");

  std::unique_ptr<TH1F> MFTTracksEta5 = std::make_unique<TH1F> ("MFT Tracks 5 MC MC eta", "-5 cm < zVertex < 5 cm", 100, etaMin, etaMax);
  MFTTracksEta5->GetXaxis()->SetTitle("\\eta ");
  std::unique_ptr<TH1F> MFTTracksEta5_10pos = std::make_unique<TH1F> ("MFT Tracks -5 -10 MC eta", "-10 cm < zVertex < -5 cm", 100, etaMin, etaMax);
  MFTTracksEta5_10pos->GetXaxis()->SetTitle("\\eta ");
  std::unique_ptr<TH1F> MFTTracksEta5_10neg = std::make_unique<TH1F> ("MFT Tracks 5 10 MC eta", "5 cm < zVertex < 10 cm", 100, etaMin, etaMax);
  MFTTracksEta5_10neg->GetXaxis()->SetTitle("\\eta ");


  std::unique_ptr<TH1F> MCTracksp5 = std::make_unique<TH1F> ("MC Tracks 5 p", "-5 cm < zVertex < 5 cm", 100, pMin, pMax);
  MCTracksp5->GetXaxis()->SetTitle("Total p");
  std::unique_ptr<TH1F> MCTracksp5_10pos = std::make_unique<TH1F> ("MC Tracks -5 -10 p", "-10 cm < zVertex < -5 cm", 100, pMin, pMax);
  MCTracksp5_10pos->GetXaxis()->SetTitle("Total p");
  std::unique_ptr<TH1F> MCTracksp5_10neg = std::make_unique<TH1F> ("MC Tracks 5 10 p", "5 cm < zVertex < 10 cm", 100, pMin, pMax);
  MCTracksp5_10neg->GetXaxis()->SetTitle("Total p");

  std::unique_ptr<TH1F> MFTTracksMCp5 = std::make_unique<TH1F> ("MFT Tracks 5 MC p", "-5 cm < zVertex < 5 cm", 100, pMin, pMax);
  MFTTracksMCp5->GetXaxis()->SetTitle("Total p");
  std::unique_ptr<TH1F> MFTTracksMCp5_10pos = std::make_unique<TH1F> ("MFT Tracks -5 -10 MC p", "-10 cm < zVertex < -5 cm", 100, pMin, pMax);
  MFTTracksMCp5_10pos->GetXaxis()->SetTitle("Total p");
  std::unique_ptr<TH1F> MFTTracksMCp5_10neg = std::make_unique<TH1F> ("MFT Tracks 5 10 MC p", "5 cm < zVertex < 10 cm", 100, pMin, pMax);
  MFTTracksMCp5_10neg->GetXaxis()->SetTitle("Total p");

  std::unique_ptr<TH1I> Trackablility = std::make_unique<TH1I> ("Trackablility", "In how many disks the tracks has hits", 6, 0, 6);
  Trackablility->GetXaxis()->SetTitle("Number of disks");

  //Histos for Missed (missed tracks that could be tracked)
  std::unique_ptr<TH1F> MissedlepT = std::make_unique<TH1F> ("Missed Tracks MC pT", "Missed Tracks pT", 100, pMin, pMax);
  std::unique_ptr<TH1F> Missedp = std::make_unique<TH1F> ("Missed Tracks MC p", "Missed Tracks p", 100, pMin, pMax);
  std::unique_ptr<TH1F> MissedEta = std::make_unique<TH1F> ("Missed Tracks MC eta", "Missed Pseudorapidity", 100, etaMin, etaMax);

  //Histos for Trackables
  std::unique_ptr<TH1F> TrackablepT = std::make_unique<TH1F> ("Trackables Tracks MC p_T", "Trackables Tracks p_T", 100, pMin, pMax);
  TrackablepT->GetXaxis()->SetTitle("Transverse p");
  std::unique_ptr<TH1F> Trackablep = std::make_unique<TH1F> ("Trackables Tracks MC p", "Trackables Tracks p", 100, pMin, pMax);
  Trackablep->GetXaxis()->SetTitle("Total p");
  std::unique_ptr<TH1F> TrackableEta = std::make_unique<TH1F> ("Trackables Tracks MC eta", "Trackables Pseudorapidity", 100, etaMin, etaMax);
  TrackableEta->GetXaxis()->SetTitle("\\eta ");

  //2D Histos
  std::unique_ptr<TH2F> MFTTrackedEtaZ = std::make_unique<TH2F> ("MFT_Tracked_MC_eta_z", "Reconstructed Tracks", 31, -15, 16, 25, etaMin, etaMax);
  MFTTrackedEtaZ->GetXaxis()->SetTitle("Vertex PosZ ~[cm]");
  std::unique_ptr<TH2F> MFTTrackablesEtaZ = std::make_unique<TH2F> ("MFT_Trackables_MC_eta_z", "MFT Trackables:", 31, -15, 16, 25, etaMin, etaMax);
  MFTTrackablesEtaZ->GetXaxis()->SetTitle("Vertex PosZ ~[cm]");
  std::unique_ptr<TH2F> MCTracksEtaZ = std::make_unique<TH2F> ("MCTracks_eta_z", "MC Tracks: Pseudorapidity vs zVertex", 31, -15, 16, 25, etaMin, etaMax);
  MCTracksEtaZ->GetXaxis()->SetTitle("Vertex PosZ ~[cm]");


  // MFT Fit Results histos
  std::unique_ptr<TH1F> MFTTracksP = std::make_unique<TH1F> ("MFT Tracks Fitted p", "Standalone MFT Tracks P", 10000, pMin, pMaxFit);
  MFTTracksP->GetXaxis()->SetTitle("p ~[GeV]");
  std::unique_ptr<TH1F> MFTTracksDeltaP = std::make_unique<TH1F> ("MFT Tracks Delta_p", "P_{Fit} - P_{MC}", 10000, -pMaxFit, pMaxFit);
  MFTTracksDeltaP->GetXaxis()->SetTitle("\\Delta p ~[GeV]");
  std::unique_ptr<TH1F> MFTTracksDeltaPt = std::make_unique<TH1F> ("MFT Tracks Delta_pt", "Pt_{Fit} - Pt_{MC}", 10000, -pMaxFit, pMaxFit);
  MFTTracksDeltaPt->GetXaxis()->SetTitle("\\Delta p_t ~[GeV]");
  std::unique_ptr<TH1F> MFTTrackDeltaEta = std::make_unique<TH1F> ("MFT Tracks Fitted Delta_eta", "\\eta_{Fit} - \\eta_{MC} ", 1000, -.1, +.1);
  MFTTrackDeltaEta->GetXaxis()->SetTitle("\\Delta \\eta");
  std::unique_ptr<TH1F> MFTTrackDeltaPhi = std::make_unique<TH1F> ("MFT Tracks Fitted Phi at Vertex", "\\phi _{Fit} - \\phi_{MC}" , 1000, phiMin, phiMax);
  MFTTrackDeltaPhi->GetXaxis()->SetTitle("\\Delta \\phi ~[rad]");
  std::unique_ptr<TH1F> MFTTrackDeltaPhiDeg = std::make_unique<TH1F> ("MFT Tracks Fitted Phi at Vertex [deg]", "\\phi _{Fit} - \\phi_{MC}" , 1000, TMath::RadToDeg()*phiMin, TMath::RadToDeg()*phiMax);
  MFTTrackDeltaPhiDeg->GetXaxis()->SetTitle("\\Delta \\phi ~[deg]");
  std::unique_ptr<TH1F> MFTTrackDeltaX = std::make_unique<TH1F> ("MFT Tracks Delta X", "Standalone MFT Tracks Delta X at Z_vertex", 1000, -.3, .3);
  MFTTrackDeltaX->GetXaxis()->SetTitle("\\Delta x ~[cm]");
  std::unique_ptr<TH1F> MFTTrackDeltaY = std::make_unique<TH1F> ("MFT Tracks Delta Y", "Standalone MFT Tracks Delta Y at Z_vertex", 1000, -.3, .3);
  MFTTrackDeltaY->GetXaxis()->SetTitle("\\Delta y ~[cm]");
  std::unique_ptr<TH1F> MFTTrackR = std::make_unique<TH1F> ("MFT Tracks Delta R", "Standalone MFT Tracks Delta R at Z_vertex", 10000, -2, +2);
  MFTTrackR->GetXaxis()->SetTitle("\\Delta r ~[cm]");
  std::unique_ptr<TH2F> MFTTrackDeltaXYVertex = std::make_unique<TH2F> ("MFT Tracks Vertex at Z_Vertex = 0", "Standalone MFT Tracks at Z_vertex", 250, -.05, .05,1000, -.05, .05);
  MFTTrackDeltaXYVertex->GetXaxis()->SetTitle("\\Delta x ~[cm]");
  MFTTrackDeltaXYVertex->GetYaxis()->SetTitle("\\Delta y ~[cm]");
  std::unique_ptr<TH1F> MFTTrackQ = std::make_unique<TH1F> ("MFT Tracks Q", "Standalone MFT Tracks Charge", 5, -2.1, 2.1);
  MFTTrackQ->GetXaxis()->SetTitle("Q");




  // MC
  TFile *o2sim_KineFileIn = new TFile(o2sim_KineFile);
  TTree *o2SimKineTree = (TTree*) o2sim_KineFileIn -> Get("o2sim");

  vector<MCTrackT<float>>* mcTr = nullptr;
  o2SimKineTree -> SetBranchAddress("MCTrack",&mcTr);
  o2::dataformats::MCEventHeader* eventHeader = nullptr;
  o2SimKineTree -> SetBranchAddress("MCEventHeader.",&eventHeader);

  // MFT Hits
  TFile *HitsMFTFileIn = new TFile(HitsMFTFile);
  TTree *o2MFTHitsTree = (TTree*) HitsMFTFileIn -> Get("o2sim");
  vector<Hit>* mfthit = nullptr;
  o2MFTHitsTree -> SetBranchAddress("MFTHit",&mfthit);

  Int_t numberOfEvents = o2SimKineTree -> GetEntries();
  Int_t numberOfMFTEvents = o2MFTHitsTree -> GetEntries();
  if (numberOfEvents == numberOfMFTEvents )
   std::cout << "numberOfEvents = " << numberOfEvents << std::endl;
  else {
    std::cout << "ERROR: Inconsistent number of entries on " << o2sim_KineFile << " and " << HitsMFTFile << std::endl;
    return -1;
  }

  // MFT Tracks
  TFile *trkFileIn = new TFile(trkFile);
  TTree *mftTrackTree = (TTree*) trkFileIn -> Get("o2sim");
  std::vector<o2::mft::TrackMFT> trackMFTVec, *trackMFTVecP = &trackMFTVec;
  mftTrackTree->SetBranchAddress("MFTTrack", &trackMFTVecP);

  Int_t nCleanTracksMFT = 0, nInvalidTracksMFT = 0, nMFTTrackable = 0;

  mftTrackTree -> GetEntry(0);
  o2SimKineTree -> GetEntry(0);
  o2MFTHitsTree -> GetEntry(0);

  for (auto event = 0 ; event < numberOfEvents ; event++) { // Resize vector to accomodate found status of all tracks in all events
    o2SimKineTree -> GetEntry(event);
    auto numberOfTracksThisEvent = eventHeader->getMCEventStats().getNKeptTracks();
    if(DEBUG_VERBOSE) std::cout << "Resizing allFoundTracks for event " << event <<  " with ntracks = " << numberOfTracksThisEvent << std::endl;
    eventFoundTracks tempFoundTracks(numberOfTracksThisEvent,false);
    allFoundTracksMFT.push_back(tempFoundTracks);// [event].resize(numberOfTracksThisEvent,false);
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
   std::cout << "Starting Part 1: Reconstructed MFT Tracks!" << std::endl;
  // TracksMFT - Identify reconstructed tracks
  for (auto &trackMFT : trackMFTVec) {
  auto thisTrackMCCompLabels = trackMFT.getMCCompLabels();
  std::map<Int_t, Int_t> trkIDs;
  for (auto iLabel = 0; iLabel < trackMFT.getNPoints(); iLabel++) { // Count trackIDs of this track
    auto id = thisTrackMCCompLabels[iLabel].getTrackID();
    trkIDs[id]=trkIDs[id]+1;
    //std::cout << thisTrackMCCompLabels[iLabel].getTrackID() << " ";
  }
  //std::cout << std::endl;

  Int_t thisTrkID = -1, thisTEventIDLabel = -1;
  //Decide if track was successfully tracked
  for (auto iLabel = 0; iLabel < trackMFT.getNPoints(); iLabel++) {
    auto id = thisTrackMCCompLabels[iLabel].getTrackID();
    //std::cout << "ID ; count ; NPoints ; count/NPoints : "<< id << " " << trkIDs[id] << " " << trackMFT.getNPoints()  << " " << 1.0*trkIDs[id]/trackMFT.getNPoints() << std::endl;
    if(1.0 * trkIDs[id] / trackMFT.getNPoints() >= 0.8 ) {
      thisTrkID = id;
      thisTEventIDLabel = iLabel;
    }

  }

  auto eventID =  thisTrackMCCompLabels[thisTEventIDLabel].getEventID();
  //std::cout << "This TrackMFT ID = " << thisTrkID << " in eventID = " << eventID << std::endl;
  if( (thisTrkID >= 0 & thisTrkID != 0x7FFFFFF) & (eventID <= numberOfEvents) ) { // If is good match and not noise ...
   allFoundTracksMFT[eventID][thisTrkID]=true;
   nCleanTracksMFT++;
   }
  else {
    //std::cout << "Noise or Mixed TrackMFT!" << std::endl;
    nInvalidTracksMFT++;
  }

  if( (thisTrkID >= 0 & thisTrkID != 0x7FFFFFF) & (eventID <= numberOfEvents) )
      if(allFoundTracksMFT[eventID][thisTrkID] == true) {
        o2SimKineTree -> GetEntry(eventID);
        MCTrackT<float>* thisTrack =  &(*mcTr).at(thisTrkID);
        auto vx_MC = thisTrack->GetStartVertexCoordinatesX();
        auto vy_MC = thisTrack->GetStartVertexCoordinatesY();
        auto vz_MC = thisTrack->GetStartVertexCoordinatesZ();
        auto pt_MC = thisTrack->GetPt();
        auto p_MC = thisTrack->GetP();
        auto phi_MC = TMath::ATan2(thisTrack->Py(),thisTrack->Px());
        auto pdgcode_MC = thisTrack->GetPdgCode();
        //std::cout << "pdgcode_MC = " <<  pdgcode_MC;
        int q_MC;
        if (TDatabasePDG::Instance()->GetParticle(pdgcode_MC)) {
         q_MC = TDatabasePDG::Instance()->GetParticle(pdgcode_MC)->Charge()/3;
         //std::cout << " => " <<  TDatabasePDG::Instance()->GetParticle(pdgcode_MC)->GetName() << " ; q = " << q_MC <<  "\n";
         }

         else {
           q_MC = 0;
           std::cout << " => pdgcode ERROR " << q_MC <<  "\n";
         }
        auto eta_MC = atanh (thisTrack->GetStartVertexMomentumZ()/p_MC); // eta;
        extrapMFTTrackHelixToZ(trackMFT, vz_MC, -5.0); // propagate track to vertex Z
        auto dx = trackMFT.getX() - vx_MC;
        auto dy = trackMFT.getY() - vy_MC;
        auto d_eta = trackMFT.getEta() - eta_MC;
        auto d_P = trackMFT.getP() - p_MC;
        auto d_Pt = trackMFT.getPt() - pt_MC;
        auto d_Phi = trackMFT.getPhi() - phi_MC;
        MFTTracksP->Fill(trackMFT.getP());
        MFTTracksDeltaP->Fill(d_P);
        MFTTracksDeltaPt->Fill(d_Pt);
        MFTTrackDeltaEta->Fill(d_eta);
        MFTTrackDeltaPhi->Fill(d_Phi);
        MFTTrackDeltaPhiDeg->Fill(TMath::RadToDeg()*d_Phi);
        MFTTrackDeltaX->Fill(dx);
        MFTTrackDeltaY->Fill(dy);
        MFTTrackDeltaXYVertex->Fill(dx,dy);
        MFTTrackR->Fill(sqrt(dx*dx+dy*dy));
        MFTTrackQ->Fill(trackMFT.getCharge()-q_MC);
        }
  } // Loop on TracksMFT

  TFile outFile("MFTFitterTrackerCheck.root","RECREATE");

  // Part 2: MC hits and tracks
  if(!DISABLE_PART2)  {
  std::cout << "Starting Part 2: MC hits and tracks!" << std::endl;
  for (Int_t event=0; event<numberOfEvents ; event++) { // Loop over events in o2sim
    //std::cout << "Loop over events in o2sim. Event = " << event << std::endl;
    //if(event % 100 == 0)  std::cout << event << " ...\n";
    o2MFTHitsTree -> GetEntry(event);
    o2SimKineTree -> GetEntry(event);
    Int_t nMFTHits = mfthit->size(); // Number of mft hits in this event
    //std::cout << "Event " << event << " has " << eventHeader->getMCEventStats().getNKeptTracks() << " tracks and " << nMFTHits << " hits\n";

    std::vector<trackHasHitsinMFTDisks> mcTrackHasHitsInMFTDisks(eventHeader->getMCEventStats().getNKeptTracks(),{0,0,0,0,0}); //

    if(DEBUG_VERBOSE) std::cout << "Loop over " << nMFTHits << " mfthits to identify trackable MFT tracks in event " <<  event << std::endl;
    for (Int_t n_hit=0 ; n_hit < nMFTHits; n_hit++) { // Loop over mfthits to identify trackable tracks
      Hit* hitp = &(*mfthit).at(n_hit);
      Int_t trID = hitp->GetTrackID(); // ID of the tracks having given the hit
      //std::cout << "n_hit = " << n_hit << " ** trID = " << trID << std::endl;

      //Float_t z = hitp->GetZ(); // Z position of the hit => Identify MFT disk
      mcTrackHasHitsInMFTDisks[trID][mftChipMapper.chip2Layer(hitp->GetDetectorID())/2] = true;
      }

    for (Int_t trID=0 ; trID < eventHeader->getMCEventStats().getNKeptTracks(); trID++) { // Loop on MC tracks
      //std::cout << "Loop on tracks to build histos. Track " << trID << " at event " << event << " -> " ;

      //fill MC histograms
      MCTrackT<float>* thisTrack =  &(*mcTr).at(trID);
      auto z = thisTrack->GetStartVertexCoordinatesZ();
      auto pt = thisTrack->GetPt();
      auto p = thisTrack->GetP();
      auto eta = atanh (thisTrack->GetStartVertexMomentumZ()/p); // eta;
      MCTrackspT->Fill(pt);
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

        nMFTTrackable++;
        MFTTrackablesEtaZ->Fill(z,eta);
        TrackablepT->Fill(pt);
        Trackablep->Fill(p);
        TrackableEta->Fill(eta);


        bool wasFound = allFoundTracksMFT[event][trID];

        if(wasFound) {
          MFTTracksMCpT->Fill(pt);
          MFTTracksMCp->Fill(p);
          MFTTrackMCEta->Fill(eta);


          if( (z>-5) & (z<5) ) {
            MFTTracksEta5->Fill(eta);
            MFTTracksMCp5->Fill(p);
          }
          if( (z>5) & (z<10) ) {
            MFTTracksEta5_10pos->Fill(eta);
            MFTTracksMCp5_10pos->Fill(p);
          }
          if( (z>-10) & (z<-5) ) {
            MFTTracksEta5_10neg->Fill(eta);
            MFTTracksMCp5_10neg->Fill(p);
          }


          if(allFoundTracksMFT[event][trID]) {
            MFTTracksMCpT->Fill(thisTrack->GetPt());
            MFTTracksMCp->Fill(thisTrack->GetP());
            MFTTrackMCEta->Fill(eta);
          }
        }
      } else {  // Fill histograms for Missed Tracks
        MissedlepT->Fill(thisTrack->GetPt());
        Missedp->Fill(thisTrack->GetP());
        MissedEta->Fill(eta);
      }
    //std::cout << " Finished Track " << trID << std::endl;

    } // end Loop on tracks
    //std::cout << "Finished event " << event << std::endl;
    } // end loop over events
    // Part 3: Calculate Efficiencies
    //std::cout << "Building efficiencies histos..." << std::endl;
    TH1F MFTEfficiencypT = (*MFTTracksMCpT)/ (*TrackablepT);
    TH1F MFTTEfficiencyp = (*MFTTracksMCp) / (*Trackablep);
    TH1F MFTEfficiencyEta = (*MFTTrackMCEta) / (*TrackableEta);
    TH2F MFTTrackerEfficiency = (*MFTTrackedEtaZ) / (*MFTTrackablesEtaZ);
    TH2F MFTEfficiency2D = (*MFTTrackedEtaZ) / (*MCTracksEtaZ);
    TH2F MFTAcceptance = (*MFTTrackablesEtaZ) / (*MCTracksEtaZ);
    TH1F MFTEffsEta5 = (*MFTTracksEta5)/(*MCTracksEta5);
    TH1F MFTEffEta5_10pos = (*MFTTracksEta5_10pos)/(*MCTracksEta5_10pos);
    TH1F MFTEffEta5_10neg = (*MFTTracksEta5_10neg)/(*MCTracksEta5_10neg);


    TH1F MFTEffsp5 = (*MFTTracksMCp5)/(*MCTracksp5);
    TH1F MFTEffp5_10pos = (*MFTTracksMCp5_10pos)/(*MCTracksp5_10pos);
    TH1F MFTEffp5_10neg = (*MFTTracksMCp5_10neg)/(*MCTracksp5_10neg);

    // Stacks
    THStack mftEtaStack("PStack","MFT Tracks");
    MFTEffsEta5.SetLineColor(kBlack);
    mftEtaStack.Add(&MFTEffsEta5,"nostack");
    MFTEffEta5_10pos.SetLineColor(kBlue);
    mftEtaStack.Add(&MFTEffEta5_10pos,"nostack");
    MFTEffEta5_10neg.SetLineColor(kRed);
    mftEtaStack.Add(&MFTEffEta5_10neg,"nostack");


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
    MFTEffsEta5.GetXaxis()->SetTitle("\\eta ");
    MFTEffEta5_10pos.SetNameTitle("MFT Eta Efficiency_5_10pos", "5 cm < z < 10 cm");
    MFTEffEta5_10pos.GetXaxis()->SetTitle("\\eta ");
    MFTEffEta5_10neg.SetNameTitle("MFT Eta Efficiency_10_5neg", "-10 cm < z < -5 cm");
    MFTEffEta5_10neg.GetXaxis()->SetTitle("\\eta ");

    MFTEffsp5.SetNameTitle("MFT P Efficiency5_5", "-5 cm < z < 5 cm");
    MFTEffsp5.GetXaxis()->SetTitle("P (GeV)");
    MFTEffp5_10pos.SetNameTitle("MFT P Efficiency_5_10pos", "5 cm < z < 10 cm");
    MFTEffp5_10pos.GetXaxis()->SetTitle("P (GeV)");
    MFTEffp5_10neg.SetNameTitle("MFT P Efficiency_10_5neg", "-10 cm < z < -5 cm");
    MFTEffp5_10neg.GetXaxis()->SetTitle("P (GeV)");

    MFTEffsEta5.SetNameTitle("MFT Eta Efficiency5_5", "-5 cm < z < 5 cm");
    MFTEffsEta5.GetXaxis()->SetTitle("\\eta ");
    MFTEffEta5_10pos.SetNameTitle("MFT Eta Efficiency_5_10pos", "5 cm < z < 10 cm");
    MFTEffEta5_10pos.GetXaxis()->SetTitle("\\eta ");
    MFTEffEta5_10neg.SetNameTitle("MFT Eta Efficiency_10_5neg", "-10 cm < z < -5 cm");
    MFTEffEta5_10neg.GetXaxis()->SetTitle("\\eta ");








    MFTTrackablesEtaZ->SetOption("CONT4");
    MFTTrackablesEtaZ->Write();

    MCTrackspT->Write();
    MCTracksp->Write();
    MCTrackEta->Write();

    MFTTracksMCpT->Write();
    MFTTracksMCp->Write();
    MFTTrackMCEta->Write();

    MFTTracksMCpT->Write();
    MFTTracksMCp->Write();
    MFTTrackMCEta->Write();

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


    MFTTracksMCp5->Write();
    MFTTracksMCp5_10pos->Write();
    MFTTracksMCp5_10neg->Write();

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



    MFTTrackedEtaZ->SetOption("CONT4");
    MFTTrackedEtaZ->Write();


    MFTTrackerEfficiency.SetOption("CONT4");
    MFTTrackerEfficiency.Write();

    MFTEfficiency2D.SetOption("CONT4");
    MFTEfficiency2D.Write();

    MFTAcceptance.SetOption("CONT4");
    MFTAcceptance.Write();


}








// Write histograms to file
//std::cout << "Writting histograms to file..." << std::endl;




// MFT Fitted Results
MFTTracksP->Write();
MFTTracksDeltaP->Write();
MFTTracksDeltaPt->Write();
MFTTrackDeltaEta->Write();
MFTTrackDeltaPhi->Write();
MFTTrackDeltaPhiDeg->Write();
MFTTrackDeltaX->Write();
MFTTrackDeltaY->Write();
MFTTrackR->Write();
MFTTrackDeltaXYVertex->Write();
MFTTrackQ->Write();

outFile.Close();

Int_t totalRecoMFTTracks = nCleanTracksMFT + nInvalidTracksMFT;
std::cout << std::endl;
std::cout << "---------------------------------------------------" << std::endl;
std::cout << "-----------   Track finding Summary   -------------" << std::endl;
std::cout << "---------------------------------------------------" << std::endl;
std::cout << "Number of MFT trackables = " << nMFTTrackable << std::endl;
std::cout << "Number of reconstructed MFT Tracks = " << totalRecoMFTTracks << std::endl;
std::cout << "Number of clean MFT Tracks = " << nCleanTracksMFT << std::endl;
std::cout << "Number of mixed MFT Tracks = " << nInvalidTracksMFT << std::endl;
std::cout << "---------------------------------------------------" << std::endl;
std::cout << "nCleanTracksMFT = " << nCleanTracksMFT << std::endl;
std::cout << "nInvalidTracksMFT = " << nInvalidTracksMFT << " (" << 100.f*nInvalidTracksMFT/(nCleanTracksMFT+nInvalidTracksMFT) << " %)" << std::endl;
std::cout << "---------------------------------------------------" << std::endl;
std::cout << std::endl;

std::cout << "---------------------------------------------------" << std::endl;
std::cout << "-------------   Fitting Summary   -----------------" << std::endl;
std::cout << "---------------------------------------------------" << std::endl;
std::cout << " P_mean = " << MFTTracksP->GetMean() << std::endl;
std::cout << " P_StdDev = " << MFTTracksP->GetStdDev() << std::endl;
std::cout << " DeltaP_mean = " << MFTTracksP->GetMean() << std::endl;
std::cout << " DeltaP_StdDev = " << MFTTracksP->GetStdDev() << std::endl;
std::cout << " DeltaPt_mean = " << MFTTracksDeltaPt->GetMean() << std::endl;
std::cout << " DeltaPt_StdDev = " << MFTTracksDeltaPt->GetStdDev() << std::endl;
std::cout << " Eta_mean = " << MFTTrackDeltaEta->GetMean() << std::endl;
std::cout << " Eta_StdDev = " << MFTTrackDeltaEta->GetStdDev() << std::endl;
std::cout << " Phi_mean = " << MFTTrackDeltaPhi->GetMean() << std::endl;
std::cout << " Phi_StdDev = " << MFTTrackDeltaPhi->GetStdDev() << std::endl;
std::cout << " Phi_mean = " << MFTTrackDeltaPhiDeg->GetMean() << std::endl;
std::cout << " Phi_StdDev = " << MFTTrackDeltaPhiDeg->GetStdDev() << std::endl;
std::cout << " X_mean = " << MFTTrackDeltaX->GetMean() << std::endl;
std::cout << " X_StdDev = " << MFTTrackDeltaX->GetStdDev() << std::endl;
std::cout << " Y_mean = " << MFTTrackDeltaY->GetMean() << std::endl;
std::cout << " Y_StdDev = " << MFTTrackDeltaY->GetStdDev() << std::endl;
std::cout << " R_mean = " << MFTTrackR->GetMean() << std::endl;
std::cout << " R_StdDev = " << MFTTrackR->GetStdDev() << std::endl;
std::cout << " Charge_mean = " << MFTTrackDeltaY->GetMean() << std::endl;
std::cout << "---------------------------------------------------" << std::endl;

}
