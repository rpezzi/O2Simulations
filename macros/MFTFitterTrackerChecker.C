#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <TMath.h>
#include "CommonConstants/MathConstants.h"
#include "SimulationDataFormat/MCTrack.h"
#include "ITSMFTSimulation/Hit.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include <TGeoGlobalMagField.h>
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include "SimulationDataFormat/MCEventHeader.h"
#include <TStyle.h>

//constexpr Double_t MFTLayerZ[] = {-45.3, -46.7, -48.6, -50.0, -52.4, -53.8, -67.7, -69.1, -76.1, -77.5};
using o2::itsmft::Hit;
using o2::MCTrackT;
using o2::mft::TrackMFT;
o2::itsmft::ChipMappingMFT mftChipMapper;
using eventFoundTracks = std::vector<bool>;
vector<eventFoundTracks> allFoundTracksMFT; // True for reconstructed tracks - one vector of bool per event

using trackHasHitsinMFTDisks = std::array<bool,5>; // Disks with hits from a MFT track

bool DEBUG_VERBOSE = false;

bool EXPORT_HISTOS_IMAGES = false;


//_________________________________________________________________________________________________
double getZField(double x, double y, double z) {
const auto grp = o2::parameters::GRPObject::loadFrom("o2sim_grp.root");
if (grp) {
  o2::base::Propagator::initFieldFromGRP(grp);
  auto field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());

  double position[3] = {x,y,z}; // Field at center of MFT
  return field->getBz(position);

} else {
  LOG(ERROR) << "Cannot retrieve GRP from file !";
  return 0;
}
}

//_________________________________________________________________________________________________
void extrapMFTTrackHelixToZ(o2::mft::TrackMFT& track, double zEnd, double Field)
{
   using TrackMFT = o2::mft::TrackMFT;

  /// Track extrapolated to the plane at "zEnd" considering a helix

  if (track.getZ() == zEnd) {
    return; // nothing to be done if same z
  }

  // Compute track parameters
  double dZ = (zEnd - track.getZ());
  double x0 = track.getX();
  double y0 = track.getY();
  double phi0 = track.getPhi();
  double cosphi0 = TMath::Cos(phi0);
  double sinphi0 = TMath::Sin(phi0);
  double invtanl0 = 1.0 / track.getTanl();
  double invqpt0 = track.getInvQPt();

  double k = 3e-4 * TMath::Abs(Field);
  double n = dZ * invtanl0;
  double theta = - invqpt0 * dZ * k * invtanl0;
  auto Hz = std::copysign(1.0,Field);
  double deltax = n * cosphi0 - 0.5 * n * theta * Hz * sinphi0;
  double deltay = n * sinphi0 + 0.5 * n * theta * Hz * cosphi0;

  double x = x0 + deltax;
  double y = y0 + deltay;
  double phi = phi0 + theta;

  track.setX(x);
  track.setY(y);
  track.setZ(zEnd);
  track.setPhi(phi);
}

//_________________________________________________________________________________________________
template <typename H>
void exportHisto(const H& histo)
{
   //gStyle->SetImageScaling(3.);
   TCanvas *c = new TCanvas;
   c->SetBatch();
   std::string imgpath{"images/"};
   gSystem->MakeDirectory(imgpath.c_str());
   H *h = new H(histo);
   h->Draw();
   gSystem->ProcessEvents();
   for (std::string type: {".pdf", ".png"}) c->Print((imgpath+std::string(h->GetName()) + type).c_str());
}

//_________________________________________________________________________________________________
int MFTFitterTrackerChecker( const Char_t *trkFile = "mfttracks.root",
                              const Char_t *o2sim_KineFile = "o2sim_Kine.root",
                              const Char_t *HitsMFTFile = "o2sim_HitsMFT.root",
                              Double_t pMin = 0.0,
                              Double_t pMax = 100.0,
                              Double_t deltaetaMin = -.1,
                              Double_t deltaetaMax = +.1,
                              Double_t etaMin = -3.4,
                              Double_t etaMax = -2.4,
                              Double_t deltaphiMin = -.2, //-3.15,
                              Double_t deltaphiMax = .2 //+3.15,
                            ) {


  // histos

  enum TH2HistosCodes {
    kMFTTrackDeltaXYVertex,
    kMFTrackQPRec_MC,
    kMFTrackPtResolution,
    kMCTracksEtaZ
  };

  std::map<int,const char *> TH2Names {
    {kMFTTrackDeltaXYVertex, "MFT Tracks Vertex at Z = 0"},
    {kMFTrackQPRec_MC, "MFT Track QP FITxMC"},
    {kMFTrackPtResolution, "MFT Track Pt Resolution"},
    {kMCTracksEtaZ, "MCTracks_eta_z"}
  };

  std::map<int,const char *> TH2Titles {
    {kMFTTrackDeltaXYVertex, "Standalone MFT Tracks at Z_vertex"},
    {kMFTrackQPRec_MC, "Charged Momentum: Reconstructed vs MC"},
    {kMFTrackPtResolution, "Pt Resolution"},
    {kMCTracksEtaZ, "MC Tracks: Pseudorapidity vs zVertex"}
  };

  std::map<int, std::array<double,6>> TH2Binning {
    {kMFTTrackDeltaXYVertex, {300, -.05, .05, 300, -.05, .05} },
    {kMFTrackQPRec_MC, {150, -10, 10, 150, -10, 10} },
    {kMFTrackPtResolution, {200, 0, 5, 200, 0, 15} },
    {kMCTracksEtaZ, {31, -15, 16, 25, etaMin, etaMax} }
  };


  std::map<int,const char *> TH2XaxisTitles {
    {kMFTTrackDeltaXYVertex, "\\Delta x ~[cm]"},
    {kMFTrackQPRec_MC, "(q.p)_{MC} [GeV]"},
    {kMFTrackPtResolution, "pt_{MC} [GeV]"},
    {kMCTracksEtaZ, "Vertex PosZ [cm]"}
  };

  std::map<int,const char *> TH2YaxisTitles {
    {kMFTTrackDeltaXYVertex, "\\Delta y ~[cm]"},
    {kMFTrackQPRec_MC, "(q.p)_{fit} [GeV]"},
    {kMFTrackPtResolution, "(p_t)_{fit} / (p_t)_{MC}"},
    {kMCTracksEtaZ, "\\eta"}
  };


  enum TH1HistosCodes {
    kMFTTracksP,
    kMFTTracksP_res,
    kMFTTracksPt_res,
    kMFTTrackDeltaEta,
    kMFTTrackDeltaPhi,
    kMFTTrackDeltaPhiDeg,
    kMFTTrackDeltaX,
    kMFTTrackDeltaY,
    kMFTTrackR,
    kMFTTrackQ,
    kMCTrackspT,
    kMCTracksp,
    kMCTrackEta
  };

  std::map<int,const char *> TH1Names {
    {kMFTTracksP, "MFT Tracks Fitted p"},
    {kMFTTracksP_res, "MFT Tracks P resolution"},
    {kMFTTracksPt_res, "MFT Tracks P_t resolution"},
    {kMFTTrackDeltaEta, "MFT Tracks Fitted Delta_eta"},
    {kMFTTrackDeltaPhi, "MFT Tracks Fitted Phi at Vertex"},
    {kMFTTrackDeltaPhiDeg,"MFT Tracks Fitted Phi at Vertex [deg]"},
    {kMFTTrackDeltaX, "MFT Tracks Delta X"},
    {kMFTTrackDeltaY, "MFT Tracks Delta Y"},
    {kMFTTrackR, "MFT Tracks Delta R"},
    {kMFTTrackQ, "MFT Tracks Q"},
    {kMCTrackspT, "MC Tracks p_T"},
    {kMCTracksp, "MC Tracks p"},
    {kMCTrackEta, "MC Tracks eta"}
  };

  std::map<int,const char *> TH1Titles {
    {kMFTTracksP, "Standalone MFT Tracks P"},
    {kMFTTracksP_res, "P_{Fit}/P_{MC}"},
    {kMFTTracksPt_res,"Pt_{Fit}/Pt_{MC}" },
    {kMFTTrackDeltaEta, "\\eta_{Fit} - \\eta_{MC} "},
    {kMFTTrackDeltaPhi, "\\phi _{Fit} - \\phi_{MC}"},
    {kMFTTrackDeltaPhiDeg, "\\phi _{Fit} - \\phi_{MC}"},
    {kMFTTrackDeltaX, "Standalone MFT Tracks Delta X at Z_vertex"},
    {kMFTTrackDeltaY, "Standalone MFT Tracks Delta Y at Z_vertex"},
    {kMFTTrackR, "Standalone MFT Tracks Delta R at Z_vertex"},
    {kMFTTrackQ, "Standalone MFT Tracks Charge"},
    {kMCTrackspT, "MC Tracks p_T"},
    {kMCTracksp, "MC Tracks p"},
    {kMCTrackEta, "MC Tracks Pseudorapidity"}
  };

  std::map<int, std::array<double,3>> TH1Binning {
    {kMFTTracksP, {500, pMin, pMax} },
    {kMFTTracksP_res,  {500, 0, 50}},
    {kMFTTracksPt_res,  {300, 0, 10}},
    {kMFTTrackDeltaEta, {1000, deltaetaMin, deltaetaMax}},
    {kMFTTrackDeltaPhi, {1000, deltaphiMin, deltaphiMax}},
    {kMFTTrackDeltaPhiDeg, {1000, TMath::RadToDeg()*deltaphiMin, TMath::RadToDeg()*deltaphiMax}},
    {kMFTTrackDeltaX, {1000, -.3, .3}},
    {kMFTTrackDeltaY, {1000, -.3, .3}},
    {kMFTTrackR, {250, 0, 0.5}},
    {kMFTTrackQ, {5, -2.1, 2.1}},
    {kMCTrackspT, {5000, 0, 50}},
    {kMCTracksp, {1000, pMin, pMax}},
    {kMCTrackEta, {1000, etaMin, etaMax}}
  };

  std::map<int,const char *> TH1XaxisTitles {
    {kMFTTracksP, "p [GeV]"},
    {kMFTTracksP_res, "P_{Fit}/P_{MC}"},
    {kMFTTracksPt_res, "Pt_{Fit}/Pt_{MC}"},
    {kMFTTrackDeltaEta, "\\Delta \\eta"},
    {kMFTTrackDeltaPhi, "\\Delta \\phi ~[rad]"},
    {kMFTTrackDeltaPhiDeg, "\\Delta \\phi ~[deg]"},
    {kMFTTrackDeltaX, "\\Delta x ~[cm]"},
    {kMFTTrackDeltaY, "\\Delta y ~[cm]"},
    {kMFTTrackR, "\\Delta r ~[cm]"},
    {kMFTTrackQ, "Q"},
    {kMCTrackspT, "p_t [GeV]"},
    {kMCTracksp, "p [GeV]"},
    {kMCTrackEta, " \\eta"}
  };



  const int nTH1Histos = TH1Names.size();
  std::vector<std::unique_ptr<TH1F>> TH1Histos(nTH1Histos);
  auto nHisto = 0;
  for (auto& h : TH1Histos) {
    h = std::make_unique<TH1F> (TH1Names[nHisto], TH1Titles[nHisto], (int)TH1Binning[nHisto][0], TH1Binning[nHisto][1], TH1Binning[nHisto][2]);
    h->GetXaxis()->SetTitle(TH1XaxisTitles[nHisto]);

    ++nHisto;
    }

  const int nTH2Histos = TH2Names.size();
  std::vector<std::unique_ptr<TH2F>> TH2Histos(nTH2Histos);
  auto n2Histo = 0;
  for (auto& h : TH2Histos) {
    h = std::make_unique<TH2F> (TH2Names[n2Histo], TH2Titles[n2Histo], (int)TH2Binning[n2Histo][0], TH2Binning[n2Histo][1], TH2Binning[n2Histo][2], (int)TH2Binning[n2Histo][3], TH2Binning[n2Histo][4], TH2Binning[n2Histo][5]);
    h->GetXaxis()->SetTitle(TH2XaxisTitles[n2Histo]);
    h->GetYaxis()->SetTitle(TH2YaxisTitles[n2Histo]);
    h->SetOption("COLZ");
    ++n2Histo;
    }

  Int_t nChargeMatch = 0;
  Int_t nChargeMiss = 0;

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

  auto field_z = getZField(0, 0, -61.4); // Get field at Center of MFT

  //char [100];
  std::string outfilename = "Fittercheck_" + std::string(trkFile);
  //strcat(outfilename,trkFile);
  //strcat(outfilename,"");
  TFile outFile(outfilename.c_str(),"RECREATE");


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


  // Reconstructed MFT Tracks
   std::cout << "Loop over reconstructed MFT Tracks!" << std::endl;
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
        if (thisTrack->getMotherTrackId() == -1) { // Only primaries
          auto vx_MC = thisTrack->GetStartVertexCoordinatesX();
          auto vy_MC = thisTrack->GetStartVertexCoordinatesY();
          auto vz_MC = thisTrack->GetStartVertexCoordinatesZ();
          auto Pt_MC = thisTrack->GetPt();
          auto P_MC = thisTrack->GetP();
          auto phi_MC = TMath::ATan2(thisTrack->Py(),thisTrack->Px());
          auto eta_MC = atanh (thisTrack->GetStartVertexMomentumZ()/P_MC); // eta;
          auto pdgcode_MC = thisTrack->GetPdgCode();
          //std::cout << "pdgcode_MC = " <<  pdgcode_MC;
          int Q_MC;
          if (TDatabasePDG::Instance()->GetParticle(pdgcode_MC)) {
           Q_MC = TDatabasePDG::Instance()->GetParticle(pdgcode_MC)->Charge()/3;
           //std::cout << " => " <<  TDatabasePDG::Instance()->GetParticle(pdgcode_MC)->GetName() << " ; q = " << Q_MC <<  "\n";
           }

           else {
             Q_MC = 0;
             std::cout << " => pdgcode ERROR " << Q_MC <<  "\n";
           }

          //auto forceP = 0.5;
          //auto forcePt = forceP/sqrt(1+trackMFT.getTanl()*trackMFT.getTanl());
          //trackMFT.setInvQPt(1.0/forcePt);
          //trackMFT.setCharge(-1);

          extrapMFTTrackHelixToZ(trackMFT, vz_MC, field_z); // propagate track to vertex Z
          auto dx = trackMFT.getX() - vx_MC;
          auto dy = trackMFT.getY() - vy_MC;
          auto d_eta = trackMFT.getEta() - eta_MC;
          auto Pt_fit = trackMFT.getPt();
          auto P_fit = trackMFT.getP();
          auto Q_fit = trackMFT.getCharge();
          auto P_res = P_fit / P_MC;
          auto Pt_res = Pt_fit / Pt_MC;
          auto d_Phi = trackMFT.getPhi() - phi_MC;
          auto d_Charge = Q_fit-Q_MC;

          TH1Histos[kMFTTracksP]->Fill(trackMFT.getP());
          TH1Histos[kMFTTracksP_res]->Fill(P_res);
          TH1Histos[kMFTTracksPt_res]->Fill(Pt_res);
          TH1Histos[kMFTTrackDeltaEta]->Fill(d_eta);
          TH1Histos[kMFTTrackDeltaPhi]->Fill(d_Phi);
          TH1Histos[kMFTTrackDeltaPhiDeg]->Fill(TMath::RadToDeg()*d_Phi);
          TH1Histos[kMFTTrackDeltaX]->Fill(dx);
          TH1Histos[kMFTTrackDeltaY]->Fill(dy);
          TH1Histos[kMFTTrackR]->Fill(sqrt(dx*dx+dy*dy));
          TH1Histos[kMFTTrackQ]->Fill(d_Charge);
          TH2Histos[kMFTTrackDeltaXYVertex]->Fill(dx,dy);
          TH2Histos[kMFTrackQPRec_MC]->Fill(P_MC*Q_MC,P_fit*Q_fit);
          TH2Histos[kMFTrackPtResolution]->Fill(Pt_MC,Pt_fit/Pt_MC);

          TH1Histos[kMCTrackspT]->Fill(Pt_MC);
          TH1Histos[kMCTracksp]->Fill(P_MC);
          TH1Histos[kMCTrackEta]->Fill(eta_MC);
          TH2Histos[kMCTracksEtaZ]->Fill(vz_MC,eta_MC);

          if (d_Charge == 0) {
            nChargeMatch++;
          }
          else
          nChargeMiss++;
          }

        }
  } // Loop on TracksMFT


// Customize histograms
gStyle->SetStatFormat("4.3g");
gStyle->SetStatW(.38);
gStyle->SetStatH(.26);

TH1Histos[kMFTTrackQ]->SetTitle(Form("nChargeMatch = %d (%.2f%%)", nChargeMatch, 100.*nChargeMatch/(nChargeMiss+nChargeMatch)));

// Write histograms to file and export images
for (auto& h : TH2Histos) {
  h->Write();
  if (EXPORT_HISTOS_IMAGES)
    exportHisto(*h);
  }

for (auto& h : TH1Histos) {
  h->Write();
  if (EXPORT_HISTOS_IMAGES)
    exportHisto(*h);
  }
outFile.Close();


Int_t totalRecoMFTTracks = nCleanTracksMFT + nInvalidTracksMFT;
std::cout << std::endl;
std::cout << "---------------------------------------------------" << std::endl;
std::cout << "-----------   Track finding Summary   -------------" << std::endl;
std::cout << "---------------------------------------------------" << std::endl;
std::cout << "Number of reconstructed MFT Tracks = " << totalRecoMFTTracks << std::endl;
std::cout << "Number of clean MFT Tracks = " << nCleanTracksMFT << std::endl;
std::cout << "Number of invalid MFT Tracks = " << nInvalidTracksMFT << std::endl;
std::cout << "---------------------------------------------------" << std::endl;
std::cout << std::endl;

std::cout << "---------------------------------------------------" << std::endl;
std::cout << "-------------   Fitting Summary   -----------------" << std::endl;
std::cout << "---------------------------------------------------" << std::endl;
std::cout << " P_mean = " << TH1Histos[kMFTTracksP]->GetMean() << std::endl;
std::cout << " P_StdDev = " << TH1Histos[kMFTTracksP]->GetStdDev() << std::endl;
std::cout << " P_Res_mean = " << TH1Histos[kMFTTracksP_res]->GetMean() << std::endl;
std::cout << " P_Res_StdDev = " << TH1Histos[kMFTTracksP_res]->GetStdDev() << std::endl;
std::cout << " Pt_Res_mean = " << TH1Histos[kMFTTracksPt_res]->GetMean() << std::endl;
std::cout << " Pt_Res_StdDev = " << TH1Histos[kMFTTracksPt_res]->GetStdDev() << std::endl;
std::cout << " Eta_mean = " << TH1Histos[kMFTTrackDeltaEta]->GetMean() << std::endl;
std::cout << " Eta_StdDev = " << TH1Histos[kMFTTrackDeltaEta]->GetStdDev() << std::endl;
std::cout << " Phi_mean = " << TH1Histos[kMFTTrackDeltaPhi]->GetMean() << std::endl;
std::cout << " Phi_StdDev = " << TH1Histos[kMFTTrackDeltaPhi]->GetStdDev() << std::endl;
std::cout << " Phi_mean = " << TH1Histos[kMFTTrackDeltaPhiDeg]->GetMean() << std::endl;
std::cout << " Phi_StdDev = " << TH1Histos[kMFTTrackDeltaPhiDeg]->GetStdDev() << std::endl;
std::cout << " DeltaX_mean = " << TH1Histos[kMFTTrackDeltaX]->GetMean() << std::endl;
std::cout << " DeltaX_StdDev = " << TH1Histos[kMFTTrackDeltaX]->GetStdDev() << std::endl;
std::cout << " DeltaY_mean = " << TH1Histos[kMFTTrackDeltaY]->GetMean() << std::endl;
std::cout << " DeltaY_StdDev = " << TH1Histos[kMFTTrackDeltaY]->GetStdDev() << std::endl;
std::cout << " R_mean = " << TH1Histos[kMFTTrackR]->GetMean() << std::endl;
std::cout << " R_StdDev = " << TH1Histos[kMFTTrackR]->GetStdDev() << std::endl;
std::cout << " Charge_mean = " << TH1Histos[kMFTTrackDeltaY]->GetMean() << std::endl;
std::cout << " nChargeMatch = " << nChargeMatch << " (" << 100.*nChargeMatch/(nChargeMiss+nChargeMatch) << "%)" << std::endl;
std::cout << "---------------------------------------------------" << std::endl;

return 0;
}
