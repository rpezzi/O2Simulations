#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TLatex.h"
#include "TEfficiency.h"
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
//#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include "SimulationDataFormat/MCEventHeader.h"
#include <TStyle.h>
#include <TProfile.h>
#include <TGraph.h>

// MFTTools

#include "mfttools/MagField.C"
#include "mfttools/MFTTrackExtrap.C"
#include "mfttools/HistosHelpers.C"



using o2::itsmft::Hit;
using o2::MCTrackT;
using o2::mft::TrackMFT;
o2::itsmft::ChipMappingMFT mftChipMapper;
using eventFoundTracks = std::vector<bool>;
using std::vector;
vector<eventFoundTracks> allFoundTracksMFT; // True for reconstructed tracks - one vector of bool per event

using trackHasHitsinMFTDisks = std::array<bool,5>; // Disks with hits from a MFT track

bool DEBUG_VERBOSE = false;
bool EXPORT_HISTOS_IMAGES = false;




//_________________________________________________________________________________________________
int MFTSmootherVsExtrap( const Char_t *trkFile = "mfttracks.root",
                              const Char_t *o2sim_KineFile = "o2sim_Kine.root",
                              const Char_t *HitsMFTFile = "o2sim_HitsMFT.root",
                              Double_t pMin = 0.0,
                              Double_t pMax = 100.0,
                              Double_t deltaetaMin = -.1/5,
                              Double_t deltaetaMax = +.1/5,
                              Double_t etaMin = -3.4,
                              Double_t etaMax = -2.4,
                              Double_t deltaphiMin = -.005, //-3.15,
                              Double_t deltaphiMax = .005 //+3.15,
                            ) {


    std::cout << "########################################################" << std::endl;
    std::cout << "#########  WARNING! WARNING! WARNING! WARNING! #########" << std::endl;
    std::cout << "########################################################" << std::endl;
    std::cout << "# This macro has been modified to propagate MFT tracks #" << std::endl;
    std::cout << "# to the last cluster, closest to the Hadron absorber. #" << std::endl;
    std::cout << "# Intended to compare smoothed and propagated track    #" << std::endl;
    std::cout << "# parameters. Please ignore references to vertex       #" << std::endl;
    std::cout << "########################################################" << std::endl;


  // Seed configuration
  std::string seed_cfg{trkFile};
  std::string trk_start{"mfttracks_"};
  std::string trk_ext{".root"};
  std::string trk_trk{"mfttracks"};
  if (seed_cfg.find(trk_start) < seed_cfg.length()) seed_cfg.replace(seed_cfg.find(trk_start),trk_start.length(),"");
  if (seed_cfg.find(trk_ext) < seed_cfg.length()) seed_cfg.replace(seed_cfg.find(trk_ext),trk_ext.length(),"");
  if (seed_cfg.find(trk_trk) < seed_cfg.length()) seed_cfg.replace(seed_cfg.find(trk_trk),trk_trk.length(),"");
  std::cout << seed_cfg << std::endl;


  // histos
  //gROOT->SetStyle("Bold");
  //gStyle->SetStatW(.38);
  //gStyle->SetStatH(.26);
  gStyle->SetPalette(1,0);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(4);
  gStyle->SetFrameFillColor(10);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  //gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(3);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(3);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(3);
  gStyle->SetLabelSize(0.06,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  //gStyle->SetLabelColor(kBlue,"xy");
  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetTitleSize(0.08,"o");
  gStyle->SetTitleOffset(.95,"Y");
  gStyle->SetTitleFillColor(10);
  //gStyle->SetTitleTextColor(kNlacBlue);
  gStyle->SetStatColor(10);



  enum TH2HistosCodes {
    kMFTTrackDeltaXYVertex,
    kMFTTrackDeltaXYVertex0_1,
    kMFTTrackDeltaXYVertex1_4,
    kMFTTrackDeltaXYVertex4plus,
    kMFTrackQPRec_MC,
    kMFTrackPtResolution,
    kMFTrackInvPtResolution,
    kMCTracksEtaZ
  };

  std::map<int,const char *> TH2Names {
    {kMFTTrackDeltaXYVertex, "MFT Tracks Vertex at Z = 0"},
    {kMFTTrackDeltaXYVertex0_1, "MFT Tracks Vertex at Z = 0 P0_1"},
    {kMFTTrackDeltaXYVertex1_4, "MFT Tracks Vertex at Z = 0 P1_4"},
    {kMFTTrackDeltaXYVertex4plus, "MFT Tracks Vertex at Z = 0 P4plus"},
    {kMFTrackQPRec_MC, "MFT Track QP ExtrapxSmoothed"},
    {kMFTrackPtResolution, "MFT Track Pt Resolution"},
    {kMFTrackInvPtResolution, "MFT Track InvPt Resolution"},
    {kMCTracksEtaZ, "MCTracks_eta_z"}
  };

  std::map<int,const char *> TH2Titles {
    {kMFTTrackDeltaXYVertex, "MFT Tracks at Z_lastCluster"},
    {kMFTTrackDeltaXYVertex0_1, "MFT Tracks at Z_lastCluster (p < 1)"},
    {kMFTTrackDeltaXYVertex1_4, "MFT Tracks at Z_lastCluster (1 < p < 4)"},
    {kMFTTrackDeltaXYVertex4plus, "MFT Tracks at Z_lastCluster (p > 4)"},
    {kMFTrackQPRec_MC, "Charged Momentum: Extrap vs Smoothed"},
    {kMFTrackPtResolution, "Pt Resolution"},
    {kMFTrackInvPtResolution, "InvPt Resolution"},
    {kMCTracksEtaZ, "Pseudorapidity vs z"}
  };

  std::map<int, std::array<double,6>> TH2Binning {
    {kMFTTrackDeltaXYVertex, {300, -.05, .05, 300, -.05, .05} },
    {kMFTTrackDeltaXYVertex0_1, {300, -.05, .05, 300, -.05, .05} },
    {kMFTTrackDeltaXYVertex1_4, {300, -.05, .05, 300, -.05, .05} },
    {kMFTTrackDeltaXYVertex4plus, {300, -.05, .05, 300, -.05, .05} },
    {kMFTrackQPRec_MC, {100, -10, 10, 100, -10, 10} },
    {kMFTrackPtResolution, {14, 0, 7, 250, 0, 25} },
    {kMFTrackInvPtResolution, {14, 0, 7, 300, -2, 2} },
//    {kMFTrackPtResolution, {100, 0, 5, 100, 0, 15} },
//    {kMFTrackInvPtResolution, {100, 0, 5, 100, -15, 15} },
    {kMCTracksEtaZ, {31, -15, 16, 25, etaMin, etaMax} }
  };


  std::map<int,const char *> TH2XaxisTitles {
    {kMFTTrackDeltaXYVertex, "\\Delta x ~[cm]"},
    {kMFTTrackDeltaXYVertex0_1, "\\Delta x ~[cm]"},
    {kMFTTrackDeltaXYVertex1_4, "\\Delta x ~[cm]"},
    {kMFTTrackDeltaXYVertex4plus, "\\Delta x ~[cm]"},
    {kMFTrackQPRec_MC, "(q.p)_{Smoothed} [GeV]"},
    {kMFTrackPtResolution, "pt_{Smoothed} [GeV]"},
    {kMFTrackInvPtResolution, "pt_{Smoothed} [GeV]"},
    {kMCTracksEtaZ, "PosZ [cm]"}
  };

  std::map<int,const char *> TH2YaxisTitles {
    {kMFTTrackDeltaXYVertex, "\\Delta y ~[cm]"},
    {kMFTTrackDeltaXYVertex0_1, "\\Delta y ~[cm]"},
    {kMFTTrackDeltaXYVertex1_4, "\\Delta y ~[cm]"},
    {kMFTTrackDeltaXYVertex4plus, "\\Delta y ~[cm]"},
    {kMFTrackQPRec_MC, "(q.p)_{Extrap} [GeV]"},
    {kMFTrackPtResolution, "pt_{Extrap} / pt_{Smoothed}"},
    {kMFTrackInvPtResolution, "(1/(p_t)_{Extrap} - 1/(p_t)_{Smoothed})*(p_t)_{Smoothed}"},
    {kMCTracksEtaZ, "\\eta"}
  };


  enum TH1HistosCodes {
    kMFTTracksP,
    //kMFTTracksP_res,
    //kMFTTracksPt_res,
    kMFTTrackDeltaEta,
    kMFTTrackDeltaEta0_1,
    kMFTTrackDeltaEta1_4,
    kMFTTrackDeltaEta4plus,
    kMFTTrackDeltaPhi,
    kMFTTrackDeltaPhi0_1,
    kMFTTrackDeltaPhi1_4,
    kMFTTrackDeltaPhi4plus,
    kMFTTrackDeltaPhiDeg,
    kMFTTrackDeltaPhiDeg0_1,
    kMFTTrackDeltaPhiDeg1_4,
    kMFTTrackDeltaPhiDeg4plus,
    kMFTTrackDeltaX,
    kMFTTrackDeltaX0_1,
    kMFTTrackDeltaX1_4,
    kMFTTrackDeltaX4plus,
    kMFTTrackDeltaY,
    kMFTTrackR,
    kMFTTrackQ,
    kMFTTrackQ0_1,
    kMFTTrackQ1_4,
    kMFTTrackQ4plus,
    kMCTrackspT,
    kMCTracksp,
    kMCTrackEta
  };



  std::map<int,const char *> TH1Names {
    {kMFTTracksP, "MFT Tracks Fitted p"},
    //{kMFTTracksP_res, "MFT Tracks P resolution"},
    //{kMFTTracksPt_res, "MFT Tracks P_t resolution"},
    {kMFTTrackDeltaEta, "MFT Tracks Fitted Delta_eta"},
    {kMFTTrackDeltaEta0_1, "MFT Tracks eta (p < 1)"},
    {kMFTTrackDeltaEta1_4, "MFT Tracks eta (1 < p < 4)"},
    {kMFTTrackDeltaEta4plus, "MFT Tracks eta (p > 4)"},
    {kMFTTrackDeltaPhi, "MFT Tracks Fitted Phi at Vertex"},
    {kMFTTrackDeltaPhi0_1,"MFT Tracks Fitted Phi at Vertex [rad] (p < 1)"},
    {kMFTTrackDeltaPhi1_4,"MFT Tracks Fitted Phi at Vertex [rad] (1 < p < 4)"},
    {kMFTTrackDeltaPhi4plus,"MFT Tracks Fitted Phi at Vertex [rad] (p > 4)"},
    {kMFTTrackDeltaPhiDeg,"MFT Tracks Fitted Phi at Vertex [deg]"},
    {kMFTTrackDeltaPhiDeg0_1,"MFT Tracks Fitted Phi at Vertex [deg] (p < 1)"},
    {kMFTTrackDeltaPhiDeg1_4,"MFT Tracks Fitted Phi at Vertex [deg] (1 < p < 4)"},
    {kMFTTrackDeltaPhiDeg4plus,"MFT Tracks Fitted Phi at Vertex [deg] (p > 4)"},
    {kMFTTrackDeltaX, "MFT Tracks Delta X"},
    {kMFTTrackDeltaX0_1, "MFT Tracks Delta X (p < 1)"},
    {kMFTTrackDeltaX1_4, "MFT Tracks Delta X (1 < p < 4)"},
    {kMFTTrackDeltaX4plus, "MFT Tracks Delta X (p > 4)"},
    {kMFTTrackDeltaY, "MFT Tracks Delta Y"},
    {kMFTTrackR, "MFT Tracks Delta R"},
    {kMFTTrackQ, "Charge Match"},
    {kMFTTrackQ0_1, "Charge Match (p < 1)"},
    {kMFTTrackQ1_4, "Charge Match (1 < p < 4)"},
    {kMFTTrackQ4plus, "Charge Match (p > 4)"},
    {kMCTrackspT, "Smoothed Tracks p_T"},
    {kMCTracksp, "Smoothed Tracks p"},
    {kMCTrackEta, "Smoothed Tracks eta"}
  };

  std::map<int,const char *> TH1Titles {
    {kMFTTracksP, "Standalone MFT Tracks P"},
    //{kMFTTracksP_res, "P_{Extrap}/P_{Smoothed}"},
    //{kMFTTracksPt_res,"Pt_{Extrap}/Pt_{Smoothed}" },
    {kMFTTrackDeltaEta, "\\eta_{Extrap} - \\eta_{Smoothed} "},
    {kMFTTrackDeltaEta0_1, "\\eta_{Extrap} - \\eta_{Smoothed} "},
    {kMFTTrackDeltaEta1_4, "\\eta_{Extrap} - \\eta_{Smoothed} "},
    {kMFTTrackDeltaEta4plus, "\\eta_{Extrap} - \\eta_{Smoothed} "},
    {kMFTTrackDeltaPhi, "\\phi _{Extrap} - \\phi_{Smoothed}"},
    {kMFTTrackDeltaPhi0_1, "\\phi _{Extrap} - \\phi_{Smoothed}"},
    {kMFTTrackDeltaPhi1_4, "\\phi _{Extrap} - \\phi_{Smoothed}"},
    {kMFTTrackDeltaPhi4plus, "\\phi _{Extrap} - \\phi_{Smoothed}"},
    {kMFTTrackDeltaPhiDeg, "\\phi _{Extrap} - \\phi_{Smoothed}"},
    {kMFTTrackDeltaPhiDeg0_1, "\\phi _{Extrap} - \\phi_{Smoothed}"},
    {kMFTTrackDeltaPhiDeg1_4, "\\phi _{Extrap} - \\phi_{Smoothed}"},
    {kMFTTrackDeltaPhiDeg4plus, "\\phi _{Extrap} - \\phi_{Smoothed}"},
    {kMFTTrackDeltaX, "MFT Tracks Delta X at Z_lastCluster"},
    {kMFTTrackDeltaX0_1, "MFT Tracks Delta X at Z_lastCluster"},
    {kMFTTrackDeltaX1_4, "MFT Tracks Delta X at Z_lastCluster"},
    {kMFTTrackDeltaX4plus, "MFT Tracks Delta X at Z_lastCluster"},
    {kMFTTrackDeltaY, "MFT Tracks Delta Y at Z_lastCluster"},
    {kMFTTrackR, "MFT Tracks Delta R at Z_lastCluster"},
    {kMFTTrackQ, "MFT Tracks Charge Match"},
    {kMFTTrackQ0_1, "MFT Tracks Charge Match (p < 1)"},
    {kMFTTrackQ1_4, "MFT Tracks Charge Match (1 < p < 4)"},
    {kMFTTrackQ4plus, "MFT Tracks Charge Match (p > 4)"},
    {kMCTrackspT, "Smoothed Tracks p_T"},
    {kMCTracksp, "Smoothed Tracks p"},
    {kMCTrackEta, "Smoothed Tracks Pseudorapidity"}
  };

  std::map<int, std::array<double,3>> TH1Binning {
    {kMFTTracksP, {500, pMin, pMax} },
    //{kMFTTracksP_res,  {500, 0, 50}},
    //{kMFTTracksPt_res,  {300, 0, 10}},
    {kMFTTrackDeltaEta, {1000, deltaetaMin, deltaetaMax}},
    {kMFTTrackDeltaEta0_1, {1000, deltaetaMin, deltaetaMax}},
    {kMFTTrackDeltaEta1_4, {1000, deltaetaMin, deltaetaMax}},
    {kMFTTrackDeltaEta4plus, {1000, deltaetaMin, deltaetaMax}},
    {kMFTTrackDeltaPhi, {1000, deltaphiMin, deltaphiMax}},
    {kMFTTrackDeltaPhi0_1, {1000, deltaphiMin, deltaphiMax}},
    {kMFTTrackDeltaPhi1_4, {1000, deltaphiMin, deltaphiMax}},
    {kMFTTrackDeltaPhi4plus, {1000, deltaphiMin, deltaphiMax}},
    {kMFTTrackDeltaPhiDeg, {1000, TMath::RadToDeg()*deltaphiMin, TMath::RadToDeg()*deltaphiMax}},
    {kMFTTrackDeltaPhiDeg0_1, {1000, TMath::RadToDeg()*deltaphiMin, TMath::RadToDeg()*deltaphiMax}},
    {kMFTTrackDeltaPhiDeg1_4, {1000, TMath::RadToDeg()*deltaphiMin, TMath::RadToDeg()*deltaphiMax}},
    {kMFTTrackDeltaPhiDeg4plus, {1000, TMath::RadToDeg()*deltaphiMin, TMath::RadToDeg()*deltaphiMax}},
    {kMFTTrackDeltaX, {1000, -.5, .5}},
    {kMFTTrackDeltaX0_1, {1000, -.5, .5}},
    {kMFTTrackDeltaX1_4, {1000, -.5, .5}},
    {kMFTTrackDeltaX4plus, {1000, -.5, .5}},
    {kMFTTrackDeltaY, {1000, -.5, .5}},
    {kMFTTrackR, {250, 0, 0.5}},
    {kMFTTrackQ, {5, -2.1, 2.1}},
    {kMFTTrackQ0_1, {5, -2.1, 2.1}},
    {kMFTTrackQ1_4, {5, -2.1, 2.1}},
    {kMFTTrackQ4plus, {5, -2.1, 2.1}},
    {kMCTrackspT, {5000, 0, 50}},
    {kMCTracksp, {1000, pMin, pMax}},
    {kMCTrackEta, {1000, etaMin, etaMax}}
  };

  std::map<int,const char *> TH1XaxisTitles {
    {kMFTTracksP, "p [GeV]"},
    //{kMFTTracksP_res, "P_{Extrap}/P_{Smoothed}"},
    //{kMFTTracksPt_res, "Pt_{Extrap}/Pt_{Smoothed}"},
    {kMFTTrackDeltaEta, "\\Delta \\eta"},
    {kMFTTrackDeltaEta0_1, "\\Delta \\eta"},
    {kMFTTrackDeltaEta1_4, "\\Delta \\eta"},
    {kMFTTrackDeltaEta4plus, "\\Delta \\eta"},
    {kMFTTrackDeltaPhi, "\\Delta \\phi ~[rad]"},
    {kMFTTrackDeltaPhi0_1, "\\Delta \\phi ~[rad]"},
    {kMFTTrackDeltaPhi1_4, "\\Delta \\phi ~[rad]"},
    {kMFTTrackDeltaPhi4plus, "\\Delta \\phi ~[rad]"},
    {kMFTTrackDeltaPhiDeg, "\\Delta \\phi ~[deg]"},
    {kMFTTrackDeltaPhiDeg0_1, "\\Delta \\phi ~[deg]"},
    {kMFTTrackDeltaPhiDeg1_4, "\\Delta \\phi ~[deg]"},
    {kMFTTrackDeltaPhiDeg4plus, "\\Delta \\phi ~[deg]"},
    {kMFTTrackDeltaX, "\\Delta x ~[cm]"},
    {kMFTTrackDeltaX0_1, "\\Delta x ~[cm]"},
    {kMFTTrackDeltaX1_4, "\\Delta x ~[cm]"},
    {kMFTTrackDeltaX4plus, "\\Delta x ~[cm]"},
    {kMFTTrackDeltaY, "\\Delta y ~[cm]"},
    {kMFTTrackR, "\\Delta r ~[cm]"},
    {kMFTTrackQ, "q_{Extrap}-q_{Smoothed}"},
    {kMFTTrackQ0_1, "q_{Extrap}-q_{Smoothed}"},
    {kMFTTrackQ1_4, "q_{Extrap}-q_{Smoothed}"},
    {kMFTTrackQ4plus, "q_{Extrap}-q_{Smoothed}"},
    {kMCTrackspT, "p_t [GeV]"},
    {kMCTracksp, "p [GeV]"},
    {kMCTrackEta, " \\eta"}
  };


 
  //Create histograms
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
    //gStyle->SetLineWidth(4);
    //gROOT->ForceStyle();
    h->GetXaxis()->SetTitle(TH2XaxisTitles[n2Histo]);
    h->GetYaxis()->SetTitle(TH2YaxisTitles[n2Histo]);
    //h->GetXaxis()->SetLabelSize(0.05);
    //h->GetXaxis()->SetTitleSize(0.05);
    //h->GetYaxis()->SetLabelSize(0.06);
    //h->GetYaxis()->SetTitleSize(0.06);
    h->SetOption("COLZ");
    ++n2Histo;
    }

  // Profiles histograms
  auto PRes_Profile = new TProfile("Pt_res_prof","Profile of p{Extrap}/p{Smoothed}",40,0,20,0,20, "s");
  PRes_Profile->GetXaxis()->SetTitle("p_{Smoothed}");
  PRes_Profile->GetYaxis()->SetTitle("mean(P_{Extrap}/P_{Smoothed})");

  auto DeltaX_Profile  = new TProfile("DeltaX_prof","Position resolution",40,0,20,-10000.,10000., "s");
  DeltaX_Profile->GetXaxis()->SetTitle("p_{Smoothed} [GeV]");
  DeltaX_Profile->GetYaxis()->SetTitle("\\sigma_x ~[\\mu m]");

  // TEfficiency histogram
  TEfficiency* qMatchEff = new TEfficiency("QMatchEff","Charge Match;p [GeV];#epsilon",10,0,10);
  //qMatchEff->GetPaintedHistogram()->GetXaxis()->SetLabelSize(0.06);
  //qMatchEff->GetPaintedHistogram()->GetYaxis()->SetLabelSize(0.06);
  //qMatchEff->GetPaintedHistogram()->GetXaxis()->SetTitleSize(0.06);
  //qMatchEff->GetPaintedHistogram()->GetYaxis()->SetTitleSize(0.06);

  // Counters
  Int_t nChargeMatch = 0;
  Int_t nChargeMiss = 0;
  Int_t nChargeMatch0_1 = 0;
  Int_t nChargeMiss0_1 = 0;
  Int_t nChargeMatch1_4 = 0;
  Int_t nChargeMiss1_4 = 0;
  Int_t nChargeMatch4plus = 0;
  Int_t nChargeMiss4plus = 0;
  Int_t nCleanTracksMFT = 0, nInvalidTracksMFT = 0, nMFTTrackable = 0;



  // Files & Trees
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
	//std::cout << " MotherID = " << thisTrack->getMotherTrackId() << std::endl; 
	// if (thisTrack->getMotherTrackId() == -1) { // Only primaries
          auto vx_MC = trackMFT.getXLast(); // MC variables storing smoothed parameters at last cluster 
          auto vy_MC = trackMFT.getYLast();
          auto vz_MC = trackMFT.getYLast();
          auto Pt_MC = trackMFT.getPtLast();
          auto P_MC = trackMFT.getPLast();
          auto phi_MC = trackMFT.getPhiLast();
          auto eta_MC = trackMFT.getEtaLast();
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
          Q_MC = TMath::Sign(1., trackMFT.getInvQPtLast());

          extrapMFTTrackHelixToZ(trackMFT, trackMFT.getZLast(), field_z); // propagate track to last cluster

          auto dx = trackMFT.getX() - vx_MC;
          auto dy = trackMFT.getY() - vy_MC;
          auto d_eta = trackMFT.getEta() - eta_MC;
          auto Pt_fit = trackMFT.getPt();
          auto Pt_fitLast = trackMFT.getPtLast();
          auto P_fit = trackMFT.getP();
          auto Q_fit = trackMFT.getCharge();
          auto P_res = P_fit / P_MC;
          auto Pt_res = Pt_fit / Pt_MC;
          auto d_Phi = trackMFT.getPhi() - phi_MC;
          auto d_Charge = Q_fit-Q_MC;

          TH1Histos[kMFTTracksP]->Fill(trackMFT.getP());
          //TH1Histos[kMFTTracksP_res]->Fill(P_res);
          //TH1Histos[kMFTTracksPt_res]->Fill(Pt_res);
          TH1Histos[kMFTTrackDeltaEta]->Fill(d_eta);
          TH1Histos[kMFTTrackDeltaPhi]->Fill(d_Phi);
          TH1Histos[kMFTTrackDeltaPhiDeg]->Fill(TMath::RadToDeg()*d_Phi);
          TH1Histos[kMFTTrackDeltaX]->Fill(dx);
          DeltaX_Profile->Fill(P_MC,dx*1e4);
          TH1Histos[kMFTTrackDeltaY]->Fill(dy);
          TH1Histos[kMFTTrackR]->Fill(sqrt(dx*dx+dy*dy));
          TH1Histos[kMFTTrackQ]->Fill(d_Charge);
          TH2Histos[kMFTTrackDeltaXYVertex]->Fill(dx,dy);
          TH2Histos[kMFTrackQPRec_MC]->Fill(P_MC*Q_MC,P_fit*Q_fit);
          TH2Histos[kMFTrackPtResolution]->Fill(Pt_MC,Pt_fit/Pt_MC);
          PRes_Profile->Fill(P_MC,P_fit/P_MC);
          TH2Histos[kMFTrackInvPtResolution]->Fill(Pt_MC,(1.0/Pt_fit-1.0/Pt_MC)*Pt_MC);

          // MC histos
          TH1Histos[kMCTrackspT]->Fill(Pt_MC);
          TH1Histos[kMCTracksp]->Fill(P_MC);
          TH1Histos[kMCTrackEta]->Fill(eta_MC);
          TH2Histos[kMCTracksEtaZ]->Fill(vz_MC,eta_MC);

          // Differential histos
          if (P_MC <= 1.0) {
               TH2Histos[kMFTTrackDeltaXYVertex0_1]->Fill(dx,dy);
               TH1Histos[kMFTTrackDeltaEta0_1]->Fill(d_eta);
               TH1Histos[kMFTTrackDeltaPhi0_1]->Fill(d_Phi);
               TH1Histos[kMFTTrackDeltaPhiDeg0_1]->Fill(TMath::RadToDeg()*d_Phi);
               TH1Histos[kMFTTrackDeltaX0_1]->Fill(dx);
               TH1Histos[kMFTTrackQ0_1]->Fill(d_Charge);
               d_Charge ? nChargeMiss0_1++ : nChargeMatch0_1++;
             }
          if (P_MC > 1.0 and P_MC <= 4 ) {
               TH2Histos[kMFTTrackDeltaXYVertex1_4]->Fill(dx,dy);
               TH1Histos[kMFTTrackDeltaEta1_4]->Fill(d_eta);
               TH1Histos[kMFTTrackDeltaPhi1_4]->Fill(d_Phi);
               TH1Histos[kMFTTrackDeltaPhiDeg1_4]->Fill(TMath::RadToDeg()*d_Phi);
               TH1Histos[kMFTTrackDeltaX1_4]->Fill(dx);
               TH1Histos[kMFTTrackQ1_4]->Fill(d_Charge);
               d_Charge ? nChargeMiss1_4++ : nChargeMatch1_4++;
             }
          if (P_MC > 4.0) {
             TH2Histos[kMFTTrackDeltaXYVertex4plus]->Fill(dx,dy);
             TH1Histos[kMFTTrackDeltaEta4plus]->Fill(d_eta);
             TH1Histos[kMFTTrackDeltaPhi4plus]->Fill(d_Phi);
             TH1Histos[kMFTTrackDeltaPhiDeg4plus]->Fill(TMath::RadToDeg()*d_Phi);
             TH1Histos[kMFTTrackDeltaX4plus]->Fill(dx);
             TH1Histos[kMFTTrackQ4plus]->Fill(d_Charge);
             d_Charge ? nChargeMiss4plus++ : nChargeMatch4plus++;
           }

           d_Charge ? nChargeMiss++ : nChargeMatch++;
	   qMatchEff->Fill(!d_Charge,P_MC);
	   //      }

        }
  } // Loop on TracksMFT


// Customize histograms
TH1Histos[kMFTTrackQ]->SetTitle(Form("nChargeMatch = %d (%.2f%%)", nChargeMatch, 100.*nChargeMatch/(nChargeMiss+nChargeMatch)));
TH1Histos[kMFTTrackQ0_1]->SetTitle(Form("nChargeMatch = %d (%.2f%%)", nChargeMatch0_1, 100.*nChargeMatch0_1/(nChargeMiss0_1+nChargeMatch0_1)));
TH1Histos[kMFTTrackQ1_4]->SetTitle(Form("nChargeMatch = %d (%.2f%%)", nChargeMatch1_4, 100.*nChargeMatch1_4/(nChargeMiss1_4+nChargeMatch1_4)));
TH1Histos[kMFTTrackQ4plus]->SetTitle(Form("nChargeMatch = %d (%.2f%%)", nChargeMatch4plus, 100.*nChargeMatch4plus/(nChargeMiss4plus+nChargeMatch4plus)));

qMatchEff->SetTitle(Form("Charge match = %.2f%%", 100.*nChargeMatch/(nChargeMiss+nChargeMatch)));


//Remove stat boxes
TH2Histos[kMFTrackQPRec_MC]->SetStats(0);
TH2Histos[kMFTrackPtResolution]->SetStats(0);
TH2Histos[kMFTrackInvPtResolution]->SetStats(0);
TH2Histos[kMCTracksEtaZ]->SetStats(0);
PRes_Profile->SetStats(0);
DeltaX_Profile->SetStats(0);
TH1Histos[kMFTTrackQ]->SetStats(0);


//Fit Slices: Pt resolution
FitSlicesy(*TH2Histos[kMFTrackInvPtResolution],*TH2Histos[kMFTrackQPRec_MC]);
FitSlicesy(*TH2Histos[kMFTrackPtResolution],*TH2Histos[kMFTrackQPRec_MC]);


// sigmaX resultion Profile
TH1D* DeltaX_Error = new  TH1D();
DeltaX_Error =  DeltaX_Profile->ProjectionX("DeltaX_Error", "C=E");
DeltaX_Error->Write();

// Summary Canvases
auto pt_resolution = summary_report(*TH2Histos[kMFTrackPtResolution],
				    *TH2Histos[kMFTrackQPRec_MC],
				    *PRes_Profile,
				    *qMatchEff,
				    "Pt Summary",
				    seed_cfg,
				    0, 0, 0, 0,
				    Form("%.2f%%", 100.0*TH2Histos[kMFTrackPtResolution]->Integral()/TH2Histos[kMFTrackPtResolution]->GetEntries()),
				    Form("%.2f%%", 100.0*TH2Histos[kMFTrackQPRec_MC]->Integral()/TH2Histos[kMFTrackQPRec_MC]->GetEntries())
				    );

auto invpt_resolution = summary_report(*TH2Histos[kMFTrackInvPtResolution],
				       *TH2Histos[kMFTrackQPRec_MC],
				       *(TH1F*)gDirectory->Get((std::string(TH2Histos[kMFTrackInvPtResolution]->GetName()) + std::string("_1")).c_str()),
				       *(TH1F*)gDirectory->Get((std::string(TH2Histos[kMFTrackInvPtResolution]->GetName()) + std::string("_2")).c_str()),
				       "InvPt Summary", seed_cfg,
				       0, 0, 0, 0,
				       Form("%.2f%%", 100.0*TH2Histos[kMFTrackInvPtResolution]->Integral()/TH2Histos[kMFTrackInvPtResolution]->GetEntries()),
				       Form("%.2f%%", 100.0*TH2Histos[kMFTrackQPRec_MC]->Integral()/TH2Histos[kMFTrackQPRec_MC]->GetEntries())
				       );

auto vertexing_resolution = summary_report(*TH2Histos[kMFTTrackDeltaXYVertex],
					   *TH1Histos[kMFTTrackDeltaX],
					   *DeltaX_Error,
					   *TH1Histos[kMFTTrackDeltaPhiDeg],
					   "Position Summary",
					   seed_cfg,
					   0, 1, 0, 1,
					   Form("%.2f%%", 100.0*TH2Histos[kMFTTrackDeltaXYVertex]->Integral()/TH2Histos[kMFTTrackDeltaXYVertex]->GetEntries()),
					   Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaX]->Integral()/TH1Histos[kMFTTrackDeltaX]->GetEntries()),
					   Form(" "),
					   Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaPhiDeg]->Integral()/TH1Histos[kMFTTrackDeltaPhiDeg]->GetEntries())
					   );
 
auto vertexing_resolution0_1 = summary_report(*TH2Histos[kMFTTrackDeltaXYVertex0_1],
					      *TH1Histos[kMFTTrackDeltaX0_1],
					      *TH1Histos[kMFTTrackDeltaEta0_1],
					      *TH1Histos[kMFTTrackDeltaPhiDeg0_1],
					      "Position Summary pt < 1",
					      seed_cfg,
					      0, 1, 1, 1,
					      Form("%.2f%%", 100.0*TH2Histos[kMFTTrackDeltaXYVertex0_1]->Integral()/TH2Histos[kMFTTrackDeltaXYVertex0_1]->GetEntries()),
					      Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaX0_1]->Integral()/TH1Histos[kMFTTrackDeltaX0_1]->GetEntries()),
					      Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaEta0_1]->Integral()/TH1Histos[kMFTTrackDeltaEta0_1]->GetEntries()),
					      Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaPhiDeg0_1]->Integral()/TH1Histos[kMFTTrackDeltaPhiDeg0_1]->GetEntries())
					      );

auto vertexing_resolution1_4 = summary_report(*TH2Histos[kMFTTrackDeltaXYVertex1_4],
					      *TH1Histos[kMFTTrackDeltaX1_4],
					      *TH1Histos[kMFTTrackDeltaEta1_4],
					      *TH1Histos[kMFTTrackDeltaPhiDeg1_4],
					      "Position Summary 1 < p_t < 4",
					      seed_cfg,
					      0, 1, 1, 1,
					      Form("%.2f%%", 100.0*TH2Histos[kMFTTrackDeltaXYVertex1_4]->Integral()/TH2Histos[kMFTTrackDeltaXYVertex1_4]->GetEntries()),
					      Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaX1_4]->Integral()/TH1Histos[kMFTTrackDeltaX1_4]->GetEntries()),
					      Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaEta1_4]->Integral()/TH1Histos[kMFTTrackDeltaEta1_4]->GetEntries()),
					      Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaPhiDeg1_4]->Integral()/TH1Histos[kMFTTrackDeltaPhiDeg1_4]->GetEntries())
					      );

auto vertexing_resolution4plus = summary_report(*TH2Histos[kMFTTrackDeltaXYVertex4plus],
						*TH1Histos[kMFTTrackDeltaX4plus],
						*TH1Histos[kMFTTrackDeltaEta4plus],
						*TH1Histos[kMFTTrackDeltaPhiDeg4plus],
						"Position Summary p_t > 4",
						seed_cfg,
						0, 1, 1, 1,
                        		        Form("%.2f%%", 100.0*TH2Histos[kMFTTrackDeltaXYVertex4plus]->Integral()/TH2Histos[kMFTTrackDeltaXYVertex4plus]->GetEntries()),
					        Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaX4plus]->Integral()/TH1Histos[kMFTTrackDeltaX4plus]->GetEntries()),
					        Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaEta4plus]->Integral()/TH1Histos[kMFTTrackDeltaEta4plus]->GetEntries()),
					        Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaPhiDeg4plus]->Integral()/TH1Histos[kMFTTrackDeltaPhiDeg4plus]->GetEntries())
						);

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

PRes_Profile->Write();
DeltaX_Profile->Write();
qMatchEff->Write(); 
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
std::cout << " Eta_mean = " << TH1Histos[kMFTTrackDeltaEta]->GetMean() << std::endl;
std::cout << " Eta_StdDev = " << TH1Histos[kMFTTrackDeltaEta]->GetStdDev() << std::endl;
std::cout << " Eta_StdDev(pt<1) = " << TH1Histos[kMFTTrackDeltaEta0_1]->GetStdDev() << std::endl;
std::cout << " Eta_StdDev(1<pt<4) = " << TH1Histos[kMFTTrackDeltaEta1_4]->GetStdDev() << std::endl;
std::cout << " Eta_StdDev(pt>4) = " << TH1Histos[kMFTTrackDeltaEta4plus]->GetStdDev() << std::endl;
std::cout << " Phi_mean = " << TH1Histos[kMFTTrackDeltaPhi]->GetMean() << std::endl;
std::cout << " Phi_StdDev = " << TH1Histos[kMFTTrackDeltaPhi]->GetStdDev() << std::endl;
std::cout << " Phi_StdDev(pt<1) = " << TH1Histos[kMFTTrackDeltaPhi0_1]->GetStdDev() << std::endl;
std::cout << " Phi_StdDev(1<pt<4) = " << TH1Histos[kMFTTrackDeltaPhi1_4]->GetStdDev() << std::endl;
std::cout << " Phi_StdDev(pt>4) = " << TH1Histos[kMFTTrackDeltaPhi4plus]->GetStdDev() << std::endl;
std::cout << " Phi_meanDeg = " << TH1Histos[kMFTTrackDeltaPhiDeg]->GetMean() << std::endl;
std::cout << " Phi_StdDevDeg = " << TH1Histos[kMFTTrackDeltaPhiDeg]->GetStdDev() << std::endl;
std::cout << " Phi_StdDevDeg(pt<1) = " << TH1Histos[kMFTTrackDeltaPhiDeg0_1]->GetStdDev() << std::endl;
std::cout << " Phi_StdDevDeg(1<pt<4) = " << TH1Histos[kMFTTrackDeltaPhiDeg1_4]->GetStdDev() << std::endl;
std::cout << " Phi_StdDevDeg(pt>4) = " << TH1Histos[kMFTTrackDeltaPhiDeg4plus]->GetStdDev() << std::endl;
std::cout << " DeltaX_mean = " << TH1Histos[kMFTTrackDeltaX]->GetMean() << std::endl;
std::cout << " DeltaX_StdDev = " << TH1Histos[kMFTTrackDeltaX]->GetStdDev() << std::endl;
std::cout << " DeltaX_StdDev(pt<1) = " << TH1Histos[kMFTTrackDeltaX0_1]->GetStdDev() << std::endl;
std::cout << " DeltaX_StdDev(1<pt<4) = " << TH1Histos[kMFTTrackDeltaX1_4]->GetStdDev() << std::endl;
std::cout << " DeltaX_StdDev(pt>4) = " << TH1Histos[kMFTTrackDeltaX4plus]->GetStdDev() << std::endl;
std::cout << " DeltaY_mean = " << TH1Histos[kMFTTrackDeltaY]->GetMean() << std::endl;
std::cout << " DeltaY_StdDev = " << TH1Histos[kMFTTrackDeltaY]->GetStdDev() << std::endl;
std::cout << " R_mean = " << TH1Histos[kMFTTrackR]->GetMean() << std::endl;
std::cout << " R_StdDev = " << TH1Histos[kMFTTrackR]->GetStdDev() << std::endl;
std::cout << " Charge_mean = " << TH1Histos[kMFTTrackDeltaY]->GetMean() << std::endl;
std::cout << " nChargeMatch = " << nChargeMatch << " (" << 100.*nChargeMatch/(nChargeMiss+nChargeMatch) << "%)" << std::endl;
std::cout << "---------------------------------------------------" << std::endl; 

return 0;
}




