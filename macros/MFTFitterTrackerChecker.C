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
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include "SimulationDataFormat/MCEventHeader.h"
#include <TStyle.h>
#include <TProfile.h>
#include <TGraph.h>

//constexpr Double_t MFTLayerZ[] = {-45.3, -46.7, -48.6, -50.0, -52.4, -53.8, -67.7, -69.1, -76.1, -77.5};
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
  double px0 = track.getPx();
  double py0 = track.getPy();
  double invtanl0 = 1.0 / track.getTanl();;
  double invqpt0 = track.getInvQPt();
  auto q = track.getCharge();
  auto Hz =  std::copysign(1, Field); 
  double k = TMath::Abs(3e-4 * Field);
  auto invk = 1.0 / k;
  double theta = -invqpt0 * dZ * k * invtanl0;
  double costheta, sintheta;
  o2::utils::sincos(theta, sintheta, costheta);
  double deltax = Hz * py0 * invk * (1.0 - costheta) - px0 * q * invk * sintheta;
  double deltay = -Hz * px0 * invk * (1.0 - costheta) - py0 * q * invk * sintheta;

  double x = x0 + deltax;
  double y = y0 + deltay;
  double phi = track.getPhi() + Hz * theta;

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
   c->Update();
   for (std::string type: {".pdf", ".png"}) c->Print((imgpath+std::string(h->GetName()) + type).c_str());
}


//_________________________________________________________________________________________________
void FitSlicesy( TH2F& histo1,  TH2F& histo2)
{
// ## FitSlicesy
TH2F *h1 = &histo1;
TH2F *h2 = &histo2;


h1->FitSlicesY(0,0,-1,1);
// Create a canvas and divide it
auto CanvasName = std::string(h1->GetName()) + std::string("FitSlicesY");
TCanvas *c1 = new TCanvas(CanvasName.c_str(),CanvasName.c_str(),700,500);
c1->Divide(2,1);
TPad *leftPad = (TPad*)c1->cd(1);;
leftPad->Divide(1,2);

// Draw 2-d original histogram histo1
leftPad->cd(1);
gPad->SetTopMargin(0.05);
h1->Draw();

// Draw histo2
leftPad->cd(2);
h2->Draw();
TPad *rightPad = (TPad*)c1->cd(2);
rightPad->Divide(1,2);
rightPad->cd(1);
TH2F *h1_Fit1 = (TH2F*)gDirectory->Get((std::string(h1->GetName()) + std::string("_1")).c_str());
h1_Fit1->SetStats(0);
h1_Fit1->SetTitle("Mean");
//h1_Fit1->GetXaxis()->SetLabelSize(0.05);
//h1_Fit1->GetXaxis()->SetTitleSize(0.05);
//h1_Fit1->GetYaxis()->SetLabelSize(0.05);
//h1_Fit1->GetYaxis()->SetTitleSize(0.05);
h1_Fit1->Draw();

// Show fitted "sigma" for each slice
rightPad->cd(2);
gPad->SetTopMargin(0.05);
gPad->SetLeftMargin(0.05);
TH2F *h1_Fit2 = (TH2F*)gDirectory->Get((std::string(h1->GetName()) + std::string("_2")).c_str());
h1_Fit2->SetStats(0);
h1_Fit2->SetTitle("Sigma");
//h1_Fit2->GetXaxis()->SetLabelSize(0.05);
//h1_Fit2->GetXaxis()->SetTitleSize(0.05);
//h1_Fit2->GetYaxis()->SetLabelSize(0.05);
//h1_Fit2->GetYaxis()->SetTitleSize(0.05);
h1_Fit2->Draw();


//h2_InvPtResolution_0->Write();
h1_Fit1->Write();
h1_Fit2->Write();
c1->Write();
}


//_________________________________________________________________________________________________
template <typename H1, typename H2, typename H3, typename H4>
TCanvas summary_report( H1& histo1,
			H2& histo2,
			H3& histo3,
			H4& histo4,
			std::string CanvasName,
			std::string tlt="Summary",
			int h1log=0,
			int h2log=0,
			int h3log=0,
			int h4log=0,
			std::string h1_foot="",
			std::string h2_foot="",
			std::string h3_foot="",
			std::string h4_foot="")
{
H1 *h1 = &histo1;
H2 *h2 = &histo2;
H3 *h3 = &histo3;
H4 *h4 = &histo4;


h1->FitSlicesY(0,0,-1,1);
// Create a canvas and divide it
TCanvas *c1 = new TCanvas(CanvasName.c_str(),CanvasName.c_str(),1920,1080);
c1->UseCurrentStyle();
//c1->SetCanvasSize(1920,1080);
//gROOT->SetStyle("Bold");

TLatex *Title = new TLatex() ;
Title->SetTextSize(0.035);
Title->SetTextAlign(23);
Title->DrawLatex(.5 ,.995, tlt.c_str());
Title->Draw();

c1->Divide(2,1);
TPad *leftPad = (TPad*)c1->cd(1);;
leftPad->SetPad(0.000,0.000,0.5,.96);
leftPad->Divide(1,2);

leftPad->cd(1);

gPad->SetBottomMargin(0.15);
gPad->SetRightMargin(0.15);
gPad->SetLogy(h1log);

 
h1->Draw();
h1->SetMarkerColor(4);
h1->SetMarkerStyle(21);
h1->SetMarkerSize(1); 
//h1->GetXaxis()->SetLabelSize(0.06);
//h1->GetXaxis()->SetTitleSize(0.05);
//h1->GetYaxis()->SetLabelSize(0.06);
//h1->GetYaxis()->SetTitleSize(0.05);
TLatex *h1_entries = new TLatex() ;
h1_entries->SetNDC();
h1_entries->SetTextSize(0.035);
h1_entries->SetTextAlign(23);
h1_entries->DrawLatex(.08 ,.08, h1_foot.c_str()) ;
h1_entries->Draw();



leftPad->cd(2);
gPad->SetBottomMargin(0.15);
gPad->SetTopMargin(0.10);
gPad->SetLogy(h2log);


h2->Draw();
h2->SetMarkerColor(4);
h2->SetMarkerStyle(21);
h2->SetMarkerSize(1);
TLatex *h2_entries = new TLatex() ;
h2_entries->SetNDC();
h2_entries->SetTextSize(0.035);
h2_entries->SetTextAlign(23);
h2_entries->DrawLatex(.08 ,.08, h2_foot.c_str()) ;
h2_entries->Draw();


TPad *rightPad = (TPad*)c1->cd(2);
rightPad->SetPad(0.5,0,1,0.97);
rightPad->Divide(1,2);
rightPad->cd(1);
gPad->SetBottomMargin(0.15);
gPad->SetTopMargin(0.10);
gPad->SetLogy(h3log);

h3->Draw();
h3->SetMarkerColor(4);
h3->SetMarkerStyle(21);
h3->SetMarkerSize(1);
TLatex *h3_entries = new TLatex() ;
h3_entries->SetNDC();
h3_entries->SetTextSize(0.035);
h3_entries->SetTextAlign(23);
h3_entries->DrawLatex(.08 ,.08, h3_foot.c_str()) ;
h3_entries->Draw();

rightPad->cd(2);
gPad->SetBottomMargin(0.15);
gPad->SetLogy(h4log);


 
h4->Draw();
h4->SetMarkerColor(4);
h4->SetMarkerStyle(21);
h4->SetMarkerSize(1);
TLatex *h4_entries = new TLatex() ;
h4_entries->SetNDC();
h4_entries->SetTextSize(0.035);
h4_entries->SetTextAlign(23);
h4_entries->DrawLatex(.08 ,.08, h4_foot.c_str()) ;
h4_entries->Draw();

c1->cd();




c1->Update();
c1->Write();
return c1;
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
    {kMFTTrackDeltaXYVertex0_1, "MFT Tracks Vertex at Z = 0 Pt0_1"},
    {kMFTTrackDeltaXYVertex1_4, "MFT Tracks Vertex at Z = 0 Pt1_4"},
    {kMFTTrackDeltaXYVertex4plus, "MFT Tracks Vertex at Z = 0 Pt4plus"},
    {kMFTrackQPRec_MC, "MFT Track QP FITxMC"},
    {kMFTrackPtResolution, "MFT Track Pt Resolution"},
    {kMFTrackInvPtResolution, "MFT Track InvPt Resolution"},
    {kMCTracksEtaZ, "MCTracks_eta_z"}
  };

  std::map<int,const char *> TH2Titles {
    {kMFTTrackDeltaXYVertex, "MFT Tracks at Z_vertex"},
    {kMFTTrackDeltaXYVertex0_1, "MFT Tracks at Z_vertex (pt < 1)"},
    {kMFTTrackDeltaXYVertex1_4, "MFT Tracks at Z_vertex (1 < pt < 4)"},
    {kMFTTrackDeltaXYVertex4plus, "MFT Tracks at Z_vertex (pt > 4)"},
    {kMFTrackQPRec_MC, "Charged Momentum: Reconstructed vs MC"},
    {kMFTrackPtResolution, "Pt Resolution"},
    {kMFTrackInvPtResolution, "InvPt Resolution"},
    {kMCTracksEtaZ, "MC Tracks: Pseudorapidity vs zVertex"}
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
    {kMFTrackQPRec_MC, "(q.p)_{MC} [GeV]"},
    {kMFTrackPtResolution, "pt_{MC} [GeV]"},
    {kMFTrackInvPtResolution, "pt_{MC} [GeV]"},
    {kMCTracksEtaZ, "Vertex PosZ [cm]"}
  };

  std::map<int,const char *> TH2YaxisTitles {
    {kMFTTrackDeltaXYVertex, "\\Delta y ~[cm]"},
    {kMFTTrackDeltaXYVertex0_1, "\\Delta y ~[cm]"},
    {kMFTTrackDeltaXYVertex1_4, "\\Delta y ~[cm]"},
    {kMFTTrackDeltaXYVertex4plus, "\\Delta y ~[cm]"},
    {kMFTrackQPRec_MC, "(q.p)_{fit} [GeV]"},
    {kMFTrackPtResolution, "pt_{fit} / pt_{MC}"},
    {kMFTrackInvPtResolution, "(1/(p_t)_{fit} - 1/(p_t)_{MC})*(p_t)_{MC}"},
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
    {kMFTTrackDeltaEta0_1, "MFT Tracks eta (pt < 1)"},
    {kMFTTrackDeltaEta1_4, "MFT Tracks eta (1 < pt < 4)"},
    {kMFTTrackDeltaEta4plus, "MFT Tracks eta (pt > 4)"},
    {kMFTTrackDeltaPhi, "MFT Tracks Fitted Phi at Vertex"},
    {kMFTTrackDeltaPhi0_1,"MFT Tracks Fitted Phi at Vertex [rad] (pt < 1)"},
    {kMFTTrackDeltaPhi1_4,"MFT Tracks Fitted Phi at Vertex [rad] (1 < pt < 4)"},
    {kMFTTrackDeltaPhi4plus,"MFT Tracks Fitted Phi at Vertex [rad] (pt > 4)"},
    {kMFTTrackDeltaPhiDeg,"MFT Tracks Fitted Phi at Vertex [deg]"},
    {kMFTTrackDeltaPhiDeg0_1,"MFT Tracks Fitted Phi at Vertex [deg] (pt < 1)"},
    {kMFTTrackDeltaPhiDeg1_4,"MFT Tracks Fitted Phi at Vertex [deg] (1 < pt < 4)"},
    {kMFTTrackDeltaPhiDeg4plus,"MFT Tracks Fitted Phi at Vertex [deg] (pt > 4)"},
    {kMFTTrackDeltaX, "MFT Tracks Delta X"},
    {kMFTTrackDeltaX0_1, "MFT Tracks Delta X (pt < 1)"},
    {kMFTTrackDeltaX1_4, "MFT Tracks Delta X (1 < pt < 4)"},
    {kMFTTrackDeltaX4plus, "MFT Tracks Delta X (pt > 4)"},
    {kMFTTrackDeltaY, "MFT Tracks Delta Y"},
    {kMFTTrackR, "MFT Tracks Delta R"},
    {kMFTTrackQ, "Charge Match"},
    {kMFTTrackQ0_1, "Charge Match (pt < 1)"},
    {kMFTTrackQ1_4, "Charge Match (1 < pt < 4)"},
    {kMFTTrackQ4plus, "Charge Match (pt > 4)"},
    {kMCTrackspT, "MC Tracks p_T"},
    {kMCTracksp, "MC Tracks p"},
    {kMCTrackEta, "MC Tracks eta"}
  };

  std::map<int,const char *> TH1Titles {
    {kMFTTracksP, "Standalone MFT Tracks P"},
    //{kMFTTracksP_res, "P_{Fit}/P_{MC}"},
    //{kMFTTracksPt_res,"Pt_{Fit}/Pt_{MC}" },
    {kMFTTrackDeltaEta, "\\eta_{Fit} - \\eta_{MC} "},
    {kMFTTrackDeltaEta0_1, "\\eta_{Fit} - \\eta_{MC} "},
    {kMFTTrackDeltaEta1_4, "\\eta_{Fit} - \\eta_{MC} "},
    {kMFTTrackDeltaEta4plus, "\\eta_{Fit} - \\eta_{MC} "},
    {kMFTTrackDeltaPhi, "\\phi _{Fit} - \\phi_{MC}"},
    {kMFTTrackDeltaPhi0_1, "\\phi _{Fit} - \\phi_{MC}"},
    {kMFTTrackDeltaPhi1_4, "\\phi _{Fit} - \\phi_{MC}"},
    {kMFTTrackDeltaPhi4plus, "\\phi _{Fit} - \\phi_{MC}"},
    {kMFTTrackDeltaPhiDeg, "\\phi _{Fit} - \\phi_{MC}"},
    {kMFTTrackDeltaPhiDeg0_1, "\\phi _{Fit} - \\phi_{MC}"},
    {kMFTTrackDeltaPhiDeg1_4, "\\phi _{Fit} - \\phi_{MC}"},
    {kMFTTrackDeltaPhiDeg4plus, "\\phi _{Fit} - \\phi_{MC}"},
    {kMFTTrackDeltaX, "MFT Tracks Delta X at Z_vertex"},
    {kMFTTrackDeltaX0_1, "MFT Tracks Delta X at Z_vertex"},
    {kMFTTrackDeltaX1_4, "MFT Tracks Delta X at Z_vertex"},
    {kMFTTrackDeltaX4plus, "MFT Tracks Delta X at Z_vertex"},
    {kMFTTrackDeltaY, "MFT Tracks Delta Y at Z_vertex"},
    {kMFTTrackR, "MFT Tracks Delta R at Z_vertex"},
    {kMFTTrackQ, "MFT Tracks Charge Match"},
    {kMFTTrackQ0_1, "MFT Tracks Charge Match (pt < 1)"},
    {kMFTTrackQ1_4, "MFT Tracks Charge Match (1 < pt < 4)"},
    {kMFTTrackQ4plus, "MFT Tracks Charge Match (pt > 4)"},
    {kMCTrackspT, "MC Tracks p_T"},
    {kMCTracksp, "MC Tracks p"},
    {kMCTrackEta, "MC Tracks Pseudorapidity"}
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
    //{kMFTTracksP_res, "P_{Fit}/P_{MC}"},
    //{kMFTTracksPt_res, "Pt_{Fit}/Pt_{MC}"},
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
    {kMFTTrackQ, "q_{fit}-q_{MC}"},
    {kMFTTrackQ0_1, "q_{fit}-q_{MC}"},
    {kMFTTrackQ1_4, "q_{fit}-q_{MC}"},
    {kMFTTrackQ4plus, "q_{fit}-q_{MC}"},
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
  auto PtRes_Profile  = new TProfile("Pt_res_prof","Profile of pt{fit}/pt{MC}",14,0,7,0,20, "s");
  PtRes_Profile->GetXaxis()->SetTitle("pt_{MC}");
  PtRes_Profile->GetYaxis()->SetTitle("mean(Pt_{Fit}/Pt_{MC})");

  auto DeltaX_Profile  = new TProfile("DeltaX_prof","Vertexing resolution",14,0,7,-10000.,10000., "s");
  DeltaX_Profile->GetXaxis()->SetTitle("pt_{MC} [GeV]");
  DeltaX_Profile->GetYaxis()->SetTitle("\\sigma_x ~[\\mu m]");

  // TEfficiency histogram
  TEfficiency* qMatchEff = new TEfficiency("QMatchEff","Charge Match;p_t [GeV];#epsilon",10,0,10);
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
          //TH1Histos[kMFTTracksP_res]->Fill(P_res);
          //TH1Histos[kMFTTracksPt_res]->Fill(Pt_res);
          TH1Histos[kMFTTrackDeltaEta]->Fill(d_eta);
          TH1Histos[kMFTTrackDeltaPhi]->Fill(d_Phi);
          TH1Histos[kMFTTrackDeltaPhiDeg]->Fill(TMath::RadToDeg()*d_Phi);
          TH1Histos[kMFTTrackDeltaX]->Fill(dx);
          DeltaX_Profile->Fill(Pt_MC,dx*1e4);
          TH1Histos[kMFTTrackDeltaY]->Fill(dy);
          TH1Histos[kMFTTrackR]->Fill(sqrt(dx*dx+dy*dy));
          TH1Histos[kMFTTrackQ]->Fill(d_Charge);
          TH2Histos[kMFTTrackDeltaXYVertex]->Fill(dx,dy);
          TH2Histos[kMFTrackQPRec_MC]->Fill(P_MC*Q_MC,P_fit*Q_fit);
          TH2Histos[kMFTrackPtResolution]->Fill(Pt_MC,Pt_fit/Pt_MC);
          PtRes_Profile->Fill(Pt_MC,Pt_fit/Pt_MC);
          TH2Histos[kMFTrackInvPtResolution]->Fill(Pt_MC,(1.0/Pt_fit-1.0/Pt_MC)*Pt_MC);

          // MC histos
          TH1Histos[kMCTrackspT]->Fill(Pt_MC);
          TH1Histos[kMCTracksp]->Fill(P_MC);
          TH1Histos[kMCTrackEta]->Fill(eta_MC);
          TH2Histos[kMCTracksEtaZ]->Fill(vz_MC,eta_MC);

          // Differential histos
          if (Pt_MC <= 1.0) {
               TH2Histos[kMFTTrackDeltaXYVertex0_1]->Fill(dx,dy);
               TH1Histos[kMFTTrackDeltaEta0_1]->Fill(d_eta);
               TH1Histos[kMFTTrackDeltaPhi0_1]->Fill(d_Phi);
               TH1Histos[kMFTTrackDeltaPhiDeg0_1]->Fill(TMath::RadToDeg()*d_Phi);
               TH1Histos[kMFTTrackDeltaX0_1]->Fill(dx);
               TH1Histos[kMFTTrackQ0_1]->Fill(d_Charge);
               d_Charge ? nChargeMiss0_1++ : nChargeMatch0_1++;
             }
          if (Pt_MC > 1.0 and Pt_MC <= 4 ) {
               TH2Histos[kMFTTrackDeltaXYVertex1_4]->Fill(dx,dy);
               TH1Histos[kMFTTrackDeltaEta1_4]->Fill(d_eta);
               TH1Histos[kMFTTrackDeltaPhi1_4]->Fill(d_Phi);
               TH1Histos[kMFTTrackDeltaPhiDeg1_4]->Fill(TMath::RadToDeg()*d_Phi);
               TH1Histos[kMFTTrackDeltaX1_4]->Fill(dx);
               TH1Histos[kMFTTrackQ1_4]->Fill(d_Charge);
               d_Charge ? nChargeMiss1_4++ : nChargeMatch1_4++;
             }
          if (Pt_MC > 4.0) {
             TH2Histos[kMFTTrackDeltaXYVertex4plus]->Fill(dx,dy);
             TH1Histos[kMFTTrackDeltaEta4plus]->Fill(d_eta);
             TH1Histos[kMFTTrackDeltaPhi4plus]->Fill(d_Phi);
             TH1Histos[kMFTTrackDeltaPhiDeg4plus]->Fill(TMath::RadToDeg()*d_Phi);
             TH1Histos[kMFTTrackDeltaX4plus]->Fill(dx);
             TH1Histos[kMFTTrackQ4plus]->Fill(d_Charge);
             d_Charge ? nChargeMiss4plus++ : nChargeMatch4plus++;
           }

           d_Charge ? nChargeMiss++ : nChargeMatch++;
	   qMatchEff->Fill(!d_Charge,Pt_MC);
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
PtRes_Profile->SetStats(0);
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
				    *PtRes_Profile,
				    *qMatchEff,
				    "Pt Summary",
				    seed_cfg,
				    0, 0, 0, 0,
				    Form("%.2f%%", 100.0*TH2Histos[kMFTrackPtResolution]->Integral()/TH2Histos[kMFTrackPtResolution]->GetEntries()),
				    Form("%.2f%%", 100.0*TH2Histos[kMFTrackQPRec_MC]->Integral()/TH2Histos[kMFTrackQPRec_MC]->GetEntries())
				    );
auto errcnt = 0;
std::cout << errcnt++ << std::endl;
auto invpt_resolution = summary_report(*TH2Histos[kMFTrackInvPtResolution],
				       *TH2Histos[kMFTrackQPRec_MC],
				       *(TH1F*)gDirectory->Get((std::string(TH2Histos[kMFTrackInvPtResolution]->GetName()) + std::string("_1")).c_str()),
				       *(TH1F*)gDirectory->Get((std::string(TH2Histos[kMFTrackInvPtResolution]->GetName()) + std::string("_2")).c_str()),
				       "InvPt Summary", seed_cfg,
				       0, 0, 0, 0,
				       Form("%.2f%%", 100.0*TH2Histos[kMFTrackInvPtResolution]->Integral()/TH2Histos[kMFTrackInvPtResolution]->GetEntries()),
				       Form("%.2f%%", 100.0*TH2Histos[kMFTrackQPRec_MC]->Integral()/TH2Histos[kMFTrackQPRec_MC]->GetEntries())
				       );

 std::cout << errcnt++ << std::endl;
auto vertexing_resolution = summary_report(*TH2Histos[kMFTTrackDeltaXYVertex],
					   *TH1Histos[kMFTTrackDeltaX],
					   *DeltaX_Error,
					   *TH1Histos[kMFTTrackDeltaPhiDeg],
					   "Vertexing Summary",
					   seed_cfg,
					   0, 1, 0, 1,
					   Form("%.2f%%", 100.0*TH2Histos[kMFTTrackDeltaXYVertex]->Integral()/TH2Histos[kMFTTrackDeltaXYVertex]->GetEntries()),
					   Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaX]->Integral()/TH1Histos[kMFTTrackDeltaX]->GetEntries()),
					   Form(" "),
					   Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaPhiDeg]->Integral()/TH1Histos[kMFTTrackDeltaPhiDeg]->GetEntries())
					   );
 
 std::cout << errcnt++ << std::endl;
auto vertexing_resolution0_1 = summary_report(*TH2Histos[kMFTTrackDeltaXYVertex0_1],
					      *TH1Histos[kMFTTrackDeltaX0_1],
					      *TH1Histos[kMFTTrackDeltaEta0_1],
					      *TH1Histos[kMFTTrackDeltaPhiDeg0_1],
					      "Vertexing Summary pt < 1",
					      seed_cfg,
					      0, 1, 1, 1,
					      Form("%.2f%%", 100.0*TH2Histos[kMFTTrackDeltaXYVertex0_1]->Integral()/TH2Histos[kMFTTrackDeltaXYVertex0_1]->GetEntries()),
					      Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaX0_1]->Integral()/TH1Histos[kMFTTrackDeltaX0_1]->GetEntries()),
					      Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaEta0_1]->Integral()/TH1Histos[kMFTTrackDeltaEta0_1]->GetEntries()),
					      Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaPhiDeg0_1]->Integral()/TH1Histos[kMFTTrackDeltaPhiDeg0_1]->GetEntries())
					      );


 std::cout << errcnt++ << std::endl;
auto vertexing_resolution1_4 = summary_report(*TH2Histos[kMFTTrackDeltaXYVertex1_4],
					      *TH1Histos[kMFTTrackDeltaX1_4],
					      *TH1Histos[kMFTTrackDeltaEta1_4],
					      *TH1Histos[kMFTTrackDeltaPhiDeg1_4],
					      "Vertexing Summary 1 < p_t < 4",
					      seed_cfg,
					      0, 1, 1, 1,
					      Form("%.2f%%", 100.0*TH2Histos[kMFTTrackDeltaXYVertex1_4]->Integral()/TH2Histos[kMFTTrackDeltaXYVertex1_4]->GetEntries()),
					      Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaX1_4]->Integral()/TH1Histos[kMFTTrackDeltaX1_4]->GetEntries()),
					      Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaEta1_4]->Integral()/TH1Histos[kMFTTrackDeltaEta1_4]->GetEntries()),
					      Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaPhiDeg1_4]->Integral()/TH1Histos[kMFTTrackDeltaPhiDeg1_4]->GetEntries())
					      );

 std::cout << errcnt++ << std::endl;
auto vertexing_resolution4plus = summary_report(*TH2Histos[kMFTTrackDeltaXYVertex4plus],
						*TH1Histos[kMFTTrackDeltaX4plus],
						*TH1Histos[kMFTTrackDeltaEta4plus],
						*TH1Histos[kMFTTrackDeltaPhiDeg4plus],
						"Vertexing Summary p_t > 4",
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

PtRes_Profile->Write();
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




