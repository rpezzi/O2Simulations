#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TLatex.h"
#include "TEfficiency.h"
#include <TMath.h>
#include "CommonConstants/MathConstants.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "ITSMFTSimulation/Hit.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "MFTTracking/TrackCA.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include <TGeoGlobalMagField.h>
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include <TH1F.h>
#include <TH2F.h>
#include "SimulationDataFormat/MCEventHeader.h"
#include <TStyle.h>
#include <TProfile.h>
#include <TGraph.h>

#endif

// MFTTools

#include "mfttools/MagField.C"
#include "mfttools/HistosHelpers.C"



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
int MFTFitterTrackerChecker( const Char_t *trkFile = "mfttracks.root",
                              const Char_t *o2sim_KineFile = "o2sim_Kine.root",
                              const Char_t *HitsMFTFile = "o2sim_HitsMFT.root"
                            )
{

// Histos parameters
Double_t pMin = 0.0;
Double_t pMax = 100.0;
Double_t deltaetaMin = -.1;
Double_t deltaetaMax = +.1;
Double_t etaMin = -3.4;
Double_t etaMax = -2.4;
Double_t deltaphiMin = -.2; //-3.15,
Double_t deltaphiMax = .2; //+3.15,
Double_t deltatanlMin = -2.0;
Double_t deltatanlMax = 2.0;

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
 gStyle->SetOptStat("emr");
 gStyle->SetStatW(.28);
 gStyle->SetStatH(.26);
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
  kMFTTrackChi2vsFitChi2,
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
	   {kMFTTrackChi2vsFitChi2, "MFT TracksChi2vsFitChi2"},
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
	   {kMFTTrackChi2vsFitChi2, "Tracks Chi2 vs FitChi2"},
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
	  {kMFTTrackChi2vsFitChi2, {500, 0, 1000, 250, 0., 500.}},
	    {kMFTrackQPRec_MC, {100, -10, 10, 100, -10, 10} },
	      {kMFTrackPtResolution, {14, 0, 7, 250, 0, 25} },
		{kMFTrackInvPtResolution, {14, 0, 7, 300, -2, 2} },
		  {kMCTracksEtaZ, {31, -15, 16, 25, etaMin, etaMax} }
};


 std::map<int,const char *> TH2XaxisTitles {
   {kMFTTrackDeltaXYVertex, "\\Delta x ~[cm]"},
     {kMFTTrackDeltaXYVertex0_1, "\\Delta x ~[cm]"},
       {kMFTTrackDeltaXYVertex1_4, "\\Delta x ~[cm]"},
	 {kMFTTrackDeltaXYVertex4plus, "\\Delta x ~[cm]"},
	   {kMFTTrackChi2vsFitChi2, "Fit ~ \\chi^2"},
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
	   {kMFTTrackChi2vsFitChi2, "Track ~ \\chi^2"},
	     {kMFTrackQPRec_MC, "(q.p)_{fit} [GeV]"},
	       {kMFTrackPtResolution, "pt_{fit} / pt_{MC}"},
		 {kMFTrackInvPtResolution, "(1/(p_t)_{fit} - 1/(p_t)_{MC})*(p_t)_{MC}"},
		   {kMCTracksEtaZ, "\\eta"}
 };


 enum TH1HistosCodes {
   kMFTTrackDeltaXErr,
   kMFTTrackDeltaYErr,
   kMFTTrackDeltaPhiErr,
   kMFTTrackDeltaTanLErr,
   kMFTTrackDeltainvQPtErr,
   kMFTTrackXChi2,
   kMFTTrackYChi2,
   kMFTTrackPhiChi2,
   kMFTTrackTanlChi2,
   kMFTTrackinvQPtChi2,
   kFitChi2,
   kMFTTracksP,
   kMFTTrackDeltaTanl,
   kMFTTrackDeltaTanl0_1,
   kMFTTrackDeltaTanl1_4,
   kMFTTrackDeltaTanl4plus,
   kMFTTrackDeltaPhi,
   kMFTTrackDeltaPhi0_1,
   kMFTTrackDeltaPhi1_4,
   kMFTTrackDeltaPhi4plus,
   kMFTTrackDeltaPhiDeg,
   kMFTTrackDeltaPhiDeg0_1,
   kMFTTrackDeltaPhiDeg1_4,
   kMFTTrackDeltaPhiDeg4plus,
   kMFTTrackDeltaInvQPt,
   kMFTTrackDeltaInvQPtSeed,
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
   kMFTTrackChi2,
   kMCTrackspT,
   kMCTracksp,
   kMCTrackEta
 };


 

 std::map<int,const char *> TH1Names {
   {kMFTTracksP, "MFT Tracks Fitted p"},
     {kMFTTrackDeltaXErr, "Delta X / SigmaX"},
       {kMFTTrackDeltaYErr, "Delta Y / SigmaY"},
	 {kMFTTrackDeltaPhiErr, "Delta Phi at Vertex / SigmaPhi"},
	   {kMFTTrackDeltaTanLErr, "Delta_Tanl / SigmaTanl"},
	     {kMFTTrackDeltainvQPtErr, "Delta_InvQPt / Sigma_{q/pt}"},
	       {kMFTTrackDeltaTanl, "MFT Tracks Fitted Delta_tanl"},
		 {kMFTTrackXChi2, "X Chi2"},
		   {kMFTTrackYChi2, "Y Chi2"},
		     {kMFTTrackPhiChi2, "Phi chi2"},
		       {kMFTTrackTanlChi2, "Tanl Chi2"},
			 {kMFTTrackinvQPtChi2, "InvQPt Chi2"},
			   {kFitChi2, "Fit Chi2"},
			     {kMFTTrackDeltaTanl0_1, "MFT Tracks tanl (pt < 1)"},
			       {kMFTTrackDeltaTanl1_4, "MFT Tracks tanl (1 < pt < 4)"},
				 {kMFTTrackDeltaTanl4plus, "MFT Tracks tanl (pt > 4)"},
				   {kMFTTrackDeltaPhi, "MFT Tracks Fitted Phi at Vertex"},
				     {kMFTTrackDeltaPhi0_1,"MFT Tracks Fitted Phi at Vertex [rad] (pt < 1)"},
				       {kMFTTrackDeltaPhi1_4,"MFT Tracks Fitted Phi at Vertex [rad] (1 < pt < 4)"},
					 {kMFTTrackDeltaPhi4plus,"MFT Tracks Fitted Phi at Vertex [rad] (pt > 4)"},
					   {kMFTTrackDeltaPhiDeg,"MFT Tracks Fitted Phi at Vertex [deg]"},
					     {kMFTTrackDeltaPhiDeg0_1,"MFT Tracks Fitted Phi at Vertex [deg] (pt < 1)"},
					       {kMFTTrackDeltaPhiDeg1_4,"MFT Tracks Fitted Phi at Vertex [deg] (1 < pt < 4)"},
						 {kMFTTrackDeltaPhiDeg4plus,"MFT Tracks Fitted Phi at Vertex [deg] (pt > 4)"},
						   {kMFTTrackDeltaInvQPt, "MFT Tracks invQPt"},
						     {kMFTTrackDeltaInvQPtSeed, "MFT Tracks invQPt Seed"},
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
									   {kMFTTrackChi2, "Tracks Chi2"},
									     {kMCTrackspT, "MC Tracks p_T"},
									       {kMCTracksp, "MC Tracks p"},
										 {kMCTrackEta, "MC Tracks eta"}
 };

 std::map<int,const char *> TH1Titles {
   {kMFTTracksP, "Standalone MFT Tracks P"},
				       {kMFTTrackDeltaXErr, "\\Delta X / \\sigma_X"},
					 {kMFTTrackDeltaYErr, "\\Delta Y / \\sigma_Y"},
					   {kMFTTrackDeltaPhiErr, "(\\Delta \\phi / \\sigma_\\phi"},
					     {kMFTTrackDeltaTanLErr, "\\Delta TanL / \\sigma_{TanL} "},
					       {kMFTTrackDeltainvQPtErr, "\\Delta(q/Pt) / \\sigma_{q/pt}"},
						 {kMFTTrackXChi2, "\\chi^2(x)"},
						   {kMFTTrackYChi2, "\\chi^2(y)"},
						     {kMFTTrackPhiChi2, "\\chi^2(\\phi)"},
						       {kMFTTrackTanlChi2, "\\chi^2(TanL)"},
							 {kMFTTrackinvQPtChi2, "\\chi^2(InvQP_t)"},
							   {kFitChi2, "Fit Chi2"},
							     {kMFTTrackDeltaTanl, "tanl_{Fit} - tanl_{MC} "},
							       {kMFTTrackDeltaTanl0_1, "tanl_{Fit} - tanl_{MC} (pt < 1)"},
								 {kMFTTrackDeltaTanl1_4, "tanl_{Fit} - tanl_{MC} (1 < p_t < 4)"},
								   {kMFTTrackDeltaTanl4plus, "tanl_{Fit} - tanl_{MC} (p_t > 4)"},
								     {kMFTTrackDeltaPhi, "\\phi _{Fit} - \\phi_{MC}"},
								       {kMFTTrackDeltaPhi0_1, "\\phi _{Fit} - \\phi_{MC}"},
									 {kMFTTrackDeltaPhi1_4, "\\phi _{Fit} - \\phi_{MC}"},
									   {kMFTTrackDeltaPhi4plus, "\\phi _{Fit} - \\phi_{MC}"},
									     {kMFTTrackDeltaPhiDeg, "\\phi _{Fit} - \\phi_{MC}"},
									       {kMFTTrackDeltaPhiDeg0_1, "\\phi _{Fit} - \\phi_{MC}"},
										 {kMFTTrackDeltaPhiDeg1_4, "\\phi _{Fit} - \\phi_{MC}"},
										   {kMFTTrackDeltaPhiDeg4plus, "\\phi _{Fit} - \\phi_{MC}"},
										     {kMFTTrackDeltaInvQPt, "MFT Tracks \\Delta invQPt"},
										       {kMFTTrackDeltaInvQPtSeed, "MFT Tracks \\Delta invQPt Seed"},
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
													     {kMFTTrackChi2, "MFT Tracks ~ \\chi^2"},
													       {kMCTrackspT, "MC Tracks p_T"},
														 {kMCTracksp, "MC Tracks p"},
														   {kMCTrackEta, "MC Tracks Pseudorapidity"}
 };

 std::map<int, std::array<double,3>> TH1Binning {
   {kMFTTracksP, {500, pMin, pMax} },
     {kMFTTrackDeltaXErr, {500, -10, 10}},
       {kMFTTrackDeltaYErr, {500, -10, 10}},
	 {kMFTTrackDeltaPhiErr, {500, -10, +10}},
	   {kMFTTrackDeltaTanLErr, {500, -10, +10}},
	     {kMFTTrackDeltainvQPtErr, {500, -50, +50}},
	       {kMFTTrackXChi2, {500, 0, 100}},
		 {kMFTTrackYChi2, {500, 0, 100}},
		   {kMFTTrackPhiChi2, {500, 0, 100}},
		     {kMFTTrackTanlChi2, {500, 0, 100}},
		       {kMFTTrackinvQPtChi2, {500, 0, 100}},
			 {kFitChi2, {500, 0, 50}},
			   {kMFTTrackDeltaTanl, {1000, deltatanlMin, deltatanlMax}},
			     {kMFTTrackDeltaTanl0_1, {1000, deltatanlMin, deltatanlMax}},
			       {kMFTTrackDeltaTanl1_4, {1000, deltatanlMin, deltatanlMax}},
				 {kMFTTrackDeltaTanl4plus, {1000, deltatanlMin, deltatanlMax}},
				   {kMFTTrackDeltaPhi, {1000, deltaphiMin, deltaphiMax}},
				     {kMFTTrackDeltaPhi0_1, {1000, deltaphiMin, deltaphiMax}},
				       {kMFTTrackDeltaPhi1_4, {1000, deltaphiMin, deltaphiMax}},
					 {kMFTTrackDeltaPhi4plus, {1000, deltaphiMin, deltaphiMax}},
					   {kMFTTrackDeltaPhiDeg, {1000, TMath::RadToDeg()*deltaphiMin, TMath::RadToDeg()*deltaphiMax}},
					     {kMFTTrackDeltaPhiDeg0_1, {1000, TMath::RadToDeg()*deltaphiMin, TMath::RadToDeg()*deltaphiMax}},
					       {kMFTTrackDeltaPhiDeg1_4, {1000, TMath::RadToDeg()*deltaphiMin, TMath::RadToDeg()*deltaphiMax}},
						 {kMFTTrackDeltaPhiDeg4plus, {1000, TMath::RadToDeg()*deltaphiMin, TMath::RadToDeg()*deltaphiMax}},
						   {kMFTTrackDeltaInvQPt, {1000, -10., 10.}},
						     {kMFTTrackDeltaInvQPtSeed, {1000, -0.5, 0.5}},
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
									   {kMFTTrackChi2, {10000, 0, 1000}},
									     {kMCTrackspT, {5000, 0, 50}},
									       {kMCTracksp, {1000, pMin, pMax}},
										 {kMCTrackEta, {1000, etaMin, etaMax}}
 };

 std::map<int,const char *> TH1XaxisTitles {
   {kMFTTracksP, "p [GeV]"},
     {kMFTTrackDeltaXErr, "\\Delta x  /\\sigma_{x}"},
       {kMFTTrackDeltaYErr, "\\Delta y  /\\sigma_{y}"},
	 {kMFTTrackDeltaPhiErr, "\\Delta \\phi  /\\sigma_{\\phi}"},
	   {kMFTTrackDeltaTanLErr, "\\Delta tanl /\\sigma_{tanl}"},
	     {kMFTTrackDeltainvQPtErr, "\\Delta (q/p_t)/\\sigma_{q/Pt}"},
	       {kMFTTrackDeltaTanl, "\\Delta tanl"},
		 {kMFTTrackXChi2, "\\chi^2"},
		   {kMFTTrackYChi2, "\\chi^2"},
		     {kMFTTrackPhiChi2, "\\chi^2"},
		       {kMFTTrackTanlChi2, "\\chi^2"},
			 {kMFTTrackinvQPtChi2, "\\chi^2"},
			   {kFitChi2, "\\chi^2"},
			     {kMFTTrackDeltaTanl0_1, "\\Delta tanl"},
			       {kMFTTrackDeltaTanl1_4, "\\Delta tanl"},
				 {kMFTTrackDeltaTanl4plus, "\\Delta tanl"},
				   {kMFTTrackDeltaPhi, "\\Delta \\phi ~[rad]"},
				     {kMFTTrackDeltaPhi0_1, "\\Delta \\phi ~[rad]"},
				       {kMFTTrackDeltaPhi1_4, "\\Delta \\phi ~[rad]"},
					 {kMFTTrackDeltaPhi4plus, "\\Delta \\phi ~[rad]"},
					   {kMFTTrackDeltaPhiDeg, "\\Delta \\phi ~[deg]"},
					     {kMFTTrackDeltaPhiDeg0_1, "\\Delta \\phi ~[deg]"},
					       {kMFTTrackDeltaPhiDeg1_4, "\\Delta \\phi ~[deg]"},
						 {kMFTTrackDeltaPhiDeg4plus, "\\Delta \\phi ~[deg]"},
						   {kMFTTrackDeltaInvQPt, "\\Delta invQPt"},
						     {kMFTTrackDeltaInvQPtSeed, "\\Delta invQPt Seed"},
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
									   {kMFTTrackChi2, "\\chi^2"},
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

 o2::dataformats::MCTruthContainer<o2::MCCompLabel>* mcLabels = nullptr;
 mftTrackTree -> SetBranchAddress("MFTTrackMCTruth",&mcLabels);

 mftTrackTree -> GetEntry(0);
 o2SimKineTree -> GetEntry(0);

  
 auto field_z = getZField(0, 0, -61.4); // Get field at Center of MFT

 //char [100];
 std::string outfilename = "Fittercheck_" + std::string(trkFile);
 //strcat(outfilename,trkFile);
 //strcat(outfilename,"");
 TFile outFile(outfilename.c_str(),"RECREATE");

 // Reconstructed MFT Tracks
 std::cout << "Loop over events and reconstructed MFT Tracks!" << std::endl;
 // TracksMFT - Identify reconstructed tracks
 auto totalTracks = 0;
 for (int iEvent = 0 ; iEvent < numberOfEvents ; iEvent++ ) {
   auto  iTrack = 0;
   if(DEBUG_VERBOSE) std::cout << "Event = " << iEvent << std::endl;
   o2SimKineTree -> GetEntry(iEvent);
   for (auto &trackMFT : trackMFTVec) {
     auto label = mcLabels->getLabels(iTrack);
     //label[0].print();
     if(iEvent==label[0].getEventID()) {
       if(DEBUG_VERBOSE) std::cout << "  Track #" << iTrack << ":  trackID = " << label[0].getTrackID() <<
			   " ; EventID = " << label[0].getEventID() <<
			   " ; SourceID = " << label[0].getSourceID() <<
			   " ; isFake = " << label[0].isFake() << std::endl;

       if(label[0].isCorrect()) { // Good track: add to histograms
	 nCleanTracksMFT++;
	 auto eventID = label[0].getEventID();
	 auto thisTrkID = label[0].getTrackID();
 	 MCTrackT<float>* thisTrack =  &(*mcTr).at(thisTrkID);
	 auto vx_MC = thisTrack->GetStartVertexCoordinatesX();
	 auto vy_MC = thisTrack->GetStartVertexCoordinatesY();
	 auto vz_MC = thisTrack->GetStartVertexCoordinatesZ();
	 auto Pt_MC = thisTrack->GetPt();
	 auto P_MC = thisTrack->GetP();
	 auto phi_MC = TMath::ATan2(thisTrack->Py(),thisTrack->Px());
	 auto eta_MC = atanh (thisTrack->GetStartVertexMomentumZ()/P_MC); // eta;
	 auto tanl_MC = thisTrack->Pz()/thisTrack->GetPt();
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

	 trackMFT.propagateToZhelix(vz_MC,field_z);
	 //trackMFT.propagateToZquadratic(vz_MC,field_z);
	 //trackMFT.propagateToZlinear(vz_MC,field_z);

	 auto Q_fit = trackMFT.getCharge();
	 auto dx = trackMFT.getX() - vx_MC;
	 auto dy = trackMFT.getY() - vy_MC;
	 auto d_eta = trackMFT.getEta() - eta_MC;
	 auto d_tanl = trackMFT.getTanl() - tanl_MC;
	 auto Pt_fit = trackMFT.getPt();
	 auto invQPt_seed = trackMFT.getInvQPtSeed();
	 auto d_invQPt = Q_fit/Pt_fit-Q_MC/Pt_MC;
	 auto d_invQPtSeed = invQPt_seed-Q_fit/Pt_fit;//Q_MC/Pt_MC;
	 auto P_fit = trackMFT.getP();
	 auto P_res = P_fit / P_MC;
	 auto Pt_res = Pt_fit / Pt_MC;
	 auto d_Phi = trackMFT.getPhi() - phi_MC;
	 auto d_Charge = Q_fit-Q_MC;
	 auto xChi2 = dx*dx/trackMFT.getCovariances()(0,0);
	 auto yChi2 = dy*dy/trackMFT.getCovariances()(1,1);
	 auto phiChi2 = d_Phi*d_Phi/trackMFT.getCovariances()(2,2);
	 auto tanlChi2 = d_tanl*d_tanl/trackMFT.getCovariances()(3,3);
	 auto invQPtChi2 = d_invQPt*d_invQPt/sqrt(trackMFT.getCovariances()(4,4));
         auto fitChi2 = xChi2 + yChi2 + phiChi2 + tanlChi2;// + invQPtChi2;
	 auto trackChi2 = trackMFT.getTrackChi2();
	 TH1Histos[kMFTTracksP]->Fill(trackMFT.getP());
	 TH1Histos[kMFTTrackDeltaTanl]->Fill(d_tanl);
	 TH1Histos[kMFTTrackDeltaPhi]->Fill(d_Phi);
	 TH1Histos[kMFTTrackDeltaInvQPt]->Fill(d_invQPt);
	 TH1Histos[kMFTTrackDeltaInvQPtSeed]->Fill(d_invQPtSeed);
	 TH1Histos[kMFTTrackDeltaPhiDeg]->Fill(TMath::RadToDeg()*d_Phi);
	 TH1Histos[kMFTTrackDeltaX]->Fill(dx);
	 
	 //std::cout << "DeltaX / sigmaX = " << dx/sqrt(trackMFT.getCovariances()(0,0)) << std::endl;
	 TH1Histos[kMFTTrackDeltaXErr]->Fill(dx/sqrt(trackMFT.getCovariances()(0,0)));
	 //std::cout << "DeltaY / sigmaY = " << dy/sqrt(trackMFT.getCovariances()(1,1)) << std::endl;
	 TH1Histos[kMFTTrackDeltaYErr]->Fill(dy/sqrt(trackMFT.getCovariances()(1,1)));
	 //std::cout << "DeltaPhi / sigmaPhi = " << d_Phi/sqrt(trackMFT.getCovariances()(2,2)) << std::endl;
	 TH1Histos[kMFTTrackDeltaPhiErr]->Fill(d_Phi/sqrt(trackMFT.getCovariances()(2,2)));
	 //std::cout << "DeltaTanl / sigmaTanl = " << d_tanl/sqrt(trackMFT.getCovariances()(3,3)) << std::endl;
	 TH1Histos[kMFTTrackDeltaTanLErr]->Fill(d_tanl/sqrt(trackMFT.getCovariances()(3,3)));
	 //std::cout << "DeltaPt / sigmaPt = " << d_Pt/sqrt(trackMFT.getCovariances()(4,4)) << std::endl;
	 TH1Histos[kMFTTrackDeltainvQPtErr]->Fill(d_invQPt/sqrt(trackMFT.getCovariances()(4,4)));

	 TH1Histos[kMFTTrackXChi2]->Fill(xChi2);
	 TH1Histos[kMFTTrackYChi2]->Fill(yChi2);
	 TH1Histos[kMFTTrackPhiChi2]->Fill(phiChi2);
	 TH1Histos[kMFTTrackTanlChi2]->Fill(tanlChi2);
	 TH1Histos[kMFTTrackinvQPtChi2]->Fill(invQPtChi2);
	 TH1Histos[kFitChi2]->Fill(fitChi2);
	 TH2Histos[kMFTTrackChi2vsFitChi2]->Fill(fitChi2,trackChi2);

	 DeltaX_Profile->Fill(Pt_MC,dx*1e4);
	 TH1Histos[kMFTTrackDeltaY]->Fill(dy);
	 TH1Histos[kMFTTrackR]->Fill(sqrt(dx*dx+dy*dy));
	 TH1Histos[kMFTTrackQ]->Fill(d_Charge);
	 TH1Histos[kMFTTrackChi2]->Fill(trackChi2);
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
	   TH1Histos[kMFTTrackDeltaTanl0_1]->Fill(d_tanl);
	   TH1Histos[kMFTTrackDeltaPhi0_1]->Fill(d_Phi);
	   TH1Histos[kMFTTrackDeltaPhiDeg0_1]->Fill(TMath::RadToDeg()*d_Phi);
	   TH1Histos[kMFTTrackDeltaX0_1]->Fill(dx);
	   TH1Histos[kMFTTrackQ0_1]->Fill(d_Charge);
	   d_Charge ? nChargeMiss0_1++ : nChargeMatch0_1++;
	 }
	 if (Pt_MC > 1.0 and Pt_MC <= 4 ) {
	   TH2Histos[kMFTTrackDeltaXYVertex1_4]->Fill(dx,dy);
	   TH1Histos[kMFTTrackDeltaTanl1_4]->Fill(d_tanl);
	   TH1Histos[kMFTTrackDeltaPhi1_4]->Fill(d_Phi);
	   TH1Histos[kMFTTrackDeltaPhiDeg1_4]->Fill(TMath::RadToDeg()*d_Phi);
	   TH1Histos[kMFTTrackDeltaX1_4]->Fill(dx);
	   TH1Histos[kMFTTrackQ1_4]->Fill(d_Charge);
	   d_Charge ? nChargeMiss1_4++ : nChargeMatch1_4++;
	 }
	 if (Pt_MC > 4.0) {
	   TH2Histos[kMFTTrackDeltaXYVertex4plus]->Fill(dx,dy);
	   TH1Histos[kMFTTrackDeltaTanl4plus]->Fill(d_tanl);
	   TH1Histos[kMFTTrackDeltaPhi4plus]->Fill(d_Phi);
	   TH1Histos[kMFTTrackDeltaPhiDeg4plus]->Fill(TMath::RadToDeg()*d_Phi);
	   TH1Histos[kMFTTrackDeltaX4plus]->Fill(dx);
	   TH1Histos[kMFTTrackQ4plus]->Fill(d_Charge);
	   d_Charge ? nChargeMiss4plus++ : nChargeMatch4plus++;
	 }

	 d_Charge ? nChargeMiss++ : nChargeMatch++;
	 qMatchEff->Fill(!d_Charge,Pt_MC);
       }
       else {
	 nInvalidTracksMFT++;
       }
     }
     iTrack++;
   } // Loop on TracksMFT
   totalTracks+=iTrack;
 } // Loop over events

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


 auto param_resolution = summary_report_3x2(*TH2Histos[kMFTTrackDeltaXYVertex],
					   *TH2Histos[kMFTrackPtResolution],
				           *PtRes_Profile,
					   *DeltaX_Error,
					   *TH2Histos[kMFTrackQPRec_MC],
				           *qMatchEff,
				           "Param Summary",
				           seed_cfg,
				           0, 0, 0, 0, 0, 0,
					   Form("%.2f%%", 100.0*TH2Histos[kMFTTrackDeltaXYVertex]->Integral()/TH2Histos[kMFTTrackDeltaXYVertex]->GetEntries()),
				           Form("%.2f%%", 100.0*TH2Histos[kMFTrackPtResolution]->Integral()/TH2Histos[kMFTrackPtResolution]->GetEntries()),
					   "-",
					   "-",
				           Form("%.2f%%", 100.0*TH2Histos[kMFTrackQPRec_MC]->Integral()/TH2Histos[kMFTrackQPRec_MC]->GetEntries()),
					   "-"
				           );

 
 auto covariances_summary = summary_report_3x2(*TH1Histos[kMFTTrackDeltaXErr],
					       *TH1Histos[kMFTTrackDeltaPhiErr],
					       *TH1Histos[kMFTTrackDeltainvQPtErr],
					       *TH1Histos[kMFTTrackDeltaYErr],
					       *TH1Histos[kMFTTrackDeltaTanLErr],
					       *TH2Histos[kMFTrackQPRec_MC],
					       "Covariances Summary",
					       seed_cfg,
					       1, 1, 1, 1, 1, 0,
					       Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaXErr]->Integral()/TH1Histos[kMFTTrackDeltaXErr]->GetEntries()),
					       Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaPhiErr]->Integral()/TH1Histos[kMFTTrackDeltaPhiErr]->GetEntries()),
					       Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltainvQPtErr]->Integral()/TH1Histos[kMFTTrackDeltainvQPtErr]->GetEntries()),
					       Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaYErr]->Integral()/TH1Histos[kMFTTrackDeltaYErr]->GetEntries()),
					       Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaTanLErr]->Integral()/TH1Histos[kMFTTrackDeltaTanLErr]->GetEntries()),
					       Form("%.2f%%", 100.0*TH2Histos[kMFTrackQPRec_MC]->Integral()/TH2Histos[kMFTrackQPRec_MC]->GetEntries())
					       );


 auto long_summary = summary_report_3x3(*TH2Histos[kMFTTrackDeltaXYVertex],
					*TH1Histos[kMFTTrackDeltaXErr],
					*TH1Histos[kMFTTrackDeltaYErr],
					*DeltaX_Error,
					*TH2Histos[kMFTrackQPRec_MC],
					*TH1Histos[kMFTTrackDeltaPhiErr],
					*qMatchEff,
					*TH1Histos[kMFTTrackDeltainvQPtErr],
					*TH1Histos[kMFTTrackDeltaTanLErr],
					"Summary3x3",
					seed_cfg,
					0, 1, 1, 0, 0, 1, 0, 1, 1,
					Form("%.2f%%", 100.0*TH2Histos[kMFTTrackDeltaXYVertex]->Integral()/TH2Histos[kMFTTrackDeltaXYVertex]->GetEntries()),
					Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaXErr]->Integral()/TH1Histos[kMFTTrackDeltaXErr]->GetEntries()),
					Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaYErr]->Integral()/TH1Histos[kMFTTrackDeltaYErr]->GetEntries()),
					"-",
					Form("%.2f%%", 100.0*TH2Histos[kMFTrackQPRec_MC]->Integral()/TH2Histos[kMFTrackQPRec_MC]->GetEntries()),
					Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaPhiErr]->Integral()/TH1Histos[kMFTTrackDeltaPhiErr]->GetEntries()),
					"-",
					Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltainvQPtErr]->Integral()/TH1Histos[kMFTTrackDeltainvQPtErr]->GetEntries()),
					Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaTanLErr]->Integral()/TH1Histos[kMFTTrackDeltaTanLErr]->GetEntries())
					);


 auto param_summary_diff_pt = summary_report_3x3(*TH1Histos[kMFTTrackDeltaX0_1],
						 *TH1Histos[kMFTTrackDeltaTanl0_1],
						 *TH1Histos[kMFTTrackDeltaPhiDeg0_1],
						 *TH1Histos[kMFTTrackDeltaX1_4],
						 *TH1Histos[kMFTTrackDeltaTanl1_4],
						 *TH1Histos[kMFTTrackDeltaPhiDeg1_4],
						 *TH1Histos[kMFTTrackDeltaX4plus],
						 *TH1Histos[kMFTTrackDeltaTanl4plus],
						 *TH1Histos[kMFTTrackDeltaPhiDeg4plus],
						 "ParamSummaryVsPt",
						 seed_cfg,
 						 1, 1, 1, 1, 1, 1, 1, 1, 1,
						 Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaX0_1]->Integral()/TH1Histos[kMFTTrackDeltaX0_1]->GetEntries()),
						 Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaTanl0_1]->Integral()/TH1Histos[kMFTTrackDeltaTanl0_1]->GetEntries()),
						 Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaPhiDeg0_1]->Integral()/TH1Histos[kMFTTrackDeltaPhiDeg0_1]->GetEntries()),
						 Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaX1_4]->Integral()/TH1Histos[kMFTTrackDeltaX1_4]->GetEntries()),
						 Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaTanl1_4]->Integral()/TH1Histos[kMFTTrackDeltaTanl1_4]->GetEntries()),
						 Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaPhiDeg1_4]->Integral()/TH1Histos[kMFTTrackDeltaPhiDeg1_4]->GetEntries()),
						 Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaX4plus]->Integral()/TH1Histos[kMFTTrackDeltaX4plus]->GetEntries()),
						 Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaTanl4plus]->Integral()/TH1Histos[kMFTTrackDeltaTanl4plus]->GetEntries()),
						 Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaPhiDeg4plus]->Integral()/TH1Histos[kMFTTrackDeltaPhiDeg4plus]->GetEntries())
						 );

 
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
					    "Vertexing Summary",
					    seed_cfg,
					    0, 1, 0, 1,
					    Form("%.2f%%", 100.0*TH2Histos[kMFTTrackDeltaXYVertex]->Integral()/TH2Histos[kMFTTrackDeltaXYVertex]->GetEntries()),
					    Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaX]->Integral()/TH1Histos[kMFTTrackDeltaX]->GetEntries()),
					    Form("-"),
					    Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaPhiDeg]->Integral()/TH1Histos[kMFTTrackDeltaPhiDeg]->GetEntries())
					    );
 
 auto vertexing_resolution0_1 = summary_report(*TH2Histos[kMFTTrackDeltaXYVertex0_1],
					       *TH1Histos[kMFTTrackDeltaX0_1],
					       *TH1Histos[kMFTTrackDeltaTanl0_1],
					       *TH1Histos[kMFTTrackDeltaPhiDeg0_1],
					       "Vertexing Summary pt < 1",
					       seed_cfg,
					       0, 1, 1, 1,
					       Form("%.2f%%", 100.0*TH2Histos[kMFTTrackDeltaXYVertex0_1]->Integral()/TH2Histos[kMFTTrackDeltaXYVertex0_1]->GetEntries()),
					       Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaX0_1]->Integral()/TH1Histos[kMFTTrackDeltaX0_1]->GetEntries()),
					       Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaTanl0_1]->Integral()/TH1Histos[kMFTTrackDeltaTanl0_1]->GetEntries()),
					       Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaPhiDeg0_1]->Integral()/TH1Histos[kMFTTrackDeltaPhiDeg0_1]->GetEntries())
					       );

 auto vertexing_resolution1_4 = summary_report(*TH2Histos[kMFTTrackDeltaXYVertex1_4],
					       *TH1Histos[kMFTTrackDeltaX1_4],
					       *TH1Histos[kMFTTrackDeltaTanl1_4],
					       *TH1Histos[kMFTTrackDeltaPhiDeg1_4],
					       "Vertexing Summary 1 < p_t < 4",
					       seed_cfg,
					       0, 1, 1, 1,
					       Form("%.2f%%", 100.0*TH2Histos[kMFTTrackDeltaXYVertex1_4]->Integral()/TH2Histos[kMFTTrackDeltaXYVertex1_4]->GetEntries()),
					       Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaX1_4]->Integral()/TH1Histos[kMFTTrackDeltaX1_4]->GetEntries()),
					       Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaTanl1_4]->Integral()/TH1Histos[kMFTTrackDeltaTanl1_4]->GetEntries()),
					       Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaPhiDeg1_4]->Integral()/TH1Histos[kMFTTrackDeltaPhiDeg1_4]->GetEntries())
					       );

 auto vertexing_resolution4plus = summary_report(*TH2Histos[kMFTTrackDeltaXYVertex4plus],
						 *TH1Histos[kMFTTrackDeltaX4plus],
						 *TH1Histos[kMFTTrackDeltaTanl4plus],
						 *TH1Histos[kMFTTrackDeltaPhiDeg4plus],
						 "Vertexing Summary p_t > 4",
						 seed_cfg,
						 0, 1, 1, 1,
						 Form("%.2f%%", 100.0*TH2Histos[kMFTTrackDeltaXYVertex4plus]->Integral()/TH2Histos[kMFTTrackDeltaXYVertex4plus]->GetEntries()),
						 Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaX4plus]->Integral()/TH1Histos[kMFTTrackDeltaX4plus]->GetEntries()),
						 Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaTanl4plus]->Integral()/TH1Histos[kMFTTrackDeltaTanl4plus]->GetEntries()),
						 Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaPhiDeg4plus]->Integral()/TH1Histos[kMFTTrackDeltaPhiDeg4plus]->GetEntries())
						 );
 /* 
 auto covariances_summary = summary_report(*TH1Histos[kMFTTrackDeltainvQPtErr],
					   *TH1Histos[kMFTTrackDeltaXErr],
					   *TH1Histos[kMFTTrackDeltaTanLErr],
					   *TH1Histos[kMFTTrackDeltaPhiErr],
					   "Covariances Summary",
					   seed_cfg,
					   1, 1, 1, 1,
					   Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltainvQPtErr]->Integral()/TH1Histos[kMFTTrackDeltainvQPtErr]->GetEntries()),
					   Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaXErr]->Integral()/TH1Histos[kMFTTrackDeltaXErr]->GetEntries()),
					   Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaTanLErr]->Integral()/TH1Histos[kMFTTrackDeltaTanLErr]->GetEntries()),
					   Form("%.2f%%", 100.0*TH1Histos[kMFTTrackDeltaPhiErr]->Integral()/TH1Histos[kMFTTrackDeltaPhiErr]->GetEntries())
					   ); */

 auto chi2_summary = summary_report(*TH1Histos[kMFTTrackChi2],
				    *TH1Histos[kMFTTrackXChi2],
				    *TH1Histos[kMFTTrackTanlChi2],
				    *TH1Histos[kMFTTrackPhiChi2],
				    "Chi2 Summary",
				    seed_cfg,
				    1, 1, 1, 1,				    
				    Form("%.2f%%", 100.0*TH1Histos[kMFTTrackChi2]->Integral()/TH1Histos[kMFTTrackChi2]->GetEntries()),
				    Form("%.2f%%", 100.0*TH1Histos[kMFTTrackXChi2]->Integral()/TH1Histos[kMFTTrackXChi2]->GetEntries()),
				    Form("%.2f%%", 100.0*TH1Histos[kMFTTrackTanlChi2]->Integral()/TH1Histos[kMFTTrackTanlChi2]->GetEntries()),
				    Form("%.2f%%", 100.0*TH1Histos[kMFTTrackPhiChi2]->Integral()/TH1Histos[kMFTTrackPhiChi2]->GetEntries())
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
 std::cout << " Tanl_mean = " << TH1Histos[kMFTTrackDeltaTanl]->GetMean() << std::endl;
 std::cout << " Tanl_StdDev = " << TH1Histos[kMFTTrackDeltaTanl]->GetStdDev() << std::endl;
 std::cout << " Tanl_StdDev(pt<1) = " << TH1Histos[kMFTTrackDeltaTanl0_1]->GetStdDev() << std::endl;
 std::cout << " Tanl_StdDev(1<pt<4) = " << TH1Histos[kMFTTrackDeltaTanl1_4]->GetStdDev() << std::endl;
 std::cout << " Tanl_StdDev(pt>4) = " << TH1Histos[kMFTTrackDeltaTanl4plus]->GetStdDev() << std::endl;
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
