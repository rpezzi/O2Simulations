#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "SimulationDataFormat/MCEventHeader.h"
#include "MFTAnaBaseTask.h"
#include <TEfficiency.h>
#include <TCanvas.h>

#endif

// 1D Histograms
std::vector<std::string> TH1HistoCodes{
    "TrackablesEta",
    "TrackablesPt"};

std::map<std::string, const char *> TH1Names{
    {"TrackablesEta", "MFTTrackablesEta"},
    {"TrackablesPt", "MFTTrackablesPt"}};

std::map<std::string, const char *> TH1Titles{
    {"TrackablesEta", "MC number of hits"},
    {"TrackablesPt", "Pt pf MC tracks with hits"}};

std::map<std::string, std::array<double, 3>> TH1Binning{
    {"TrackablesEta", {40, -4., -2.}},
    {"TrackablesPt", {100, 0., 10.}}};

std::map<std::string, const char *> TH1XaxisTitles{
    {"TrackablesEta", "\\eta"},
    {"TrackablesPt", "pt [GeV/c]"}};

// 2D Histograms
std::vector<std::string> TH2HistoCodes{
    "TrackablesZ_Eta",
    "TracksZ_Eta"};

std::map<std::string, const char *> TH2Names{
    {"TrackablesZ_Eta", "MFTTrackablesZ_Eta"},
    {"TracksZ_Eta", "MFTTracks_Z_Eta"}};

std::map<std::string, const char *> TH2Titles{
    {"TrackablesZ_Eta", "MFT Trackables Z-Eta"},
    {"TracksZ_Eta", "MFT Tracks Z-Eta"}};

std::map<std::string, std::array<double, 6>> TH2Binning{
    {"TrackablesZ_Eta", {160, -20, 20, 40, -4., -2.}},
    {"TracksZ_Eta", {160, -20, 20, 40, -4., -2.}}};

std::map<std::string, const char *> TH2XaxisTitles{
    {"TrackablesZ_Eta", "Vertex Z position [cm]"},
    {"TracksZ_Eta", "Vertex Z position [cm]"}};

std::map<std::string, const char *> TH2YaxisTitles{
    {"TrackablesZ_Eta", "\\eta"},
    {"TracksZ_Eta", "\\eta"}};

// 3D Histograms
std::vector<std::string> TH3HistoCodes{
    "GenTracksPt_Z_Eta",
    "TrackablesPt_Z_Eta",
    "TracksPt_Z_Eta"};

std::map<std::string, const char *> TH3Names{
    {"GenTracksPt_Z_Eta", "GenTracksPt_Z_Eta"},
    {"TrackablesPt_Z_Eta", "MFTTrackablesPt_Z_Eta"},
    {"TracksPt_Z_Eta", "MFTTracksPt_Z_Eta"}};

std::map<std::string, const char *> TH3Titles{
    {"GenTracksPt_Z_Eta", "Generated tracks Pt-Z-Eta"},
    {"TrackablesPt_Z_Eta", "MFT Trackables Pt-Z-Eta"},
    {"TracksPt_Z_Eta", "MFT Tracks Pt-Z-Eta"}};

std::map<std::string, std::array<double, 9>> TH3Binning{
    {"GenTracksPt_Z_Eta", {100, 0., 10., 160, -20, 20, 40, -4., -2.}},
    {"TrackablesPt_Z_Eta", {100, 0., 10., 160, -20, 20, 40, -4., -2.}},
    {"TracksPt_Z_Eta", {100, 0., 10., 160, -20, 20, 40, -4., -2.}}};

std::map<std::string, const char *> TH3XaxisTitles{
    {"GenTracksPt_Z_Eta", "p_t [GeV/c]"},
    {"TrackablesPt_Z_Eta", "p_t [GeV/c]"},
    {"TracksPt_Z_Eta", "p_t [GeV/c]"}};

std::map<std::string, const char *> TH3YaxisTitles{
    {"GenTracksPt_Z_Eta", "Vertex Z position [cm]"},
    {"TrackablesPt_Z_Eta", "Vertex Z position [cm]"},
    {"TracksPt_Z_Eta", "Vertex Z position [cm]"}};

std::map<std::string, const char *> TH3ZaxisTitles{
    {"GenTracksPt_Z_Eta", "\\eta"},
    {"TrackablesPt_Z_Eta", "\\eta"},
    {"TracksPt_Z_Eta", "\\eta"}};

class MFTAcceptanceTask_ : public MFTAnaBaseTask
{
public:
    void process()
    {
        // MC tracks
        int iMCTrack = 0;
        for (auto &mcTrack : anaSimTracks)
        {
            TH3Histos["GenTracksPt_Z_Eta"]->Fill(mcTrack.getPt(), mcTrack.getVertexZ(), mcTrack.getEta());
            if (mcTrack.isTrackable())
            {
                TH1Histos["TrackablesPt"]->Fill(mcTrack.getPt());
                TH1Histos["TrackablesEta"]->Fill(mcTrack.getEta());
                TH2Histos["TrackablesZ_Eta"]->Fill(mcTrack.getVertexZ(), mcTrack.getEta());
                TH3Histos["TrackablesPt_Z_Eta"]->Fill(mcTrack.getPt(), mcTrack.getVertexZ(), mcTrack.getEta());
                if (mcTrack.getNSATracks() > 0)
                {
                    TH2Histos["TracksZ_Eta"]->Fill(mcTrack.getVertexZ(), mcTrack.getEta());
                    TH3Histos["TracksPt_Z_Eta"]->Fill(mcTrack.getPt(), mcTrack.getVertexZ(), mcTrack.getEta());
                }
            }
        }

        // Efficiency 2D
        TH2F MFTEff_Z_Eta = (*TH2Histos["TracksZ_Eta"]) / (*TH2Histos["TrackablesZ_Eta"]);
        MFTEff_Z_Eta.SetNameTitle("MFT_Tracking efficiencyZ_Eta", "MFT Tracking Efficiency");
        MFTEff_Z_Eta.SetOption("COLZ");
        MFTEff_Z_Eta.Write();

        // Efficiency 3D
        TH3F MFTEff_Pt_Z_Eta = (*TH3Histos["TracksPt_Z_Eta"]) / (*TH3Histos["TrackablesPt_Z_Eta"]);
        MFTEff_Pt_Z_Eta.SetNameTitle("MFT_Tracking_EfficiencyPt_Z_Eta", "MFT Tracking Efficiency");
        MFTEff_Pt_Z_Eta.Write();

        // Acceptance * Efficiency 3D
        TH3F MFTAccptEffPt_Z_Eta = (*TH3Histos["TracksPt_Z_Eta"]) / (*TH3Histos["GenTracksPt_Z_Eta"]);
        MFTAccptEffPt_Z_Eta.SetNameTitle("MFT_Acc_Eff_Pt_Z_Eta", "MFT Acceptance * Efficiency");
        MFTAccptEffPt_Z_Eta.Write();

        // MFT Standalone tracks (SA tracks)
        int iSATrack = 0;
        for (auto &saTrack : anaSimSATracks)
        {
            //printSATrack(&saTrack, iSATrack++);
        }
    }
};

void MFTAcceptanceTask(std::string run = "")
{
    gSystem->Load("${ALIBUILD_WORK_DIR}/MFTAna/libMFTAnaSim");

    MFTAcceptanceTask_ MFTtask;
    MFTtask.setTH1HistoCodes(TH1HistoCodes);
    MFTtask.setTH1Names(TH1Names);
    MFTtask.setTH1Titles(TH1Titles);
    MFTtask.setTH1Binning(TH1Binning);
    MFTtask.setTH1XaxisTitles(TH1XaxisTitles);

    MFTtask.setTH2HistoCodes(TH2HistoCodes);
    MFTtask.setTH2Names(TH2Names);
    MFTtask.setTH2Titles(TH2Titles);
    MFTtask.setTH2Binning(TH2Binning);
    MFTtask.setTH2XaxisTitles(TH2XaxisTitles);
    MFTtask.setTH2YaxisTitles(TH2YaxisTitles);

    MFTtask.setTH3HistoCodes(TH3HistoCodes);
    MFTtask.setTH3Names(TH3Names);
    MFTtask.setTH3Titles(TH3Titles);
    MFTtask.setTH3Binning(TH3Binning);
    MFTtask.setTH3XaxisTitles(TH3XaxisTitles);
    MFTtask.setTH3YaxisTitles(TH3YaxisTitles);
    MFTtask.setTH3ZaxisTitles(TH3ZaxisTitles);
    MFTtask.run(run);
}
