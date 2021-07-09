#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "SimulationDataFormat/MCEventHeader.h"
#include "MFTAnaBaseTask.h"

#endif

// 1D Histograms
std::vector<std::string> TH1HistoCodes{
    "MCNrOfHits",
    "MCall_Pt",
    "MCinSA_Pt"};

std::map<std::string, const char *> TH1Names{
    {"MCNrOfHits", "MC number of hits"},
    {"MCall_Pt", "Pt of MC tracks with hits"},
    {"MCinSA_Pt", "Pt of MC tracks with clusters in SA tracks"}};

std::map<std::string, const char *> TH1Titles{
    {"MCNrOfHits", "MC number of hits"},
    {"MCall_Pt", "Pt of MC tracks with hits"},
    {"MCinSA_Pt", "Pt of MC tracks with clusters in SA tracks"}};

std::map<std::string, std::array<double, 3>> TH1Binning{
    {"MCNrOfHits", {50, 0., 50.}},
    {"MCall_Pt", {100, 0., 5.}},
    {"MCinSA_Pt", {100, 0., 5.}}};

std::map<std::string, const char *> TH1XaxisTitles{
    {"MCNrOfHits", "Nr of hits"},
    {"MCall_Pt", "pt [GeV/c]"},
    {"MCinSA_Pt", "pt [GeV/c]"}};

// 2D Histograms
std::vector<std::string> TH2HistoCodes{
    "TracksZ_Eta"};

std::map<std::string, const char *> TH2Names{
    {"TracksZ_Eta", "MFTTracks_Z_Eta"}};

std::map<std::string, const char *> TH2Titles{
    {"TracksZ_Eta", "MFT Tracks Z-Eta"}};

std::map<std::string, std::array<double, 6>> TH2Binning{
    {"TracksZ_Eta", {40, -10, 10, 40, -4., -2.}}};

std::map<std::string, const char *> TH2XaxisTitles{
    {"TracksZ_Eta", "\\eta"}};

std::map<std::string, const char *> TH2YaxisTitles{
    {"TracksZ_Eta", "Vertex Z position [cm]"}};

//_____________________________________________________________________________
class MFTAnaSampleTask_ : public MFTAnaBaseTask
{
public:
    void process()
    {
        // MC tracks
        int iMCTrack = 0;
        for (auto &mcTrack : anaSimTracks)
        {
            TH1Histos["MCNrOfHits"]->Fill(mcTrack.getNHits());
            TH1Histos["MCall_Pt"]->Fill(mcTrack.getPt());
            if (mcTrack.getNSATracks() > 0)
            {
                TH2Histos["TracksZ_Eta"]->Fill(mcTrack.getVertexZ(), mcTrack.getEta());
                TH1Histos["MCinSA_Pt"]->Fill(mcTrack.getPt());
            }
        }
        // SA tracks
        int iSATrack = 0;
        for (auto &saTrack : anaSimSATracks)
        {
            //printSATrack(&saTrack, iSATrack++);
        }
    }
};

//_____________________________________________________________________________
void MFTAnaSampleTask(std::string run = "")
{
    gSystem->Load("${ALIBUILD_WORK_DIR}/MFTAna/libMFTAnaSim");

    MFTAnaSampleTask_ MFTtask;
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
    MFTtask.run(run);
}
