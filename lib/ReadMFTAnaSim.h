#include <vector>

#include <TH1F.h>

enum TH1HistoCodes {
  MCNrOfHits,
  MCHasRecPt,
  NHistograms
};

std::map<int, const char *> TH1Names {
  {MCNrOfHits, "MC number of hits"},
  {MCHasRecPt, "Pt of MC tracks with clusters in SA track(s)"}
};

std::map<int, const char *> TH1Titles {
  {MCNrOfHits, "MC number of hits"},
  {MCHasRecPt, "Pt of MC tracks with clusters in SA track(s)"}
};

std::map<int, std::array<double, 3>> TH1Binning {
  {MCNrOfHits, {50, 0., 50.}},
  {MCHasRecPt, {100, 0., 10.}}
};

std::map<int, const char *> TH1XaxisTitles {
  {MCNrOfHits, "Nr of hits"},
  {MCHasRecPt, "pt [GeV/c]"}
};

std::vector<TH1F*> TH1Histos(NHistograms);

//_____________________________________________________________________________
void createHistograms() {  
  auto nHisto = 0;
  for (auto& h : TH1Histos) {
    h = new TH1F(TH1Names[nHisto], TH1Titles[nHisto], (int)TH1Binning[nHisto][0], TH1Binning[nHisto][1], TH1Binning[nHisto][2]);
    h->GetXaxis()->SetTitle(TH1XaxisTitles[nHisto]);
    ++nHisto;
  }
}
