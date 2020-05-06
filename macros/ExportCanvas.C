#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TLatex.h"
#include <TMath.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>


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
template <typename H1, typename H2, typename H3, typename H4>
TCanvas summary_report(const H1& histo1, const H2& histo2, const H3& histo3, const H4& histo4, std::string CanvasName, std::string tlt="Summary")
{
H1 *h1 = new H1(histo1);
H2 *h2 = new H2(histo2);
H3 *h3 = new H3(histo3);
H4 *h4 = new H4(histo4);


h1->FitSlicesY(0,0,-1,1);
// Create a canvas and divide it
TCanvas *c1 = new TCanvas(CanvasName.c_str(),CanvasName.c_str(),700,500);

TLatex *Title = new TLatex() ;
Title->SetTextSize(0.035);
Title->SetTextAlign(23);
Title->DrawLatex(.5 ,1, tlt.c_str());
Title->Draw();

c1->Divide(2,1);
TPad *leftPad = (TPad*)c1->cd(1);;
leftPad->SetPad(0.0086,0.010,0.49,0.97);
leftPad->Divide(1,2);

leftPad->cd(1);
gPad->SetTopMargin(.12);
h1->Draw();

leftPad->cd(2);
h2->Draw();
TPad *rightPad = (TPad*)c1->cd(2);
rightPad->SetPad(0.51,0.010,0.99,0.97);
rightPad->Divide(1,2);
rightPad->cd(1);
h3->Draw();

rightPad->cd(2);
gPad->SetTopMargin(.12);
h4->Draw();

c1->cd();
c1->Write();
return c1;
}


//_________________________________________________________________________________________________
int ExportCanvas( const Char_t *inFileName,
                  const Char_t *canvasName) {


  // Seed configuration
  std::string export_name{inFileName};
  std::string cut_1{"Fittercheck_mfttracks"};
  std::string cut_2{".root"};
  //std::string trk_trk{"Fittercheck_mfttracks"};
  if (export_name.find(cut_1) < export_name.length()) export_name.replace(export_name.find(cut_1),cut_1.length(),"");
  if (export_name.find(cut_2) < export_name.length()) export_name.replace(export_name.find(cut_2),cut_2.length(),"");
  //if (export_name.find(trk_trk) < export_name.length()) export_name.replace(export_name.find(trk_trk),trk_trk.length(),"");
  std::cout << export_name << std::endl;

  //gStyle->SetLabelFont(32,"XY");
  //gROOT->SetStyle("Bold");
  //gStyle->SetLineWidth(6);
  //gStyle->SetLabelSize(0.06,"xyz");
  //gStyle->SetTitleSize(0.06,"xyz");
  //gROOT->ForceStyle();

  //Input file

  TFile *inFile = new TFile(inFileName);
  if (!inFile) {
    std::cout << inFileName << " error." << std::endl;
    return 1;
  }


 //Get canvas and export
 TCanvas *canvas = (TCanvas*)inFile->Get(canvasName);
 canvas->SetBatch();
 canvas->Draw();
 //canvas->SetCanvasSize(1920,1080);




 //canvas->UseCurrentStyle();
 gSystem->ProcessEvents();
 canvas->Update();
 for (std::string type: {".pdf", ".png", ".svg"}) canvas->Print((std::string(canvasName) + std::string(export_name) + type).c_str());
 return 0;
}
