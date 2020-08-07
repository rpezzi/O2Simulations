#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "DataFormatsITSMFT/Digit.h"
#include "DataFormatsITSMFT/ROFRecord.h"

#include "MFTBase/GeometryTGeo.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include <stdio.h>
#include <math.h>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <sstream>
#include <vector>
#include <fstream>

#endif

using namespace std;
using o2::itsmft::Digit;


int getZone(int layer, int ladderID);

struct DigitChipInfo {
  int zone;
  int transID;
  int layer;
};


vector<DigitChipInfo> getDigitChipInfo();

void mft_sw2hw_digits(int FrameNumber = -1)
{

  // Create output streams  
  vector<fstream> output_layers;
  if(FrameNumber > 0)
  for (auto i = 0 ; i < 10 ; i++) {
    auto filename = std::string("Frame") + std::to_string(FrameNumber) + std::string("_Layer") + std::to_string(i) + std::string(".txt");
    output_layers.emplace_back(filename, ios::out);
  }
  
  auto zone_transIDs = getDigitChipInfo();
  
  // MFT Digits
  const Char_t *DigitsMFTFile = "mftdigits.root";
  TFile *DigitsMFTFileIn = new TFile(DigitsMFTFile);
  TTree *o2MFTDigitsTree = (TTree*) DigitsMFTFileIn -> Get("o2sim");
  std::vector<o2::itsmft::Digit>* mftdigit = nullptr;
  o2MFTDigitsTree -> SetBranchAddress("MFTDigit",&mftdigit);
  std::vector<o2::itsmft::ROFRecord>* mftdigitROFs = nullptr;
  o2MFTDigitsTree -> SetBranchAddress("MFTDigitROF",&mftdigitROFs);

  Int_t numberOfDigitsTreeEntries = o2MFTDigitsTree -> GetEntries();
  o2MFTDigitsTree -> GetEntry(0);  
  auto nROFs = mftdigitROFs->size();
  auto nDigits = mftdigit->size();

  std::cout << "numberOfDigitsTreeEntries = " << numberOfDigitsTreeEntries << std::endl;
  std::cout << "Number of RO Frames = " << nROFs << std::endl;
  std::cout << "Number of MFT Digits = " << nDigits << std::endl;

  if(FrameNumber < 0) {
    std::cout << "No frame specified. Printing ROFRecors.\n";
    int n = 0;
    for (auto rofRec : *mftdigitROFs ) {
      std::cout << "ROFRecord = " << n << " ; First = "<< rofRec.getFirstEntry() << " ; nDigits = " << rofRec.getNEntries() << std::endl ;
      n++;
    }
  }
  else {
    auto firstDigit = (*mftdigitROFs)[FrameNumber].getFirstEntry();
    auto nDigits = (*mftdigitROFs)[FrameNumber].getNEntries();
    std::cout << "## FrameNumber = " << FrameNumber << " ; FirstDigit = "<< firstDigit << " ; nDigits = " << nDigits << std::endl ;
    std::cout << "# zone ; transID ; digit_row ; digit_column " << std::endl;

    for (int i_digit = firstDigit ; i_digit < firstDigit+nDigits ; i_digit++ ) {
      auto& thisdigit = (*mftdigit)[i_digit];
      auto chipID = thisdigit.getChipIndex();
      if (chipID < 936/2) continue;
      auto zone = zone_transIDs[chipID].zone;
      auto transID = zone_transIDs[chipID].transID;   
      auto layer = zone_transIDs[chipID].layer;   
      //std::cout << zone << " ; " << transID << " ; " << thisdigit.getRow() << " ; " << thisdigit.getColumn() << std::endl;
      output_layers[layer] << zone << " ; " << transID << " ; " << thisdigit.getRow() << " ; " << thisdigit.getColumn() << std::endl;
    }

    //Close output streams
    for (auto i = 0 ; i < 10 ; i++) {
      output_layers[i].close();
    }
  }
 
}

//__________________________________________________________________________
vector<DigitChipInfo> getDigitChipInfo() {
  
  std::vector<DigitChipInfo> zoneAndTransIDs;
  std::cout << "Loading Transceiver_IDs" << std::endl;
  const std::string inputGeom = "o2sim_geometry.root";
  o2::base::GeometryManager::loadGeometry(inputGeom);
  auto gm = o2::mft::GeometryTGeo::Instance(); // geometry manager
  const o2::itsmft::ChipMappingMFT map;

  int half, disk, layer, zone, module, ladder, ladderID;
  int chipIDglo, chipIDlocSW, chipIDlocHW, chipOnModule;
  int ruIDSW, ruIDHW, ruType;
  uint16_t ruOnLayer, chipOnRU, link = 0;
  uint8_t connector, cableHW, cableSW;
  int faceID;
  const int NChips = map.getNChips();

  const o2::itsmft::ChipOnRUInfo* chipOnRUInfo;
  const o2::itsmft::RUInfo* ruInfo;
  o2::itsmft::ChipInfo chipInfo;



 for (int chip = 0; chip < NChips; ++chip) {
    
    gm->getSensorID(chip, half, disk, ladder, chipIDlocSW);
    layer = gm->getLayer(chip);
    ladderID = gm->getLadderID(disk, ladder);
    zone = getZone(layer, ladderID);
    map.getChipInfoSW(chip, chipInfo);    
    chipOnRUInfo = chipInfo.chOnRU;
    cableHW = chipOnRUInfo->cableHW;
    int transID = (int)cableHW;
    auto& temp = zoneAndTransIDs.emplace_back();
    temp.zone = getZone(layer, ladderID);
    temp.transID = (int)cableHW;
    temp.layer = layer;
    //std::cout << "chipID = " << chip << " ; half = " << half << " ; disk = " << disk << " ; ladder = " << ladder << " ; chipIDlocSW = " << chipIDlocSW << " ; layer = " << layer << " ; zone = " << zone << " ; transID = " << transID << std::endl ;   

    //std::cout << "chipID = " << chip << " ; zone = " << getZone(layer, ladderID) << " ; TransID = " << transID   << std::endl ; 
 }
 return zoneAndTransIDs;


}



//__________________________________________________________________________
int getZone(int layer, int ladderID)
{
  int zone = -1;
  if (layer == 0) {
    if (ladderID >=  0 && ladderID <=  2) zone = 0;
    if (ladderID >=  3 && ladderID <=  5) zone = 1;
    if (ladderID >=  6 && ladderID <=  8) zone = 2;
    if (ladderID >=  9 && ladderID <= 11) zone = 3;
  }
  if (layer == 1) {
    if (ladderID >= 12 && ladderID <= 14) zone = 3;
    if (ladderID >= 15 && ladderID <= 17) zone = 2;
    if (ladderID >= 18 && ladderID <= 20) zone = 1;
    if (ladderID >= 21 && ladderID <= 23) zone = 0;
  }
  if (layer == 2) {
    if (ladderID >=  0 && ladderID <=  2) zone = 0;
    if (ladderID >=  3 && ladderID <=  5) zone = 1;
    if (ladderID >=  6 && ladderID <=  8) zone = 2;
    if (ladderID >=  9 && ladderID <= 11) zone = 3;
  }
  if (layer == 3) {
    if (ladderID >= 12 && ladderID <= 14) zone = 3;
    if (ladderID >= 15 && ladderID <= 17) zone = 2;
    if (ladderID >= 18 && ladderID <= 20) zone = 1;
    if (ladderID >= 21 && ladderID <= 23) zone = 0;
  }
  if (layer == 4) {
    if (ladderID >=  0 && ladderID <=  2) zone = 0;
    if (ladderID >=  3 && ladderID <=  5) zone = 1;
    if (ladderID >=  6 && ladderID <=  8) zone = 2;
    if (ladderID >=  9 && ladderID <= 12) zone = 3;
  }
  if (layer == 5) {
    if (ladderID >= 13 && ladderID <= 16) zone = 3;
    if (ladderID >= 17 && ladderID <= 19) zone = 2;
    if (ladderID >= 20 && ladderID <= 22) zone = 1;
    if (ladderID >= 23 && ladderID <= 25) zone = 0;
  }
  if (layer == 6) {
    if (ladderID >=  0 && ladderID <=  3) zone = 0;
    if (ladderID >=  4 && ladderID <=  7) zone = 1;
    if (ladderID >=  8 && ladderID <= 11) zone = 2;
    if (ladderID >= 12 && ladderID <= 15) zone = 3;
  }
  if (layer == 7) {
    if (ladderID >= 16 && ladderID <= 19) zone = 3;
    if (ladderID >= 20 && ladderID <= 23) zone = 2;
    if (ladderID >= 24 && ladderID <= 27) zone = 1;
    if (ladderID >= 28 && ladderID <= 31) zone = 0;
  }
  if (layer == 8) {
    if (ladderID >=  0 && ladderID <=  4) zone = 0;
    if (ladderID >=  5 && ladderID <=  8) zone = 1;
    if (ladderID >=  9 && ladderID <= 12) zone = 2;
    if (ladderID >= 13 && ladderID <= 16) zone = 3;
  }
  if (layer == 9) {
    if (ladderID >= 17 && ladderID <= 20) zone = 3;
    if (ladderID >= 21 && ladderID <= 24) zone = 2;
    if (ladderID >= 25 && ladderID <= 28) zone = 1;
    if (ladderID >= 29 && ladderID <= 33) zone = 0;
  }

  return zone;
}
