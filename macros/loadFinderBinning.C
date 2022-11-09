

  void loadFinderBinning()
  {

  
    std::cout << "Loading MFTFinderBinning.root" << std::endl;
    TFile* f = new TFile(Form("MFTFinderBinning.root"));

    auto* Bins = f->GetObjectChecked("mBinsS", "std::array<std::array<std::vector<std::vector<Int_t>>, 9>, 9>");
    auto* mBins = (static_cast<std::array<std::array<std::vector<std::vector<Int_t>>, 9>, 9>*>(Bins));
    auto castedBins = *mBins;
    for (Int_t layer1 = 0; layer1 < (10 - 1); ++layer1) {
        for (Int_t layer2 = (layer1 + 1); layer2 < 10; ++layer2) {
	  std::cout << "layer1 = " << layer1 << " ; layer2 = " << layer2 << std::endl;
	  int binindex = 0;
            for (auto& bins : castedBins[0][1]) {
              std::cout << "bin #" << binindex++ << ": ";
              for (auto& i : bins) {
                std::cout << i << " ";
              }
           std::cout << std::endl;
           }
	}
    }
  }
