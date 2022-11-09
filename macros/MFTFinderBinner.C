#include "MFTTracking/Constants.h"
#include "MFTTracking/TrackerConfig.h"
#include "MathUtils/Utils.h"
#include "MathUtils/Cartesian.h"

namespace o2
{
namespace mft
{

typedef std::array<std::array<std::vector<std::vector<Int_t>>, 9>, 9> MFTBinMap;

class Binner : public TrackerConfig
{
 public:
  MFTBinMap mBinsS; // Track finding bin map (track seeds)
  MFTBinMap mBins;  // Track finding bin map (intermediate layers)

  void MFTFinderBinner(Int_t RBins = 50, Int_t PhiBins = 50)
  {
    mRBins = RBins;
    mPhiBins = PhiBins;
    mRPhiBins = mRBins * mPhiBins;

    for (Int_t layer1 = 0; layer1 < (10 - 1); ++layer1) {
      for (Int_t layer2 = (layer1 + 1); layer2 < 10; ++layer2) {
        std::cout << "layer1 = " << layer1 << " ; layer2 = " << layer2 << std::endl;
        int binindex = 0;
        for (auto& bins : mBinsS[0][1]) {
          std::cout << "bin #" << binindex++ << ": ";
          for (auto& i : bins) {
            std::cout << i << " ";
          }
          std::cout << std::endl;
        }
      }
    }

    /// calculate Look-Up-Table of the R-Phi bins projection from one layer to another
    /// layer1 + global R-Phi bin index ---> layer2 + R bin index + Phi bin index

    Float_t dz, x, y, r, phi, x_proj, y_proj, r_proj, phi_proj;
    Int_t binIndex1, binIndex2, binIndex2S, binR_proj, binPhi_proj;

    for (Int_t layer1 = 0; layer1 < (constants::mft::LayersNumber - 1); ++layer1) {
      for (Int_t layer2 = (layer1 + 1); layer2 < constants::mft::LayersNumber; ++layer2) {
        mBinsS[layer1][layer2 - 1].resize(mRBins * mPhiBins);
        mBins[layer1][layer2 - 1].resize(mRBins * mPhiBins);
      }
    }

    // return;

    for (Int_t layer1 = 0; layer1 < (constants::mft::LayersNumber - 1); ++layer1) {

      for (Int_t iRBin = 0; iRBin < mRBins; ++iRBin) {

        r = (iRBin + 0.5) * mRBinSize + constants::index_table::RMin;

        for (Int_t iPhiBin = 0; iPhiBin < mPhiBins; ++iPhiBin) {

          phi = (iPhiBin + 0.5) * mPhiBinSize + constants::index_table::PhiMin;

          binIndex1 = getBinIndex(iRBin, iPhiBin);

          x = r * TMath::Cos(phi);
          y = r * TMath::Sin(phi);

          for (Int_t layer2 = (layer1 + 1); layer2 < constants::mft::LayersNumber; ++layer2) {

            dz = constants::mft::LayerZCoordinate()[layer2] - constants::mft::LayerZCoordinate()[layer1];
            x_proj = x + dz * x * constants::mft::InverseLayerZCoordinate()[layer1];
            y_proj = y + dz * y * constants::mft::InverseLayerZCoordinate()[layer1];
            auto clsPoint2D = math_utils::Point2D<Float_t>(x_proj, y_proj);
            r_proj = clsPoint2D.R();
            phi_proj = clsPoint2D.Phi();
            o2::math_utils::bringTo02PiGen(phi_proj);

            binR_proj = getRBinIndex(r_proj);
            binPhi_proj = getPhiBinIndex(phi_proj);

            int binRS, binPhiS;

            int binwRS = mLTFseed2BinWin;
            int binhwRS = binwRS / 2;

            int binwPhiS = mLTFseed2BinWin;
            int binhwPhiS = binwPhiS / 2;

            for (Int_t iR = 0; iR < binwRS; ++iR) {
              binRS = binR_proj + (iR - binhwRS);
              if (binRS < 0) {
                continue;
              }

              for (Int_t iPhi = 0; iPhi < binwPhiS; ++iPhi) {
                binPhiS = binPhi_proj + (iPhi - binhwPhiS);
                if (binPhiS < 0) {
                  continue;
                }

                binIndex2S = getBinIndex(binRS, binPhiS);
                mBinsS[layer1][layer2 - 1][binIndex1].emplace_back(binIndex2S);
              }
            }

            int binR, binPhi;

            int binwR = mLTFinterBinWin;
            int binhwR = binwR / 2;

            int binwPhi = mLTFinterBinWin;
            int binhwPhi = binwPhi / 2;

            for (Int_t iR = 0; iR < binwR; ++iR) {
              binR = binR_proj + (iR - binhwR);
              if (binR < 0) {
                continue;
              }

              for (Int_t iPhi = 0; iPhi < binwPhi; ++iPhi) {
                binPhi = binPhi_proj + (iPhi - binhwPhi);
                if (binPhi < 0) {
                  continue;
                }

                binIndex2 = getBinIndex(binR, binPhi);
                mBins[layer1][layer2 - 1][binIndex1].emplace_back(binIndex2);
              }
            }

          } // end loop layer2
        }   // end loop PhiBinIndex
      }     // end loop RBinIndex
    }       // end loop layer1

    std::cout << "Finished binning calculation!" << std::endl;

    std::cout << " mBinsS[1][1][3].size() = " << mBinsS[0][0][0].size() << std::endl;
    std::cout << " mBinsS[1][1][6].size() = " << mBinsS[0][0][1].size() << std::endl;
    std::cout << " mBinsS[1][1][22].size() = " << mBinsS[0][0][2].size() << std::endl;

    for (Int_t layer1 = 0; layer1 < (10 - 1); ++layer1) {
      for (Int_t layer2 = (layer1 + 1); layer2 < 10; ++layer2) {
        std::cout << "layer1 = " << layer1 << " ; layer2 = " << layer2 << std::endl;
        int binindex = 0;
        for (auto& bins : mBinsS[0][1]) {
          std::cout << "bin #" << binindex++ << ": ";
          for (auto& i : bins) {
            std::cout << i << " ";
          }
          std::cout << std::endl;
        }
      }
    }



    TFile* fout = new TFile("MFTFinderBinning.root", "RECREATE");
    fout->WriteObjectAny(&mBins, "std::array<std::array<std::vector<std::vector<Int_t>>, 9>, 9>", "mBins");
    fout->WriteObjectAny(&mBinsS, "std::array<std::array<std::vector<std::vector<Int_t>>, 9>, 9>", "mBinsS");
    fout->Close();
  }
};

} // namespace mft
} // namespace o2

void MFTFinderBinner(Int_t RBins = 50, Int_t PhiBins = 50)  {

o2::mft::Binner binner;
binner.MFTFinderBinner(RBins, PhiBins);

}
