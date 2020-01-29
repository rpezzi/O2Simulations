#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "DataFormatsITSMFT/ROFRecord.h"
#include "MFTTracking/TrackCA.h"
#include <FairLogger.h>
#include <TChain.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <TTree.h>
#include <vector>
#endif

void MFTTrackInspector(
    std::string outName =
        "MFTTrackInspector.root",           // name of the output binary file
    std::string inpName = "mfttracks.root", // name of the input MFT tracks
    std::string trackTreeName = "o2sim",    // name of the tracks tree
    std::string trackLTFBranchName = "MFTTrackLTF", // name of the tracks branch
    std::string trackCABranchName = "MFTTrackCA",   // name of the tracks branch

    std::string rofRecName =
        "MFTTracksROF" // name of the ROF records tree and its branch
) {
  TStopwatch swTot;
  swTot.Start();
  using ROFR = o2::itsmft::ROFRecord;
  using ROFRVEC = std::vector<o2::itsmft::ROFRecord>;

  ///-------> input
  TChain trackLTFTree(trackTreeName.c_str());
  TChain trackCATree(trackTreeName.c_str());
  TChain rofTree(rofRecName.c_str());

  trackLTFTree.AddFile(inpName.c_str());
  trackCATree.AddFile(inpName.c_str());
  rofTree.AddFile(inpName.c_str());

  // TrackLTF record entries in the track tree
  std::vector<o2::mft::TrackLTF> trackLTFVec, *trackLTFVecP = &trackLTFVec;
  if (!trackLTFTree.GetBranch(trackLTFBranchName.c_str())) {
    LOG(FATAL) << "Failed to find the branch " << trackLTFBranchName
               << " in the tree " << trackTreeName;
  }
  trackLTFTree.SetBranchAddress(trackLTFBranchName.c_str(), &trackLTFVecP);

  // TrackCA record entries in the track tree
  std::vector<o2::mft::TrackCA> trackCAVec, *trackCAVecP = &trackCAVec;
  if (!trackCATree.GetBranch(trackCABranchName.c_str())) {
    LOG(FATAL) << "Failed to find the branch " << trackCABranchName
               << " in the tree " << trackTreeName;
  }
  trackCATree.SetBranchAddress(trackCABranchName.c_str(), &trackCAVecP);

  // ROF record entries in the track tree
  ROFRVEC rofRecVec, *rofRecVecP = &rofRecVec;
  if (!rofTree.GetBranch(rofRecName.c_str())) {
    LOG(FATAL) << "Failed to find the branch " << rofRecName << " in the tree "
               << rofRecName;
  }
  rofTree.SetBranchAddress(rofRecName.c_str(), &rofRecVecP);

  ///-------< input

  // Loop on trackLTFs
  int lastTrackLTFtreeID = -1;
  long trackLTFoffs = 0;
  // LOG(INFO) << "trackLTFTree.GetEntries(): " << trackLTFTree.GetEntries();
  for (int i = 0; i < trackLTFTree.GetEntries(); i++) {
    trackLTFTree.GetEntry(i);
    if (trackLTFTree.GetTreeNumber() >
        lastTrackLTFtreeID) {       // this part is needed for chained input
      if (lastTrackLTFtreeID > 0) { // new chunk, increase the offset
        trackLTFoffs += trackLTFTree.GetTree()->GetEntries();
      }
      lastTrackLTFtreeID = trackLTFTree.GetTreeNumber();
    }

    // trackLTFTree-------------------------------------------------------------------------------<<<<
    int trackLTFNumber = 0, trackLTFmixed = 0;
    for (const auto &trackLTF : trackLTFVec) {
      // LOG(INFO) << "===== TrackLTF # " << trackLTFNumber << " nPoints = " <<
      // trackLTF.getNPoints() << " =====";
      auto thisTrackMCCompLabels = trackLTF.getMCCompLabels();
      auto firstTrackID = thisTrackMCCompLabels[0].getTrackID();
      for (auto iLabel = 0; iLabel < trackLTF.getNPoints(); iLabel++) {
        if (firstTrackID != thisTrackMCCompLabels[iLabel].getTrackID()) {
          // LOG(INFO) << "===== TrackLTF # " << trackLTFNumber << " TrackIDs
          // from different MCtracks!" ;
          trackLTFmixed++;
          break;
        }
      }
      trackLTFNumber++;
    }
    LOG(INFO) << "========== Mixed TrackLTFs = " << trackLTFmixed << " of "
              << trackLTFNumber << " reconstructed tracks ("
              << 100.0 * trackLTFmixed / trackLTFNumber << " %)";
  } // loop over multiple TrackLTF (in case of chaining)

  // Loop on trackCAs
  int lastTrackCAtreeID = -1;
  long trackCAoffs = 0;
  // LOG(INFO) << "trackCATree.GetEntries(): " << trackCATree.GetEntries();
  for (int i = 0; i < trackCATree.GetEntries(); i++) {
    trackCATree.GetEntry(i);
    if (trackCATree.GetTreeNumber() >
        lastTrackCAtreeID) {       // this part is needed for chained input
      if (lastTrackCAtreeID > 0) { // new chunk, increase the offset
        trackCAoffs += trackCATree.GetTree()->GetEntries();
      }
      lastTrackCAtreeID = trackCATree.GetTreeNumber();
    }

    // trackCATree-------------------------------------------------------------------------------<<<<
    int trackCANumber = 0, trackCAmixed = 0;
    for (const auto &trackCA : trackCAVec) {
      // LOG(INFO) << "===== TrackCA # " << trackCANumber << " nPoints = " <<
      // trackCA.getNPoints() << " =====";
      auto thisTrackMCCompLabels = trackCA.getMCCompLabels();
      auto firstTrackID = thisTrackMCCompLabels[0].getTrackID();
      for (auto iLabel = 0; iLabel < trackCA.getNPoints(); iLabel++) {
        if (firstTrackID != thisTrackMCCompLabels[iLabel].getTrackID()) {
          // LOG(INFO) << "===== TrackCA # " << trackCANumber << " TrackIDs from
          // different MCtracks!" ;
          trackCAmixed++;
          break;
        }
      }
      trackCANumber++;
    }
    LOG(INFO) << "========== Mixed TrackCAs = " << trackCAmixed << " of "
              << trackCANumber << " reconstructed tracks ("
              << 100.0 * trackCAmixed / trackCANumber << " %)";
  } // loop over multiple TrackCA (in case of chaining)

  //
  swTot.Stop();
  swTot.Print();
}
