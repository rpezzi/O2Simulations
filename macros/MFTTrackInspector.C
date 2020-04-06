#include "DataFormatsITSMFT/ROFRecord.h"
#include "MFTTracking/TrackCA.h"
#include <FairLogger.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <TTree.h>
#include <vector>

bool DEBUG_VERBOSE = false;


void MFTTrackInspector(std::string path = "./",
    std::string outName =
        "MFTTrackInspector.root",           // name of the output binary file
    std::string inpName = "mfttracks.root", // name of the input MFT tracks
    std::string trackTreeName = "o2sim",    // name of the tracks tree
    std::string trackLTFBranchName = "MFTTrackLTF", // name of the tracks branch
    std::string trackCABranchName = "MFTTrackCA",   // name of the tracks branch

    std::string rofRecName =
        "MFTTracksROF" // name of the ROF records tree and its branch
)
{
  if (path.back() != '/') {
    path += '/';
  }

  std::unique_ptr<TFile> tracFile(TFile::Open((path + inpName).c_str()));
  if (!tracFile || tracFile->IsZombie()) {
    LOG(ERROR) << "Failed to open input tracks file " << (path + inpName);
    return;
  }

  TTree* tracTree = (TTree*)tracFile->Get("o2sim");
  if (!tracTree) {
    LOG(ERROR) << "Failed to get tracks tree";
    return;
  }

  // TrackLTF record entries in the track tree
  std::vector<o2::mft::TrackLTF> trackLTFVec, *trackLTFVecP = &trackLTFVec;
  tracTree->SetBranchAddress("MFTTrackLTF", &trackLTFVecP);

  // TrackCA record entries in the track tree
  std::vector<o2::mft::TrackCA> trackCAVec, *trackCAVecP = &trackCAVec;
  tracTree->SetBranchAddress("MFTTrackCA", &trackCAVecP);

  // ROF record entries in the track tree
  std::vector<o2::itsmft::ROFRecord>* rofRecVec = nullptr;
  tracTree->SetBranchAddress("MFTTracksROF", &rofRecVec);

  tracTree->GetEntry(0);
  int nTracksLTF =  trackLTFVecP->size();
  int nTracksCA = trackCAVecP->size();
  printf("Tracks LTF %d CA %d \n", nTracksLTF, nTracksCA);

  // Loop on trackLTFs
  int lastTrackLTFtreeID = -1;
  long trackLTFoffs = 0;
  // LOG(INFO) << "trackLTFTree.GetEntries(): " << trackLTFTree.GetEntries();
  for (int i = 0; i < nTracksLTF; i++) {

    // trackLTFTree----<<<<
    int trackLTFNumber = 0, trackLTFmixed = 0;
    for (const auto &trackLTF : trackLTFVec) {
      if(DEBUG_VERBOSE) LOG(INFO) << "===== TrackLTF # " << trackLTFNumber << " nPoints = " <<  trackLTF.getNPoints() << " =====";
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
  }

  // Loop on trackCAs
  int lastTrackCAtreeID = -1;
  long trackCAoffs = 0;
  // LOG(INFO) << "trackCATree.GetEntries(): " << trackCATree.GetEntries();
  for (int i = 0; i < nTracksCA; i++) {

    // trackCATree----<<<<
    int trackCANumber = 0, trackCAmixed = 0;
    for (const auto &trackCA : trackCAVec) {
      if(DEBUG_VERBOSE) LOG(INFO) << "===== TrackCA # " << trackCANumber << " nPoints = " <<
      trackCA.getNPoints() << " =====";
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
  }

}
