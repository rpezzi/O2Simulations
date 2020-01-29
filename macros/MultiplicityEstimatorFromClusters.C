#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

using o2::itsmft::Cluster;
using o2::itsmft::ROFRecord;
using o2::itsmft::MC2ROFRecord;
using trackHasClustersinMFTDisks = std::array<bool,5>; // Records on which MFT disks a track has left clusters
using trackMap = std::map<Int_t, trackHasClustersinMFTDisks>; // Maps tracks to records of "on which disks this track has clusters"
using eventMap = std::map<Int_t, trackMap>; // Maps tracks to event
using sourceMap = std::map<Int_t, eventMap>; // Maps eventIDs to source
sourceMap trackabilityMap;
std::map<Int_t, Int_t> srcIDs;

void MultiplicityEstimatorFromClusters(const Char_t *ClustersFile = "mftclusters.root") {
  std::unique_ptr<TH1I> Trackablility = std::make_unique<TH1I> ("Trackablility", "In how many disks each track has clusters", 6, 0, 6);
  Trackablility->GetXaxis()->SetTitle("Number of disks");

  // Clusters
  TFile *clusFileIn = new TFile(ClustersFile);
  TTree* clusTree = (TTree*)clusFileIn->Get("o2sim");
  std::vector<Cluster>* mftcluster = nullptr;
  clusTree->SetBranchAddress("MFTCluster", &mftcluster);

  // ROFrecords
  std::vector<ROFRecord> rofRecVec, *rofRecVecP = &rofRecVec;
  clusTree->SetBranchAddress("MFTClustersROF", &rofRecVecP);

  // Cluster MC labels
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* clusLabArray = nullptr;
  std::vector<MC2ROFRecord> mc2rofVec, *mc2rofVecP = &mc2rofVec;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> labels, *plabels = &labels;
  clusTree->SetBranchAddress("MFTClustersMC2ROF", &mc2rofVecP);

  // ROFrames
  clusTree->GetEntry(0);
  int nROFRec = (int)rofRecVec.size();
  std::vector<int> mcEvMin(nROFRec, 0xFFFFFFF);
  std::vector<int> mcEvMax(nROFRec, -1);

  Int_t nMFTTrackables=0;
  Int_t nMFTTrackablesQED=0;
  Int_t nMFTTrackablesGenerator=0;
  Int_t nCountedLabels=0;

  // MFT Clusters and trackability
  //   - Identify trackable tracks (clusters in at least 4 disks)
  //   - Associate tracks to sourceID: Noise, o2sim, QED...


Int_t numberOfEntries = clusTree -> GetEntries();
std::cout << "numberOfEntries = " << numberOfEntries << std::endl;
for (Int_t entry=0; entry<numberOfEntries ; entry++) { // Loop over clusTree entries
    // std::cout << "Loop over entries in o2sim. Event = " << entry << std::endl;
    clusTree -> GetEntry(entry);
    Int_t nMFTClusters = mftcluster->size(); // Number of mft clusters in this entry
    Int_t nMFTClusterMCLabels = plabels->getNElements(); // Number of MFTClustersMCLabels in this Entry
    Int_t nMFTClusterIndexedSize = plabels->getIndexedSize(); // Number of original data indexed in this entry

    std::cout << "This entry has nMFTClusters = " << nMFTClusters << " ; nMFTClusterMCLabels = " << nMFTClusterMCLabels <<  " ; nMFTClusterIndexedSize = " << nMFTClusterIndexedSize << "\n";


    //std::cout << std::endl << "Loop over ROFRecords ...\n";
  //  for (int irof = 0; irof < nROFRec; irof++) {
  //      const auto& rofRec = rofRecVec[irof];
  //      auto firstcluster =  rofRec.getEntry().getIndex();
  //      auto nClustersinROF = rofRec.getNEntries();
  //      auto lastcluster = firstcluster + nClustersinROF;
        //rofRec.print();
        //std::cout << "Loop over MFTClusters...\n";
        for (Int_t n_cluster=0 ; n_cluster < nMFTClusters; n_cluster++) { // Loop over mftclusters in ROF
          Cluster* clusterp = &(*mftcluster).at(n_cluster);

          auto index = plabels->getMCTruthHeader(n_cluster).index;
          auto nextindex = plabels->getMCTruthHeader(n_cluster+1).index;
          //auto thissize = plabels->getMCTruthHeader(n_cluster+1).index - index;

          auto thissize = (n_cluster < nMFTClusterIndexedSize - 1) ? (nextindex - index) : (nMFTClusterMCLabels - index);

          //auto thisMClabel = plabels->getElement(nlabel);
          //std::cout << "Cluster # " << n_cluster << " index: " << index << " (size = " << thissize << ")" << " is from disk " << mftChipMapper.chip2Layer(clusterp->getSensorID())/2 << std::endl;
          for (auto label = index ; label < index+thissize ; label++ ) {
            auto sourceID =  plabels->getElement(label).getSourceID();
            auto eventID =  plabels->getElement(label).getEventID();
            auto trackID =  plabels->getElement(label).getTrackID();
            static o2::itsmft::ChipMappingMFT mftChipMapper;
            auto mftDisk =  mftChipMapper.chip2Layer(clusterp->getSensorID())/2;
            if (plabels->getElement(label).isValid())  trackabilityMap[sourceID][eventID][trackID][mftDisk] = true;

            srcIDs[sourceID]=srcIDs[sourceID]+1;
            nCountedLabels++;

            //std::cout << " MCLabel # " << label << ": " << plabels->getElement(label) << " sourceID = " << sourceID << " eventID = " << eventID << " trackID = " << trackID << " MFTLayer = " << mftChipMapper.chip2Layer(clusterp->getSensorID()) << " ROF = " << irof << std::endl;
            //if (n_cluster == 1757 | n_cluster == 1758 ) std::cout << n_cluster << " -> " << clusterp->getX() << " " << clusterp->getY() << " " << clusterp->getZ() << " " << std::endl;
          }
        }
    //}

    // Evaluate trackability
    for(auto source = trackabilityMap.begin(); source != trackabilityMap.end(); ++source) {
      //std::cout << " *** SourceID = "<< source->first << " \n";
      for(auto event = source->second.begin(); event != source->second.end(); ++event) {
        //std::cout << "   *** eventID = "<< event->first << " \n";
        for(auto track = event->second.begin(); track != event->second.end(); ++track) {
          //std::cout << "     *** trackID = "<< track->first << " \n";
           auto ndisks = 0;
           for (auto mftDisk = 0 ; mftDisk < 5 ; mftDisk++)
             ndisks += track->second[mftDisk];
             Trackablility->Fill(ndisks);
           if (ndisks >=4) {
             //std::cout << "     *** isTrackable!\n";
             nMFTTrackables++;
             if (source->first == 99 ) nMFTTrackablesQED++;
             if (source->first == 0 ) nMFTTrackablesGenerator++;
           }
        }
      }
    }
    std::cout << std::endl;
  //std::cout << "Finished entry " << entry << std::endl;
} // end loop over clusTree entries

std::cout << "\n|=========================|\n";
std::cout << "| " << std::setw(10) << "SourceIDs"  << " | " << std::setw(10) << "nClusters"  << " |" << "\n";
for(auto source = srcIDs.begin(); source != srcIDs.end(); ++source)
    std::cout << "| " << std::setw(10) << source->first  << " | " << std::setw(10) << source->second  << " |" << "\n";
std::cout << "|=========================|\n\n";


// Write histograms to file
//std::cout << "Writting histograms to file..." << std::endl;
TFile outFile("MultiplicityEstimatorFromClusters.root","RECREATE");
Trackablility->Write();

std::cout << "\n|===================================|\n";
std::cout << "| " << std::setw(15) << "SourceIDs"  << " | " << std::setw(15) << "MFT Trackables"  << " |" << "\n";
std::cout << "| " << std::setw(15) << "GENERATOR"  << " | " << std::setw(15) << nMFTTrackablesGenerator  << " |" << "\n";
std::cout << "| " << std::setw(15) << "QED Background"  << " | " << std::setw(15) << nMFTTrackablesQED  << " |" << "\n";
std::cout << "| " << std::setw(15) << "TOTAL:"  << " | " << std::setw(15) << nMFTTrackables  << " |" << "\n";
std::cout << "|===================================|\n";

std::cout << "Number of nCountedLabels = " << nCountedLabels << std::endl;


//new TBrowser;
}
