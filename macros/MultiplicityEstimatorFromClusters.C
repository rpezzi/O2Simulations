#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

void MultiplicityEstimatorFromClusters(const Char_t *ClustersFile = "mftclusters.root", bool debug = false) {

  using o2::itsmft::Cluster;
  using trackHasClustersinMFTDisks = std::array<bool,5>; // Disks with hits from a MFT track
  using trackMap = std::map<Int_t, trackHasClustersinMFTDisks>;
  using eventMap = std::map<Int_t, trackMap>;
  using sourceMap = std::map<Int_t, eventMap>;
  sourceMap trackabilityMap;

  // Clusters
  TFile *clusFileIn = new TFile(ClustersFile);
  TTree* clusTree = (TTree*)clusFileIn->Get("o2sim");
  std::vector<Cluster>* mftcluster = nullptr;
  clusTree->SetBranchAddress("MFTCluster", &mftcluster);

  // Cluster MC labels
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* clusLabArray = nullptr;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> labels, *plabels = &labels;
  clusTree->SetBranchAddress("MFTClusterMCTruth", &plabels);

  // Get the only TTree entry :(
  clusTree->GetEntry(0);

  Int_t nMFTTrackables=0;
  Int_t nMFTTrackablesQED=0;
  Int_t nMFTTrackablesHIJING=0;
  Int_t nCountedLabels=0;

  //## MFT Clusters and trackability
  //   - Identify trackable tracks (clusters in at least 4 disks)
  //   - Associate tracks to sourceID: Noise, o2sim, QED...
  //###

  Int_t nMFTClusters = mftcluster->size(); // Number of mft clusters in this entry
  
  Int_t nMFTClusterMCLabels = plabels->getNElements(); // Number of MFTClustersMCLabels in this Entry
  Int_t nMFTClusterIndexedSize = plabels->getIndexedSize(); // Number of original data indexed in this entry
  
  std::cout << "This entry has nMFTClusters = " << nMFTClusters << " ; nMFTClusterMCLabels = " << nMFTClusterMCLabels <<  " ; nMFTClusterIndexedSize = " << nMFTClusterIndexedSize << "\n";

  std::map<Int_t, Int_t> srcIDs;
  std::map<Int_t, Int_t> evtIDs;

  for(vector<Cluster>::iterator itCluster = mftcluster->begin(); itCluster!= mftcluster->end(); ++itCluster){ // Loop over mftclusters in ROF
    Int_t n_cluster = std::distance(mftcluster->begin(), itCluster);
    if(debug) 
      std::cout << "Looking at cluster "<< n_cluster << std::endl; 

    auto index = plabels->getMCTruthHeader(n_cluster).index;
    auto nextindex = plabels->getMCTruthHeader(n_cluster+1).index;
    auto thissize = (n_cluster < nMFTClusterIndexedSize - 1) ? (nextindex - index) : (nMFTClusterMCLabels - index);
  
    for(auto label = index ; label < index+thissize ; label++ ) {
       auto sourceID =  plabels->getElement(label).getSourceID();
       auto eventID =  plabels->getElement(label).getEventID();
       auto trackID =  plabels->getElement(label).getTrackID();
       static o2::itsmft::ChipMappingMFT mftChipMapper;
       auto mftDisk =  mftChipMapper.chip2Layer(itCluster->getSensorID())/2;
       if (plabels->getElement(label).isValid())  trackabilityMap[sourceID][eventID][trackID][mftDisk] = true;
  
      srcIDs[sourceID]=srcIDs[sourceID]+1;
      evtIDs[eventID]=evtIDs[eventID]+1;
      nCountedLabels++;

      if(debug) {  
         std::cout << "    - MCLabel # " << label << ": " << plabels->getElement(label) << std::endl 
                   << "    - sourceID = " << sourceID << std::endl 
                   << "    - eventID = " << eventID << std::endl 
                   << "    - trackID = " << trackID << std::endl 
                   << "    - MFTLayer = " << mftChipMapper.chip2Layer(itCluster->getSensorID()) << std::endl; 
      }
    }//end of label loop
   }//end of cluster loop
  
  // Evaluate trackability
  for(auto source = trackabilityMap.begin(); source != trackabilityMap.end(); ++source) {
    if(debug) 
       std::cout << " *** SourceID = "<< source->first << " \n";
    for(auto event = source->second.begin(); event != source->second.end(); ++event) {
      if(debug) 
         std::cout << "   *** eventID = "<< event->first << " \n";
      for(auto track = event->second.begin(); track != event->second.end(); ++track) {
        if(debug) 
          std::cout << "     *** trackID = "<< track->first << " \n";
          auto ndisks = 0;
          for (auto mftDisk = 0 ; mftDisk < 5 ; mftDisk++)
            ndisks += track->second[mftDisk];
          if (ndisks >=4) {
             if(debug) std::cout << "     *** isTrackable!\n";
             nMFTTrackables++;
             if (source->first == 99 ) nMFTTrackablesQED++;
             if (source->first == 0 ) nMFTTrackablesHIJING++;
         }
      }//end of trackID loop
    }//end of eventID loop
  }//end of sourceID loop
// 
//  // Write histograms to file
//  std::cout << "Writting histograms to file..." << std::endl;
//  TFile outFile("MultiplicityEstimatorFromClusters.root","RECREATE");
//  Trackablility->Write();
//  MultiplicityDistrib->Write();
//  outFile.Close();
//  
  //Int_t totalRecoMFTTracks = nCleanTracksLTF + nCleanTracksCA + nInvalidTracksLTF + nInvalidTracksCA;
  std::cout << "Number of MFT trackables (GENERATOR) = " << nMFTTrackablesHIJING << std::endl;
  std::cout << "Number of MFT trackables (QED Background) = " << nMFTTrackablesQED << std::endl;
  std::cout << "Number of MFT trackables (Total) = " << nMFTTrackables << std::endl;
  std::cout << "Number of nCountedLabels = " << nCountedLabels << std::endl;
  
  
  //new TBrowser;
}
