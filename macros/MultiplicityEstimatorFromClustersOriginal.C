#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

void MultiplicityEstimatorFromClustersOriginal(const Char_t *ClustersFile = "mftclusters.root", Int_t nentries = -1, bool debug = false) {
  std::unique_ptr<TH1I> Trackablility = std::make_unique<TH1I> ("Trackablility", "In how many disks the tracks has clusters", 6, 0, 6);
  Trackablility->GetXaxis()->SetTitle("Number of disks");
  std::unique_ptr<TH1I> MultiplicityDistrib = std::make_unique<TH1I> ("Multiplicity", "MFT Trackables Distribution", 10000, 0, 10000);
  MultiplicityDistrib->GetXaxis()->SetTitle("Trackable Multiplicity");

  using o2::itsmft::Cluster;
  using o2::itsmft::ROFRecord;
  using o2::itsmft::MC2ROFRecord;
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

  // ROFrecords
  //## Deprecated because of the change in data format
  //std::vector<ROFRecord> rofRecVec, *rofRecVecP = &rofRecVec;
  //TTree* ROFRecTree = (TTree*) clusFileIn->Get("MFTClustersROF");
  //ROFRecTree->SetBranchAddress("MFTClustersROF", &rofRecVecP);
  //###
  std::vector<ROFRecord>* rofRecVec = nullptr;
  clusTree->SetBranchAddress("MFTClustersROF", &rofRecVec);

  // Cluster MC labels
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* clusLabArray = nullptr;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> labels, *plabels = &labels;
  clusTree->SetBranchAddress("MFTClusterMCTruth", &plabels);
  //## Deprecated because of the change in data format
  //std::vector<MC2ROFRecord> mc2rofVec, *mc2rofVecP = &mc2rofVec;
  //TTree* MC2ROFRecTree = nullptr;
  //MC2ROFRecTree = (TTree*)clusFileIn->Get("MFTClustersMC2ROF");
  //MC2ROFRecTree->SetBranchAddress("MFTClustersMC2ROF", &mc2rofVecP);
  //###
  std::vector<MC2ROFRecord>* mc2rofVec = nullptr;
  clusTree->SetBranchAddress("MFTClustersMC2ROF", &mc2rofVec);


  // ROFrames
  clusTree->GetEntry(0);
  int nROFRec = (int)rofRecVec->size();
  if(debug) 
    std::cout << "Reconstructed ROF size = " << nROFRec << std::endl;

  std::vector<int> mcEvMin(nROFRec, 0xFFFFFFF);
  std::vector<int> mcEvMax(nROFRec, -1);

  Int_t nMFTTrackables=0;
  Int_t nMFTTrackablesQED=0;
  Int_t nMFTTrackablesHIJING=0;
  Int_t nCountedLabels=0;

  //## MFT Clusters and trackability
  //   - Identify trackable tracks (clusters in at least 4 disks)
  //   - Associate tracks to sourceID: Noise, o2sim, QED...
  //###


  Int_t numberOfEntries;
  if(nentries == -1)
     numberOfEntries = clusTree -> GetEntries();
  else
     numberOfEntries = nentries;

  std::cout << "numberOfEntries = " << numberOfEntries << std::endl;

  for (Int_t entry=0; entry<numberOfEntries ; entry++) { // Loop over clusTree entries
      std::cout << "Loop over entries in o2sim. Event = " << entry << std::endl;

      clusTree -> GetEntry(entry);
      Int_t nMFTClusters = mftcluster->size(); // Number of mft clusters in this entry
  
      Int_t nMFTClusterMCLabels = plabels->getNElements(); // Number of MFTClustersMCLabels in this Entry
      Int_t nMFTClusterIndexedSize = plabels->getIndexedSize(); // Number of original data indexed in this entry
  
      if(debug)
         std::cout << "This entry has nMFTClusters = " << nMFTClusters << " ; nMFTClusterMCLabels = " << nMFTClusterMCLabels <<  " ; nMFTClusterIndexedSize = " << nMFTClusterIndexedSize << "\n";

     /*
      for (auto nlabel = 0 ; nlabel < nMFTClusterMCLabels ; nlabel++)
      {
        auto index = plabels->getMCTruthHeader(nlabel).index;
        auto thissize = plabels->getMCTruthHeader(nlabel+1).index - index;
        auto thisMClabel = plabels->getElement(nlabel);
        //std::cout << "Label # " << nlabel << " index: " << index << " (size = " << thissize << " getElement(nlabel) " <<  thisMClabel << std::endl;
      } */
  
  
      if(debug)
        std::cout << std::endl << "Loop over ROFRecords ...\n";

      for (int irof = 0; irof < nROFRec; irof++) {
          std::map<Int_t, Int_t> srcIDs;
          std::map<Int_t, Int_t> evtIDs;
          //const auto& rofRec = rofRecVec->at(irof);
          const auto& rofRec = (*rofRecVec)[irof];
          //## Deprecated due to change in data format
          //auto firstcluster =  rofRec.getROFEntry().getIndex();
          //auto nClustersinROF = rofRec.getNROFEntries();
          //auto lastcluster = firstcluster + nClustersinROF;
          //rofRec.print();
          //####
          auto firstcluster =  rofRec.getFirstEntry();
          auto nClustersinROF = rofRec.getNEntries();
          auto lastcluster = firstcluster + nClustersinROF;
          if(debug)
            std::cout << "Loop over MFTClusters in ROF " << irof << std::endl;

          for(Int_t n_cluster=firstcluster ; n_cluster < lastcluster; n_cluster++) { // Loop over mftclusters in ROF
            if(debug) 
              std::cout << "Looking at cluster "<< n_cluster << " in this ROF" << std::endl; 
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
              evtIDs[eventID]=evtIDs[eventID]+1;
              nCountedLabels++;

              if(debug) {  
                 std::cout << "    - MCLabel # " << label << ": " << plabels->getElement(label) << std::endl 
                           << "    - sourceID = " << sourceID << std::endl 
                           << "    - eventID = " << eventID << std::endl 
                           << "    - trackID = " << trackID << std::endl 
                           << "    - MFTLayer = " << mftChipMapper.chip2Layer(clusterp->getSensorID()) << std::endl 
                           << "    - ROF = " << irof << std::endl;
                 //if (n_cluster == 1757 | n_cluster == 1758 ) std::cout << n_cluster << " -> " << clusterp->getX() << " " << clusterp->getY() << " " << clusterp->getZ() << " " << std::endl;
              }
            }
          }
//  
//  
//              /*
//              std::cout << "Source IDs: " ;
//              for(map<int,int>::iterator it = srcIDs.begin(); it != srcIDs.end(); ++it) {
//                std::cout << it->first << "(" << it->second << ") ";
//              }
//              std::cout << std::endl;
//  
//              std::cout << "eventIDs: " ;
//              for(map<int,int>::iterator it = evtIDs.begin(); it != evtIDs.end(); ++it) {
//                std::cout << it->first << "(" << it->second << ") ";
//              }
//              std::cout << std::endl << std::endl;
//              */
//  
      }
  
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
          }
        }
      }
      if(debug) std::cout << "Finished entry " << entry << std::endl;
  } // end loop over clusTree entries

  
  
  // Write histograms to file
  std::cout << "Writting histograms to file..." << std::endl;
  TFile outFile("MultiplicityEstimatorFromClusters.root","RECREATE");
  Trackablility->Write();
  MultiplicityDistrib->Write();
  outFile.Close();
  
  //Int_t totalRecoMFTTracks = nCleanTracksLTF + nCleanTracksCA + nInvalidTracksLTF + nInvalidTracksCA;
  std::cout << "Number of MFT trackables (GENERATOR) = " << nMFTTrackablesHIJING << std::endl;
  std::cout << "Number of MFT trackables (QED Background) = " << nMFTTrackablesQED << std::endl;
  std::cout << "Number of MFT trackables (Total) = " << nMFTTrackables << std::endl;
  std::cout << "Number of nCountedLabels = " << nCountedLabels << std::endl;
  
  
  //new TBrowser;
}
