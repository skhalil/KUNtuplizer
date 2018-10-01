#include <string>
#include <vector>
#include <TTree.h>
#include <TLorentzVector.h>

class GenInfoTree {
public:
   float genWt;
   std::vector<float> lheWts;
   std::vector<int> lheWtIDs;
   void clearTreeVectors() {
      lheWtIDs.clear();
      lheWts.clear();
   }
   void RegisterTree(TTree* tree, std::string name="EventWeights"){
    tree->Branch((name+"_genWt").c_str(), &genWt, (name+"_genWt/F").c_str());  
    tree->Branch((name+"_lheWts").c_str(), &lheWts, (name+"_lheWts/F").c_str());
    tree->Branch((name+"_lheWtIDs").c_str(), &lheWtIDs, (name+"_lheWtIDs/F").c_str());
   } 

};
