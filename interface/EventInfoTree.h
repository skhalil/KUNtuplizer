#include <string>
#include <vector>
#include <TTree.h>
#include <TLorentzVector.h>

class EventInfoTree {
  public:
    int    countEvents ;
    int    runno ;
    int    lumisec ;
    int    evtno;
    int    npv;   
    int    npuTrue;
    int    puBX;
    int    npuInt;
    int    nGoodVtx;
    std::vector<int>   ndofVtx;
    std::vector<float> chi2Vtx;
    std::vector<float> zVtx;
    std::vector<float> rhoVtx;
    std::vector<float> ptVtx;

    int    nSV;
    int    nGoodSV;
    std::vector<float> ptSV;
    std::vector<float> etaSV;
    std::vector<float> phiSV;
    std::vector<float> massSV;
    std::vector<float> chi2SV;
    std::vector<float> dxySV;
    std::vector<float> dxyerrSV;
    std::vector<float> dxysigSV;
    std::vector<float> d3dSV;
    std::vector<float> d3derrSV;
    std::vector<float> d3dsigSV;
    std::vector<float> costhetasvpvSV;
    std::vector<int>   ntrkSV;
    std::vector<int>   ndauSV;
    std::vector<int>   ndofSV;

    void clearTreeVectors() {
       ndofVtx.clear();
       chi2Vtx.clear();
       zVtx.clear();
       rhoVtx.clear();
       ptVtx.clear();
       
       ptSV.clear();
       etaSV.clear();
       phiSV.clear();
       massSV.clear();
       chi2SV.clear();
       dxySV.clear();
       dxyerrSV.clear();
       dxysigSV.clear();
       d3dSV.clear();
       d3derrSV.clear();
       d3dsigSV.clear();       
       costhetasvpvSV.clear();
       ntrkSV.clear();
       ndauSV.clear();
       ndofSV.clear();
    }

    void RegisterTree(TTree* tree, std::string name="SelectedEvents") {
      tree->Branch((name+"_countEvents").c_str(), &countEvents, (name+"_countEvents/I").c_str());
      tree->Branch((name+"_runno").c_str(),   &runno,    (name+"_runno/I").c_str());
      tree->Branch((name+"_lumisec").c_str(), &lumisec,  (name+"_lumisec/I").c_str());
      tree->Branch((name+"_evtno").c_str(),   &evtno,    (name+"_evtno/I").c_str());
      tree->Branch((name+"_npv").c_str(),     &npv,      "npv/I");
      tree->Branch((name+"_npuTrue").c_str(), &npuTrue,  "npuTrue/I");
      tree->Branch((name+"_puBX").c_str(),    &puBX,     "puBX/I");
      tree->Branch((name+"_npuInt").c_str(),  &npuInt,   "npuInt/I");      
      tree->Branch((name+"_nGoodVtx").c_str(),&nGoodVtx, (name+"_nGoodVtx/I").c_str());
      tree->Branch((name+"_ndofVtx").c_str(), &ndofVtx,  (name+"_ndofVtx/I").c_str());
      tree->Branch((name+"_chi2Vtx").c_str(), &chi2Vtx,  (name+"_chi2Vtx/F").c_str());
      tree->Branch((name+"_zVtx").c_str(),    &zVtx,     (name+"_zVtx/F").c_str());
      tree->Branch((name+"_rhoVtx").c_str(),  &rhoVtx,   (name+"_rhoVtx/F").c_str()); 
      tree->Branch((name+"_ptVtx").c_str(),   &ptVtx,    (name+"_ptVtx/F").c_str());
      
      tree->Branch((name+"_nSV").c_str(),     &nSV,      (name+"_nSV/I").c_str());
      tree->Branch((name+"_nGoodSV").c_str(), &nGoodSV,  (name+"_nGoodSV/I").c_str());
      tree->Branch((name+"_ptSV").c_str(),    &ptSV,     (name+"_ptSV/F").c_str());
      tree->Branch((name+"_etaSV").c_str(),   &etaSV,    (name+"_etaSV/F").c_str());
      tree->Branch((name+"_phiSV").c_str(),   &phiSV,    (name+"_phiSV/F").c_str());
      tree->Branch((name+"_massSV").c_str(),  &massSV,   (name+"_massSV/F").c_str());
      tree->Branch((name+"_chi2SV").c_str(),  &chi2SV,   (name+"_chi2SV/F").c_str()); 
      tree->Branch((name+"_dxySV").c_str(),   &dxySV,    (name+"_dxySV/F").c_str()); 
      tree->Branch((name+"_dxyerrSV").c_str(),&dxyerrSV, (name+"_dxyerrSV/F").c_str()); 
      tree->Branch((name+"_dxysigSV").c_str(),&dxysigSV, (name+"_dxysigSV/F").c_str());                   
      tree->Branch((name+"_d3dSV").c_str(),   &d3dSV,    (name+"_d3dSV/F").c_str()); 
      tree->Branch((name+"_d3derrSV").c_str(),&d3derrSV, (name+"_d3derrSV/F").c_str()); 
      tree->Branch((name+"_d3dsigSV").c_str(),&d3dsigSV, (name+"_d3dsigSV/F").c_str());
      tree->Branch((name+"_costhetasvpvSV").c_str(),&costhetasvpvSV, (name+"_costhetasvpvSV/F").c_str()); 
      tree->Branch((name+"_ntrkSV").c_str(),  &ntrkSV,   (name+"_ntrkSV/I").c_str());
      tree->Branch((name+"_ndauSV").c_str(),  &ndauSV,   (name+"_ndauSV/I").c_str());
      tree->Branch((name+"_ndofSV").c_str(),  &ndofSV,   (name+"_ndofSV/I").c_str());
    }
};
