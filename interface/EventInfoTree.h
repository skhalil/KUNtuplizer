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
    std::vector<float> energySV;
    std::vector<float> massSV;
    std::vector<int>   isMatchedToJetSV; 
    std::vector<float> chi2SV;
    std::vector<float> normChi2SV;
    std::vector<float> dxySV;
    std::vector<float> dxyerrSV;
    std::vector<float> dxysigSV;
    std::vector<float> d3dSV;
    std::vector<float> d3derrSV;
    std::vector<float> d3dsigSV;
    std::vector<float> costhetaSvPvSV;
    std::vector<int>   ntrkSV;
    std::vector<int>   ndauSV;
    std::vector<int>   ndofSV;
    std::vector<int>   isSVb;
    std::vector<int>   isSVc;
    std::vector<int>   isSVgen;
    std::vector<float> matchedGenPtSV;
    std::vector<float> matchedGenEtaSV;
    std::vector<float> matchedGenMassSV;
    std::vector<float> matchedGenPhiSV;

    void clearTreeVectors() {
       ndofVtx.clear();
       chi2Vtx.clear();
       zVtx.clear();
       rhoVtx.clear();
       ptVtx.clear();
       
       ptSV.clear();
       etaSV.clear();
       phiSV.clear();
       energySV.clear();
       massSV.clear();
       isMatchedToJetSV.clear();
       chi2SV.clear();
       normChi2SV.clear();
       dxySV.clear();
       dxyerrSV.clear();
       dxysigSV.clear();
       d3dSV.clear();
       d3derrSV.clear();
       d3dsigSV.clear();       
       costhetaSvPvSV.clear();
       ntrkSV.clear();
       ndauSV.clear();
       ndofSV.clear();
       isSVb.clear();
       isSVc.clear();
       isSVgen.clear();
       matchedGenPtSV.clear();
       matchedGenEtaSV.clear();
       matchedGenMassSV.clear();
       matchedGenPhiSV.clear();
    }

    void RegisterTree(TTree* tree, std::string name="Evts") {
      tree->Branch((name+"_countEvents").c_str(), &countEvents, (name+"_countEvents/I").c_str());
      tree->Branch((name+"_runno").c_str(),   &runno,    (name+"_runno/I").c_str());
      tree->Branch((name+"_lumisec").c_str(), &lumisec,  (name+"_lumisec/I").c_str());
      tree->Branch((name+"_evtno").c_str(),   &evtno,    (name+"_evtno/I").c_str());
      tree->Branch((name+"_npv").c_str(),     &npv,      (name+"_npv/I").c_str());
      tree->Branch((name+"_npuTrue").c_str(), &npuTrue,  (name+"_npuTrue/I").c_str());
      tree->Branch((name+"_puBX").c_str(),    &puBX,     (name+"_puBX/I").c_str());
      tree->Branch((name+"_npuInt").c_str(),  &npuInt,   (name+"_npuInt/I").c_str());      
      tree->Branch((name+"_nGoodVtx").c_str(),&nGoodVtx, (name+"_nGoodVtx/I").c_str());
      tree->Branch((name+"_ndofVtx").c_str(), &ndofVtx);
      tree->Branch((name+"_chi2Vtx").c_str(), &chi2Vtx);
      tree->Branch((name+"_zVtx").c_str(),    &zVtx);
      tree->Branch((name+"_rhoVtx").c_str(),  &rhoVtx); 
      tree->Branch((name+"_ptVtx").c_str(),   &ptVtx);
      
      tree->Branch((name+"_nSV").c_str(),             &nSV,      (name+"_nSV/I").c_str());
      tree->Branch((name+"_nGoodSV").c_str(),         &nGoodSV,  (name+"_nGoodSV/I").c_str());
      tree->Branch((name+"_ptSV").c_str(),            &ptSV);
      tree->Branch((name+"_etaSV").c_str(),           &etaSV);
      tree->Branch((name+"_phiSV").c_str(),           &phiSV);
      tree->Branch((name+"_energySV").c_str(),        &energySV);
      tree->Branch((name+"_massSV").c_str(),          &massSV);
      tree->Branch((name+"_isMatchedToJetSV").c_str(),&isMatchedToJetSV);
      tree->Branch((name+"_chi2SV").c_str(),          &chi2SV); 
      tree->Branch((name+"_normChi2SV").c_str(),      &normChi2SV);
      tree->Branch((name+"_dxySV").c_str(),           &dxySV); 
      tree->Branch((name+"_dxyerrSV").c_str(),        &dxyerrSV); 
      tree->Branch((name+"_dxysigSV").c_str(),        &dxysigSV);
      tree->Branch((name+"_d3dSV").c_str(),           &d3dSV); 
      tree->Branch((name+"_d3derrSV").c_str(),        &d3derrSV); 
      tree->Branch((name+"_d3dsigSV").c_str(),        &d3dsigSV);
      tree->Branch((name+"_costhetaSvPvSV").c_str(),  &costhetaSvPvSV); 
      tree->Branch((name+"_ntrkSV").c_str(),          &ntrkSV);
      tree->Branch((name+"_ndauSV").c_str(),          &ndauSV);
      tree->Branch((name+"_ndofSV").c_str(),          &ndofSV);
      tree->Branch((name+"_isSVb").c_str(),           &isSVb);
      tree->Branch((name+"_isSVc").c_str(),           &isSVc);
      tree->Branch((name+"_isSVgen").c_str(),         &isSVgen);
      tree->Branch((name+"_matchedGenPtSV").c_str(),  &matchedGenPtSV);
      tree->Branch((name+"_matchedGenEtaSV").c_str(), &matchedGenEtaSV);
      tree->Branch((name+"_matchedGenMassSV").c_str(),&matchedGenMassSV);
      tree->Branch((name+"_matchedGenPhiSV").c_str(), &matchedGenPhiSV);
   
    }
};
