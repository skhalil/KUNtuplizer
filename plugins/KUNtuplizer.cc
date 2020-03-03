// -*- C++ -*-
//
// Package:    Framework/KUNtuplizer
// Class:      KUNtuplizer
// 
/**\class KUNtuplizer KUNtuplizer.cc Framework/KUNtuplizer/plugins/KUNtuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Sadia Khalil
//         Created:  Fri, 21 Sep 2018 22:01:49 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include "Framework/KUNtuplizer/interface/EventInfoTree.h"
#include "Framework/KUNtuplizer/interface/GenInfoTree.h"
#include "Framework/KUNtuplizer/interface/ElectronTree.h"
//#include "Framework/KUNtuplizer/interface/MuonTree.h"
//#include "Framework/KUNtuplizer/interface/METTree.h"
//#include "Framework/KUNtuplizer/interface/JetTree.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <Math/VectorUtil.h>
//
// class declaration
//

class KUNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit KUNtuplizer(const edm::ParameterSet&);
      ~KUNtuplizer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   protected:
      float catchInfsAndBound(const float& in,const float& replace_value, const float& lowerbound, const float& upperbound);
      const float& catchInfs (const float& in,const float& replace_value);     

   private:
      virtual void beginRun(edm::Run const&, edm::EventSetup const&); 
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endRun(edm::Run const&, edm::EventSetup const&); 
      virtual void endJob() override;

      // ----------member data ---------------------------  
   //edm::EDGetTokenT<std::vector<pat::PackedGenParticle>>         genPartsToken_;
   edm::EDGetTokenT<std::vector<reco::GenParticle>>              genPartsToken_;
   edm::EDGetTokenT<GenEventInfoProduct>                         genToken_;
   edm::EDGetTokenT<LHEEventProduct>                             genlheToken_;
   edm::InputTag                                                 puInfo_;
   edm::EDGetTokenT<std::vector<reco::Vertex>>                   vtxToken_;
   edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> svToken_;
   edm::EDGetTokenT<std::vector<pat::Electron>>                  elecsToken_;
   edm::EDGetTokenT<pat::PackedCandidateCollection>              packedPFCands_;
   edm::EDGetTokenT<edm::View<pat::Jet> >                        ak4jetsToken_;
   edm::EDGetTokenT< reco::ConversionCollection >                conv_;
   edm::EDGetTokenT< reco::BeamSpot >                            beamSpot_;
   edm::EDGetTokenT< double >                                    rho_;
   edm::Service<TFileService>                                    fs_;
   std::map<std::string, TH1D*>                                  h1_; 
   TTree*        tree_;    
   GenInfoTree   genevt_; 
   EventInfoTree evt_;
   ElectronTree  ele_;
   int           ndofPV_;
   double        zPV_;
   double        rhoPV_;
   double        ptMaxSV_;
   double        dlenSigSV_;
   double        d2dSV_;
   double        cosPdotVSV_;
   int           ndauSV_;
   double        jetdRSV_;
   double        genPartdRSV_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
KUNtuplizer::KUNtuplizer(const edm::ParameterSet& iConfig):
   //genPartsToken_  (consumes<std::vector<pat::PackedGenParticle>>(iConfig.getParameter<edm::InputTag>("genParts"))),
   genPartsToken_  (consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParts"))),
   genToken_       (consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
   genlheToken_    (consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer",""))),
   puInfo_         (iConfig.getParameter<edm::InputTag>("puInfo")),
   vtxToken_       (consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
   svToken_        (consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secVertices"))),
   elecsToken_     (consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
   packedPFCands_  (consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedPFCands"))), 
   ak4jetsToken_   (consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),
   conv_           (consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversion"))),
   beamSpot_       (consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))), 
   rho_            (consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
   ndofPV_         (iConfig.getParameter<int>("ndofPV")),
   zPV_            (iConfig.getParameter<double>("zPV")),
   rhoPV_          (iConfig.getParameter<double>("rhoPV")),
   ptMaxSV_        (iConfig.getParameter<double>("ptMaxSV")), 
   dlenSigSV_      (iConfig.getParameter<double>("dlenSigSV")),
   d2dSV_          (iConfig.getParameter<double>("d2dSV")),
   cosPdotVSV_     (iConfig.getParameter<double>("cosPdotVSV")),
   ndauSV_         (iConfig.getParameter<int>("ndauSV")),
   jetdRSV_        (iConfig.getParameter<double>("jetdRSV")),
   genPartdRSV_    (iConfig.getParameter<double>("genPartdRSV"))
                    
{
   consumes<std::vector<PileupSummaryInfo>>(puInfo_);
   //now do what ever initialization is needed
   usesResource("TFileService");
   
}

KUNtuplizer::~KUNtuplizer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
KUNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm; 
   genevt_.clearTreeVectors(); 
   evt_.clearTreeVectors(); 

   //Generator weights
   float genwt = 1.0;
   edm::Handle<GenEventInfoProduct> genInfo;
   iEvent.getByToken( genToken_,genInfo);
   
   if (genInfo.isValid()){
      genwt = genInfo->weight();
      genwt /= std::abs(genwt);   
   }
   genevt_.genWt = genwt;
   
   //LHE weight
   edm::Handle<LHEEventProduct> lheInfo;
   iEvent.getByToken(genlheToken_, lheInfo);
   
   if(lheInfo.isValid()) {
      int size = lheInfo->weights().size();
      double asdd = 1.;
      if (size != 0){
         //asdd = lheInfo->weights()[0].wgt;
         asdd = lheInfo->originalXWGTUP();
      }
      
      genevt_.lheWtIDs.reserve(size);
      genevt_.lheWts.reserve(size);
      for(unsigned int i=0; i<lheInfo->weights().size(); ++i) {         
         double asdde =lheInfo->weights()[i].wgt;
         int asddeID = std::stoi(lheInfo->weights()[i].id);
         genevt_.lheWtIDs.push_back(asddeID);
         genevt_.lheWts.push_back(asdde/asdd);
      }
   }

   // Basic event info
   evt_.runno = iEvent.eventAuxiliary().run();
   evt_.lumisec = iEvent.eventAuxiliary().luminosityBlock();
   evt_.evtno = iEvent.eventAuxiliary().event();
   evt_.countEvents = 1.;//test

   //Pileup
   edm::Handle<std::vector< PileupSummaryInfo > > puInfo;
   iEvent.getByLabel(puInfo_, puInfo);
   std::vector<PileupSummaryInfo>::const_iterator pvi;
   for(pvi = puInfo->begin(); pvi != puInfo->end(); ++pvi) {
      evt_.npuTrue = pvi->getTrueNumInteractions(); 
      evt_.npuInt = pvi->getBunchCrossing(); 
      evt_.puBX = pvi->getPU_NumInteractions();
   }
   
   // Vertex
   edm::Handle<std::vector<reco::Vertex>> vertices;
   iEvent.getByToken(vtxToken_, vertices);

   if (vertices->empty()) return;
   evt_.npv=0; evt_.nGoodVtx = 0;
   evt_.npv = vertices->size();
   
   int prVtx = -1;
   // store the vertex related variables after basic PV selection
   for (size_t i = 0; i < vertices->size(); i++) {
      if (vertices->at(i).isFake()) continue;
      if (vertices->at(i).ndof() <= ndofPV_) continue;
      if (fabs(vertices->at(i).position().z()) > zPV_) continue;
      if (vertices->at(i).position().rho() > rhoPV_) continue;    
      if (prVtx < 0) {prVtx = i; }
      evt_.ndofVtx.push_back(vertices->at(i).ndof());
      evt_.chi2Vtx.push_back(vertices->at(i).chi2());
      evt_.zVtx.push_back(vertices->at(i).position().z());
      evt_.rhoVtx.push_back(vertices->at(i).position().rho());
      evt_.etaVtx.push_back(vertices->at(i).position().eta());
      evt_.phiVtx.push_back(vertices->at(i).position().phi());
      evt_.nGoodVtx++;
   }
   
   if (evt_.nGoodVtx==0) return;
   //auto primaryVertex=vertices->at(prVtx);   
   const auto & pv = (*vertices)[0];

   // Jets
   edm::Handle<edm::View<pat::Jet> > ak4jets;
   iEvent.getByToken(ak4jetsToken_, ak4jets);

   // GenParticles
   //edm::Handle<std::vector<pat::PackedGenParticle>> genParts;
   edm::Handle<std::vector<reco::GenParticle>> genParts;
   iEvent.getByToken(genPartsToken_, genParts);
   
   // Secondary Vertices
   edm::Handle<reco::VertexCompositePtrCandidateCollection> secVertices;
   iEvent.getByToken(svToken_, secVertices);
   evt_.nSV = 0; evt_.nGoodSV = 0;
   evt_.nSV = secVertices->size();
    
   VertexDistance3D vdist;  
   VertexDistanceXY vdistXY;
   
   for (const auto & sv : *secVertices) {

      h1_["cutflow"] -> Fill(1, genwt);
 
      //sum pt of all tracks in SV
      if(sv.pt() > ptMaxSV_) continue;//20
      h1_["cutflow"] -> Fill(2, genwt);



      bool isSVInJetDir = false;  float dRSVJet = -100.0; float dRSVFlightJet = -100.0;    
      for (const pat::Jet & jet : *ak4jets ){
         GlobalVector svDir(sv.px(), sv.py(), sv.pz());
         GlobalVector flightDir(sv.vertex().x() - pv.x(), sv.vertex().y() - pv.y(), sv.vertex().z() - pv.z());
         GlobalVector jetDir(jet.px(),jet.py(),jet.pz());
         dRSVJet = reco::deltaR(svDir, jetDir );
         dRSVFlightJet = reco::deltaR(flightDir, jetDir);
         if(dRSVFlightJet <=  jetdRSV_) {//0.4
            isSVInJetDir = true; break;}
      }
      bool cut3 = !isSVInJetDir; 
      if(cut3)                      {h1_["cutflow"] -> Fill(3, genwt);         }
      else continue;
      
      // rest of the selections will be stored as cutflow but will be not applied to store the ntuple branches
      // ------------------------------

      //store the branches after gen particle matching 
      float genSVdR = -1000.0; float matchGenSVdR = -1000.0;
      int matchGenMomID = 0; bool matchGenIsHardProcess = false; bool matchGenIsLastCopyBeforeFSR = false;     
      float matchGenPt = -1000.0; float matchGenEta = -1000.0; float matchGenMass = -1000.0; 
      float matchGenPhi = -1000.0; float matchGenEnergy = -1000.0;
      int nbs=0; int ngs=0; int ncs=0; int nls=0; int nothers=0; 
      bool b = false; bool c = false; bool l = false; bool g = false;
      bool o = false; bool p = false; bool bg = false; bool cg = false;         
      bool isbMeson = false; bool isbBaryon = false; 
      bool iscMeson = false; bool iscBaryon = false;
      bool islMeson = false; bool islBaryon = false;      
      bool fill1stHeavyGen = false;
      bool fill1stNotHeavyGen = false;
      float mindR = 100.0; float maxdR = 0.0;
      
      for (const reco::GenParticle & gp : *genParts) {
         genSVdR = reco::deltaR(gp.p4(), sv.p4());
         if (genSVdR <= mindR) { mindR = genSVdR;}
         if (genSVdR > maxdR)  { maxdR = genSVdR;}
         
         if (genSVdR <= genPartdRSV_){ 
            p = true; 
            
            //std::cout << "SV with pt "<< sv.p4().pt() << " is matched to p = " << gp.pdgId() << std::endl;
            //std::cout << "matchGenSVdR : " << matchGenSVdR << std::endl;

            int pdgid = abs(gp.pdgId());

            isbMeson = pdgid%1000 > 500 && pdgid%1000 < 600;
            isbBaryon = pdgid > 5000 &&  pdgid < 6000;
            
            iscMeson = pdgid%1000 > 400 && pdgid%1000 < 500;
            iscBaryon = pdgid > 4000 &&  pdgid < 5000;

            islMeson = pdgid%1000 > 100 && pdgid%1000 < 400;
            islBaryon = pdgid > 1000 &&  pdgid < 4000;
                  
            if      (pdgid == 5  || isbMeson || isbBaryon){b=true;nbs++;}
            else if (pdgid == 4  || iscMeson || iscBaryon){c=true;ncs++;} 
            else if (pdgid == 21)                         {g=true; ngs++;} // 9 is excluded
            else if (pdgid == 1  || pdgid == 2 || pdgid == 3 || islMeson || islBaryon) {l=true;nls++;}
            else    {nothers++;}
            
            if (nbs+ncs == 0 && fill1stNotHeavyGen == false){// if matched to a non-b particle, store the P4 of the leading non-b particle
               //std::cout << "leading not gen P4 stored " << std::endl;
               
               matchGenSVdR   = genSVdR; 
               matchGenPt     = gp.p4().pt();               
               matchGenEta    = gp.p4().eta();
               matchGenPhi    = gp.p4().phi();
               matchGenMass   = gp.p4().mass();
               matchGenEnergy = gp.p4().energy();
               matchGenMomID               = gp.mother()->pdgId();               
               matchGenIsHardProcess       = gp.isHardProcess();              
               matchGenIsLastCopyBeforeFSR = gp.statusFlags().isLastCopyBeforeFSR();               
               fill1stNotHeavyGen = true;
               
            }
            else if (nbs+ncs != 0 && fill1stHeavyGen == false) {// if a b or c quark is found, then overwrite the P4 to the leading one
               //std::cout << "leading heavy quark P4 stored " << std::endl;
               matchGenSVdR   = genSVdR; 
               matchGenPt     = gp.p4().pt();
               matchGenEta    = gp.p4().eta();
               matchGenPhi    = gp.p4().phi();
               matchGenMass   = gp.p4().mass();
               matchGenEnergy = gp.p4().energy();
               matchGenMomID               = gp.mother()->pdgId();
               matchGenIsHardProcess       = gp.isHardProcess();
               matchGenIsLastCopyBeforeFSR = gp.statusFlags().isLastCopyBeforeFSR();
               fill1stHeavyGen = true;
            }
         }
      }  
      //std::cout << "loop end" <<std::endl;
      bool cut4=p; bool cut5 = (b || c);   
      if(p){
         h1_["cutflow"] -> Fill(4, genwt);
         if(cut5){ h1_["cutflow"] -> Fill(5, genwt);} 
         if(b && g) {bg = true;} 
         if(c && g) {cg = true;}
         //if(ncs>0)    {c = true;}
         //if(nls>0)    {l = true;} 
         if(nbs+ncs+nls == 0){o = true;}            
      }
      //else continue;

      //std::cout << "end the gen loop ----> " << std::endl;
      //std::cout << "genSVdR : " << genSVdR << ", matchGenSVdR :" << matchGenSVdR << std::endl;
      //std::cout << "mindR : " << mindR <<", maxdR : " << maxdR << std::endl;
      //std::cout << "nbs " << nbs <<", ncs " << ncs << ", nls " << nls << ", ngs " << ngs << std::endl;
      //std::cout << "b " << b << ", bg " <<bg << ", c " << c <<", l " << l << ", other " << o << std::endl;
      //SV track separation from jets     

      
      //distance in the transverse plane between the SV and PV
      Measurement1D d2d= vdistXY.distance(pv, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error()))); 
      bool cut6 = d2d.value() < d2dSV_;
      if(cut4 && cut5 && cut6)              {h1_["cutflow"] -> Fill(6, genwt);} //3
      
      //pointing angle (i.e. the angle between the sum of the momentum of the 
      //tracks in the SV and the flight direction betwen PV and SV)  
      double dx = (sv.vx() - pv.x()), dy = (sv.vy() - pv.y()), dz = (sv.vz() - pv.z());    
      double cosPdotV = (dx * sv.px() + dy*sv.py() + dz*sv.pz())/(sv.p()*std::sqrt(dx * dx + dy * dy + dz * dz));
      bool cut7 =  cosPdotV > cosPdotVSV_;    
      if(cut4 && cut5 && cut6 && cut7)      {h1_["cutflow"] -> Fill(7, genwt);}//0.98
      
      //Number of tracks
      double ndaus  = sv.numberOfDaughters();
      //double ntrk   = sv.numberOfSourceCandidatePtrs();
      bool cut8 = ndaus > ndauSV_;
      if(cut4 && cut5 && cut6 && cut7 && cut8)         {h1_["cutflow"] -> Fill(8, genwt);} //2 
      
      //significance of distance in 3d space point between the SV and PV
      Measurement1D dl= vdist.distance(pv, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
      bool cut9 = dl.significance() > dlenSigSV_;
      if(cut4 && cut5 && cut6 && cut7 && cut8 && cut9) {h1_["cutflow"] -> Fill(9, genwt);}//4
      
      //Fill the branches
      evt_.nGoodSV++;
      evt_.ptSV.push_back(sv.pt()); 
      evt_.etaSV.push_back(sv.eta());
      evt_.phiSV.push_back(sv.phi());
      evt_.energySV.push_back(sv.energy());
      evt_.massSV.push_back(sv.mass());
      evt_.isMatchedToJetSV.push_back(isSVInJetDir);
      evt_.dRFlightJetSV.push_back(dRSVFlightJet);
      evt_.dRJetSV.push_back(dRSVJet);
      evt_.chi2SV.push_back(sv.vertexChi2()); 
      evt_.normChi2SV.push_back(sv.vertexNormalizedChi2());
      evt_.dxySV.push_back(d2d.value());
      //evt_.dxyerrSV.push_back(catchInfsAndBound(d2d.error()-2,0,-2,0));
      evt_.dxyerrSV.push_back(d2d.error());
      evt_.dxysigSV.push_back(d2d.significance());  
      evt_.d3dSV.push_back(dl.value());
      evt_.d3derrSV.push_back(dl.error());
      evt_.d3dsigSV.push_back(dl.significance());
      evt_.costhetaSvPvSV.push_back(cosPdotV);
      evt_.ntrkSV.push_back(sv.numberOfSourceCandidatePtrs());
      evt_.ndauSV.push_back(ndaus);
      evt_.ndofSV.push_back(sv.vertexNdof());
      evt_.isMatchedtoGenSV.push_back(matchGenSVdR);
      evt_.isSVgen.push_back(p);      
      evt_.isSVb.push_back(b);//
      evt_.isSVc.push_back(c);//
      evt_.isSVl.push_back(l);//
      evt_.isSVo.push_back(o);//
      evt_.isSVbg.push_back(bg);//
      evt_.isSVcg.push_back(cg);
      evt_.nbsSV.push_back(nbs);
      evt_.ncsSV.push_back(ncs);
      evt_.nlsSV.push_back(nls);
      evt_.ngsSV.push_back(ngs);
      evt_.nothersSV.push_back(nothers);
      evt_.matchedGenIsHardProcess.push_back(matchGenIsHardProcess);
      evt_.matchedGenIsLastCopyBeforeFSR.push_back(matchGenIsLastCopyBeforeFSR);
      evt_.matchedGenPtSV.push_back(matchGenPt);
      evt_.matchedGenEtaSV.push_back(matchGenEta);
      evt_.matchedGenMassSV.push_back(matchGenMass);
      evt_.matchedGenPhiSV.push_back(matchGenPhi);
      evt_.matchedGenEnergySV.push_back(matchGenEnergy);
      evt_.matchedGenMomSV.push_back(matchGenMomID);
      evt_.matchedGenMindRSV.push_back(mindR);
      evt_.matchedGendRSV.push_back(matchGenSVdR);
   }
   //std::cout << "???" << std::endl; 
   tree_->Fill();
}

// ------------ method called once each run ----------------
void
KUNtuplizer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
   //edm::ESHandle<ME0Geometry> hGeom;
   //iSetup.get<MuonGeometryRecord>().get(hGeom);
   //ME0Geometry_ =( &*hGeom);   
}

// ------------ method called when ending the processing of a run  ------------
void
KUNtuplizer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called once each job just before starting event loop  ------------
void 
KUNtuplizer::beginJob()
{
   tree_ = fs_->make<TTree>("svtree", "svtree") ;
   evt_.RegisterTree(tree_, "Evt") ;
   genevt_.RegisterTree(tree_, "GenEvt") ;

   h1_["cutflow"] = fs_->make<TH1D>("cutflow SV", "cut flow", 9, 0.5, 9.5) ;
   h1_["cutflow"] -> GetXaxis() -> SetBinLabel(1, "allSV") ;
   h1_["cutflow"] -> GetXaxis() -> SetBinLabel(2, "SV p_{T} #leq 20") ;
   h1_["cutflow"] -> GetXaxis() -> SetBinLabel(3, "#DeltaR(SV,j)>0.4") ;
   h1_["cutflow"] -> GetXaxis() -> SetBinLabel(4, "#DeltaR(SV,pgen) #leq 0.05") ;
   h1_["cutflow"] -> GetXaxis() -> SetBinLabel(5, "#isHeavyFlavor") ;
   h1_["cutflow"] -> GetXaxis() -> SetBinLabel(6, "IP2D < 3") ;
   h1_["cutflow"] -> GetXaxis() -> SetBinLabel(7, "cos#theta_{PV,SV} > 0.98") ;
   h1_["cutflow"] -> GetXaxis() -> SetBinLabel(8, "ndau > 2") ; 
   h1_["cutflow"] -> GetXaxis() -> SetBinLabel(9, "SIP3D > 4") ;    
}

// ------------ method called once each job just after ending the event loop  ------------
void 
KUNtuplizer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
KUNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


const float& 
KUNtuplizer::catchInfs(const float& in,const float& replace_value) {
  if(in==in){
    if(std::isinf(in))
      return replace_value;
    else if(in < -1e32 || in > 1e32)
      return replace_value;
    return in;
  }
  return replace_value;
}

float 
KUNtuplizer::catchInfsAndBound(const float& in,const float& replace_value, const float& lowerbound, const float& upperbound) {
  float withoutinfs=catchInfs(in,replace_value);
  if(withoutinfs<lowerbound) return lowerbound;
  if(withoutinfs>upperbound) return upperbound;
  return withoutinfs;
}

//define this as a plug-in
DEFINE_FWK_MODULE(KUNtuplizer);
