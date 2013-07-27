// -*- C++ -*-
//
// Package:    ZJetsExpress
// Class:      ZJetsExpress
// 
/**\class ZJetsExpress ZJetsExpress.cc ElectroWeakAnalysis/VPlusJets/src/ZJetsExpress.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  A. Marini, K. Kousouris,  K. Theofilatos
//         Created:  Mon Oct 31 07:52:10 CDT 2011
// $Id: PATZJetsExpress.cc,v 1.1 2012/03/30 12:53:54 kkousour Exp $
//
//

// system include files
#include <memory>

#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>


#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1I.h"
#include "TH1F.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h" 
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h" 
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h" 
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRawData/interface/EcalRawDataCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

//
// class declaration
//
using namespace edm;
using namespace std;
using namespace reco;

class PATZJetsExpress : public edm::EDAnalyzer {
   public:
      explicit PATZJetsExpress(const edm::ParameterSet&);
      ~PATZJetsExpress();

   private:
      virtual void beginJob();
      virtual void beginRun(edm::Run const &, edm::EventSetup const& iSetup);
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob();
      virtual bool checkTriggerName(string,std::vector<string>); //checks if string belongs to any of the vector<string>
      // ---- method that builds the tree -------------------------------
      void buildTree();
      // ---- method that re-initializes the tree branches --------------
      void clearTree();
      // ---- structures for temporary object storage -------------------
      struct PARTICLE {
        // ---- momentum 4-vector ---------------------------------------
        TLorentzVector p4; 
        // ---- standard isolation --------------------------------------
        float iso;
        // ---- modified isolation --------------------------------------
        float isoPF; 
        // ---- modified isolation --------------------------------------
        float isoRho; 
        // ---- charge id (+/-1: electrons, +/-2: muons) ---------------- 
        int chid; 
        // ---- tight id ------------------------------------------------
        int id;
        // --- flag bit
        int bit;
        // --- float parameters
        std::vector<float> parameters;
      };
      struct GENPARTICLE {
        // ---- momentum 4-vector ---------------------------------------
        TLorentzVector p4;  
        // ---- pdgid ---------------------------------------------------
        int pdgId; 
      };
      struct JET {
        // ---- momentum 4-vector ---------------------------------------
        TLorentzVector p4; 
        // ---- tight id ------------------------------------------------
        int   id; 
        // ---- jet area (needed for pu subtraction) --------------------
        float area;
        // ---- track pt fraction associated with the PV ----------------
        float beta;
        // ---- btagger -------------------------------------------------
        float btag;
        // ---- jet energy correction factor ----------------------------
        float jec; 
        // ---- jet energy uncertainty ----------------------------------
        float unc;
        // ---- charged hadron energy fraction --------------------------
        float chf;
        // ---- neutral hadron energy fraction --------------------------
        float nhf;
        // ---- photon energy fraction ---------------------------------- 
        float phf;
        // ---- muon energy fraction ------------------------------------ 
        float muf;
        // ---- electron energy fraction --------------------------------
        float elf;
      };
      // ---- sorting rules ---------------------------------------------
      static bool lepSortingRule(PARTICLE x, PARTICLE y)                {return x.p4.Pt() > y.p4.Pt();}
      static bool lepSortingRuleGEN(GENPARTICLE x, GENPARTICLE y)       {return x.p4.Pt() > y.p4.Pt();}
      static bool phoSortingRuleGEN(GENPARTICLE x, GENPARTICLE y)       {return x.p4.Pt() > y.p4.Pt();}
      //static bool jetSortingRule(JET x, JET y)                      {return x.p4.Pt() > y.p4.Pt();}
      static bool p4SortingRule(TLorentzVector x, TLorentzVector y) {return x.Pt() > y.Pt();}
      // ---------- member data -----------------------------------------
      edm::Service<TFileService> fTFileService;	
      TTree *myTree_;
      // ---- histogram to record the number of events ------------------
      TH1I  *hRecoLeptons_,*hGenLeptons_,*hEvents_,*hWEvents_;
      TH1F  *hMuMuMass_,*hElElMass_,*hElMuMass_;
      TH1F  *hZMuMuMass_,*hZElElMass_;
      TH1F  *hMuMuMassWeighted_,*hElElMassWeighted_,*hElMuMassWeighted_;
      TH1F  *hElElEBMass_,*hLepLepMass_;
      TH1F  *hTriggerNames_,*hTriggerPass_;
      TH1F  *hGenPhotonPt_,*hGenPhotonMatchedPt_;
      TH1F  *hGenPhotonEta_,*hGenPhotonMatchedEta_;
      TH1F  *hGenMuonPt_,*hGenMuonMatchedPt_;
      TH1F  *hGenMuonEta_,*hGenMuonMatchedEta_;
      TH1F  *hGenElectronPt_,*hGenElectronMatchedPt_;
      TH1F  *hGenElectronEta_,*hGenElectronMatchedEta_;
      // ---- simulated in-time pile-up ---------------------------------
      TH1D  *mcPU_;
      // ---- flag to set the JEC uncertainty object --------------------
      //bool mIsJECuncSet;
      // ---- jet energy corrector object -------------------------------
      //const JetCorrector *mJEC;
      // ---- jet energy uncertainty object -----------------------------
      //JetCorrectionUncertainty *mJECunc;
      // ---- trigger ---------------------------------------------------
      std::string   processName_;
      std::vector<std::string> triggerNames_,triggerNamesFull_;
      std::vector<std::string> triggerFamily1_;
      std::vector<std::string> triggerFamily2_;
      std::vector<std::string> triggerFamily3_;
      std::vector<std::string> triggerFamily4_;
      std::vector<std::string> prescaleDontAsk_;
      std::vector<unsigned int> triggerIndex_;
      edm::InputTag triggerResultsTag_;
      edm::InputTag triggerEventTag_;
      edm::Handle<edm::TriggerResults>   triggerResultsHandle_;
      edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
      HLTConfigProvider hltConfig_;
      // ---- configurable parameters -----------------------------------
      bool          mIsMC;
      int           mMinNjets;
      double        mMinJetPt,mMaxJetEta,mMinLepPt,mMaxLepEta,mMaxCombRelIso03,mJetLepIsoR,mJetPhoIsoR,mMinLLMass,mMinPhoPt,mMaxPhoEta;
      //string        mJECserviceMC, mJECserviceDATA, mPayloadName;
      edm::InputTag mJetsName,mSrcRho,mSrcRho25;
      // ---- tree variables --------------------------------------------
      // ---- event number ----------------------------------------------
      ULong64_t eventNum_;
      // ---- run number ------------------------------------------------  
      int runNum_; 
      // ---- lumi section ----------------------------------------------
      int lumi_;
      // ---- flag to identify real data --------------------------------
      int isRealData_;
      // ---- trigger decisions -----------------------------------------
      std::vector<int> *fired_;
      int isTriggered_;
      // ---- L1 prescale -----------------------------------------------
      std::vector<int> *prescaleL1_;
      // ---- HLT prescale -----------------------------------------------
      std::vector<int> *prescaleHLT_;
      // ---- dilepton mass ---------------------------------------------
      float llM_,llMGEN_;
      // ---- dilepton rapidity -----------------------------------------  
      float llY_,llYGEN_;
      // ---- dilepton pseudorapidity -----------------------------------------  
      float llEta_,llEtaGEN_;
      // ---- dilepton pt -----------------------------------------------
      float llPt_,llPtGEN_;
      // ---- dilepton phi ----------------------------------------------
      float llPhi_,llPhiGEN_;
      // ---- dphi between the two leptons ------------------------------
      float llDPhi_,llDPhiGEN_; 
      // ---- dphi between jets and dilepton ----------------------------
      vector<float> *jetllDPhi_,*jetllDPhiGEN_;
      // ---- lepton kinematics -----------------------------------------
      vector<float> *lepPt_,*lepEta_,*lepPhi_,*lepE_,*lepPtGEN_,*lepEtaGEN_,*lepPhiGEN_,*lepEGEN_;
      // ---- lepton properties ----------------------------------------- 
      vector<int>   *lepChId_,*lepId_,*lepChIdGEN_,*lepMatchedGEN_;
      vector<float> *lepIso_,*lepIsoPF_,*lepIsoRho_,*lepMatchedDRGEN_;
      // ---- number of leptons -----------------------------------------
      int nLeptons_,nLeptonsGEN_;
      // ---- photon variables
      int nPhotons_;
      int isPhotonlead_;
      float photonE_;
      float photonPt_;
      float photonEta_;
      float photonPhi_;
      float photonIso_;
      float photonID_;
      int   photonBit_;
      float pfhadPhoPt_; // hadronic recoil computed for photon+jet interpretation
      float mPhotonj1_;
      float ptPhotonj1_;
      vector<float> *jetPhotonDPhi_;
      vector<float> *photonPar_;
      int nPhotonsGEN_;
      float photonEGEN_;
      float photonPtGEN_;
      float photonEtaGEN_;
      float photonPhiGEN_;
      float photonRECODRGEN_;
      // ---- FSR photon
      float FSRphotonE_;
      float FSRphotonPt_;
      float FSRphotonEta_;
      float FSRphotonPhi_;
      float FSRphotonIso_;
      float FSRphotonID_;
      float FSRphotonllM_;
      int   FSRphotonBit_;
      int   FSRphotonJet_;
      vector<float> *FSRphotonPar_;
      // ---- VBParton variables
      float VBPartonM_;
      float VBPartonE_;
      float VBPartonPt_;
      float VBPartonEta_;
      float VBPartonPhi_;
      int VBPartonDM_; // decay mode
      // ---- jet kinematics --------------------------------------------
      vector<float> *jetPt_,*jetEta_,*jetY_,*jetPhi_,*jetE_,*jetPtGEN_,*jetEtaGEN_,*jetYGEN_,*jetPhiGEN_,*jetEGEN_;
      // ---- jet composition fractions ---------------------------------
      vector<float> *jetCHF_,*jetPHF_,*jetNHF_,*jetMUF_,*jetELF_;
      // ---- other jet properties --------------------------------------
      vector<float> *jetBeta_,*jetBtag_,*jetArea_,*jetJEC_,*jetUNC_;
      // ---- tight jet id ----------------------------------------------
      vector<int>   *jetId_; 
      // ---- DR rejected jet kinematics --------------------------------------------
      vector<float> *rjetPt_,*rjetEta_,*rjetY_,*rjetPhi_,*rjetE_,*rjetPtGEN_,*rjetEtaGEN_,*rjetYGEN_,*rjetPhiGEN_,*rjetEGEN_;
      // ---- rjet composition fractions ---------------------------------
      vector<float> *rjetCHF_,*rjetPHF_,*rjetNHF_,*rjetMUF_,*rjetELF_;
      // ---- other rjet properties --------------------------------------
      vector<float> *rjetBeta_,*rjetBtag_,*rjetArea_,*rjetJEC_,*rjetUNC_;
      // ---- tight rjet id ----------------------------------------------
      vector<int>   *rjetId_; 
      // ---- number of jets --------------------------------------------
      int nJets_,nJetsGEN_;
      int nRJets_;
      // ---- flag to determine if the Z is one of the 2 leading objects-
      int isZlead_,isZleadGEN_;
      // ---- HT of the two leading objects -----------------------------
      float htLead_,htLeadGEN_;
      // ---- jet pt geometric mean -------------------------------------
      float jetPtGeMean_,jetPtGeMeanGEN_;
      // ---- jet pt arithmetic mean ------------------------------------
      float jetPtArMean_,jetPtArMeanGEN_;
      // ---- dphi between jets -----------------------------------------
      float j1j2DPhi_,j1j3DPhi_,j2j3DPhi_,j1j2DPhiGEN_,j1j3DPhiGEN_,j2j3DPhiGEN_;
      // ---- dR between jets -------------------------------------------
      float j1j2DR_,j1j3DR_,j2j3DR_,j1j2DRGEN_,j1j3DRGEN_,j2j3DRGEN_; 
      // ---- pf met ----------------------------------------------------
      float pfmet_;
      // ---- pf sumEt --------------------------------------------------
      float pfSumEt_;
      // ---- sum of pt of all jets --------------------------------------
      float HTJetSum_;
      // ---- sum of pt of all GEN jets --------------------------------------
      float HTJetSumGEN_;
      // ---- pf met phi ------------------------------------------------
      float pfmetPhi_;
      // ---- pt of the hadronic recoil ---------------------------------
      float pfhadPt_; 
      // ---- invariant mass of the Z and the leading jets -------------- 
      float mZj1_,mZj1GEN_;
      float mZj1j2_;
      float mZj1j2j3_;
      // ---- transverse momentum of Z and the leading jet -------------- 
      float ptZj1_,ptZj1GEN_;
      // ---- transverse momentum of Z and the leading jet -------------- 
      float costhetaZj1_,costhetaZj1GEN_;
      float costhetaPhotonj1_;
      // ---- invariant mass of jets (aka potatoes) ---------------------
      float mj1j2_,mj1j2GEN_;
      float mj1j2j3_;
      float mrj1rj2_;
      // ---- invariant mas of all leptons ------------------------------
      float mLep_,mLepGEN_;
      // ---- pf pt density ---------------------------------------------
      float rho_;
      float rho25_;
      // ---- reconstructed vertices' prperties -------------------------
      vector<float> *vtxZ_,*vtxNdof_;
      // ---- number of good reconstructed vertices ---------------------
      int   nVtx_;
      // ---- number of simulated pu interactions -----------------------
      int   puINT_,puOOT_;
      // ---- RECO, GEN accepted flags for MC ---------------------------
      int selRECO_,selGEN_;
      // ---- MC weight
      float mcWeight_;
};
//
// class implemetation
//
// ---- constructor -----------------------------------------------------
PATZJetsExpress::PATZJetsExpress(const ParameterSet& iConfig)
{
  mMinNjets          = iConfig.getParameter<int>                       ("minNjets");   
  mJetLepIsoR        = iConfig.getParameter<double>                    ("jetLepIsoRadius");
  mJetPhoIsoR        = iConfig.getParameter<double>                    ("jetLepIsoRadius");
  mMinJetPt          = iConfig.getParameter<double>                    ("minJetPt");
  mMaxJetEta         = iConfig.getParameter<double>                    ("maxJetEta");
  mMinLepPt          = iConfig.getParameter<double>                    ("minLepPt");
  mMaxLepEta         = iConfig.getParameter<double>                    ("maxLepEta");
  mMinPhoPt          = iConfig.getParameter<double>                    ("minPhoPt");
  mMaxPhoEta         = iConfig.getParameter<double>                    ("maxPhoEta");
  mMinLLMass         = iConfig.getParameter<double>                    ("minLLMass");
  mMaxCombRelIso03   = iConfig.getParameter<double>                    ("maxCombRelIso03");
  mJetsName          = iConfig.getParameter<edm::InputTag>             ("jets");
  mSrcRho            = iConfig.getParameter<edm::InputTag>             ("srcRho");
  mSrcRho25          = iConfig.getParameter<edm::InputTag>             ("srcRho25");
  //mJECserviceDATA    = iConfig.getParameter<std::string>               ("jecServiceDATA");
  //mJECserviceMC      = iConfig.getParameter<std::string>               ("jecServiceMC");
  //mPayloadName       = iConfig.getParameter<std::string>               ("payload");
  processName_       = iConfig.getParameter<std::string>               ("processName");
  triggerNames_      = iConfig.getParameter<std::vector<std::string> > ("triggerName");
  triggerFamily1_    = iConfig.getParameter<std::vector<std::string> > ("triggerFamily1");
  triggerFamily2_    = iConfig.getParameter<std::vector<std::string> > ("triggerFamily2");
  triggerFamily3_    = iConfig.getParameter<std::vector<std::string> > ("triggerFamily3");
  triggerFamily4_    = iConfig.getParameter<std::vector<std::string> > ("triggerFamily4");
  prescaleDontAsk_   = iConfig.getParameter<std::vector<std::string> > ("prescaleDontAsk");
  triggerResultsTag_ = iConfig.getParameter<edm::InputTag>             ("triggerResults");
  triggerEventTag_   = iConfig.getParameter<edm::InputTag>             ("triggerEvent");   
}
// ---- destructor ------------------------------------------------------
PATZJetsExpress::~PATZJetsExpress()
{
}
// ---
bool PATZJetsExpress::checkTriggerName(string aString,std::vector<string> aFamily)
{
  bool result(false);	
  for(unsigned int i=0;i<aFamily.size();i++) // checks if any of the aFamily strings contains aString
  {
    size_t found = aFamily[i].find(aString);
    if(found!=string::npos)result=true;
  }
  return result;
} 
// ---- method called once each job just before starting event loop -----
void PATZJetsExpress::beginJob()
{
  // ---- create the objects --------------------------------------------
  hRecoLeptons_          = fTFileService->make<TH1I>("RecoLeptons", "RecoLeptons",6,0,6);hRecoLeptons_->Sumw2();
  hGenLeptons_           = fTFileService->make<TH1I>("GenLeptons", "GenLeptons",6,0,6);hGenLeptons_->Sumw2();
  hEvents_               = fTFileService->make<TH1I>("Events", "Events",1,0,1);hEvents_->Sumw2();
  hMuMuMass_             = fTFileService->make<TH1F>("MuMuMass", "MuMuMass",40,71,111);hMuMuMass_->Sumw2();
  hElElMass_             = fTFileService->make<TH1F>("ElElMass", "ElElMass",40,71,111);hElElMass_->Sumw2();
  hElMuMass_             = fTFileService->make<TH1F>("ElMuMass", "ElMuMass",40,71,111);hElMuMass_->Sumw2();
  hElElEBMass_           = fTFileService->make<TH1F>("ElElEBMass", "ElElEBMass",40,71,111);hElElEBMass_->Sumw2();
  hLepLepMass_           = fTFileService->make<TH1F>("LepLepMass", "LepLepMass",980,40,2000);hLepLepMass_->Sumw2();
  hTriggerNames_         = fTFileService->make<TH1F>("TriggerNames","TriggerNames",1,0,1);
  hTriggerNames_ ->SetBit(TH1::kCanRebin);
  for(unsigned i=0;i<triggerNames_.size();i++)
    hTriggerNames_->Fill(triggerNames_[i].c_str(),1);
  hTriggerPass_          = fTFileService->make<TH1F>("TriggerPass","TriggerPass",1,0,1);
  hTriggerPass_  ->SetBit(TH1::kCanRebin);
  mcPU_                  = fTFileService->make<TH1D>("mcPU", "mcPU",40,0,40);

  hWEvents_              = fTFileService->make<TH1I>("WEvents", "Weighted Events",1,0,1);hWEvents_->Sumw2();

  hMuMuMassWeighted_     = fTFileService->make<TH1F>("MuMuMassWeighted", "MuMuMassWeighted",40,71,111);hMuMuMassWeighted_->Sumw2();
  hElElMassWeighted_     = fTFileService->make<TH1F>("ElElMassWeighted", "ElElMassWeighted",40,71,111);hElElMassWeighted_->Sumw2();
  hElMuMassWeighted_     = fTFileService->make<TH1F>("ElMuMassWeighted", "ElMuMassWeighted",40,71,111);hElMuMassWeighted_->Sumw2();

  hGenPhotonPt_          = fTFileService->make<TH1F>("hGenPhotonPt", "hGenPhotonPt;photon p_{T} [GeV];events",300,50,350);hGenPhotonPt_->Sumw2();
  hGenPhotonEta_         = fTFileService->make<TH1F>("hGenPhotonEta","hGenPhotonPt;photon p_{T} [GeV];events",300,-3.0,3.0);hGenPhotonEta_->Sumw2();
  hGenPhotonMatchedPt_   = fTFileService->make<TH1F>("hGenPhotonMatchedPt", "hGenPhotonMatchedPt;photon p_{T} [GeV];events",300,50,350);
  hGenPhotonMatchedEta_  = fTFileService->make<TH1F>("hGenPhotonMatchedEta","hGenPhotonMatchedPt;photon p_{T} [GeV];events",300,-3.0,3.0);
  hGenPhotonMatchedPt_->Sumw2();
  hGenPhotonMatchedEta_->Sumw2();

  hGenMuonPt_            = fTFileService->make<TH1F>("hGenMuonPt", "hGenMuonPt;muon p_{T} [GeV];events",300,0,300);hGenMuonPt_->Sumw2();
  hGenMuonEta_           = fTFileService->make<TH1F>("hGenMuonEta","hGenMuonPt;muon p_{T} [GeV];events",300,-3.0,3.0);hGenMuonEta_->Sumw2();
  hGenMuonMatchedPt_     = fTFileService->make<TH1F>("hGenMuonMatchedPt", "hGenMuonMatchedPt;muon p_{T} [GeV];events",300,0,300);
  hGenMuonMatchedEta_    = fTFileService->make<TH1F>("hGenMuonMatchedEta","hGenMuonMatchedPt;muon p_{T} [GeV];events",300,-3.0,3.0);
  hGenMuonMatchedPt_->Sumw2();
  hGenMuonMatchedEta_->Sumw2();

  hGenElectronPt_            = fTFileService->make<TH1F>("hGenElectronPt", "hGenElectronPt;muon p_{T} [GeV];events",300,0,300);hGenElectronPt_->Sumw2();
  hGenElectronEta_           = fTFileService->make<TH1F>("hGenElectronEta","hGenElectronPt;muon p_{T} [GeV];events",300,-3.0,3.0);hGenElectronEta_->Sumw2();
  hGenElectronMatchedPt_     = fTFileService->make<TH1F>("hGenElectronMatchedPt", "hGenElectronMatchedPt;muon p_{T} [GeV];events",300,0,300);
  hGenElectronMatchedEta_    = fTFileService->make<TH1F>("hGenElectronMatchedEta","hGenElectronMatchedPt;muon p_{T} [GeV];events",300,-3.0,3.0);
  hGenElectronMatchedPt_->Sumw2();
  hGenElectronMatchedEta_->Sumw2();

  hZMuMuMass_             = fTFileService->make<TH1F>("ZMuMuMass", "ZMuMuMass",40,71,111);hZMuMuMass_->Sumw2();
  hZElElMass_             = fTFileService->make<TH1F>("ZElElMass", "ZElElMass",40,71,111);hZElElMass_->Sumw2();


  

  myTree_                = fTFileService->make<TTree>("events", "events");
  // ---- build the tree ------------------------------------------------
  buildTree();
  // ---- set the jec uncertainty flag ----------------------------------
  //mIsJECuncSet = false; 
}
// ---- method called everytime there is a new run ----------------------
void PATZJetsExpress::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  if (triggerNames_.size() > 0) {
    bool changed(true);
    if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
      if (changed) {
        triggerNamesFull_.clear();
        // check if trigger names in (new) config
        cout<<"New trigger menu found !!!"<<endl;
        triggerIndex_.clear(); 
        const unsigned int n(hltConfig_.size());
        for(unsigned itrig=0;itrig<triggerNames_.size();itrig++) {
          bool found(false);
          for(unsigned iName=0;iName<n;iName++) {
            std::string ss(hltConfig_.triggerName(iName));
            if (int(ss.find(triggerNames_[itrig])) > -1) {
              triggerNamesFull_.push_back(ss);
              found = true;
              continue;
            } 
          }
          if (!found) {
            triggerNamesFull_.push_back("");
          }
          triggerIndex_.push_back(hltConfig_.triggerIndex(triggerNamesFull_[itrig]));
          cout<<triggerNames_[itrig]<<" "<<triggerNamesFull_[itrig]<<" "<<triggerIndex_[itrig]<<" ";  
          if (triggerIndex_[itrig] >= n)
            cout<<"does not exist in the current menu"<<endl;
          else
            cout<<"exists"<<endl;
        }// trigger names loop
      }
    } 
    else {
      cout << "ProcessedTreeProducer::analyze:"
           << " config extraction failure with process name "
           << processName_ << endl;
    }
  }
}
// ---- event loop ------------------------------------------------------
void PATZJetsExpress::analyze(const Event& iEvent, const EventSetup& iSetup)
{
  // ---- event counter -------------------------------------------------
  hEvents_->Fill(0.5);  
  // ---- initialize the tree branches ----------------------------------
  clearTree();
  isRealData_ = iEvent.isRealData() ? 1:0;


  // ----  MC truth block -----------------------------------------------
  vector<GENPARTICLE>      myGenLeptons;
  vector<TLorentzVector> myGenJets;  
  vector<GENPARTICLE> myGenPhotons;  
  TLorentzVector VBParton(0,0,0,0);
  if (!isRealData_) {
    // ---- PU ----------------------------------------------------------
    Handle<vector<PileupSummaryInfo> > pileupInfo;
    iEvent.getByLabel("addPileupInfo", pileupInfo);
    vector<PileupSummaryInfo>::const_iterator PUI;
    puOOT_ = 0;
    puINT_ = 0;
    for(PUI = pileupInfo->begin(); PUI != pileupInfo->end(); ++PUI) {
      if (PUI->getBunchCrossing() == 0)
        puINT_ += PUI->getPU_NumInteractions();     
      else
        puOOT_ += PUI->getPU_NumInteractions();
    }// PUI loop
    mcPU_->Fill(puINT_);
    // --- MC weight
    Handle<GenEventInfoProduct> geninfo;  
    iEvent.getByLabel("generator",geninfo);
    mcWeight_ = geninfo->weight();
    hWEvents_->Fill(0.5,mcWeight_);
    // --- Gen Jets
    Handle<GenJetCollection> genjets;
    iEvent.getByLabel("ak5GenJets",genjets);
    Handle<GenParticleCollection> gen;
    iEvent.getByLabel("genParticles", gen);
    GenParticleCollection::const_iterator i_gen;
    GenJetCollection::const_iterator i_genjet;
 

    int VBPartonDM=0;
    // ---- loop over the gen particles ---------------------------------
    for(i_gen = gen->begin(); i_gen != gen->end(); i_gen++) {
     
    // save MC Vector Boson partons momenta and their decay mode, when applicable
    if( (i_gen->pdgId() ==23 || i_gen->pdgId()==22) && i_gen->status()==3) {
      float Px = i_gen->p4().Px();
      float Py = i_gen->p4().Py();
      float Pz = i_gen->p4().Pz();
      float E  = i_gen->p4().E();
      VBParton.SetPxPyPzE(Px,Py,Pz,E);
      const GenParticle* gen_Dau;
      for(int kk = 0; kk < int(i_gen-> numberOfDaughters()); ++kk) {
	gen_Dau = static_cast<const GenParticle*>(i_gen->daughter(kk)); // find daughter
	if(gen_Dau->pdgId()==23) continue;
	VBPartonDM = abs(gen_Dau->pdgId());
      }
      if(VBParton.Pt()>0) {
	VBPartonDM_  = VBPartonDM;
	VBPartonE_   = VBParton.E();
	VBPartonPt_  = VBParton.Pt();
	VBPartonEta_ = VBParton.Eta();
	VBPartonPhi_ = VBParton.Phi();
       	VBPartonM_   = VBParton.M();

        if(VBPartonDM_==13)hZMuMuMass_->Fill(VBPartonM_,mcWeight_);
        if(VBPartonDM_==11)hZElElMass_->Fill(VBPartonM_,mcWeight_);
      }
    }
 
      // ---- consider only final state particles -----------------------
      if (i_gen->status() == 1) {   
        // ---- consider only electron and muon flavors -----------------
        if (abs(i_gen->pdgId()) == 11 || abs(i_gen->pdgId()) == 13) {
           // ---- apply geometric and kinematic acceptance -------------
          if ((i_gen->pt() > mMinLepPt) && (fabs(i_gen->eta())) < mMaxLepEta) {
            GENPARTICLE aGenLepton;
            TLorentzVector lepP4GEN(i_gen->p4().Px(),i_gen->p4().Py(),i_gen->p4().Pz(),i_gen->p4().E());
            aGenLepton.p4    = lepP4GEN; 
            aGenLepton.pdgId = i_gen->pdgId();
            myGenLeptons.push_back(aGenLepton);
          }
        }

        // ---- consider only electron and muon flavors -----------------
        if (abs(i_gen->pdgId()) == 22) {
           // ---- apply geometric and kinematic acceptance -------------
          if ((i_gen->pt() > mMinPhoPt) && (fabs(i_gen->eta())) < mMaxPhoEta) {
            GENPARTICLE aGenPhoton;
            TLorentzVector phoP4GEN(i_gen->p4().Px(),i_gen->p4().Py(),i_gen->p4().Pz(),i_gen->p4().E());
            aGenPhoton.pdgId = i_gen->pdgId();
            aGenPhoton.p4    = phoP4GEN;
            myGenPhotons.push_back(aGenPhoton);
          }
        }
      }
    }
    hGenLeptons_->Fill(int(myGenLeptons.size()));
    // ---- sort the genLeptons -----------------------------------------
    sort(myGenLeptons.begin(),myGenLeptons.end(),lepSortingRuleGEN);
    // ---- sort the genPhotons -----------------------------------------
    sort(myGenPhotons.begin(),myGenPhotons.end(),phoSortingRuleGEN);
    // ---- genjets -----------------------------------------------------
    for(i_genjet = genjets->begin(); i_genjet != genjets->end(); i_genjet++) {
      // ---- genlepton - genjet cross cleaning -------------------------
      bool isISO(true);
      for(unsigned l=0;l<myGenLeptons.size();l++) { 
        // ---- genjet vs 2 leading genlepton cleaning ------------------
        if (l >= 2) continue; 
        if (deltaR(i_genjet->eta(),i_genjet->phi(),myGenLeptons[l].p4.Eta(),myGenLeptons[l].p4.Phi()) < mJetLepIsoR) {
          isISO = false;
          continue;
        }
      }

      // ---- genjet vs leading leading photon cleaning ------------------
      if (myGenPhotons.size()>0) {
       if (deltaR(i_genjet->eta(),i_genjet->phi(),myGenPhotons[0].p4.Eta(),myGenPhotons[0].p4.Phi()) < mJetPhoIsoR) { 
          isISO = false;
          continue;
        }
      }

      if (!isISO) continue;
      // ---- preselection on genjets -----------------------------------
      if ((i_genjet->pt() < mMinJetPt) || (fabs(i_genjet->eta()) > mMaxJetEta)) continue;
      TLorentzVector aGenJet(i_genjet->p4().Px(),i_genjet->p4().Py(),i_genjet->p4().Pz(),i_genjet->p4().E());
      myGenJets.push_back(aGenJet);  
    }// genjet loop
  }// if MC
  //---- Rho ------------------------------------------------------------
  Handle<double> rho;
  iEvent.getByLabel(mSrcRho,rho);
  //---- Rho25 ------------------------------------------------------------
  Handle<double> rho25;
  iEvent.getByLabel(mSrcRho25,rho25);
  //---- reco vertices block --------------------------------------------
  edm::Handle<VertexCollection> vertices_;
  iEvent.getByLabel("offlinePrimaryVertices", vertices_);
  const reco::Vertex *primVtx = &(*(vertices_.product()))[0];
  for(VertexCollection::const_iterator i_vtx = vertices_->begin(); i_vtx != vertices_->end(); ++i_vtx) {  
    if (!i_vtx->isFake() && (fabs(i_vtx->z()) < 24) && (i_vtx->ndof() >= 4)) {
      vtxZ_   ->push_back(i_vtx->z());
      vtxNdof_->push_back(i_vtx->ndof());
    }
  }
  //---- leptons block --------------------------------------------------
  Handle<View<GsfElectron> > electrons_;
  iEvent.getByLabel("gsfElectrons", electrons_);
  Handle<View<Muon> > muons_;
  iEvent.getByLabel("muons",muons_);
  vector<PARTICLE> myLeptons;
  // ---- loop over muons -----------------------------------------------
  for(View<Muon>::const_iterator i_mu = muons_->begin(); i_mu != muons_->end(); ++i_mu) {
    //---- apply kinematic and geometric acceptance 
    if ((i_mu->pt() < mMinLepPt) || (fabs(i_mu->eta()) > mMaxLepEta))  continue;

    //---- apply VBTF-like id (GlobalMuonPromptTight)
    //---- https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
    if (!muon::isGoodMuon(*i_mu,muon::GlobalMuonPromptTight))          continue; // should be redundant wrt the cuts below
    if (!(i_mu->isGlobalMuon()))                                       continue;
    if (fabs(i_mu->globalTrack()->normalizedChi2()) >= 10)             continue;
    if (i_mu->globalTrack()->hitPattern().numberOfValidMuonHits()==0)  continue;
    if (i_mu->numberOfMatchedStations() <=1)                           continue; 
    if (fabs(i_mu->innerTrack()->dxy(primVtx->position())) >= 0.2)     continue;
    if (i_mu->track()->hitPattern().numberOfValidPixelHits() == 0)     continue; 
    if (i_mu->track()->hitPattern().numberOfValidTrackerHits() <= 10)  continue; 

    
    //---- subdetector isolation rho corrected, for efficiency studies look at
    //---- https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs
    //---- https://indico.cern.ch/getFile.py/access?contribId=2&resId=0&materialId=slides&confId=128982 
    float muonIso = (i_mu->isolationR03().sumPt + i_mu->isolationR03().emEt + i_mu->isolationR03().hadEt)/i_mu->pt();


    float muEta = i_mu->eta(); // essentially track direction at Vtx (recommended prescription)
    float Aecal=0.041; // initiallize with EE value
    float Ahcal=0.032; // initiallize with HE value
    if (fabs(muEta)<1.48) { 
      Aecal = 0.074;   // substitute EB value 
      Ahcal = 0.023;   // substitute EE value
    }
    float muonIsoRho = (i_mu->isolationR03().sumPt + max(0.,i_mu->isolationR03().emEt -Aecal*(*rho)) 
                       + max(0.,i_mu->isolationR03().hadEt-Ahcal*(*rho)))/i_mu->pt();

    if (muonIsoRho > mMaxCombRelIso03)                                 continue;
    if (i_mu->innerTrack()->hitPattern().numberOfValidHits() < 11)     continue;
    if (i_mu->innerTrack()->hitPattern().numberOfValidPixelHits() < 1) continue;
    if (i_mu->numberOfMatches() < 2)                                   continue;
    if (i_mu->innerTrack()->ptError()/i_mu->pt() > 0.1)                continue;
      
    PARTICLE aLepton;
    TLorentzVector lepP4(i_mu->p4().Px(),i_mu->p4().Py(),i_mu->p4().Pz(),i_mu->p4().E());
    aLepton.p4   = lepP4;
    aLepton.chid = 2*i_mu->charge();
    aLepton.id   = 1; // all muons are tight
    aLepton.iso  = muonIso;
    aLepton.isoPF = -999;
    aLepton.isoRho = muonIsoRho;
    myLeptons.push_back(aLepton);
  }// muon loop
  // ---- loop over electrons -------------------------------------------
  for(View<GsfElectron>::const_iterator i_el = electrons_->begin();i_el != electrons_->end(); ++i_el) {
    float elPt                           = i_el->p4().Pt();
    float elEta                          = i_el->p4().Eta();
    float trackIso                      = i_el->dr03TkSumPt();
    float ecalIso                        = i_el->dr03EcalRecHitSumEt();
    float hcalIso                        = i_el->dr03HcalTowerSumEt();
    float sigmaIetaIeta                  = i_el->sigmaIetaIeta();
    float hadronicOverEm                 = i_el->hadronicOverEm();
    float deltaPhiSuperClusterTrackAtVtx = i_el->deltaPhiSuperClusterTrackAtVtx();
    float deltaEtaSuperClusterTrackAtVtx = i_el->deltaEtaSuperClusterTrackAtVtx();
    bool  isTight(false);
    float combinedIso03 = (i_el->dr03TkSumPt()+max(0.,i_el->dr03EcalRecHitSumEt()-1.)+i_el->dr03HcalTowerSumEt())/i_el->p4().Pt();
    float combinedIso03Rho =0 ;

    float Aecal[2] = {0.101,0.046};
    float Ahcal[2] = {0.021,0.040};
    enum detID {barrel=0,endcap};

    if ((elPt < mMinLepPt) || (fabs(elEta) > mMaxLepEta)) continue;
    // ---- use WP90 as default preselection, store also WP80 subset (https://twiki.cern.ch/twiki/bin/view/CMS/VbtfEleID2011)
    if (i_el->isEB()) {
      combinedIso03Rho = (trackIso + max(0. ,ecalIso - Aecal[barrel]*(*rho)) + max(0.,hcalIso - Ahcal[barrel]*(*rho)) )/elPt; 
      if (combinedIso03Rho>mMaxCombRelIso03)              continue;
      if (sigmaIetaIeta > 0.01)                           continue;
      if (deltaPhiSuperClusterTrackAtVtx > 0.8)           continue;
      if (deltaEtaSuperClusterTrackAtVtx > 0.007)         continue;
      if (hadronicOverEm > 0.15)                          continue; 
      if (sigmaIetaIeta < 0.01) {  // WP80 subset
	if (deltaPhiSuperClusterTrackAtVtx < 0.06) 
	  if (deltaEtaSuperClusterTrackAtVtx < 0.004) 
	    if (hadronicOverEm < 0.04) 
              isTight = true;
      }
    }// if EE
    if (i_el->isEE()) {
      combinedIso03Rho = (trackIso + max(0. ,ecalIso - Aecal[endcap]*(*rho)) + max(0.,hcalIso - Ahcal[endcap]*(*rho)) )/elPt; 
      if (combinedIso03Rho>mMaxCombRelIso03)              continue;
      if (sigmaIetaIeta > 0.03)                           continue;
      if (deltaPhiSuperClusterTrackAtVtx > 0.7)           continue;
      if (deltaEtaSuperClusterTrackAtVtx > 0.009)          continue;
      if (hadronicOverEm > 0.15)                          continue;
      if (sigmaIetaIeta<0.03) {  // WP80 subset
	if (deltaPhiSuperClusterTrackAtVtx < 0.03)
	  if (deltaEtaSuperClusterTrackAtVtx < 0.007)
	    if (hadronicOverEm < 0.15) 
              isTight = true;
      } 
    }
    PARTICLE aLepton;
    TLorentzVector lepP4(i_el->p4().Px(),i_el->p4().Py(),i_el->p4().Pz(),i_el->p4().E());
    aLepton.p4   = lepP4;
    aLepton.chid = i_el->charge();
    aLepton.id   = 0;
    if (isTight) {
      aLepton.id = 1;
    }
    aLepton.iso   = combinedIso03;
    aLepton.isoPF = -999;
    aLepton.isoRho = combinedIso03Rho;
    myLeptons.push_back(aLepton);
  } // electrons loop
  hRecoLeptons_->Fill(int(myLeptons.size()));
  // ---- sort the leptons according to their pt ------------------------
  sort(myLeptons.begin(),myLeptons.end(),lepSortingRule); 
  // ---- PF isolation for leptons --------------------------------------
  Handle<View<PFCandidate> > pfCandidates;
  iEvent.getByLabel("particleFlow", pfCandidates);
  for(unsigned il=0;il<myLeptons.size();il++) {
    float sumPt(0.0);
    for(unsigned ipf=0;ipf<pfCandidates->size();ipf++) {
      float dR = deltaR(myLeptons[il].p4.Eta(),myLeptons[il].p4.Phi(),(*pfCandidates)[ipf].eta(),(*pfCandidates)[ipf].phi());
      if (dR < 0.3) {
        sumPt += (*pfCandidates)[ipf].pt();
      }
    }
    float isoPF = (sumPt-myLeptons[il].p4.Pt()-(*rho)*0.2827434)/myLeptons[il].p4.Pt();
    myLeptons[il].isoPF = isoPF;
  }
  //---- photons block --------------------------------------------------
  vector<PARTICLE> myPhotons;
  vector<PARTICLE> myFSRphotons;

  if(true) 
  {
    Handle<reco::PhotonCollection> photons_;
    iEvent.getByLabel("photons",photons_);

    std::vector<std::string> photonIDCollectionTags_;
    photonIDCollectionTags_.push_back("PhotonCutBasedIDLoose");
    photonIDCollectionTags_.push_back("PhotonCutBasedIDLooseEM");
    photonIDCollectionTags_.push_back("PhotonCutBasedIDTight");

    const int nPhoIDC = photonIDCollectionTags_.size();
    std::vector< const edm::ValueMap<Bool_t>* > phoIds;

    for(int j=0; j<nPhoIDC; j++) {
      edm::Handle<edm::ValueMap<Bool_t> > phoIDValueMap_;
      iEvent.getByLabel("PhotonIDProd",photonIDCollectionTags_[j], phoIDValueMap_);
      phoIds.push_back( phoIDValueMap_.product() );
    }

  //---- loop over the photon collection ---------------------------------
    int ipho = 0;
    for(reco::PhotonCollection::const_iterator it = photons_->begin();it != photons_->end(); it++) {

  //---- don't bother if it has pt less than 15 GeV ----------------------
      if(it->pt() < 15) continue;
      if(abs(it->eta()) > mMaxPhoEta) continue;
      if(it->isEBEEGap())continue;
      if(it->isEB() && it->hadronicOverEm()>0.15)continue; // on-line requirement 
      if(it->isEB() && it->sigmaIetaIeta()>0.024)continue; // on-line requirement 
      if(it->isEE() && it->hadronicOverEm()>0.10)continue; // on-line requirement 
      if(it->isEE() && it->sigmaIetaIeta()>0.040)continue; // on-line requirement 

      reco::PhotonRef phoRef(photons_,ipho++);
      int photonID=0;

      TLorentzVector aPhoton(it->p4().Px(),it->p4().Py(),it->p4().Pz(),it->p4().E());
//   photonBit |= (it->isEB()          << 0);

      std::map<TString,UChar_t> idPairs;
      for(int k=0; k<nPhoIDC; k++) {
        idPairs[ TString(photonIDCollectionTags_[k].c_str()) ] = (*phoIds[k])[phoRef];
        if(photonIDCollectionTags_[k] == "PhotonCutBasedIDLoose")photonID |= (*phoIds[k])[phoRef] <<1;
        if(photonIDCollectionTags_[k] == "PhotonCutBasedIDTight")photonID |= (*phoIds[k])[phoRef] <<2;
      }// for id

      float hcalTowerSumEtConeDR03            = it->hcalTowerSumEtConeDR03(); // hcalTowerSumEtConeDR03
      float ecalRecHitSumEtConeDR03           = it->ecalRecHitSumEtConeDR03(); // ecalRecHitSumEtConeDR03
      float trkSumPtHollowConeDR03            = it->trkSumPtHollowConeDR03();


      float hcalTowerSumEtConeDR04            = it->hcalTowerSumEtConeDR04(); // hcalTowerSumEtConeDR04
      float ecalRecHitSumEtConeDR04           = it->ecalRecHitSumEtConeDR04(); // ecalRecHitSumEtConeDR04
      float trkSumPtHollowConeDR04            = it->trkSumPtHollowConeDR04();
      float nTrkSolidConeDR04                 = it->nTrkSolidConeDR04();
      float nTrkHollowConeDR04                = it->nTrkHollowConeDR04();
      float trkSumPtSolidConeDR04             = it->trkSumPtSolidConeDR04();


      float sigmaIetaIeta                     = it->sigmaIetaIeta();
      float phoHasConvTrks                    = it->hasConversionTracks();
      float r9                                = it->r9();
      float hadronicOverEm                    = it->hadronicOverEm();
      float sigmaIphiIphi                     = -1; // to be computed later

      float gammaPt = aPhoton.Pt();
      bool  isTriggerISO = false;

      // --- get iphi-iphi
      EcalClusterLazyTools lazyTool(iEvent, iSetup, InputTag("reducedEcalRecHitsEB"), InputTag("reducedEcalRecHitsEE"));

       // Next few lines get sigma-iphi-iphi
       const reco::CaloClusterPtr  seed = it->superCluster()->seed();
       if (seed.isAvailable()) {
         // Get sigma-iphi-iphi
         std::vector<float> vCov = lazyTool.covariances(*seed);
         sigmaIphiIphi  = sqrt(vCov[2]);   // This is Sqrt(covPhiPhi)
       }


      // --- https://twiki.cern.ch/twiki/bin/viewauth/CMS/QCDPhotonsTrigger2011
      if(ecalRecHitSumEtConeDR03 < 6.0 + 0.012*gammaPt) // requirement of _IsoVL_ type photon triggers 
      if(hcalTowerSumEtConeDR03  < 4.0 + 0.005*gammaPt)
      if(trkSumPtHollowConeDR03  < 4.0 + 0.002*gammaPt)isTriggerISO=true;

      // --- https://twiki.cern.ch/twiki/bin/viewauth/CMS/Vgamma2011PhotonID
      // --- recommended photon isolation + id in one step
      bool isVgamma2011 = false;
      float Rho25 = *rho25;
      float ATr = 0.0167 ;if(it->isEE()) ATr = 0.032;
      float AEc = 0.183  ;if(it->isEE()) AEc = 0.090;
      float AHc = 0.062  ;if(it->isEE()) AHc = 0.180;
      float sigmaIetaIetaMax = 0.011 ;if(it->isEE()) sigmaIetaIetaMax = 0.03;
      if(trkSumPtHollowConeDR04  < 2.0 + 0.001*gammaPt + ATr*Rho25)
      if(ecalRecHitSumEtConeDR04 < 4.2 + 0.006*gammaPt + AEc*Rho25) 
      if(hcalTowerSumEtConeDR04  < 2.2 + 0.0025*gammaPt + AHc*Rho25)
      if(sigmaIetaIeta<sigmaIetaIetaMax)
      if(!it->hasPixelSeed())
      if(it->isEE() || (it->isEB() && sigmaIetaIeta > 0.001 && sigmaIphiIphi > 0.001)) // additional EB spike cleaning
      if(hadronicOverEm<0.05) isVgamma2011 = true;
      photonID |= isVgamma2011 << 3;

      // --- online isolation + Vgamma2011 id
      if(isTriggerISO) 
      if(sigmaIetaIeta<sigmaIetaIetaMax)
      if(!it->hasPixelSeed())
      if(it->isEE() || (it->isEB() && sigmaIetaIeta > 0.001 && sigmaIphiIphi > 0.001)) // additional EB spike cleaning
      if(hadronicOverEm<0.05) 
      photonID |= 1 << 4;

      // --- Vgamma2011 photon id w/o isolation
      if(sigmaIetaIeta<sigmaIetaIetaMax) 
      if(!it->hasPixelSeed())
      if(it->isEE() || (it->isEB() && sigmaIetaIeta > 0.001 && sigmaIphiIphi > 0.001)) // additional EB spike cleaning
      if(hadronicOverEm<0.05) 
      photonID |= 1 << 5;
      

      // photon near masked region
      float gammaEta = aPhoton.Eta();
      float gammaPhi = aPhoton.Phi();
      bool mask_0  = ( gammaPhi>=-2.72 && gammaPhi<=-2.61 && gammaEta>=-1.33 && gammaEta<=-1.25 );
      bool mask_1  = ( gammaPhi<=-3.05 && gammaEta<=-1.40 );
      bool mask_2  = ( gammaPhi<=-3.05 && gammaEta>=1.04 && gammaEta<=1.15 );
      bool mask_3  = ( gammaPhi>=-2.64 && gammaPhi<=-2.52 && gammaEta>=-0.20 && gammaEta<=-0.08 );
      bool mask_4  = ( gammaPhi>=-2.64 && gammaPhi<=-2.52 && gammaEta>=1.38 && gammaEta<=1.46 );
      bool mask_5  = ( gammaPhi>=-2.03 && gammaPhi<=-1.91 && gammaEta>=-0.47 && gammaEta<=-0.3 );
      bool mask_6  = ( gammaPhi>=-1.25 && gammaPhi<=-1.12 && gammaEta>=-1.25 && gammaEta<=-1.13 );
      bool mask_7  = ( gammaPhi>=-0.81 && gammaPhi<=-0.69 && gammaEta>=-0.82 && gammaEta<=-0.69 );
      bool mask_8  = ( gammaPhi>=-0.55 && gammaPhi<=-0.42 && gammaEta>=1.21 && gammaEta<=1.33 );
      bool mask_9  = ( gammaPhi>=-0.29 && gammaPhi<=-0.16  && gammaEta>=1.38 );
      bool mask_10 = ( gammaPhi>=0.07 && gammaPhi<=0.19 && gammaEta>=-0.29 && gammaEta<=-0.16 );
      bool mask_11 = ( gammaPhi>=0.15 && gammaPhi<=0.28 && gammaEta>=-0.37 && gammaEta<=-0.25 );
      bool mask_12 = ( gammaPhi>=0.69 && gammaPhi<=0.79 && gammaEta>=-0.20 && gammaEta<=-0.07 );
      bool mask_13 = ( gammaPhi>=0.86 && gammaPhi<=0.97 && gammaEta>=-0.11 && gammaEta<=0.02 );
      bool mask_14 = ( gammaPhi>=0.60 && gammaPhi<=0.70 && gammaEta>=0.96 && gammaEta<=1.06 );
      bool mask_15 = ( gammaPhi>=1.74 && gammaPhi<=1.84 && gammaEta>=0.08 && gammaEta<=0.19 );
      bool mask_16 = ( gammaPhi>=1.65 && gammaPhi<=1.75 && gammaEta>=0.87 && gammaEta<=0.97 );
      bool mask_17 = ( gammaPhi>=1.99 && gammaPhi<=2.10 && gammaEta>=-0.97 && gammaEta<=-0.87 );
      bool mask_18 = ( gammaPhi>=2.95 && gammaPhi<=3.05 && gammaEta>=-0.98 && gammaEta<=-0.87 );
      bool mask_19 = ( gammaPhi>=2.78 && gammaPhi<=2.89 && gammaEta>=0.86 && gammaEta<=0.98);
      bool mask_20 = ( gammaPhi>=2.69 && gammaPhi<=2.81 && gammaEta>= 1.39);
 
      bool isMasked = mask_0 || mask_1 || mask_2 || mask_3 || mask_4 || mask_5 || mask_6 || mask_7 || mask_8 || mask_9 || mask_10 || mask_11;
      isMasked      = isMasked || mask_12 || mask_13 || mask_14 || mask_15 || mask_16 || mask_17 || mask_18 || mask_19 || mask_20;      




      int photonBit=0;
      photonBit |= (it->isEB()          << 0);
      photonBit |= (it->isEE()          << 1);
      photonBit |= (it->isEBEtaGap()    << 2);
      photonBit |= (it->isEBPhiGap()    << 3);
      photonBit |= (it->isEERingGap()   << 4);
      photonBit |= (it->isEEDeeGap()    << 5);
      photonBit |= (it->isEBEEGap()     << 6);
      photonBit |= (it->hasPixelSeed()  << 7);
      photonBit |= (isTriggerISO        << 8);
      photonBit |= (isMasked            << 9);




      PARTICLE gamma; // define pseudo lepton half of the p4 of the real photon
      gamma.p4    = aPhoton;
      gamma.chid  = 0;
      gamma.id    = photonID;
      gamma.iso   = isTriggerISO;    // this bool 0/1
      gamma.isoPF = 0; 
      gamma.bit   = photonBit;
      // ok, now close your eyes what follows is a scandal
      gamma.parameters.push_back(hcalTowerSumEtConeDR04);     // 0   need to remember the ordering offline
      gamma.parameters.push_back(ecalRecHitSumEtConeDR04);    // 1  
      gamma.parameters.push_back(nTrkSolidConeDR04);          // 2  
      gamma.parameters.push_back(trkSumPtSolidConeDR04);      // 3 
      gamma.parameters.push_back(nTrkHollowConeDR04);         // 4 
      gamma.parameters.push_back(trkSumPtHollowConeDR04);     // 5
      gamma.parameters.push_back(sigmaIetaIeta);              // 6  
      gamma.parameters.push_back(phoHasConvTrks);             // 7  
      gamma.parameters.push_back(r9);                         // 8  
      gamma.parameters.push_back(hadronicOverEm);             // 9  

      // --- save FSR photon candidates and gamma+jets in seperate paths
      
      myFSRphotons.push_back(gamma);                          // FSR photons are *not* used in the photon+jet DR cone rejection

      //---- consider single event interpretation, exclusive di-lepton/photon interpretation (myLeptons.size()==0) 
      if(it->pt() > mMinPhoPt && myLeptons.size()==0 && isTriggerISO) myPhotons.push_back(gamma);    //note: hard photons imply later a DR cone rejection wrt the jets
    }
  }
  sort(myPhotons.begin(),myPhotons.end(),lepSortingRule);
  sort(myFSRphotons.begin(),myFSRphotons.end(),lepSortingRule);
  //---- jets block -----------------------------------------------------
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(mJetsName,jets);
  edm::View<pat::Jet> pat_jets = *jets;

  //Handle<PFJetCollection> jets_;
  //iEvent.getByLabel(mJetsName,jets_);
  vector<JET> myJets;
  vector<JET> myRJets;
  //---- get the jet energy corrector -----------------------------------
  //if(isRealData_)mJEC = JetCorrector::getJetCorrector(mJECserviceDATA,iSetup);
  //if(!isRealData_)mJEC = JetCorrector::getJetCorrector(mJECserviceMC,iSetup);
  //---- the JEC uncertainty only needs to be set-up once ---------------
  //if (!mIsJECuncSet) {
    //edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
    //iSetup.get<JetCorrectionsRecord>().get(mPayloadName,JetCorParColl); 
    //JetCorrectorParameters const& JetCorPar = (*JetCorParColl)["Uncertainty"];
    //mJECunc = new JetCorrectionUncertainty(JetCorPar);
    //mIsJECuncSet = true;
  //}
  // ---- jets loop -----------------------------------------------------
  for(edm::View<pat::Jet>::const_iterator i_jet = pat_jets.begin();i_jet != pat_jets.end(); ++i_jet) {
    TLorentzVector jetP4(i_jet->px(),i_jet->py(),i_jet->pz(),i_jet->energy());
    bool jetIsDuplicate(false);
    bool jetIsInAcceptance(true);
    bool jetIsIDed(true);

    //----- remove the leptons ------------------------------------------
    for(unsigned int i_lep = 0; i_lep < myLeptons.size(); i_lep++) {
      float DR = myLeptons[i_lep].p4.DeltaR(jetP4);
      if (DR < mJetLepIsoR) {
        jetIsDuplicate = true; 
      }
    }// lepton loop 

    //----- remove the leading photon ------------------------------------ (reminder nPhotons>0 only IF nLeptons==0)
    for(unsigned int i_pho = 0; i_pho < myPhotons.size(); i_pho++) {
      float DR = myPhotons[i_pho].p4.DeltaR(jetP4);
      if (DR < mJetPhoIsoR) {
        jetIsDuplicate = true;
      }
    }// photon loop

    // ---- get the jec and the uncertainty -----------------------------    
    //int index = i_jet - jets_->begin();
    //edm::RefToBase<reco::Jet> jetRef(edm::Ref<PFJetCollection>(jets_,index));
    //double jec = mJEC->correction(*i_jet,jetRef,iEvent,iSetup);
    // ---- only keep jets within the kinematic acceptance --------------
    if ((i_jet->pt() < 15) || (fabs(i_jet->eta()) > 3.0)) continue; // apply first a hardcode preselection
    if ((i_jet->pt() < mMinJetPt) || (fabs(i_jet->eta()) > mMaxJetEta)) jetIsInAcceptance = false;
    //mJECunc->setJetEta(i_jet->eta());
    // ---- the unc is a function of the corrected pt -------------------
    //mJECunc->setJetPt(jec * i_jet->pt());
    //double unc = mJECunc->getUncertainty(true);
    float jec  = 1./i_jet->jecFactor(0);
    float unc  = i_jet->userFloat("jecUnc");
    float btag = i_jet->bDiscriminator("combinedSecondaryVertexBJetTags");
    float beta = i_jet->userFloat("beta");
    // ---- keep only jets that pass the tight id -----------------------
    float chf = i_jet->chargedHadronEnergyFraction();
    float nhf = i_jet->neutralHadronEnergyFraction() + i_jet->HFHadronEnergyFraction();
    float phf = i_jet->photonEnergy()/(i_jet->jecFactor(0) * i_jet->energy());
    float elf = i_jet->electronEnergy()/(i_jet->jecFactor(0) * i_jet->energy());
    float muf = i_jet->muonEnergy()/(i_jet->jecFactor(0) * i_jet->energy());
    int chm   = i_jet->chargedHadronMultiplicity();
    int npr   = i_jet->chargedMultiplicity() + i_jet->neutralMultiplicity();
    bool id = (npr>1 && phf<0.99 && nhf<0.99 && ((fabs(i_jet->eta())<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && chf>0 && chm>0) || fabs(i_jet->eta())>2.4));
    if (!id) jetIsIDed = false;
    // ---- jet vertex association --------------------------------------
    // ---- get the vector of tracks ------------------------------------ 
    //reco::TrackRefVector vTrks(i_jet->getTrackRefs());
    //float sumTrkPt(0.0),sumTrkPtBeta(0.0),sumTrkPtBetaStar(0.0),beta(-1.0),betaStar(-1.0);
    // ---- loop over the tracks of the jet -----------------------------
    /*
    for(reco::TrackRefVector::const_iterator i_trk = vTrks.begin(); i_trk != vTrks.end(); i_trk++) {
      sumTrkPt += (*i_trk)->pt();
      // ---- loop over all vertices ------------------------------------
      for(unsigned i_vtx = 0;i_vtx < vertices_->size();i_vtx++) {
        // ---- loop over the tracks associated with the vertex ---------
        if ((*vertices_)[i_vtx].isFake() || (*vertices_)[i_vtx].ndof() < 4) continue; 
        for(reco::Vertex::trackRef_iterator i_vtxTrk = (*vertices_)[i_vtx].tracks_begin(); i_vtxTrk != (*vertices_)[i_vtx].tracks_end(); ++i_vtxTrk) {
          // ---- match the jet track to the track from the vertex ------
          reco::TrackRef trkRef(i_vtxTrk->castTo<reco::TrackRef>());
          // ---- check if the tracks match -----------------------------
          if (trkRef == (*i_trk)) {
            if (i_vtx == 0) {
              sumTrkPtBeta += (*i_trk)->pt();
            }
            else {
              sumTrkPtBetaStar += (*i_trk)->pt();
            }   
            continue;
          }
        }
      }// vertices loop 
    }// jet tracks loop
    if (sumTrkPt > 0) {
      beta     = sumTrkPtBeta/sumTrkPt;
      betaStar = sumTrkPtBetaStar/sumTrkPt;
    }
    */
    JET aJet; 
    aJet.p4       = jetP4;
    aJet.jec      = jec;
    aJet.unc      = unc;
    aJet.area     = i_jet->jetArea();
    aJet.chf      = chf;
    aJet.nhf      = nhf;
    aJet.phf      = phf;
    aJet.elf      = elf;
    aJet.muf      = muf;
    aJet.id       = 0;
    if (id) {
      aJet.id     = 1;
    }
    aJet.beta     = beta;
    aJet.btag     = btag;
    if(!jetIsDuplicate && jetIsInAcceptance && jetIsIDed)myJets.push_back(aJet);
    if(jetIsDuplicate){aJet.p4       = jetP4; myRJets.push_back(aJet);}  // store the uncorrected jet (this is virtually the matched lepton in DR) 
  }// jet loop

  // ---- sort jets according to their corrected pt ---------------------
  //sort(myJets.begin(),myJets.end(),jetSortingRule);    
  //sort(myRJets.begin(),myRJets.end(),jetSortingRule);    
  // ---- MET block -----------------------------------------------------
  Handle<View<PFMET> > pfmetCol_;
  iEvent.getByLabel("pfMet", pfmetCol_);
  float pfmetPx = (pfmetCol_->front()).px();
  float pfmetPy = (pfmetCol_->front()).py();
  // ---- define di-lepton pair (if exist)
  TLorentzVector llP4(0,0,0,0);
  if(myLeptons.size()>1) llP4 = myLeptons[0].p4 + myLeptons[1].p4;
  TLorentzVector llP4GEN(0,0,0,0); 
  if(myGenLeptons.size()>1) llP4GEN = myGenLeptons[0].p4 + myGenLeptons[1].p4;  
  // ---- counters ------------------------------------------------------
  nVtx_        = int(vtxZ_->size());
  nLeptons_    = int(myLeptons.size()); 
  nPhotons_    = int(myPhotons.size());
  nPhotonsGEN_ = int(myGenPhotons.size());
  nJets_       = int(myJets.size());
  nRJets_      = int(myRJets.size());
  nLeptonsGEN_ = int(myGenLeptons.size()); 
  nJetsGEN_    = int(myGenJets.size()); 
  // ---- Gen To Reco Matching for leptons ------------------------------------------------------
  for(unsigned gg=0; gg < myGenLeptons.size(); gg++) {                               
    lepMatchedDRGEN_->push_back(999);
    lepMatchedGEN_->push_back(999);
    for(unsigned rr = 0; rr < myLeptons.size(); rr++) {
	float DR = myGenLeptons[gg].p4.DeltaR(myLeptons[rr].p4);
	bool isSameFlavor = (abs(myGenLeptons[gg].pdgId) == 11 && abs(myLeptons[rr].chid)==1);
	isSameFlavor = isSameFlavor || (abs(myGenLeptons[gg].pdgId) == 13 && abs(myLeptons[rr].chid)==2);
	if(DR<lepMatchedDRGEN_->at(gg) && isSameFlavor) {
	    lepMatchedGEN_->at(gg) = rr;                              //store recolepton matched index
	    lepMatchedDRGEN_->at(gg) = DR;                           //store DR with matched GEN lepton
	}
      }
  }
  // ---- plot some inclusive di-lepton spectra (prior jet requirement)
  if(nLeptons_>1) {
    int dileptonId = myLeptons[0].chid*myLeptons[1].chid;
    hLepLepMass_->Fill(llP4.M());
    if(dileptonId==-4)hMuMuMass_->Fill(llP4.M());
    if(dileptonId==-1)hElElMass_->Fill(llP4.M());
    if(dileptonId==-2)hElMuMass_->Fill(llP4.M());
    if(dileptonId==-1 && abs(myLeptons[0].p4.Eta())<1.4 && abs(myLeptons[1].p4.Eta())<1.4)hElElEBMass_->Fill(llP4.M());
    if(dileptonId==-4)hMuMuMassWeighted_->Fill(llP4.M(),mcWeight_);
    if(dileptonId==-1)hElElMassWeighted_->Fill(llP4.M(),mcWeight_);
    if(dileptonId==-2)hElMuMassWeighted_->Fill(llP4.M(),mcWeight_);
  }
  if(nPhotonsGEN_>0) {
   float pt = myGenPhotons[0].p4.Pt();
   float eta = myGenPhotons[0].p4.Eta();
   float DR=999;
   if(myPhotons.size()>0)DR = myPhotons[0].p4.DeltaR(myGenPhotons[0].p4);
     hGenPhotonPt_->Fill(pt);
     hGenPhotonEta_->Fill(eta);
     if(photonRECODRGEN_<0.2) {
       hGenPhotonMatchedPt_->Fill(pt);
       hGenPhotonMatchedEta_->Fill(eta);
     }  
  }

  for(int rr=0; rr< int(lepMatchedDRGEN_->size()); rr++) {
    if(rr>=2)continue; // do this only for first 2 leptons
    float pt = myGenLeptons[rr].p4.Pt();
    float eta = myGenLeptons[rr].p4.Eta();
    int pdgId =  myGenLeptons[rr].pdgId;
    float DR = lepMatchedDRGEN_->at(rr);
    if(abs(pdgId)==13) {
      hGenMuonPt_->Fill(pt);
      hGenMuonEta_->Fill(eta);
    }
    if(abs(pdgId)==13 && DR<0.2) {
      hGenMuonMatchedPt_->Fill(pt);
      hGenMuonMatchedEta_->Fill(eta);
    }
    if(abs(pdgId)==11) {
      hGenElectronPt_->Fill(pt);
      hGenElectronEta_->Fill(eta);
    }
    if(abs(pdgId)==11 && DR<0.2) {
      hGenElectronMatchedPt_->Fill(pt);
      hGenElectronMatchedEta_->Fill(eta);
    }
  }
  // ---- keep only selected events -------------------------------------
  bool selectionRECO = ((nVtx_ > 0) && (nLeptons_ > 1) && (nJets_ >= mMinNjets) && llP4.M()>mMinLLMass);
  selectionRECO = selectionRECO || ((nVtx_ > 0) &&  nPhotons_>0 && (nJets_ >= mMinNjets)); // add photon logic for RECO
  bool selection(selectionRECO);
  bool selectionGEN(false);
  selRECO_ = 0;
  if (!isRealData_) {
    selectionGEN = ((nLeptonsGEN_ > 1) && (nJetsGEN_ >= mMinNjets) && llP4GEN.M()>mMinLLMass); 
    selectionGEN = selectionGEN || ((nPhotonsGEN_ > 0) && (nJetsGEN_ >= mMinNjets));      // add photon logic for GEN 
    selection +=  selectionGEN;
    selGEN_ = 0;
  }
  if (selection) {
    eventNum_   = iEvent.id().event();
    runNum_     = iEvent.id().run();
    lumi_       = iEvent.luminosityBlock();
    isRealData_ = isRealData_; // just pass through 
    rho_        = *rho; 
    rho25_      = *rho25; 
    pfmet_      = (pfmetCol_->front()).pt();
    pfmetPhi_   = (pfmetCol_->front()).phi();
    pfSumEt_    = (pfmetCol_->front()).sumEt();
    if (selectionRECO) {
      selRECO_ = 1;
      // ---- hadronic recoil vector --------------------------------------
      TLorentzVector pfmetP4(pfmetPx,pfmetPy,0,sqrt(pfmetPx * pfmetPx + pfmetPy * pfmetPy));
      TLorentzVector pfhadP4(0,0,0,0);
      if(myLeptons.size()>1) { 
	llM_                   = llP4.M();
	llPt_                  = llP4.Pt();
	llPhi_                 = llP4.Phi();
	if(llPt_>0)llY_        = llP4.Rapidity();
	if(llPt_>0)llEta_      = llP4.Eta();
	llDPhi_     = fabs(myLeptons[0].p4.DeltaPhi(myLeptons[1].p4));
        pfhadP4     = -pfmetP4 - llP4;      
      }
      TLorentzVector lepP4(0,0,0,0); 
      for(unsigned l = 0; l < myLeptons.size(); l++) {
        lepP4 += myLeptons[l].p4;
        lepPt_     ->push_back(myLeptons[l].p4.Pt());
        lepEta_    ->push_back(myLeptons[l].p4.Eta());
        lepPhi_    ->push_back(myLeptons[l].p4.Phi());
        lepE_      ->push_back(myLeptons[l].p4.Energy());
        lepIso_    ->push_back(myLeptons[l].iso);
        lepIsoPF_  ->push_back(myLeptons[l].isoPF);
        lepIsoRho_ ->push_back(myLeptons[l].isoRho);
        lepId_     ->push_back(myLeptons[l].id);
        lepChId_   ->push_back(myLeptons[l].chid);
      }      
      pfhadPt_    = pfhadP4.Pt();
      mLep_       = lepP4.M(); 
      // ---- store photon variables ------------------------------------
      nPhotons_                 = myPhotons.size();
      if(nPhotons_>0) {
        photonE_   = myPhotons[0].p4.Energy();
        photonPt_  = myPhotons[0].p4.Pt();
        photonEta_ = myPhotons[0].p4.Eta();
        photonPhi_ = myPhotons[0].p4.Phi();
        photonIso_ = myPhotons[0].iso;
        photonID_  = myPhotons[0].id;
        photonBit_ = myPhotons[0].bit;
        *photonPar_= myPhotons[0].parameters;
        pfhadPhoPt_=  (-pfmetP4 - myPhotons[0].p4).Pt();
      }

      if(myFSRphotons.size()>0 && myLeptons.size()>1) { // save FSR photons for di-lepton only events

	bool writeFSR = true;
        if(myLeptons[0].p4.DeltaR(myFSRphotons[0].p4)<0.2 && abs(myLeptons[0].chid)==1)writeFSR=false; //don't consider FSR when very close with electrons
        if(myLeptons[1].p4.DeltaR(myFSRphotons[0].p4)<0.2 && abs(myLeptons[1].chid)==1)writeFSR=false;
 
	if(writeFSR) {
          FSRphotonE_   = myFSRphotons[0].p4.Energy();
          FSRphotonPt_  = myFSRphotons[0].p4.Pt();
          FSRphotonEta_ = myFSRphotons[0].p4.Eta();
          FSRphotonPhi_ = myFSRphotons[0].p4.Phi();
          FSRphotonIso_ = myFSRphotons[0].iso;
          FSRphotonID_  = myFSRphotons[0].id;
          FSRphotonBit_ = myFSRphotons[0].bit;
          *FSRphotonPar_= myFSRphotons[0].parameters;
          FSRphotonllM_  = (llP4 + myFSRphotons[0].p4).M();
  
	  bool isFSRphotonCountedAsJet = false;
          for(unsigned j = 0; j < myJets.size(); j++) {
  	  float DR = myJets[j].p4.DeltaR(myFSRphotons[0].p4);
  	  if(DR<mJetPhoIsoR)isFSRphotonCountedAsJet=true;
	  }

          FSRphotonJet_ = isFSRphotonCountedAsJet;
	}
      }


      vector<TLorentzVector> allP4;
      if(myLeptons.size()>1)allP4.push_back(llP4);
      if(myPhotons.size()>0)allP4.push_back(myPhotons[0].p4);
      double prod(1.0),sum(0.0);
      for(unsigned j = 0; j < myJets.size(); j++) {
        prod *= myJets[j].p4.Pt();
        sum  += myJets[j].p4.Pt();
        allP4.push_back(myJets[j].p4); 
        if(nLeptons_ > 1) jetllDPhi_     ->push_back(fabs(llP4.DeltaPhi(myJets[j].p4)));
        if(nPhotons_ > 0) jetPhotonDPhi_ ->push_back(fabs(myPhotons[0].p4.DeltaPhi(myJets[j].p4)));
        jetPt_       ->push_back(myJets[j].p4.Pt()); 
        jetEta_      ->push_back(myJets[j].p4.Eta()); 
        jetPhi_      ->push_back(myJets[j].p4.Phi()); 
        jetE_        ->push_back(myJets[j].p4.Energy()); 
        jetY_        ->push_back(myJets[j].p4.Rapidity()); 
        jetArea_     ->push_back(myJets[j].area);
        jetBeta_     ->push_back(myJets[j].beta);
        jetBtag_     ->push_back(myJets[j].btag);
        jetJEC_      ->push_back(myJets[j].jec);
        jetUNC_      ->push_back(myJets[j].unc);
        jetCHF_      ->push_back(myJets[j].chf);
        jetPHF_      ->push_back(myJets[j].phf);
        jetNHF_      ->push_back(myJets[j].nhf);
        jetMUF_      ->push_back(myJets[j].muf);
        jetELF_      ->push_back(myJets[j].elf);
        jetId_       ->push_back(myJets[j].id);  
      }
      for(unsigned j = 0; j < myRJets.size(); j++) {
        rjetPt_       ->push_back(myRJets[j].p4.Pt()); 
        rjetEta_      ->push_back(myRJets[j].p4.Eta()); 
        rjetPhi_      ->push_back(myRJets[j].p4.Phi()); 
        rjetE_        ->push_back(myRJets[j].p4.Energy()); 
        rjetY_        ->push_back(myRJets[j].p4.Rapidity()); 
        rjetArea_     ->push_back(myRJets[j].area);
        rjetBeta_     ->push_back(myRJets[j].beta);
        rjetBtag_     ->push_back(myRJets[j].btag);
        rjetJEC_      ->push_back(myRJets[j].jec);
        rjetUNC_      ->push_back(myRJets[j].unc);
        rjetCHF_      ->push_back(myRJets[j].chf);
        rjetPHF_      ->push_back(myRJets[j].phf);
        rjetNHF_      ->push_back(myRJets[j].nhf);
        rjetMUF_      ->push_back(myRJets[j].muf);
        rjetELF_      ->push_back(myRJets[j].elf);
        rjetId_       ->push_back(myRJets[j].id);  
      }

      if(nLeptons_>1) isZlead_ = 1;
      if(nPhotons_>0) isPhotonlead_ = 1;
      if (nRJets_ > 1) mrj1rj2_  = (myRJets[0].p4 + myRJets[1].p4).M();

      if (nJets_ > 0) {
        if(nLeptons_ > 1) {
	  mZj1_        = (llP4+myJets[0].p4).M();
	  ptZj1_       = (llP4+myJets[0].p4).Pt();
	  costhetaZj1_ = tanh(0.5*(llP4.Rapidity() - myJets[0].p4.Rapidity()));
	}

        if(nPhotons_>0) {
	  mPhotonj1_        = (myJets[0].p4 + myPhotons[0].p4).M();
	  ptPhotonj1_       = (myJets[0].p4 + myPhotons[0].p4).Pt();
	  costhetaPhotonj1_ = tanh(0.5*(myPhotons[0].p4.Rapidity() - myJets[0].p4.Rapidity()));
	}

        jetPtGeMean_ = pow(prod,1./nJets_);
        jetPtArMean_ = sum/nJets_;
	HTJetSum_    = sum;

 	  
        if (nJets_ > 1) {
	  if(nLeptons_ > 1) mZj1j2_   = (llP4 + myJets[0].p4 + myJets[1].p4).M();
          mj1j2_    = (myJets[0].p4 + myJets[1].p4).M();
          j1j2DPhi_ = fabs(myJets[0].p4.DeltaPhi(myJets[1].p4));
          j1j2DR_   = myJets[0].p4.DeltaR(myJets[1].p4);

          if(nLeptons_>1) {
  	    if (llPt_ < myJets[1].p4.Pt()) {
  	      isZlead_ = 0;
	    }
	  }



          if (nJets_ > 2) {
	    if(nLeptons_ > 1) mZj1j2j3_ = (llP4+myJets[0].p4 + myJets[1].p4 + myJets[2].p4).M();
	    mj1j2j3_  = (myJets[0].p4 + myJets[1].p4 + myJets[2].p4).M();
            j1j3DPhi_ = fabs(myJets[0].p4.DeltaPhi(myJets[2].p4));
            j2j3DPhi_ = fabs(myJets[1].p4.DeltaPhi(myJets[2].p4));
            j1j3DR_   = myJets[0].p4.DeltaR(myJets[2].p4);
            j2j3DR_   = myJets[1].p4.DeltaR(myJets[2].p4);
          }  
        }
      }
      sort(allP4.begin(),allP4.end(),p4SortingRule);
      htLead_ = 0.0;
      for(int i=0;i<TMath::Min(int(allP4.size()),2);i++) {
        htLead_ += allP4[i].Pt();
      }

      // ----  Trigger block: Bother for trigger info only if event is selected (saves time)-------------
      iEvent.getByLabel(triggerResultsTag_,triggerResultsHandle_);
      if (!triggerResultsHandle_.isValid()) {
        cout << "ProcessedTreeProducer::analyze: Error in getting TriggerResults product from Event!" << endl;
        return;
      }
      iEvent.getByLabel(triggerEventTag_,triggerEventHandle_);
      if (!triggerEventHandle_.isValid()) {
        cout << "ProcessedTreeProducer::analyze: Error in getting TriggerEvent product from Event!" << endl;
        return;
      }
      // sanity check
      assert(triggerResultsHandle_->size() == hltConfig_.size());
      //------ loop over all trigger names ---------
      for(unsigned itrig=0;itrig<triggerNames_.size();itrig++) {
        bool accept(false);
        int preL1(-1);
        int preHLT(-1);
        int tmpFired(-1); 
        
    
        if (triggerIndex_[itrig] < hltConfig_.size()) {
          accept = triggerResultsHandle_->accept(triggerIndex_[itrig]);
          // --- check if your trigger bit is in the list which we don't ask for prescale (emu paths)
	  bool doCheckForPrescale = true;
          string reducedTriggerName = "";
	  int arraySize = int(triggerNamesFull_[itrig].size());
	  if(arraySize-1>0)reducedTriggerName=triggerNamesFull_[itrig].substr(0,triggerNamesFull_[itrig].size()-1); // remove last char from the str
	  for(int nn = 0; nn<int(prescaleDontAsk_.size()); nn++) {
	    if(reducedTriggerName==prescaleDontAsk_[nn])doCheckForPrescale=false;
	  }
	  //if(!doCheckForPrescale)cout << "skipping to check = " << triggerNamesFull_[itrig] << endl;

          if (triggerNamesFull_[itrig] != "" && doCheckForPrescale ) {
            const std::pair<int,int> prescales(hltConfig_.prescaleValues(iEvent,iSetup,triggerNamesFull_[itrig]));
            preL1  = prescales.first;
            preHLT = prescales.second;
          }  
          if (!accept) {
            tmpFired = 0;
          }
          else { 
            std::string ss(triggerNames_[itrig]); 
            hTriggerPass_->Fill((ss.erase(ss.find("v")-1,ss.find("v"))).c_str(),1);
            tmpFired = 1;
            // save trigger bit (0001) if family1 has fired, (0100) if family 3 has triggered
            isTriggered_ |= checkTriggerName(triggerNames_[itrig],triggerFamily1_) << 0; // if true 0001
            isTriggered_ |= checkTriggerName(triggerNames_[itrig],triggerFamily2_) << 1; // if true 0010
            isTriggered_ |= checkTriggerName(triggerNames_[itrig],triggerFamily3_) << 2; // if true 0100
            isTriggered_ |= checkTriggerName(triggerNames_[itrig],triggerFamily4_) << 3; // if true 1000
//          std::cout << "f1 " << checkTriggerName(triggerNames_[itrig],triggerFamily1_) << " " <<triggerNames_[itrig] << " " << isTriggered_ << std::endl;
//          std::cout << "f2 " << checkTriggerName(triggerNames_[itrig],triggerFamily2_) << " " <<triggerNames_[itrig] << " " << isTriggered_ << std::endl;
//          std::cout << "f3 " << checkTriggerName(triggerNames_[itrig],triggerFamily3_) << " " <<triggerNames_[itrig] << " " << isTriggered_ << std::endl;
//          std::cout << "f4 " << checkTriggerName(triggerNames_[itrig],triggerFamily4_) << " " <<triggerNames_[itrig] << " " << isTriggered_ << std::endl;
          }
        }
        
        fired_      ->push_back(tmpFired);
        prescaleL1_ ->push_back(preL1);
        prescaleHLT_->push_back(preHLT);
      }
    
    }// if selection GEN 
    if (selectionGEN) {
      selGEN_ = 1;

      if(llP4GEN.Pt()>0) {
       	llMGEN_                      = llP4GEN.M();
        llPtGEN_                     = llP4GEN.Pt();
        llPhiGEN_                    = llP4GEN.Phi();
        if(llPtGEN_>0)llYGEN_        = llP4GEN.Rapidity();
        if(llPtGEN_>0)llEtaGEN_      = llP4GEN.Eta();
        llDPhiGEN_                   = fabs(myGenLeptons[0].p4.DeltaPhi(myGenLeptons[1].p4));
      }
      TLorentzVector lepP4GEN(0,0,0,0); 

      for(unsigned l = 0; l < myGenLeptons.size(); l++) {
        lepP4GEN += myGenLeptons[l].p4;
        lepPtGEN_     ->push_back(myGenLeptons[l].p4.Pt());
        lepEtaGEN_    ->push_back(myGenLeptons[l].p4.Eta());
        lepPhiGEN_    ->push_back(myGenLeptons[l].p4.Phi());
        lepEGEN_      ->push_back(myGenLeptons[l].p4.Energy());
        lepChIdGEN_   ->push_back(myGenLeptons[l].pdgId);
      }   

      if(myGenPhotons.size()>0) {       
        photonPtGEN_     = myGenPhotons[0].p4.Pt();
        photonEtaGEN_    = myGenPhotons[0].p4.Eta();
        photonPhiGEN_    = myGenPhotons[0].p4.Phi();
        photonEGEN_      = myGenPhotons[0].p4.Energy();
        if(myPhotons.size()>0)photonRECODRGEN_ = myPhotons[0].p4.DeltaR(myGenPhotons[0].p4); // GEN TO RECO matching
      }   

      mLepGEN_ = lepP4GEN.M(); 
      vector<TLorentzVector> allP4GEN;
      allP4GEN.push_back(llP4GEN);
      double prod(1.0),sum(0.0);
      for(unsigned j = 0; j < myGenJets.size(); j++) {
        prod *= myGenJets[j].Pt();
        sum  += myGenJets[j].Pt();
        allP4GEN.push_back(myGenJets[j]); 
        jetllDPhiGEN_   ->push_back(fabs(llP4GEN.DeltaPhi(myGenJets[j])));
        jetPtGEN_       ->push_back(myGenJets[j].Pt()); 
        jetEtaGEN_      ->push_back(myGenJets[j].Eta()); 
        jetPhiGEN_      ->push_back(myGenJets[j].Phi()); 
        jetEGEN_        ->push_back(myGenJets[j].Energy()); 
        jetYGEN_        ->push_back(myGenJets[j].Rapidity());  
      }
      isZleadGEN_ = 1;
      if (nJetsGEN_ > 0) {
        mZj1GEN_        = (llP4GEN + myGenJets[0]).M();
        ptZj1GEN_       = (llP4GEN + myGenJets[0]).Pt();
        costhetaZj1GEN_ = tanh(0.5*(llP4GEN.Rapidity() - myGenJets[0].Rapidity()));
        jetPtGeMeanGEN_ = pow(prod,1./nJetsGEN_);
        jetPtArMeanGEN_ = sum/nJetsGEN_;
	HTJetSumGEN_    = sum;
        if (nJetsGEN_ > 1) {
          mj1j2GEN_    = (myGenJets[0] + myGenJets[1]).M();
          j1j2DPhiGEN_ = fabs(myGenJets[0].DeltaPhi(myGenJets[1]));
          j1j2DRGEN_   = myGenJets[0].DeltaR(myGenJets[1]);
          if (llPtGEN_ < myGenJets[1].Pt()) {
            isZleadGEN_ = 0;
          }
          if (nJetsGEN_ > 2) {
            j1j3DPhiGEN_ = fabs(myGenJets[0].DeltaPhi(myGenJets[2]));
            j2j3DPhiGEN_ = fabs(myGenJets[1].DeltaPhi(myGenJets[2]));
            j1j3DRGEN_   = myGenJets[0].DeltaR(myGenJets[2]);
            j2j3DRGEN_   = myGenJets[1].DeltaR(myGenJets[2]);
          }  
        }
      }
      sort(allP4GEN.begin(),allP4GEN.end(),p4SortingRule);
      htLeadGEN_ = 0.0;
      for(int i=0;i<TMath::Min(int(allP4GEN.size()),2);i++) {
        htLeadGEN_ += allP4GEN[i].Pt();
      }
    }// if selection GEN 
    myTree_->Fill();
  }// if selectionRECO || selectionGEN
}
// ---- method called once each job just after ending the event loop  ---
void PATZJetsExpress::endJob() 
{
}
// ---- method for tree building ----------------------------------------
void PATZJetsExpress::buildTree()
{
  fired_             = new std::vector<int>();
  prescaleL1_        = new std::vector<int>();
  prescaleHLT_       = new std::vector<int>();
  lepPt_             = new std::vector<float>();
  lepEta_            = new std::vector<float>();
  lepPhi_            = new std::vector<float>();
  lepE_              = new std::vector<float>(); 
  lepIso_            = new std::vector<float>();
  lepIsoPF_          = new std::vector<float>();
  lepIsoRho_         = new std::vector<float>();
  lepChId_           = new std::vector<int>();
  lepId_             = new std::vector<int>();
  jetPt_             = new std::vector<float>(); 
  jetEta_            = new std::vector<float>();
  jetPhi_            = new std::vector<float>();
  jetE_              = new std::vector<float>();
  jetY_              = new std::vector<float>();
  jetArea_           = new std::vector<float>();
  jetBeta_           = new std::vector<float>();
  jetBtag_           = new std::vector<float>();
  jetJEC_            = new std::vector<float>();
  jetUNC_            = new std::vector<float>();
  jetCHF_            = new std::vector<float>();
  jetPHF_            = new std::vector<float>();
  jetNHF_            = new std::vector<float>();
  jetMUF_            = new std::vector<float>();
  jetELF_            = new std::vector<float>();
  jetId_             = new std::vector<int>();
  rjetPt_            = new std::vector<float>(); 
  rjetEta_           = new std::vector<float>();
  rjetPhi_           = new std::vector<float>();
  rjetE_             = new std::vector<float>();
  rjetY_             = new std::vector<float>();
  rjetArea_          = new std::vector<float>();
  rjetBeta_          = new std::vector<float>();
  rjetBtag_          = new std::vector<float>();
  rjetJEC_           = new std::vector<float>();
  rjetUNC_           = new std::vector<float>();
  rjetCHF_           = new std::vector<float>();
  rjetPHF_           = new std::vector<float>();
  rjetNHF_           = new std::vector<float>();
  rjetMUF_           = new std::vector<float>();
  rjetELF_           = new std::vector<float>();
  rjetId_            = new std::vector<int>();
  vtxZ_              = new std::vector<float>();
  vtxNdof_           = new std::vector<float>();
  lepPtGEN_          = new std::vector<float>();
  lepEtaGEN_         = new std::vector<float>();
  lepPhiGEN_         = new std::vector<float>();
  lepEGEN_           = new std::vector<float>(); 
  lepChIdGEN_        = new std::vector<int>();
  lepMatchedGEN_     = new std::vector<int>();
  lepMatchedDRGEN_   = new std::vector<float>();
  jetPtGEN_          = new std::vector<float>(); 
  jetEtaGEN_         = new std::vector<float>();
  jetPhiGEN_         = new std::vector<float>();
  jetEGEN_           = new std::vector<float>();
  jetYGEN_           = new std::vector<float>();
  jetllDPhiGEN_      = new std::vector<float>();
  jetllDPhi_         = new std::vector<float>();
  jetPhotonDPhi_     = new std::vector<float>();
  photonPar_         = new std::vector<float>();
  FSRphotonPar_      = new std::vector<float>();
  // ---- global event variables ----------------------------------------
  myTree_->Branch("isRealData"       ,&isRealData_        ,"isRealData/I");
  myTree_->Branch("selRECO"          ,&selRECO_           ,"selRECO/I");
  myTree_->Branch("eventNum"         ,&eventNum_          ,"eventNum/l"); // 64 bit unsigned integer (ULong64_t)
  myTree_->Branch("runNum"           ,&runNum_            ,"runNum/I");
  myTree_->Branch("lumi"             ,&lumi_              ,"lumi/I");
  myTree_->Branch("nVtx"             ,&nVtx_              ,"nVtx/I");
  myTree_->Branch("nLeptons"         ,&nLeptons_          ,"nLeptons/I");
  myTree_->Branch("nPhotons"         ,&nPhotons_          ,"nPhotons/I");
  myTree_->Branch("nJets"            ,&nJets_             ,"nJets/I"); 
  myTree_->Branch("nRJets"           ,&nRJets_            ,"nRJets/I"); 
  myTree_->Branch("isZlead"          ,&isZlead_           ,"isZlead/I");
  myTree_->Branch("isPhotonlead"     ,&isPhotonlead_      ,"isPhotonlead/I");
  myTree_->Branch("rho"              ,&rho_               ,"rho/F");
  myTree_->Branch("rho25"            ,&rho25_             ,"rho25/F");
  myTree_->Branch("mZj1"             ,&mZj1_              ,"mZj1/F");
  myTree_->Branch("mZj1j2"           ,&mZj1j2_            ,"mZj1j2/F");
  myTree_->Branch("mZj1j2j3"         ,&mZj1j2j3_          ,"mZj1j2j3/F");
  myTree_->Branch("ptZj1"            ,&ptZj1_             ,"ptZj1/F");
  myTree_->Branch("costhetaZj1"      ,&costhetaZj1_       ,"costhetaZj1/F");
  myTree_->Branch("costhetaPhotonj1" ,&costhetaPhotonj1_  ,"costhetaPhotonj1/F");
  myTree_->Branch("mj1j2"            ,&mj1j2_             ,"mj1j2/F");
  myTree_->Branch("mj1j2j3"          ,&mj1j2j3_           ,"mj1j2j3/F");
  myTree_->Branch("mrj1rj2"          ,&mrj1rj2_           ,"mrj1rj2/F");
  myTree_->Branch("mLep"             ,&mLep_              ,"mLep/F");
  myTree_->Branch("htLead"           ,&htLead_            ,"htLead/F");
  myTree_->Branch("j1j2DPhi"         ,&j1j2DPhi_          ,"j1j2DPhi/F");
  myTree_->Branch("j1j3DPhi"         ,&j1j3DPhi_          ,"j1j3DPhi/F");
  myTree_->Branch("j2j3DPhi"         ,&j2j3DPhi_          ,"j2j3DPhi/F");
  myTree_->Branch("j1j2DR"           ,&j1j2DR_            ,"j1j2DR/F");
  myTree_->Branch("j1j3DR"           ,&j1j3DR_            ,"j1j3DR/F");
  myTree_->Branch("j2j3DR"           ,&j2j3DR_            ,"j2j3DR/F");
  myTree_->Branch("jetPtGeMean"      ,&jetPtGeMean_       ,"jetPtGeMean/F");
  myTree_->Branch("jetPtArMean"      ,&jetPtArMean_       ,"jetPtArMean/F");
  // ---- met variables -------------------------------------------------
  myTree_->Branch("pfmet"            ,&pfmet_             ,"pfmet/F");
  myTree_->Branch("pfmetPhi"         ,&pfmetPhi_          ,"pfmetPhi/F");
  myTree_->Branch("pfhadPt"          ,&pfhadPt_           ,"pfhadPt/F");
  myTree_->Branch("pfSumEt"          ,&pfSumEt_           ,"pfSumEt/F");  
  myTree_->Branch("HTJetSum"         ,&HTJetSum_          ,"HTJetSum/F");  
  // ---- dilepton variables --------------------------------------------
  myTree_->Branch("llM"              ,&llM_               ,"llM/F");
  myTree_->Branch("llPt"             ,&llPt_              ,"llPt/F");
  myTree_->Branch("llPhi"            ,&llPhi_             ,"llPhi/F");
  myTree_->Branch("llDPhi"           ,&llDPhi_            ,"llDPhi/F");
  myTree_->Branch("llY"              ,&llY_               ,"llY/F");
  myTree_->Branch("llEta"            ,&llEta_             ,"llEta/F");
  // ---- photon variables ----------------------------------------------
  myTree_->Branch("photonPt"         ,&photonPt_          ,"photonPt/F");
  myTree_->Branch("photonE"          ,&photonE_           ,"photonE/F");
  myTree_->Branch("photonEta"        ,&photonEta_         ,"photonEta/F");
  myTree_->Branch("photonPhi"        ,&photonPhi_         ,"photonPhi/F");
  myTree_->Branch("photonIso"        ,&photonIso_         ,"photonIso/F");
  myTree_->Branch("photonID"         ,&photonID_          ,"photonID/F");
  myTree_->Branch("photonBit"        ,&photonBit_         ,"photonBit/I");
  myTree_->Branch("pfhadPhoPt"       ,&pfhadPhoPt_        ,"pfhadPhoPt/F");
  myTree_->Branch("mPhotonj1"        ,&mPhotonj1_         ,"mPhotonj1/F");
  myTree_->Branch("ptPhotonj1"       ,&ptPhotonj1_        ,"ptPhotonj1/F");
  myTree_->Branch("jetPhotonDPhi"    ,"vector<float>"     ,&jetPhotonDPhi_);
  myTree_->Branch("photonPar"        ,"vector<float>"     ,&photonPar_);
  // ---- FSRphoton variables ----------------------------------------------
  myTree_->Branch("FSRphotonPt"      ,&FSRphotonPt_       ,"FSRphotonPt/F");
  myTree_->Branch("FSRphotonE"       ,&FSRphotonE_        ,"FSRphotonE/F");
  myTree_->Branch("FSRphotonEta"     ,&FSRphotonEta_      ,"FSRphotonEta/F");
  myTree_->Branch("FSRphotonPhi"     ,&FSRphotonPhi_      ,"FSRphotonPhi/F");
  myTree_->Branch("FSRphotonIso"     ,&FSRphotonIso_      ,"FSRphotonIso/F");
  myTree_->Branch("FSRphotonID"      ,&FSRphotonID_       ,"FSRphotonID/F");
  myTree_->Branch("FSRphotonBit"     ,&FSRphotonBit_      ,"FSRphotonBit/I");
  myTree_->Branch("FSRphotonJet"     ,&FSRphotonJet_      ,"FSRphotonJet/I");
  myTree_->Branch("FSRphotonPar"     ,"vector<float>"     ,&FSRphotonPar_);
  myTree_->Branch("FSRphotonllM"     ,&FSRphotonllM_      ,"FSRphotonllM/F");
  // ---- trigger variables ---------------------------------------------
  myTree_->Branch("fired"            ,"vector<int>"       ,&fired_);
  myTree_->Branch("prescaleL1"       ,"vector<int>"       ,&prescaleL1_);
  myTree_->Branch("prescaleHLT"      ,"vector<int>"       ,&prescaleHLT_);
  myTree_->Branch("isTriggered"      ,&isTriggered_        ,"isTriggered/I");
  // ---- lepton variables ----------------------------------------------
  myTree_->Branch("lepPt"            ,"vector<float>"     ,&lepPt_);
  myTree_->Branch("lepEta"           ,"vector<float>"     ,&lepEta_);
  myTree_->Branch("lepPhi"           ,"vector<float>"     ,&lepPhi_);
  myTree_->Branch("lepE"             ,"vector<float>"     ,&lepE_);
  myTree_->Branch("lepIso"           ,"vector<float>"     ,&lepIso_);
  myTree_->Branch("lepIsoPF"         ,"vector<float>"     ,&lepIsoPF_);
  myTree_->Branch("lepIsoRho"        ,"vector<float>"     ,&lepIsoRho_);
  myTree_->Branch("lepChId"          ,"vector<int>"       ,&lepChId_);
  myTree_->Branch("lepId"            ,"vector<int>"       ,&lepId_);
  // ---- jet variables -------------------------------------------------
  myTree_->Branch("jetPt"            ,"vector<float>"     ,&jetPt_);
  myTree_->Branch("jetEta"           ,"vector<float>"     ,&jetEta_);
  myTree_->Branch("jetPhi"           ,"vector<float>"     ,&jetPhi_);
  myTree_->Branch("jetE"             ,"vector<float>"     ,&jetE_);
  myTree_->Branch("jetY"             ,"vector<float>"     ,&jetY_);
  myTree_->Branch("jetArea"          ,"vector<float>"     ,&jetArea_);
  myTree_->Branch("jetBeta"          ,"vector<float>"     ,&jetBeta_);
  myTree_->Branch("jetBtag"          ,"vector<float>"     ,&jetBtag_);
  myTree_->Branch("jetJEC"           ,"vector<float>"     ,&jetJEC_);
  myTree_->Branch("jetUNC"           ,"vector<float>"     ,&jetUNC_);
  myTree_->Branch("jetllDPhi"        ,"vector<float>"     ,&jetllDPhi_);
  myTree_->Branch("jetCHF"           ,"vector<float>"     ,&jetCHF_);
  myTree_->Branch("jetPHF"           ,"vector<float>"     ,&jetPHF_);
  myTree_->Branch("jetNHF"           ,"vector<float>"     ,&jetNHF_);
  myTree_->Branch("jetMUF"           ,"vector<float>"     ,&jetMUF_);
  myTree_->Branch("jetELF"           ,"vector<float>"     ,&jetELF_);
  myTree_->Branch("jetId"            ,"vector<int>"       ,&jetId_);
  // ---- DR rejected jet variables ------------------------------------
  myTree_->Branch("rjetPt"           ,"vector<float>"     ,&rjetPt_);
  myTree_->Branch("rjetEta"          ,"vector<float>"     ,&rjetEta_);
  myTree_->Branch("rjetPhi"          ,"vector<float>"     ,&rjetPhi_);
  myTree_->Branch("rjetE"            ,"vector<float>"     ,&rjetE_);
  myTree_->Branch("rjetY"            ,"vector<float>"     ,&rjetY_);
  myTree_->Branch("rjetArea"         ,"vector<float>"     ,&rjetArea_);
  myTree_->Branch("rjetBeta"         ,"vector<float>"     ,&rjetBeta_);
  myTree_->Branch("rjetBtag"         ,"vector<float>"     ,&rjetBtag_);
  myTree_->Branch("rjetJEC"          ,"vector<float>"     ,&rjetJEC_);
  myTree_->Branch("rjetUNC"          ,"vector<float>"     ,&rjetUNC_);
  myTree_->Branch("rjetCHF"          ,"vector<float>"     ,&rjetCHF_);
  myTree_->Branch("rjetPHF"          ,"vector<float>"     ,&rjetPHF_);
  myTree_->Branch("rjetNHF"          ,"vector<float>"     ,&rjetNHF_);
  myTree_->Branch("rjetMUF"          ,"vector<float>"     ,&rjetMUF_);
  myTree_->Branch("rjetELF"          ,"vector<float>"     ,&rjetELF_);
  myTree_->Branch("rjetId"           ,"vector<int>"       ,&rjetId_);
  // ---- vertex variables ----------------------------------------------
  myTree_->Branch("vtxZ"             ,"vector<float>"     ,&vtxZ_);
  myTree_->Branch("vtxNdof"          ,"vector<float>"     ,&vtxNdof_);
  // ---- gen variables ----------------------------------------------
  myTree_->Branch("selGEN"           ,&selGEN_            ,"selGEN/I");
  myTree_->Branch("puINT"            ,&puINT_             ,"puINT/I");
  myTree_->Branch("puOOT"            ,&puOOT_             ,"puOOT/I");
  myTree_->Branch("nLeptonsGEN"      ,&nLeptonsGEN_       ,"nLeptonsGEN/I");
  myTree_->Branch("nJetsGEN"         ,&nJetsGEN_          ,"nJetsGEN/I");
  myTree_->Branch("isZleadGEN"       ,&isZleadGEN_        ,"isZleadGEN/I");
  myTree_->Branch("mZj1GEN"          ,&mZj1GEN_           ,"mZj1GEN/F");
  myTree_->Branch("ptZj1GEN"         ,&ptZj1GEN_          ,"ptZj1GEN/F");
  myTree_->Branch("costhetaZj1GEN"   ,&costhetaZj1GEN_    ,"costhetaZj1GEN/F");
  myTree_->Branch("mj1j2GEN"         ,&mj1j2GEN_          ,"mj1j2GEN/F");
  myTree_->Branch("mLepGEN"          ,&mLepGEN_           ,"mLepGEN/F");
  myTree_->Branch("htLeadGEN"        ,&htLeadGEN_         ,"htLeadGEN/F");
  myTree_->Branch("j1j2DPhiGEN"      ,&j1j2DPhiGEN_       ,"j1j2DPhiGEN/F");
  myTree_->Branch("j1j3DPhiGEN"      ,&j1j3DPhiGEN_       ,"j1j3DPhiGEN/F");
  myTree_->Branch("j2j3DPhiGEN"      ,&j2j3DPhiGEN_       ,"j2j3DPhiGEN/F");
  myTree_->Branch("j1j2DRGEN"        ,&j1j2DRGEN_         ,"j1j2DRGEN/F");
  myTree_->Branch("j1j3DRGEN"        ,&j1j3DRGEN_         ,"j1j3DRGEN/F");
  myTree_->Branch("j2j3DRGEN"        ,&j2j3DRGEN_         ,"j2j3DRGEN/F");
  myTree_->Branch("jetPtGeMeanGEN"   ,&jetPtGeMeanGEN_    ,"jetPtGeMeanGEN/F");
  myTree_->Branch("jetPtArMeanGEN"   ,&jetPtArMeanGEN_    ,"jetPtArMeanGEN/F");
  myTree_->Branch("llMGEN"           ,&llMGEN_            ,"llMGEN/F");
  myTree_->Branch("llPtGEN"          ,&llPtGEN_           ,"llPtGEN/F");
  myTree_->Branch("llPhiGEN"         ,&llPhiGEN_          ,"llPhiGEN/F");
  myTree_->Branch("llDPhiGEN"        ,&llDPhiGEN_         ,"llDPhiGEN/F");
  myTree_->Branch("llYGEN"           ,&llYGEN_            ,"llYGEN/F");
  myTree_->Branch("llEtaGEN"         ,&llEtaGEN_          ,"llEtaGEN/F");
  myTree_->Branch("lepPtGEN"         ,"vector<float>"     ,&lepPtGEN_);
  myTree_->Branch("lepEtaGEN"        ,"vector<float>"     ,&lepEtaGEN_);
  myTree_->Branch("lepPhiGEN"        ,"vector<float>"     ,&lepPhiGEN_);
  myTree_->Branch("lepEGEN"          ,"vector<float>"     ,&lepEGEN_);
  myTree_->Branch("lepChIdGEN"       ,"vector<int>"       ,&lepChIdGEN_);
  myTree_->Branch("lepMatchedDRGEN"  ,"vector<float>"     ,&lepMatchedDRGEN_);
  myTree_->Branch("lepMatchedGEN"    ,"vector<int>"       ,&lepMatchedGEN_);
  myTree_->Branch("jetPtGEN"         ,"vector<float>"     ,&jetPtGEN_);
  myTree_->Branch("jetEtaGEN"        ,"vector<float>"     ,&jetEtaGEN_);
  myTree_->Branch("jetPhiGEN"        ,"vector<float>"     ,&jetPhiGEN_);
  myTree_->Branch("jetEGEN"          ,"vector<float>"     ,&jetEGEN_);
  myTree_->Branch("jetllDPhiGEN"     ,"vector<float>"     ,&jetllDPhiGEN_);
  myTree_->Branch("HTJetSumGEN"      ,&HTJetSumGEN_       ,"HTJetSumGEN/F");  
  myTree_->Branch("mcWeight"         ,&mcWeight_          ,"mcWeight/F");
  myTree_->Branch("nPhotonsGEN"      ,&nPhotonsGEN_       ,"nPhotonsGEN/I");
  myTree_->Branch("photonPtGEN"      ,&photonPtGEN_       ,"photonPtGEN/F");
  myTree_->Branch("photonEGEN"       ,&photonEGEN_        ,"photonEGEN/F");
  myTree_->Branch("photonEtaGEN"     ,&photonEtaGEN_      ,"photonEtaGEN/F");
  myTree_->Branch("photonPhiGEN"     ,&photonPhiGEN_      ,"photonPhiGEN/F");
  myTree_->Branch("photonRECODRGEN"  ,&photonRECODRGEN_   ,"photonRECODRGEN/F");
  myTree_->Branch("VBPartonDM"       ,&VBPartonDM_        ,"VBPartonDM/I");
  myTree_->Branch("VBPartonM"        ,&VBPartonM_         ,"VBPartonM/F");
  myTree_->Branch("VBPartonE"        ,&VBPartonE_         ,"VBPartonE/F");
  myTree_->Branch("VBPartonPt"       ,&VBPartonPt_        ,"VBPartonPt/F");
  myTree_->Branch("VBPartonEta"      ,&VBPartonEta_       ,"VBPartonEta/F");
  myTree_->Branch("VBPartonPhi"      ,&VBPartonPhi_       ,"VBPartonPhi/F");
}
// ---- method for tree initialization ----------------------------------
void PATZJetsExpress::clearTree()
{
  selRECO_           = -999;
  isRealData_        = -999;
  eventNum_          = 0; //ULong64_t can't be negative
  runNum_            = -999;
  lumi_              = -999;
  nVtx_              = -999;
  nLeptons_          = -999;
  nPhotonsGEN_       = -999;
  nPhotons_          = -999;
  nJets_             = -999;
  nRJets_            = -999;
  isZlead_           = -999;
  isPhotonlead_      = -999;
  rho_               = -999;
  rho25_             = -999;
  pfmet_             = -999;
  pfmetPhi_          = -999;
  pfhadPt_           = -999;
  pfSumEt_           = -999;
  HTJetSum_          = -999;
  HTJetSumGEN_       = -999;
  mZj1_              = -999; 
  mZj1j2_            = -999; 
  mZj1j2j3_          = -999; 
  ptZj1_             = -999; 
  costhetaZj1_       = -999; 
  costhetaPhotonj1_  = -999; 
  mj1j2_             = -999;
  mj1j2j3_           = -999;
  mrj1rj2_           = -999;
  mLep_              = -999;
  htLead_            = -999;
  jetPtGeMean_       = -999;
  jetPtArMean_       = -999; 
  j1j2DPhi_          = -999;
  j1j3DPhi_          = -999; 
  j2j3DPhi_          = -999;
  j1j2DR_            = -999;
  j1j3DR_            = -999; 
  j2j3DR_            = -999;
  llM_               = -999;
  llPt_              = -999; 
  llPhi_             = -999;
  llDPhi_            = -999;
  llY_               = -999;
  llEta_             = -999;
  photonE_           = -999;
  photonPt_          = -999;
  photonEta_         = -999;
  photonPhi_         = -999;
  photonIso_         = -999;
  photonID_          = -999;
  photonBit_         =    0; // please keep this 0
  FSRphotonE_        = -999;
  FSRphotonPt_       = -999;
  FSRphotonllM_      = -999;
  FSRphotonEta_      = -999;
  FSRphotonPhi_      = -999;
  FSRphotonIso_      = -999;
  FSRphotonID_       = -999;
  FSRphotonBit_      =    0; // please keep this 0
  FSRphotonJet_      =    0; // please keep this 0
  pfhadPhoPt_        = -999;
  mPhotonj1_         = -999;
  ptPhotonj1_        = -999;
  isTriggered_       =    0; // please keep this 0
  jetPhotonDPhi_     ->clear();
  photonPar_         ->clear();
  FSRphotonPar_      ->clear();
  fired_             ->clear();
  prescaleL1_        ->clear();
  prescaleHLT_       ->clear();
  lepPt_             ->clear();
  lepEta_            ->clear();
  lepPhi_            ->clear();
  lepE_              ->clear();
  lepIso_            ->clear();
  lepIsoPF_          ->clear();
  lepIsoRho_         ->clear();
  lepChId_           ->clear();
  lepMatchedDRGEN_   ->clear();
  lepMatchedGEN_     ->clear();
  lepId_             ->clear();
  jetPt_             ->clear();
  jetEta_            ->clear();
  jetPhi_            ->clear();
  jetE_              ->clear();
  jetY_              ->clear();
  jetArea_           ->clear();
  jetBeta_           ->clear();
  jetBtag_           ->clear();
  jetJEC_            ->clear();
  jetUNC_            ->clear();
  jetllDPhi_         ->clear();
  jetCHF_            ->clear();
  jetPHF_            ->clear();
  jetNHF_            ->clear();
  jetMUF_            ->clear();
  jetELF_            ->clear();
  jetId_             ->clear();
  rjetPt_            ->clear();
  rjetEta_           ->clear();
  rjetPhi_           ->clear();
  rjetE_             ->clear();
  rjetY_             ->clear();
  rjetArea_          ->clear();
  rjetBeta_          ->clear();
  rjetBtag_          ->clear();
  rjetJEC_           ->clear();
  rjetUNC_           ->clear();
  rjetCHF_           ->clear();
  rjetPHF_           ->clear();
  rjetNHF_           ->clear();
  rjetMUF_           ->clear();
  rjetELF_           ->clear();
  rjetId_            ->clear();
  vtxZ_              ->clear();
  vtxNdof_           ->clear();
  selGEN_            = -999;
  puINT_             = -999;
  puOOT_             = -999;
  isRealData_        = -999;
  nLeptonsGEN_       = -999;
  nJetsGEN_          = -999;
  isZleadGEN_        = -999;
  mZj1GEN_           = -999; 
  ptZj1GEN_          = -999; 
  costhetaZj1GEN_    = -999; 
  mj1j2GEN_          = -999;
  mLepGEN_           = -999;
  htLeadGEN_         = -999;
  jetPtGeMeanGEN_    = -999;
  jetPtArMeanGEN_    = -999; 
  j1j2DPhiGEN_       = -999;
  j1j3DPhiGEN_       = -999; 
  j2j3DPhiGEN_       = -999;
  j1j2DRGEN_         = -999;
  j1j3DRGEN_         = -999; 
  j2j3DRGEN_         = -999;
  llMGEN_            = -999;
  llPtGEN_           = -999; 
  llPhiGEN_          = -999;
  llDPhiGEN_         = -999;
  llYGEN_            = -999;
  llEtaGEN_          = -999;
  lepPtGEN_          ->clear();
  lepEtaGEN_         ->clear();
  lepPhiGEN_         ->clear();
  lepEGEN_           ->clear();
  lepChIdGEN_        ->clear();
  jetPtGEN_          ->clear();
  jetEtaGEN_         ->clear();
  jetPhiGEN_         ->clear();
  jetEGEN_           ->clear();
  jetYGEN_           ->clear();
  jetllDPhiGEN_      ->clear();
  photonEGEN_        = -999;
  photonPtGEN_       = -999;
  photonEtaGEN_      = -999;
  photonPhiGEN_      = -999;
  photonRECODRGEN_   = +999; // please keep this positive (will cut offline to <0.2 for matched)
  VBPartonDM_        = -999;
  VBPartonM_         = -999;
  VBPartonE_         = -999;
  VBPartonPt_        = -999;
  VBPartonEta_       = -999;
  VBPartonPhi_       = -999;
  mcWeight_          = -999;
}
// ---- define this as a plug-in ----------------------------------------
DEFINE_FWK_MODULE(PATZJetsExpress);
