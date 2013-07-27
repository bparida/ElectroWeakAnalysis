/******************************************************************************
 * Project: CMS detector at the CERN
 *
 * Package: ElectroWeakAnalysis/VPlusJets
 *
 *
 * Authors:
 *
 *   Kalanand Mishra, Fermilab - kalanand@fnal.gov
 *
 * Description:
 *   To fill jet related quantities into a specified TTree
 *   Can work with CaloJet, GenJet, JPT jet, PF jet.
 *   Can work with jets in RECO/AOD/PAT data formats.
 * History:
 *   
 *
 * Copyright (C) 2010 FNAL 
 *****************************************************************************/

#ifndef ElectroWeakAnalysis_VPlusJets_JetTreeFiller_h
#define ElectroWeakAnalysis_VPlusJets_JetTreeFiller_h

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h" 
#include "TMath.h" 
#include <TLorentzVector.h>

#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "ElectroWeakAnalysis/VPlusJets/interface/QGLikelihoodCalculator.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "CMGTools/External/interface/PileupJetIdentifier.h"

namespace ewk {

  class JetTreeFiller {
  public:

    /// specify the name of the TTree, and the configuration for it
    JetTreeFiller(const char *name, TTree* tree, 
		  const std::string jetType,
		  const edm::ParameterSet iConfig );


    /// default constructor
    JetTreeFiller() {};


    /// Destructor, does nothing 
      ~JetTreeFiller() {delete qglikeli;};




    /// To be called once per event to fill the values for jets
     void fill(const edm::Event &iEvent);

    static const int NUM_JET_MAX = 8;

  protected:

    /// To be called once per event, to initialize variable values to default
    void init() const ;
    /// Helper function for main constructor
    void SetBranches(); 
    void SetBranch( float* x, std::string name);
    void SetBranch( int* x, std::string name);
    void SetBranch( bool* x, std::string name);
    void SetBranchSingle( float* x, std::string name);
    void SetBranchSingle( int* x, std::string name);
    void SetBranchSingle( bool* x, std::string name);

    void FillBranches() const;
    void init();


    template <typename T1> 
      void fillBasicJetQuantities(int iJet, const T1& pfjet, 
				const reco::Candidate* Vboson, 
				const reco::Candidate* Vboson2, 
				const reco::MET met); 
    template<typename T1> 
      void fillEnergyFractionsPFjets(const T1& pfjet, int iJet);

    void fillQGLH(int iJet, float fastjet_rho, 
		  std::vector<reco::PFCandidatePtr> pfCandidates);

    template<typename T1>
      void fillPileUpJetID ( const edm::Handle<edm::View<T1> > &);

    void fillBtagInfoRECO(int iJet, 
		      edm::Handle<reco::SecondaryVertexTagInfoCollection> svTagInfos,
		      const reco::JetTagCollection  &  bTagsSSVHE,
		      const reco::JetTagCollection  &  bTagsTCHE,
		      const reco::JetTagCollection  &  bTagsCSV,
		      const reco::JetTagCollection  &  bTagsJP,
                      const reco::JetTagCollection  &  bTagsSSVHP,
		      const reco::JetTagCollection  &  bTagsTCHP);

    void fillBtagInfoPAT(int iJet, const pat::Jet* pjet);
    void fillHelicityIn4bodyFrame(TLorentzVector& p4lepton1, 
				  TLorentzVector& p4lepton2);

    void fillInvariantMasses(TLorentzVector& p4lepton1, 
			     TLorentzVector& p4lepton2);

    void computeLeptonAndNu4Vectors(const reco::Candidate *Vboson, 
				    const reco::MET met, 
				    TLorentzVector& p4lepton1, 
				    TLorentzVector& p4lepton2); 

/*     void getConstituents(const reco::Jet* jet,  */
/* 			 std::vector<reco::PFCandidatePtr> &jetpfc); */
/*     void getConstituents(const pat::Jet* jet,  */
/* 			 std::vector<reco::PFCandidatePtr> &jetpfc); */

/*     template <typename T>  */
/*       std::vector<reco::PFCandidatePtr> getConstituents(const T&  jet, bool isReco); */
    
    TTree* tree_;
    std::string jetType_;
    std::string Vtype_;
    std::string LeptonType_;

    edm::InputTag mInputJets;
	edm::InputTag mInputMet;
	edm::InputTag mInputMetMVA;
    edm::InputTag mInputBoson;
    edm::InputTag sourceByValue;
  	bool runoverAOD;
	std::string bTagger;
    // 'mutable' because we will fill it from a 'const' method
    mutable std::vector<std::string> bnames;
    QGLikelihoodCalculator *qglikeli;


  private:
    // private data members
    
    int NumJets; 
    float Et[NUM_JET_MAX];
    float Pt[NUM_JET_MAX];
    float Eta[NUM_JET_MAX];
    float Phi[NUM_JET_MAX];
    float Theta[NUM_JET_MAX];
    float E[NUM_JET_MAX];
    float Y[NUM_JET_MAX];
    float Mass[NUM_JET_MAX];
    float etaetaMoment[NUM_JET_MAX];  
    float phiphiMoment[NUM_JET_MAX];      
    float etaphiMoment[NUM_JET_MAX];      
    float maxDistance[NUM_JET_MAX];
    int   nConstituents[NUM_JET_MAX];
    float Area[NUM_JET_MAX];

    float Px[NUM_JET_MAX];
    float Py[NUM_JET_MAX];
    float Pz[NUM_JET_MAX];

    float VjetMass[NUM_JET_MAX];
    float Dphi[NUM_JET_MAX];
    float Deta[NUM_JET_MAX];
    float DR[NUM_JET_MAX];
    float DphiMET[NUM_JET_MAX];
    float Response[NUM_JET_MAX];
    float bDiscriminator[NUM_JET_MAX];
    float bDiscriminatorSSVHE[NUM_JET_MAX];
    float bDiscriminatorTCHE[NUM_JET_MAX];
    float bDiscriminatorCSV[NUM_JET_MAX];
    float bDiscriminatorJP[NUM_JET_MAX];
    float bDiscriminatorSSVHP[NUM_JET_MAX];
    float bDiscriminatorTCHP[NUM_JET_MAX];

    float secVertexMass[NUM_JET_MAX];
    int numBTags;
    float VjetMass2[NUM_JET_MAX];
    float DR2[NUM_JET_MAX];
    float Dphi2[NUM_JET_MAX];
    float Deta2[NUM_JET_MAX];
    float Response2[NUM_JET_MAX];


    float PFChargedHadronEnergy[NUM_JET_MAX];
    float PFChargedHadronEnergyFraction[NUM_JET_MAX];
    float PFNeutralHadronEnergy[NUM_JET_MAX];
    float PFNeutralHadronEnergyFraction[NUM_JET_MAX];
    float PFChargedEmEnergy[NUM_JET_MAX];
    float PFChargedEmEnergyFraction[NUM_JET_MAX];
    float PFChargedMuEnergy[NUM_JET_MAX];
    float PFChargedMuEnergyFraction[NUM_JET_MAX];
    float PFNeutralEmEnergy[NUM_JET_MAX];
    float PFNeutralEmEnergyFraction[NUM_JET_MAX];
    int   PFChargedMultiplicity[NUM_JET_MAX];
    int   PFNeutralMultiplicity[NUM_JET_MAX];
    int   PFMuonMultiplicity[NUM_JET_MAX];
    float PFPhotonEnergy[NUM_JET_MAX];
    float PFPhotonEnergyFraction[NUM_JET_MAX];
    float PFElectronEnergy[NUM_JET_MAX];
    float PFElectronEnergyFraction[NUM_JET_MAX];
    float PFMuonEnergy[NUM_JET_MAX];
    float PFMuonEnergyFraction[NUM_JET_MAX];
    float PFHFHadronEnergy[NUM_JET_MAX];
    float PFHFHadronEnergyFraction[NUM_JET_MAX];
    float PFHFEMEnergy[NUM_JET_MAX];
    float PFHFEMEnergyFraction[NUM_JET_MAX];	 
    int   PFChargedHadronMultiplicity[NUM_JET_MAX];
    int   PFNeutralHadronMultiplicity[NUM_JET_MAX];
    int   PFPhotonMultiplicity[NUM_JET_MAX];
    int   PFElectronMultiplicity[NUM_JET_MAX];
    int   PFHFHadronMultiplicity[NUM_JET_MAX];
    int   PFHFEMMultiplicity[NUM_JET_MAX];    
    float PFsumPtCands[NUM_JET_MAX];
    float PFsumPt2Cands[NUM_JET_MAX];
    float PFrmsCands[NUM_JET_MAX];
    float PFptD[NUM_JET_MAX];
    float PFqgLikelihood[NUM_JET_MAX];

    float V2jMassMVAMET;
    float V2jMass;
    float V3jMass;
    float V4jMass;
    float V5jMass;
    float V6jMass;
    float c2jMass;
    float c3jMass;
    float c4jMass;
    float c5jMass;
    float c6jMass;

    float V2jCosJacksonAngle;
    float c2jCosJacksonAngle;
    float V3jCosJacksonAngle;
    float c3jCosJacksonAngle12;
    float c3jCosJacksonAngle23;
    float c3jCosJacksonAngle31;

    float cosphiDecayPlane; 
    float cosThetaLnu; 
    float cosThetaJJ;
  
    float colorCorr01;
    float colorCorr02;
    float colorCorr12;
    float colorCorr03;
    float colorCorr13;
    float colorCorr23;
    float colorCorr04;
    float colorCorr14;
    float colorCorr24;
    float colorCorr34;
    float colorCorr05;
    float colorCorr15;
    float colorCorr25;
    float colorCorr35;
    float colorCorr45;

    float j1Hel_HiggsCM;
    float j2Hel_HiggsCM;
    float l1Hel_HiggsCM;
    float l2Hel_HiggsCM;
    float b1Hel_HiggsCM;
    float b2Hel_HiggsCM;

    bool isPileUpJetLoose[NUM_JET_MAX];
    bool isPileUpJetMedium[NUM_JET_MAX];
    bool isPileUpJetTight[NUM_JET_MAX];

  };

} //namespace

#endif


