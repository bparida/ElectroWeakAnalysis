/*****************************************************************************
 * Project: CMS detector at the CERN
 *
 * Package: ElectroWeakAnalysis/VPlusJets
 * Class:   VtoElectronTreeFiller
 *
 * Authors:
 *
 *   Kalanand Mishra, Fermilab - kalanand@fnal.gov
 *
 * Description:
 *   To fill W--> enu or Z-->ee related quantities into a specified TTree
 *   Can work with jets in RECO/AOD/PAT data formats.
 * History:
 *   
 *
 * Copyright (C) 2010 FNAL 
 *****************************************************************************/

#ifndef ElectroWeakAnalysis_VPlusJets_VtoElectronTreeFiller_h
#define ElectroWeakAnalysis_VPlusJets_VtoElectronTreeFiller_h
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h" 
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"


namespace ewk {

  class VtoElectronTreeFiller {
  public:
    /// specify the name of the TTree, and the configuration for it
    VtoElectronTreeFiller(const char *name, TTree* tree, 
			  const edm::ParameterSet iConfig );
    
    /// default constructor
    VtoElectronTreeFiller() {};
    
    
    /// Destructor, does nothing 
    ~VtoElectronTreeFiller() {};
    
    
    /// To be called once per event to fill the values for jets
    void fill(const edm::Event &iEvent, int vecBosonIndex);

    
  protected:

    /// To be called once per event, to initialize variable values to default
    void init();
    /// Helper function for main constructor
    void SetBranches(); 
    void SetBranch( float* x, std::string name );
    void SetBranch( int* x, std::string name );
    void SetBranch( bool* x, std::string name );
    bool isTightElectron(const edm::Event& iEvent, const reco::GsfElectron& ele);
    bool isLooseElectron(const edm::Event& iEvent, const reco::GsfElectron& ele);

    TTree* tree_;
    const char *  name_;
    std::string Vtype_;
    std::string LeptonType_;
    edm::InputTag mInputBoson;
    edm::InputTag mInputMet;
    edm::InputTag mInputMetMVA;
	edm::InputTag mInputElectrons;
	edm::InputTag mInputBeamSpot;
    bool runoverAOD;	

  private:
    // private data members

    int nLooseElectron;
    int nTightElectron;
    
    float V_mass;
    float V_mt;
    float V_mtMVA;
    float V_px;
    float V_py;
    float V_pz;
    float V_E;
    float V_Pt;
    float V_Et;
    float V_Eta;
    float V_Phi;
    float V_Vx;
    float V_Vy;
    float V_Vz;
    float V_Y;
    float V_pzNu1;
    float V_pzNu2;

    int e1Classification; 
    int e1Charge;
    int e2Classification;
    int e2Charge;
    int e1_missingHits;
    int e2_missingHits;

    bool ise1WP95;
    bool ise1WP80;
    bool ise1WP70;
    bool ise2WP95;
    bool ise2WP80;
    bool ise2WP70;

    float e1px;
    float e1py;
    float e1pz;
    float e1E;
    float e1Et;
    float e1Pt;
    float e1Eta;
    float e1Theta;
    float e1Phi;
    float e1Vx;
    float e1Vy;
    float e1Vz;
    float e1Y;
    float e1_trackiso;
    float e1_hcaliso;
    float e1_ecaliso;
    float e1_sc_x;
    float e1_sc_y;
    float e1_sc_z;
    float e1_sc_Theta;
    float e1_sc_Eta;    
    float e1_sc_Phi;
    float e1_sc_E;
    float e1_sc_px;
    float e1_sc_py;
    float e1_sc_pz;
    float e1_sc_Pt;
    float e1_sc_Et;	  
    float e1_EoverPout;
    float e1_EoverPin;
    float e1_numberOfBrems;
    float e1_BremFraction;
    float e1_DeltaEtaIn;
    float e1_DeltaPhiIn;
    float e1_DeltaPhiOut;
    float e1_DeltaEtaOut;
    float e1_Trackmom_calo;
    float e1_Trackmom_vtx;	  
    float e1_HoverE;	  	  	  
    float e1_E9overE25;
    float e1_SigmaEtaEta;
    float e1_SigmaIetaIeta;
    float e1_escale;
    float e1_dist;
    float e1_dcot;
    float e1_convradius;
    float e1_d0bsp;
    float e1_dz000; 
    float e1_mvaTrigV0;
    float e1_mvaNonTrigV0;
  
    float e1_pfiso_chargedHadronIso;
    float e1_pfiso_photonIso;
    float e1_pfiso_neutralHadronIso;    
    float e1_pfIsoEA;    
    float e1_EffArea;    

    ///////////////////

    float e2px;
    float e2py;
    float e2pz;
    float e2E;
    float e2Pt;
    float e2Et;
    float e2Eta;
    float e2Theta;
    float e2Phi;
    float e2Vx;
    float e2Vy;
    float e2Vz;
    float e2Y;
    float e2_trackiso;
    float e2_hcaliso;
    float e2_ecaliso;
    float e2_sc_x;
    float e2_sc_y;
    float e2_sc_z;
    float e2_sc_Theta;
    float e2_sc_Eta;    
    float e2_sc_Phi;
    float e2_sc_E;
    float e2_sc_px;
    float e2_sc_py;
    float e2_sc_pz;
    float e2_sc_Pt;
    float e2_sc_Et;	  
    float e2_EoverPout;
    float e2_EoverPin;
    float e2_numberOfBrems;
    float e2_BremFraction;
    float e2_DeltaEtaIn;
    float e2_DeltaPhiIn;
    float e2_DeltaPhiOut;
    float e2_DeltaEtaOut;
    float e2_Trackmom_calo;
    float e2_Trackmom_vtx;	  
    float e2_HoverE;	  	  	 	  
    float e2_E9overE25;
    float e2_SigmaEtaEta;
    float e2_SigmaIetaIeta;
    float e2_escale;
    float e2_dist;
    float e2_dcot;
    float e2_convradius;
  };

} //namespace

#endif


