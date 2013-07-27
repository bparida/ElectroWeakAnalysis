/*****************************************************************************
 * Project: CMS detector at the CERN
 *
 * Package: ElectroWeakAnalysis/VPlusJets
 * Class:   VtoMuonTreeFiller
 *
 * Authors:
 *
 *   Kalanand Mishra, Fermilab - kalanand@fnal.gov
*   Anil Singh, Punjab University - anil79@fnal.gov
 *
 * Description:
 *   To fill W--> munu or Z-->mumu related quantities into a specified TTree
 *   Can work with jets in RECO/AOD/PAT data formats.
 * History:
 *   
 *
 * Copyright (C) 2010 FNAL 
 *****************************************************************************/

#ifndef ElectroWeakAnalysis_VPlusJets_VtoMuonTreeFiller_h
#define ElectroWeakAnalysis_VPlusJets_VtoMuonTreeFiller_h

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h" 
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"


namespace ewk {

  class VtoMuonTreeFiller {
  public:
    /// specify the name of the TTree, and the configuration for it
    VtoMuonTreeFiller(const char *name, TTree* tree, 
			  const edm::ParameterSet iConfig );

    /// default constructor
    VtoMuonTreeFiller() {};


    /// Destructor, does nothing 
      ~VtoMuonTreeFiller() {};


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

    TTree* tree_;
    const char *  name_;
    std::string Vtype_;
    std::string LeptonType_;
    edm::InputTag mInputBoson;
    edm::InputTag mInputMet;
    edm::InputTag mInputMetMVA;
	edm::InputTag mInputMuons;
	edm::InputTag mInputBeamSpot;
	bool runoverAOD;
  private:
    // private data members
    
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

    /// muon type - type of the algorithm that reconstructed this muon
    /// multiple algorithms can reconstruct the same muon
    /// GlobalMuon = 1; TrackerMuon = 2; StandAloneMuon = 3; CaloMuon = 4;
    int mu1Type; 
    int mu1Charge;
    /// number of chambers
    int mu1_numberOfChambers;
   /// get number of chambers with matched segments
    int mu1_numberOfMatches;

    int mu2Type;
    int mu2Charge;
    int mu2_numberOfChambers;
    int mu2_numberOfMatches;

    float mu1px;
    float mu1py;
    float mu1pz;
    float mu1E;
    float mu1Et;
    float mu1Pt;
    float mu1Eta;
    float mu1Theta;
    float mu1Phi;
    float mu1Vx;
    float mu1Vy;
    float mu1Vz;
    float mu1Y;
    float mu1_trackiso;
    float mu1_hcaliso;
    float mu1_ecaliso;
    
    float mu1_d0bsp;
    float mu1_dz000;
    float mu1_dzPV;
    float mu1_IP3D;
  
    float mu1_pfiso_chargedHadronIso;
    float mu1_pfiso_photonIso;
    float mu1_pfiso_neutralHadronIso;    

    float mu1pfiso_sumChargedHadronPt;
    float mu1pfiso_sumChargedParticlePt;
    float mu1pfiso_sumNeutralHadronEt;
    float mu1pfiso_sumPhotonEt;
    float mu1pfiso_sumPUPt;

    ///////////////////
    float mu2px;
    float mu2py;
    float mu2pz;
    float mu2E;
    float mu2Pt;
    float mu2Et;
    float mu2Eta;
    float mu2Theta;
    float mu2Phi;
    float mu2Vx;
    float mu2Vy;
    float mu2Vz;
    float mu2Y;
    float mu2_trackiso;
    float mu2_hcaliso;
    float mu2_ecaliso;
   
    float mu2_d0bsp;
    float mu2_dz000;
    float mu2_dzPV;
    float mu2_IP3D;
  
    float mu2_pfiso_chargedHadronIso;
    float mu2_pfiso_photonIso;
    float mu2_pfiso_neutralHadronIso;    

    float mu2pfiso_sumChargedHadronPt;
    float mu2pfiso_sumChargedParticlePt;
    float mu2pfiso_sumNeutralHadronEt;
    float mu2pfiso_sumPhotonEt;
    float mu2pfiso_sumPUPt;


  };

} //namespace

#endif


