/*****************************************************************************
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
 *   To fill lepton, jet, MET information into a specified TTree
 * History:
 *   
 *
 * Copyright (C) 2013 FNAL 
 *****************************************************************************/

#ifndef ElectroWeakAnalysis_VPlusJets_SnowmassTreeProducer_h
#define ElectroWeakAnalysis_VPlusJets_SnowmassTreeProducer_h

#include <memory>
#include <string>
#include <iostream>
#include <map>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TTree.h" 

#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "EGamma/EGammaAnalysisTools/interface/PFIsolationEstimator.h"


namespace ewk {

  class SnowmassTreeProducer : public edm::EDAnalyzer {
  public:
    explicit SnowmassTreeProducer(const edm::ParameterSet iConfig );
    ~SnowmassTreeProducer();


    // To be called once per job
    virtual void beginJob();

    /// To be called once per event
    virtual void analyze(const edm::Event &iEvent, const edm::EventSetup& iSetup);

    // To be called once per job
    virtual void endJob() ;

  protected:

    /// Helper function for main constructor 
    void SetBranch( float* x, std::string brName );
    void SetBranch( int* x, std::string brName );
    void SetBranch( bool* x, std::string brName );
    void SetBranch( float x, std::string brName );
    float EAch( float x);
    float EAnh( float x);
    float EApho( float x);

  private:
    // private data members

    /// output ROOT file for the tree and histograms
    std::string fOutputFileName ;
    TFile*  hOutputFile ;

    /// output ROOT file for the tree and histograms
    TTree*  tree_;
    edm::InputTag mInputJets;
    PFIsolationEstimator isolator;
    typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;

    int run;
    int event; 
    int lumi; 
    int bunch; 
    int nPV; 
    float fastJetRho;
    
    int l1size;
    int l2size;
    int jet_size;
    int photon_size;
    int genPart_size;
    int genjet_size;
    int track_size;
    int gsftrack_size;
    int superCluster_size;
    int superCluster5x5_size;
    int caloTower_size;
    int tau_size;

    int l1Charge[999];
    int l2Charge[999];

    float l1px[999];
    float l1py[999];
    float l1pz[999];
    float l1E[999];
    float l1Pt[999];
    float l1Eta[999];
    float l1Phi[999];
    float l1Vx[999];
    float l1Vy[999];
    float l1Vz[999];
    float l1Y[999];
    float l1trackiso[999];
    float l1ecaliso[999];
    float l1hcaliso[999];
    int   l1Type[999];
    float l1_numberOfChambers[999];      
    float l1_numberOfMatches[999];
    float l1pfiso_sumChargedHadronPt[999];
    float l1pfiso_sumChargedParticlePt[999];
    float l1pfiso_sumNeutralHadronEt[999];
    float l1pfiso_sumPhotonEt[999];
    float l1pfiso_sumPUPt[999];
    float l1_d0bsp[999];
    float l1_dz000[999];
    float l1_IP3D[999];
    float l1_dzPV[999];
    float l1_globalChi2[999];
    float l1_innerChi2[999];
    float l1_nPixelHits[999];
    float l1_nTrackerHits[999];
    int l1_isPF[999];
    int l1_isGlobal[999];
    int l1_isTracker[999];
    int l1_hasMuonSegment[999];

    ///////////////////
    float l2px[999];
    float l2py[999];
    float l2pz[999];
    float l2E[999];
    float l2Pt[999];
    float l2Eta[999];
    float l2Phi[999];
    float l2Vx[999];
    float l2Vy[999];
    float l2Vz[999];
    float l2Y[999];  
    float l2trackiso[999];
    float l2ecaliso[999];
    float l2hcaliso[999];
    int   l2_classification[999];
    float l2_HoverE[999]; 
    float l2_EoverP[999];
    float l2_DeltaEta[999];
    float l2_DeltaPhi[999];
    int   l2_numberOfBrems[999];      
    float l2_BremFraction[999];
    float l2_SigmaIetaIeta[999];
    int   l2_missingHits[999];
    float l2_dist[999];
    float l2_dcot[999];
    float l2_convradius[999];
    float l2pfiso_chargedHadronIso[999];
    float l2pfiso_photonIso[999];
    float l2pfiso_neutralHadronIso[999];
    float l2pfiso_EffAreaPU[999];
    float l2pfiso_pfIsoEA[999];
    float l2_d0bsp[999];
    float l2_dz000[999];
    float l2_IP3D[999];
    float l2_dzPV[999];


    ///////////////////
    float jet_px[999];
    float jet_py[999];
    float jet_pz[999];
    float jet_E[999];
    float jet_Pt[999];
    float jet_Eta[999];
    float jet_Phi[999];
    float jet_Y[999];
    float jet_area[999];
    float jet_bDiscriminatorSSVHE[999];
    float jet_bDiscriminatorTCHE[999];
    float jet_bDiscriminatorCSV[999];
    float jet_bDiscriminatorJP[999];
    float jet_bDiscriminatorSSVHP[999];
    float jet_bDiscriminatorTCHP[999];
 
    ///////////////////
    float Met_px;
    float Met_py;
    float Met_Et;
    float Met_Phi;
    float Met_SumET;

    ///////////////////
    float photon_px[999];
    float photon_py[999];
    float photon_pz[999];
    float photon_E[999];
    float photon_pt[999];
    float photon_eta[999];
    float photon_phi[999];
    float photon_vx[999];
    float photon_vy[999];
    float photon_vz[999];
    float photon_pfiso_charged[999];
    float photon_pfiso_photon[999];
    float photon_pfiso_neutral[999];
    float photon_trackiso[999];
    float photon_ecaliso[999];
    float photon_hcaliso[999];
    float photon_HoverE[999];
    float photon_SigmaIetaIeta[999];
    int photon_hasPixelSeed[999];
    int photon_passElecVeto[999];

    ///////////////////
    float track_px[999];
    float track_py[999];
    float track_pz[999];
    float track_Vx[999];
    float track_Vy[999];
    float track_Vz[999];
    float track_Pt[999];
    float track_Eta[999];
    float track_Phi[999];

    ///////////////////
    float gsftrack_px[999];
    float gsftrack_py[999];
    float gsftrack_pz[999];
    float gsftrack_Vx[999];
    float gsftrack_Vy[999];
    float gsftrack_Vz[999];
    float gsftrack_Pt[999];
    float gsftrack_Eta[999];
    float gsftrack_Phi[999];

    ///////////////////
    float superCluster_E    [999];
    float superCluster_rawE [999];
    float superCluster_x    [999];
    float superCluster_y    [999];
    float superCluster_z    [999];
    float superCluster_Eta  [999];
    float superCluster_Phi  [999];
    float superCluster_nHits[999];

    ///////////////////
    float superCluster5x5_E    [999];
    float superCluster5x5_rawE [999];
    float superCluster5x5_x    [999];
    float superCluster5x5_y    [999];
    float superCluster5x5_z    [999];
    float superCluster5x5_Eta  [999];
    float superCluster5x5_Phi  [999];
    float superCluster5x5_nHits[999];

    ///////////////////
    float caloTower_hadE [999];
    float caloTower_emE  [999];
    float caloTower_hadEt[999];
    float caloTower_emEt [999];
    float caloTower_Eta  [999];
    float caloTower_Phi  [999];

    ///////////////////
    float tau_px[999];
    float tau_py[999];
    float tau_pz[999];
    float tau_Pt[999];
    float tau_Eta[999];
    float tau_Phi[999];

    ///////////////////
    int genPart_Charge[9999];
    float genPart_px[9999];
    float genPart_py[9999];
    float genPart_pz[9999];
    float genPart_E[9999];
    float genPart_Pt[9999];
    float genPart_Eta[9999];
    float genPart_Phi[9999];
    float genPart_Vx[9999];
    float genPart_Vy[9999];
    float genPart_Vz[9999];
    float genPart_Y[9999];

    ///////////////////
    float genjet_px[999];
    float genjet_py[999];
    float genjet_pz[999];
    float genjet_E[999];
    float genjet_Pt[999];
    float genjet_Eta[999];
    float genjet_Phi[999];
    float genjet_Y[999];
    float genjet_area[999];
    int genPart_Status[999];
    int genPart_pdgId[999];

  };

} //namespace

#endif

