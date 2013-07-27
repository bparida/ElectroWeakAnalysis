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

// CMS includes
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TMath.h" 
// Header file
#include "ElectroWeakAnalysis/VPlusJets/interface/SnowmassTreeProducer.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidatePhotonExtraFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidatePhotonExtra.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/TauReco/interface/PFTau.h"

ewk::SnowmassTreeProducer::SnowmassTreeProducer(const edm::ParameterSet iConfig) :

  fOutputFileName ( iConfig.getParameter<std::string>("HistOutFile") ),
  hOutputFile ( new TFile( fOutputFileName.c_str(), "RECREATE" ) ), 
  tree_ ( new TTree("tree","tree") )
{
  mInputJets = iConfig.getParameter<edm::InputTag>("srcJets");
}


ewk::SnowmassTreeProducer::~SnowmassTreeProducer() {}



void ewk::SnowmassTreeProducer::beginJob()
{

  // Declare all the branches of the tree
  tree_->Branch("event_runNo",  &run,   "event_runNo/I");
  tree_->Branch("event_evtNo",  &event, "event_evtNo/I");
  tree_->Branch("event_lumi",   &lumi,  "event_lumi/I"); 
  tree_->Branch("event_bunch",  &bunch, "event_bunch/I"); 
  tree_->Branch("event_nPV",    &nPV,   "event_nPV/I"); 
  tree_->Branch("event_met",    &Met_Et,  "event_met/F"); 
  tree_->Branch("event_met_px", &Met_px,  "event_met_px/F"); 
  tree_->Branch("event_met_py", &Met_py,  "event_met_py/F"); 
  tree_->Branch("event_sumet",  &Met_SumET,"event_sumet/F"); 
tree_->Branch("event_met_phi",&Met_Phi,  "event_met_phi/F"); 
  tree_->Branch("event_fastJetRho",&fastJetRho,  "event_fastJetRho/F"); 

  // ----------------------- Declare branches -----------------------
  ///////////////////////////////////////////////
  SetBranch( &l1size,           "muon_size" );
  SetBranch( l1px,             "muon_px[muon_size]" );
  SetBranch( l1py,             "muon_py[muon_size]" );
  SetBranch( l1pz,             "muon_pz[muon_size]" );
  SetBranch( l1E,              "muon_e[muon_size]" );
  SetBranch( l1Pt,             "muon_pt[muon_size]" );
  SetBranch( l1Eta,            "muon_eta[muon_size]" ); 
  SetBranch( l1Phi,            "muon_phi[muon_size]" );
  SetBranch( l1Charge,         "muon_charge[muon_size]" );
  SetBranch( l1Vx,             "muon_vx[muon_size]" );
  SetBranch( l1Vy,             "muon_vy[muon_size]" );
  SetBranch( l1Vz,             "muon_vz[muon_size]" );
  SetBranch( l1Y,              "muon_y[muon_size]" );
  SetBranch( l1trackiso,              "muon_trackiso[muon_size]" );
  SetBranch( l1ecaliso,              "muon_ecaliso[muon_size]" );
  SetBranch( l1hcaliso,              "muon_hcaliso[muon_size]" );
  SetBranch( l1Type,              "muon_Type[muon_size]" );
  SetBranch( l1_numberOfChambers,              "muon_numberOfChambers[muon_size]" );
  SetBranch( l1_numberOfMatches,              "muon_numberOfMatches[muon_size]" );
  SetBranch( l1pfiso_sumChargedHadronPt,              "muon_pfiso_sumChargedHadronPt[muon_size]" );
  SetBranch( l1pfiso_sumChargedParticlePt,              "muon_pfiso_sumChargedParticlePt[muon_size]" );
  SetBranch( l1pfiso_sumNeutralHadronEt,              "muon_pfiso_sumNeutralHadronEt[muon_size]" );
  SetBranch( l1pfiso_sumPhotonEt,              "muon_pfiso_sumPhotonEt[muon_size]" );
  SetBranch( l1pfiso_sumPUPt,              "muon_pfiso_sumPUPt[muon_size]" );
  SetBranch( l1_d0bsp,              "muon_d0bsp[muon_size]" );
  SetBranch( l1_dz000,              "muon_dz000[muon_size]" );
  SetBranch( l1_IP3D,              "muon_IP3D[muon_size]" );
  SetBranch( l1_dzPV,              "muon_dzPV[muon_size]" );
  SetBranch( l1_globalChi2,          "muon_globalChi2[muon_size]");
  SetBranch( l1_innerChi2,          "muon_innerChi2[muon_size]");
  SetBranch( l1_nPixelHits,        "muon_nPixelHits[muon_size]");
  SetBranch( l1_nTrackerHits,      "muon_nTrackerHits[muon_size]");
  SetBranch( l1_isPF,        "muon_isPF[muon_size]");
  SetBranch( l1_isGlobal,      "muon_isGlobal[muon_size]");
  SetBranch( l1_isTracker,      "muon_isTracker[muon_size]");
  SetBranch( l1_hasMuonSegment, "muon_hasMuonSegment[muon_size]");

  ////////////////////////////////////////////////////////
  SetBranch( &l2size,           "electron_size" );
  SetBranch( l2px,             "electron_px[electron_size]" );
  SetBranch( l2py,             "electron_py[electron_size]" );
  SetBranch( l2pz,             "electron_pz[electron_size]" );
  SetBranch( l2E,              "electron_e[electron_size]" );
  SetBranch( l2Pt,             "electron_pt[electron_size]" );
  SetBranch( l2Eta,            "electron_eta[electron_size]" ); 
  SetBranch( l2Phi,            "electron_phi[electron_size]" );
  SetBranch( l2Charge,         "electron_charge[electron_size]" );
  SetBranch( l2Vx,             "electron_vx[electron_size]" );
  SetBranch( l2Vy,             "electron_vy[electron_size]" );
  SetBranch( l2Vz,             "electron_vz[electron_size]" );
  SetBranch( l2Y,              "electron_y[electron_size]" );
  SetBranch( l2trackiso,              "electron_trackiso[electron_size]" );
  SetBranch( l2ecaliso,              "electron_ecaliso[electron_size]" );
  SetBranch( l2hcaliso,              "electron_hcaliso[electron_size]" );
  SetBranch( l2_classification,              "electron_classification[electron_size]" );
  SetBranch( l2_HoverE,              "electron_HoverE[electron_size]" );
  SetBranch( l2_EoverP,              "electron_EoverP[electron_size]" );
  SetBranch( l2_DeltaEta,              "electron_DeltaEta[electron_size]" );
  SetBranch( l2_DeltaPhi,              "electron_DeltaPhi[electron_size]" );
  SetBranch( l2_numberOfBrems,              "electron_numberOfBrems[electron_size]" );
  SetBranch( l2_BremFraction,              "electron_BremFraction[electron_size]" );
  SetBranch( l2_SigmaIetaIeta,              "electron_SigmaIetaIeta[electron_size]" );
  SetBranch( l2_missingHits,              "electron_missingHits[electron_size]" );
  SetBranch( l2_dist,              "electron_convDist[electron_size]" );
  SetBranch( l2_dcot,              "electron_convDcot[electron_size]" );
  SetBranch( l2_convradius,              "electron_convRadius[electron_size]" );
  SetBranch( l2pfiso_chargedHadronIso,              "electron_pfiso_chargedHadron[electron_size]" );
  SetBranch( l2pfiso_photonIso,              "electron_pfiso_photon[electron_size]" );
  SetBranch( l2pfiso_neutralHadronIso,              "electron_pfiso_neutralHadron[electron_size]" );
  SetBranch( l2pfiso_EffAreaPU,              "electron_pfiso_EffAreaPU[electron_size]" );
  SetBranch( l2pfiso_pfIsoEA,              "electron_pfiso_pfIsoEA[electron_size]" );
  SetBranch( l2_d0bsp,              "electron_d0bsp[electron_size]" );
  SetBranch( l2_dz000,              "electron_dz000[electron_size]" );
  SetBranch( l2_IP3D,              "electron_IP3D[electron_size]" );
  SetBranch( l2_dzPV,              "electron_dzPV[electron_size]" );

  ////////////////////////////////////////////////////////
  SetBranch( &jet_size,          "jet_size" );
  SetBranch( jet_px,             "jet_px[jet_size]" );
  SetBranch( jet_py,             "jet_py[jet_size]" );
  SetBranch( jet_pz,             "jet_pz[jet_size]" );
  SetBranch( jet_E,              "jet_E[jet_size]" );
  SetBranch( jet_Pt,             "jet_pt[jet_size]" );
  SetBranch( jet_Eta,            "jet_eta[jet_size]" ); 
  SetBranch( jet_Phi,            "jet_phi[jet_size]" );
  SetBranch( jet_Y,              "jet_y[jet_size]" );
  SetBranch( jet_area,           "jet_area[jet_size]" );
  SetBranch( jet_bDiscriminatorSSVHE,           "jet_bDiscriminatorSSVHE[jet_size]" );
  SetBranch( jet_bDiscriminatorTCHE,           "jet_bDiscriminatorTCHE[jet_size]" );
  SetBranch( jet_bDiscriminatorCSV,           "jet_bDiscriminatorCSV[jet_size]" );
  SetBranch( jet_bDiscriminatorJP,           "jet_bDiscriminatorJP[jet_size]" );
  SetBranch( jet_bDiscriminatorSSVHP,           "jet_bDiscriminatorSSVHP[jet_size]" );
  SetBranch( jet_bDiscriminatorTCHP,           "jet_bDiscriminatorTCHP[jet_size]" );

  ////////////////////////////////////////////////////////
  SetBranch( &photon_size,          "photon_size" );
  SetBranch( photon_px, "photon_px[photon_size]");
  SetBranch( photon_py, "photon_py[photon_size]");
  SetBranch( photon_pz, "photon_pz[photon_size]");
  SetBranch( photon_E, "photon_E[photon_size]");
  SetBranch( photon_pt, "photon_pt[photon_size]");
  SetBranch( photon_eta, "photon_eta[photon_size]");
  SetBranch( photon_phi, "photon_phi[photon_size]");
  SetBranch( photon_vx, "photon_vx[photon_size]");
  SetBranch( photon_vy, "photon_vy[photon_size]");
  SetBranch( photon_vz, "photon_vz[photon_size]");
//   SetBranch( photon_pfiso_charged, "photon_pfiso_charged[photon_size]");
//   SetBranch( photon_pfiso_photon, "photon_pfiso_photon[photon_size]");
//   SetBranch( photon_pfiso_neutral, "photon_pfiso_neutral[photon_size]");
  SetBranch( photon_trackiso, "photon_trackiso[photon_size]");
  SetBranch( photon_ecaliso, "photon_ecaliso[photon_size]");
  SetBranch( photon_hcaliso, "photon_hcaliso[photon_size]");
  SetBranch( photon_HoverE, "photon_HoverE[photon_size]");
  SetBranch( photon_SigmaIetaIeta, "photon_SigmaIetaIeta[photon_size]");
  SetBranch( photon_hasPixelSeed, "photon_hasPixelSeed[photon_size]");
  SetBranch( photon_passElecVeto, "photon_passElecVeto[photon_size]");

  ////////////////////////////////////////////////////////  
  SetBranch( &track_size,            "track_size" );
  SetBranch( track_px,               "track_px[track_size]" );
  SetBranch( track_py,               "track_py[track_size]" );
  SetBranch( track_pz,               "track_pz[track_size]" );
  SetBranch( track_Vx,               "track_Vx[track_size]" );
  SetBranch( track_Vy,               "track_Vy[track_size]" );
  SetBranch( track_Vz,               "track_Vz[track_size]" );
  SetBranch( track_Pt,               "track_Pt[track_size]" );
  SetBranch( track_Eta,              "track_Eta[track_size]" );
  SetBranch( track_Phi,              "track_Phi[track_size]" );
  
  ////////////////////////////////////////////////////////
  SetBranch( &gsftrack_size,            "gsftrack_size" );
  SetBranch( gsftrack_px,               "gsftrack_px[gsftrack_size]" );
  SetBranch( gsftrack_py,               "gsftrack_py[gsftrack_size]" );
  SetBranch( gsftrack_pz,               "gsftrack_pz[gsftrack_size]" );
  SetBranch( gsftrack_Vx,               "gsftrack_Vx[gsftrack_size]" );
  SetBranch( gsftrack_Vy,               "gsftrack_Vy[gsftrack_size]" );
  SetBranch( gsftrack_Vz,               "gsftrack_Vz[gsftrack_size]" );
  SetBranch( gsftrack_Pt,               "gsftrack_Pt[gsftrack_size]" );
  SetBranch( gsftrack_Eta,              "gsftrack_Eta[gsftrack_size]" );
  SetBranch( gsftrack_Phi,              "gsftrack_Phi[gsftrack_size]" );

  ////////////////////////////////////////////////////////  
  SetBranch( &superCluster_size,              "superCluster_size"  );
  SetBranch( superCluster_E    ,              "superCluster_E[superCluster_size]"     );
  SetBranch( superCluster_rawE ,              "superCluster_rawE[superCluster_size]"  );
  SetBranch( superCluster_x    ,              "superCluster_x[superCluster_size]"     );
  SetBranch( superCluster_y    ,              "superCluster_y[superCluster_size]"     );
  SetBranch( superCluster_z    ,              "superCluster_z[superCluster_size]"     );
  SetBranch( superCluster_Eta  ,              "superCluster_Eta[superCluster_size]"   );
  SetBranch( superCluster_Phi  ,              "superCluster_Phi[superCluster_size]"   );
  SetBranch( superCluster_nHits,              "superCluster_nHits[superCluster_size]" );

  ////////////////////////////////////////////////////////  
  SetBranch( &superCluster5x5_size,              "superCluster5x5_size"  );
  SetBranch( superCluster5x5_E    ,              "superCluster5x5_E[superCluster5x5_size]"     );
  SetBranch( superCluster5x5_rawE ,              "superCluster5x5_rawE[superCluster5x5_size]"  );
  SetBranch( superCluster5x5_x    ,              "superCluster5x5_x[superCluster5x5_size]"     );
  SetBranch( superCluster5x5_y    ,              "superCluster5x5_y[superCluster5x5_size]"     );
  SetBranch( superCluster5x5_z    ,              "superCluster5x5_z[superCluster5x5_size]"     );
  SetBranch( superCluster5x5_Eta  ,              "superCluster5x5_Eta[superCluster5x5_size]"   );
  SetBranch( superCluster5x5_Phi  ,              "superCluster5x5_Phi[superCluster5x5_size]"   );
  SetBranch( superCluster5x5_nHits,              "superCluster5x5_nHits[superCluster5x5_size]" );

  ////////////////////////////////////////////////////////  
  SetBranch( &caloTower_size     ,               "caloTower_size"                 );
  SetBranch( caloTower_hadE      ,               "caloTower_hadE[caloTower_size]" );
  SetBranch( caloTower_emE       ,               "caloTower_emE[caloTower_size]"  );
  SetBranch( caloTower_hadEt     ,               "caloTower_hadEt[caloTower_size]");
  SetBranch( caloTower_emEt      ,               "caloTower_emEt[caloTower_size]" );
  SetBranch( caloTower_Eta       ,               "caloTower_Eta[caloTower_size]"  );
  SetBranch( caloTower_Phi       ,               "caloTower_Phi[caloTower_size]"  );

  ////////////////////////////////////////////////////////  
  SetBranch( &tau_size,            "tau_size" );
  SetBranch( tau_px,               "tau_px[tau_size]" );
  SetBranch( tau_py,               "tau_py[tau_size]" );
  SetBranch( tau_pz,               "tau_pz[tau_size]" );
  SetBranch( tau_Pt,               "tau_Pt[tau_size]" );
  SetBranch( tau_Eta,              "tau_Eta[tau_size]" );
  SetBranch( tau_Phi,              "tau_Phi[tau_size]" );

  ///////////////////////////////////////////////
  SetBranch( &genPart_size,          "genParticle_size" );
  SetBranch( genPart_px,             "genParticle_px[genParticle_size]" );
  SetBranch( genPart_py,             "genParticle_py[genParticle_size]" );
  SetBranch( genPart_pz,             "genParticle_pz[genParticle_size]" );
  SetBranch( genPart_E,              "genParticle_e[genParticle_size]" );
  SetBranch( genPart_Pt,             "genParticle_pt[genParticle_size]" );
  SetBranch( genPart_Eta,            "genParticle_eta[genParticle_size]" ); 
  SetBranch( genPart_Phi,            "genParticle_phi[genParticle_size]" );
  SetBranch( genPart_Charge,         "genParticle_charge[genParticle_size]" );
  SetBranch( genPart_Vx,             "genParticle_vx[genParticle_size]" );
  SetBranch( genPart_Vy,             "genParticle_vy[genParticle_size]" );
  SetBranch( genPart_Vz,             "genParticle_vz[genParticle_size]" );
  SetBranch( genPart_Y,              "genParticle_y[genParticle_size]" );
  SetBranch( genPart_Status,         "genPart_Status[genParticle_size]" );
  SetBranch( genPart_pdgId,          "genPart_pdgId[genParticle_size]" );

  ////////////////////////////////////////////////////////
  SetBranch( &genjet_size,          "genjet_size" );
  SetBranch( genjet_px,             "genjet_px[genjet_size]" );
  SetBranch( genjet_py,             "genjet_py[genjet_size]" );
  SetBranch( genjet_pz,             "genjet_pz[genjet_size]" );
  SetBranch( genjet_E,              "genjet_E[genjet_size]" );
  SetBranch( genjet_Pt,             "genjet_pt[genjet_size]" );
  SetBranch( genjet_Eta,            "genjet_eta[genjet_size]" ); 
  SetBranch( genjet_Phi,            "genjet_phi[genjet_size]" );
  SetBranch( genjet_Y,              "genjet_y[genjet_size]" );
  SetBranch( genjet_area,           "genjet_area[genjet_size]" );

}
/////////////////////////////////////////////////////////////////////////







void ewk::SnowmassTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // first initialize to the default values
  genPart_size       = 0;
  l1size             = 0; 
  l2size             = 0;
  jet_size           = 0;
  photon_size        = 0;
  genjet_size        = 0;
  track_size         = 0;
  gsftrack_size      = 0;
  superCluster_size  = 0;
  superCluster5x5_size = 0;
  caloTower_size     = 0;
  tau_size         = 0;

  Met_px              = -99999.;
  Met_py              = -99999.;
  Met_Et              = -99999.;
  Met_Phi             = -99999.;
  Met_SumET           = -99999.;

  for ( int i =0; i<999; i++){

    l1Charge[i]           = -99999;
    l1px[i]               = -99999.;
    l1py[i]               = -99999.;
    l1pz[i]               = -99999.;
    l1E[i]                = -99999.;
    l1Pt[i]               = -99999.;
    l1Eta[i]              = -99999.;
    l1Phi[i]              = -99999.;
    l1Vx[i]               = -99999.;
    l1Vy[i]               = -99999.;
    l1Vz[i]               = -99999.;
    l1Y[i]                = -99999.;
    l1trackiso[i]= -99999.;
    l1ecaliso[i]= -99999.;
    l1hcaliso[i]= -99999.;
    l1Type[i]= -99999;
    l1_numberOfChambers[i]= -99999.;      
    l1_numberOfMatches[i]= -99999.;
    l1pfiso_sumChargedHadronPt[i]= -99999.;
    l1pfiso_sumChargedParticlePt[i]= -99999.;
    l1pfiso_sumNeutralHadronEt[i]= -99999.;
    l1pfiso_sumPhotonEt[i]= -99999.;
    l1pfiso_sumPUPt[i]= -99999.;
    l1_d0bsp[i]= -99999.;
    l1_dz000[i]= -99999.;
    l1_IP3D[i]= -99999.;
    l1_dzPV[i]= -99999.;
    l1_globalChi2[i]=-99999.;
    l1_innerChi2[i]=-99999.;
    l1_nPixelHits[i]=-99999.;
    l1_nTrackerHits[i]=-99999.;
    l1_isPF[i]=-99999;
    l1_isGlobal[i]=-99999;
    l1_isTracker[i]=-99999;
    l1_hasMuonSegment[i]=-99999;

    l2Charge[i]           = -99999;  
    l2px[i]               = -99999.;
    l2py[i]               = -99999.;
    l2pz[i]               = -99999.;
    l2E[i]                = -99999.;
    l2Pt[i]               = -99999.;
    l2Eta[i]              = -99999.;
    l2Phi[i]              = -99999.;
    l2Vx[i]               = -99999.;
    l2Vy[i]               = -99999.;
    l2Vz[i]               = -99999.;
    l2Y[i]                = -99999.;
    l2trackiso[i]       = -99999.;
    l2ecaliso[i]        = -99999.;
    l2hcaliso[i]        = -99999.;
    l2_classification[i]  = -99999;
    l2_HoverE[i] = -99999.; 
    l2_EoverP[i]       = -99999.;
    l2_DeltaEta[i]     = -99999.;
    l2_DeltaPhi[i]     = -99999.;
    l2_numberOfBrems[i]  = -99999;      
    l2_BremFraction[i]   = -99999.;
    l2_SigmaIetaIeta[i] = -99999.;
    l2_missingHits[i]   = -99999;
    l2_dist[i]          = -99999.;
    l2_dcot[i]          = -99999.;
    l2_convradius[i]    = -99999.;
    l2pfiso_chargedHadronIso[i]   = -99999.;
    l2pfiso_photonIso[i] = -99999.;
    l2pfiso_neutralHadronIso[i]   = -99999.;
    l2pfiso_EffAreaPU[i] = -99999.;
    l2pfiso_pfIsoEA[i] = -99999.;
    l2_d0bsp[i] = -99999.;
    l2_dz000[i] = -99999.;
    l2_IP3D[i] = -99999.;
    l2_dzPV[i] = -99999.;


    jet_px[i]             = -99999.;
    jet_py[i]             = -99999.;
    jet_pz[i]             = -99999.;
    jet_E[i]              = -99999.;
    jet_Pt[i]             = -99999.;
    jet_Eta[i]            = -99999.;
    jet_Phi[i]            = -99999.;
    jet_Y[i]              = -99999.;
    jet_area[i]           = -99999.;
    jet_bDiscriminatorSSVHE[i] = -99999.;
    jet_bDiscriminatorTCHE[i]  = -99999.;
    jet_bDiscriminatorCSV[i]   = -99999.;
    jet_bDiscriminatorJP[i]    = -99999.;
    jet_bDiscriminatorSSVHP[i] = -99999.;
    jet_bDiscriminatorTCHP[i]  = -99999.;
    
    photon_px[i]  = -99999.;
    photon_py[i]  = -99999.;
    photon_pz[i]  = -99999.;
    photon_E[i]  = -99999.;
    photon_pt[i]  = -99999.;
    photon_eta[i]  = -99999.;
    photon_phi[i]  = -99999.;
    photon_vx[i]  = -99999.;
    photon_vy[i]  = -99999.;
    photon_vz[i]  = -99999.;
    photon_pfiso_charged[i]  = -99999.;
    photon_pfiso_photon[i]  = -99999.;
    photon_pfiso_neutral[i]  = -99999.;
    photon_trackiso[i]  = -99999.;
    photon_ecaliso[i]  = -99999.;
    photon_hcaliso[i]  = -99999.;
    photon_HoverE[i]  = -99999.;
    photon_SigmaIetaIeta[i]  = -99999.;
    photon_hasPixelSeed[i]  = 0;
    photon_passElecVeto[i]  = 0;

    track_px[i]                 = -99999.;
    track_py[i]                 = -99999.;
    track_pz[i]                 = -99999.;
    track_Vx[i]                 = -99999.;
    track_Vy[i]                 = -99999.;
    track_Vz[i]                 = -99999.;
    track_Pt[i]                 = -99999.;
    track_Eta[i]                = -99999.;
    track_Phi[i]                = -99999.;

    gsftrack_px[i]                 = -99999.;
    gsftrack_py[i]                 = -99999.;
    gsftrack_pz[i]                 = -99999.;
    gsftrack_Vx[i]                 = -99999.;
    gsftrack_Vy[i]                 = -99999.;
    gsftrack_Vz[i]                 = -99999.;
    gsftrack_Pt[i]                 = -99999.;
    gsftrack_Eta[i]                = -99999.;
    gsftrack_Phi[i]                = -99999.;

    superCluster_E[i]      = -99999.;
    superCluster_rawE[i]   = -99999.;
    superCluster_x[i]      = -99999.;
    superCluster_y[i]      = -99999.;
    superCluster_z[i]      = -99999.;
    superCluster_Eta[i]    = -99999.;
    superCluster_Phi[i]    = -99999.;
    superCluster_nHits[i]  = -99999.;

    superCluster5x5_E[i]      = -99999.;
    superCluster5x5_rawE[i]   = -99999.;
    superCluster5x5_x[i]      = -99999.;
    superCluster5x5_y[i]      = -99999.;
    superCluster5x5_z[i]      = -99999.;
    superCluster5x5_Eta[i]    = -99999.;
    superCluster5x5_Phi[i]    = -99999.;
    superCluster5x5_nHits[i]  = -99999.;

    caloTower_hadE[i]      = -99999.;
    caloTower_emE[i]       = -99999.;
    caloTower_hadEt[i]     = -99999.;
    caloTower_emEt[i]      = -99999.;
    caloTower_Eta[i]       = -99999.;
    caloTower_Phi[i]       = -99999.;

    tau_px[i]                 = -99999.;
    tau_py[i]                 = -99999.;
    tau_pz[i]                 = -99999.;
    tau_Pt[i]                 = -99999.;
    tau_Eta[i]                = -99999.;
    tau_Phi[i]                = -99999.;

    genjet_px[i]             = -99999.;
    genjet_py[i]             = -99999.;
    genjet_pz[i]             = -99999.;
    genjet_E[i]              = -99999.;
    genjet_Pt[i]             = -99999.;
    genjet_Eta[i]            = -99999.;
    genjet_Phi[i]            = -99999.;
    genjet_Y[i]              = -99999.;
    genjet_area[i]           = -99999.;

    genPart_Charge[i]           = -99999;
    genPart_px[i]               = -99999.;
    genPart_py[i]               = -99999.;
    genPart_pz[i]               = -99999.;
    genPart_E[i]                = -99999.;
    genPart_Pt[i]               = -99999.;
    genPart_Eta[i]              = -99999.;
    genPart_Phi[i]              = -99999.;
    genPart_Vx[i]               = -99999.;
    genPart_Vy[i]               = -99999.;
    genPart_Vz[i]               = -99999.;
    genPart_Y[i]                = -99999.;  
    genPart_Status[i]           = -99999;
    genPart_pdgId[i]            = -99999;
  }
  // initialization done



  // write event information: run, event, bunch crossing, ....
  run   = iEvent.id().run();
  event = iEvent.id().event();
  lumi  = iEvent.luminosityBlock();
  bunch = iEvent.bunchCrossing();



  // primary/secondary vertices
  // edm::Handle<reco::VertexCollection > recVtxs;
  edm::Handle <edm::View<reco::Vertex> > recVtxs;
  iEvent.getByLabel( "goodOfflinePrimaryVertices", recVtxs);
  nPV = recVtxs->size();




  /////// --- Pileup density "rho" in the event from fastJet pileup calculation --------
  edm::Handle<double> rho;
  const edm::InputTag eventrho("kt6PFJetsPFlow", "rho");
  iEvent.getByLabel(eventrho,rho);
  fastJetRho = *rho;



  //------------------ Met filling ------------------
  edm::Handle<edm::View<reco::MET> > pfmet;
  iEvent.getByLabel("patType1CorrectedPFMet", pfmet);
  if (pfmet->size() > 0) {
    Met_Phi             = (*pfmet)[0].phi();
    Met_px              = (*pfmet)[0].px();
    Met_py              = (*pfmet)[0].py();
    Met_Et              = (*pfmet)[0].et();
    Met_SumET           = (*pfmet)[0].sumEt(); 
  }




  //------------------ Muon filling ------------------
  typedef edm::View<reco::Muon> MuonView;
  edm::Handle<MuonView> muons;
  iEvent.getByLabel( "selectedPatMuonsPFlow", muons);
  l1size  = (int) (muons->size()); 
  if( l1size > 0 ) {
    int iMuon = 0;
    edm::View<reco::Muon>::const_iterator muon, endpmuons = muons->end(); 
    for (muon = muons->begin();  muon != endpmuons;  ++muon, ++iMuon) {

      l1Charge[iMuon]           = (*muon).charge();
      l1px[iMuon]               = (*muon).px();
      l1py[iMuon]               = (*muon).py();
      l1pz[iMuon]               = (*muon).pz();
      l1E[iMuon]                = (*muon).energy();
      l1Pt[iMuon]               = (*muon).pt();
      l1Eta[iMuon]              = (*muon).eta();
      l1Phi[iMuon]              = (*muon).phi();
      l1Vx[iMuon]               = (*muon).vx();
      l1Vy[iMuon]               = (*muon).vy();
      l1Vz[iMuon]               = (*muon).vz();
      l1Y[iMuon]                = (*muon).rapidity();    
      /// detector isolation 
      l1trackiso[iMuon]       = (*muon).isolationR03().sumPt;
      l1ecaliso[iMuon]        = (*muon).isolationR03().emEt;
      l1hcaliso[iMuon]        = (*muon).isolationR03().hadEt;
      /// ID
      l1Type[iMuon]  = (*muon).type();
      l1_numberOfChambers[iMuon]  = (*muon).numberOfChambers();      
      l1_numberOfMatches[iMuon]   = (*muon).numberOfMatches();
      /// PF isolation 
      l1pfiso_sumChargedHadronPt[iMuon]   = (*muon).pfIsolationR04().sumChargedHadronPt;
      l1pfiso_sumChargedParticlePt[iMuon] = (*muon).pfIsolationR04().sumChargedParticlePt;
      l1pfiso_sumNeutralHadronEt[iMuon]   = (*muon).pfIsolationR04().sumNeutralHadronEt;
      l1pfiso_sumPhotonEt[iMuon]          = (*muon).pfIsolationR04().sumPhotonEt;
      l1pfiso_sumPUPt[iMuon]              = (*muon).pfIsolationR04().sumPUPt;
      // vertex 
      const pat::Muon* patmuon1 = dynamic_cast<const pat::Muon *>( &*muon);
      l1_d0bsp[iMuon] = patmuon1->dB(pat::Muon::BS2D) ;
      l1_dz000[iMuon] = patmuon1->dB(pat::Muon::PV2D);
      l1_IP3D[iMuon] = patmuon1->dB(pat::Muon::PV3D);
      if(fabs(l1_IP3D[iMuon])>fabs(l1_dz000[iMuon])&&l1_IP3D[iMuon]<1000) 
      l1_dzPV[iMuon] = sqrt(l1_IP3D[iMuon]*l1_IP3D[iMuon]-l1_dz000[iMuon]*l1_dz000[iMuon]);

      if((*muon).isGlobalMuon()) l1_globalChi2[iMuon] = (*muon).globalTrack()->normalizedChi2();
      l1_innerChi2[iMuon] = (*muon).innerTrack()->normalizedChi2();
      l1_nPixelHits[iMuon] = (*muon).innerTrack()->hitPattern().numberOfValidPixelHits();
      l1_nTrackerHits[iMuon] = (*muon).track()->hitPattern().trackerLayersWithMeasurement();
      l1_isPF[iMuon] = (*muon).isPFMuon();
      l1_isGlobal[iMuon] = (*muon).isGlobalMuon();
      l1_isTracker[iMuon] = (*muon).isTrackerMuon();
      l1_hasMuonSegment[iMuon] = muon::isGoodMuon((*muon),muon::TMOneStationTight);

    }
  }


  //------------------ Electron filling ------------------
  typedef edm::View<reco::GsfElectron> ElectronView;
  edm::Handle<ElectronView> electrons;
  iEvent.getByLabel( "selectedPatElectronsPFlow", electrons);
  l2size  = (int) (electrons->size()); 
  if( l2size > 0 ) {
    int iElectron = 0;
    edm::View<reco::GsfElectron>::const_iterator electron, endpelectrons = electrons->end(); 
    for (electron = electrons->begin();  electron != endpelectrons;  ++electron, ++iElectron) {

      l2Charge[iElectron]           = (*electron).charge();
      l2px[iElectron]               = (*electron).px();
      l2py[iElectron]               = (*electron).py();
      l2pz[iElectron]               = (*electron).pz();
      l2E[iElectron]                = (*electron).energy();
      l2Pt[iElectron]               = (*electron).pt();
      l2Eta[iElectron]              = (*electron).eta();
      l2Phi[iElectron]              = (*electron).phi();
      l2Vx[iElectron]               = (*electron).vx();
      l2Vy[iElectron]               = (*electron).vy();
      l2Vz[iElectron]               = (*electron).vz();
      l2Y[iElectron]                = (*electron).rapidity();    
      /// detector isolation 
      l2trackiso[iElectron]       = (*electron).dr03TkSumPt();
      l2ecaliso[iElectron]        = (*electron).dr03EcalRecHitSumEt();
      l2hcaliso[iElectron]        = (*electron).dr03HcalTowerSumEt();
      /// ID
      l2_classification[iElectron]  = (*electron).classification();
      l2_HoverE[iElectron] = (*electron).hadronicOverEm(); 
      l2_EoverP[iElectron]       = (*electron).eSuperClusterOverP();
      l2_DeltaEta[iElectron]     = (*electron).deltaEtaSuperClusterTrackAtVtx();
      l2_DeltaPhi[iElectron]     = (*electron).deltaPhiSuperClusterTrackAtVtx();
      l2_numberOfBrems[iElectron]  = (*electron).numberOfBrems();      
      l2_BremFraction[iElectron]   = (*electron).fbrem();
      l2_SigmaIetaIeta[iElectron] = (*electron).sigmaIetaIeta();
      l2_missingHits[iElectron]   = (*electron).gsfTrack()->trackerExpectedHitsInner().numberOfHits();
      l2_dist[iElectron]          = (*electron).convDist();
      l2_dcot[iElectron]          = (*electron).convDcot();
      l2_convradius[iElectron]    = (*electron).convRadius();
      /// PF isolation 
      l2pfiso_chargedHadronIso[iElectron]   = (*electron).pfIsolationVariables().chargedHadronIso;
      l2pfiso_photonIso[iElectron] = (*electron).pfIsolationVariables().photonIso;
      l2pfiso_neutralHadronIso[iElectron]   = (*electron).pfIsolationVariables().neutralHadronIso;
      l2pfiso_EffAreaPU[iElectron] = ElectronEffectiveArea::GetElectronEffectiveArea( ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03 ,  
										      (*electron).superCluster()->eta(), ElectronEffectiveArea::kEleEAData2011);
      l2pfiso_pfIsoEA[iElectron] = (l2pfiso_chargedHadronIso[iElectron] +
							        max((float)(0.), l2pfiso_neutralHadronIso[iElectron]+
								    l2pfiso_photonIso[iElectron] -
								    l2pfiso_EffAreaPU[iElectron]*fastJetRho)) / l2Pt[iElectron];
      // vertex 
      const pat::Electron* patelectron1 = dynamic_cast<const pat::Electron *>( &*electron);
      l2_d0bsp[iElectron] = patelectron1->dB(pat::Electron::BS2D) ;
      l2_dz000[iElectron] = patelectron1->dB(pat::Electron::PV2D);
      l2_IP3D[iElectron] = patelectron1->dB(pat::Electron::PV3D);
      if(fabs(l2_IP3D[iElectron])>fabs(l2_dz000[iElectron])&&l2_IP3D[iElectron]<1000) 
      l2_dzPV[iElectron] = sqrt(l2_IP3D[iElectron]*l2_IP3D[iElectron]-l2_dz000[iElectron]*l2_dz000[iElectron]);
    }
  }



  //--------------- jet filling ----------------------------
  edm::Handle<edm::View<reco::Jet> > jets;
  iEvent.getByLabel( mInputJets, jets ); 
  jet_size = (int) (jets->size());
  if( jet_size > 0) {
    int iJet = 0;
    edm::View<reco::Jet>::const_iterator jet, endpjets = jets->end(); 
    for (jet = jets->begin();  jet != endpjets;  ++jet, ++iJet) {

      jet_Y[iJet]               = (*jet).rapidity();
      jet_Eta[iJet]             = (*jet).eta();
      jet_Phi[iJet]             = (*jet).phi();
      jet_E[iJet]               = (*jet).energy();
      jet_px[iJet]              = (*jet).px();
      jet_py[iJet]              = (*jet).py();
      jet_pz[iJet]              = (*jet).pz();
      jet_Pt[iJet]              = (*jet).pt();
      jet_area[iJet]            = (*jet).jetArea();

      edm::Ptr<reco::Jet> ptrJet = jets->ptrAt( jet - jets->begin() );  
      if ( ptrJet.isNonnull() && ptrJet.isAvailable() ) {
      const pat::Jet* pjet = dynamic_cast<const pat::Jet *>(ptrJet.get()) ;
      jet_bDiscriminatorSSVHE[iJet] = (*pjet).bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
      jet_bDiscriminatorTCHE[iJet] = (*pjet).bDiscriminator("trackCountingHighEffBJetTags");
      jet_bDiscriminatorCSV[iJet] = (*pjet).bDiscriminator("combinedSecondaryVertexBJetTags");
      jet_bDiscriminatorJP[iJet] = (*pjet).bDiscriminator("jetProbabilityBJetTags");
      jet_bDiscriminatorSSVHP[iJet] = (*pjet).bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
      jet_bDiscriminatorTCHP[iJet] = (*pjet).bDiscriminator("trackCountingHighPurBJetTags");
      }
    }
  }


  //--------------- photon filling ----------------------------
  edm::Handle<reco::PhotonCollection> photons;
  iEvent.getByLabel("photons",photons);
  edm::Handle<reco::BeamSpot> bsHandle;
  iEvent.getByLabel("offlineBeamSpot", bsHandle);
  const reco::BeamSpot &beamspot = *bsHandle.product();

  edm::Handle<reco::ConversionCollection> hConversions;
  iEvent.getByLabel("allConversions", hConversions);

  edm::Handle<reco::GsfElectronCollection> hElectrons;
  iEvent.getByLabel("gsfElectrons", hElectrons);

//   IsoDepositVals photonIsoValPFId(3);
//   iEvent.getByLabel("phoPFIso:chIsoForGsfEle", photonIsoValPFId[0]);
//   iEvent.getByLabel("phoPFIso:phIsoForGsfEle", photonIsoValPFId[1]);
//   iEvent.getByLabel("phoPFIso:nhIsoForGsfEle", photonIsoValPFId[2]);

//   const IsoDepositVals * photonIsoVals = &photonIsoValPFId;

  photon_size = (int) (photons->size());
  if( photon_size > 0) {
    int iPhoton = 0;
    for(unsigned ipho=0; ipho < photons->size(); ++ipho) {
      reco::PhotonRef myPhotonRef(photons,ipho);
      photon_px[iPhoton]  = myPhotonRef->px();
      photon_py[iPhoton]  = myPhotonRef->py();
      photon_pz[iPhoton]  = myPhotonRef->pz();
      photon_E[iPhoton]   = myPhotonRef->energy();
      photon_pt[iPhoton]  = myPhotonRef->et();
      photon_eta[iPhoton] = myPhotonRef->eta();
      photon_phi[iPhoton] = myPhotonRef->phi();
      photon_vx[iPhoton]  = myPhotonRef->vx();
      photon_vy[iPhoton]  = myPhotonRef->vy();
      photon_vz[iPhoton]  = myPhotonRef->vz();

//       photon_pfiso_charged[iPhoton]  = max((*(*photonIsoVals)[0])[myPhotonRef] - fastJetRho*EAch(fabs(myPhotonRef->eta())),0.);
//       photon_pfiso_photon[iPhoton]  = max((*(*photonIsoVals)[1])[myPhotonRef] - fastJetRho*EAnh(fabs(myPhotonRef->eta())),0.);
//       photon_pfiso_neutral[iPhoton]  = max((*(*photonIsoVals)[2])[myPhotonRef] - fastJetRho*EApho(fabs(myPhotonRef->eta())),0.);

      photon_trackiso[iPhoton]  = myPhotonRef->trkSumPtHollowConeDR04();
      photon_ecaliso[iPhoton]   = myPhotonRef->ecalRecHitSumEtConeDR04();
      photon_hcaliso[iPhoton]   = myPhotonRef->hcalTowerSumEtConeDR04();
      photon_HoverE[iPhoton]    = myPhotonRef->hadTowOverEm();
      photon_SigmaIetaIeta[iPhoton]  = myPhotonRef->sigmaIetaIeta();

      if(myPhotonRef->hasPixelSeed()) photon_hasPixelSeed[iPhoton] = 1;
      else photon_hasPixelSeed[iPhoton]   = 0;

      bool passConv = !(ConversionTools::hasMatchedPromptElectron(myPhotonRef->superCluster(), hElectrons, hConversions, beamspot.position()));
      if(passConv) photon_passElecVeto[iPhoton]   = 1;
      else photon_passElecVeto[iPhoton]   = 0;
    }
  }


  //  ************* Load Tracks ********************
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel("generalTracks", tracks);
  reco::TrackCollection::const_iterator trk;
  int iTracks = 0;
  if(tracks->size() > 0) {
    for ( trk = tracks->begin(); trk != tracks->end(); ++trk){
      if (fabs(trk->pt()) < 5.) continue;
      track_px[iTracks]  = trk->px();
      track_py[iTracks]  = trk->py();
      track_pz[iTracks]  = trk->pz();
      track_Vx[iTracks]  = trk->px();
      track_Vy[iTracks]  = trk->py();
      track_Vz[iTracks]  = trk->pz();
      track_Pt[iTracks]  = trk->pt();
      track_Eta[iTracks] = trk->eta();
      track_Phi[iTracks] = trk->phi();
      iTracks++;
    } 
  } 
  track_size =       iTracks;


  //  ************* Load Gsf Tracks ********************
  edm::Handle<reco::GsfTrackCollection> gsftracks;
  iEvent.getByLabel("electronGsfTracks", gsftracks);
  reco::GsfTrackCollection::const_iterator gsftrk;
  int igsf = 0;
  if(gsftracks->size() > 0) {
    for ( gsftrk = gsftracks->begin(); gsftrk != gsftracks->end(); ++gsftrk){
      if (fabs(gsftrk->pt()) < 5.) continue;
      gsftrack_px[igsf]  = gsftrk->px();
      gsftrack_py[igsf]  = gsftrk->py();
      gsftrack_pz[igsf]  = gsftrk->pz();
      gsftrack_Vx[igsf]  = gsftrk->px();
      gsftrack_Vy[igsf]  = gsftrk->py();
      gsftrack_Vz[igsf]  = gsftrk->pz();
      gsftrack_Pt[igsf]  = gsftrk->pt();
      gsftrack_Eta[igsf] = gsftrk->eta();
      gsftrack_Phi[igsf] = gsftrk->phi();
      igsf++;
    } 
  } 
  gsftrack_size =       igsf;



  //  ************* Load SuperClusters  ********************
  edm::Handle<reco::SuperClusterCollection> superClusters;
  iEvent.getByLabel("correctedHybridSuperClusters", superClusters);
  reco::SuperClusterCollection::const_iterator sc;
  int iSC = 0;
  if(superClusters->size() > 0) {
    for ( sc = superClusters->begin(); sc != superClusters->end(); ++sc){
      if (fabs(sc->energy()) < 5.) continue;
      superCluster_E[iSC]      = sc->energy();
      superCluster_rawE[iSC]   = sc->rawEnergy();
      superCluster_x[iSC]      = sc->x();
      superCluster_y[iSC]      = sc->y();
      superCluster_z[iSC]      = sc->z();
      superCluster_Eta[iSC]    = sc->eta();
      superCluster_Phi[iSC]    = sc->phi();
      superCluster_nHits[iSC]  = sc->size();
      iSC++;
    } 
  } 
  superCluster_size = iSC;


  //  ************* Load endcap superclusters  ********************
  edm::Handle<reco::SuperClusterCollection> superClusters5x5;
  iEvent.getByLabel("correctedMulti5x5SuperClustersWithPreshower", superClusters5x5);
  reco::SuperClusterCollection::const_iterator sc5x5;
  int iSC5x5 = 0;
  if(superClusters5x5->size() > 0) {
    for ( sc5x5 = superClusters5x5->begin(); sc5x5 != superClusters5x5->end(); ++sc5x5){
      if (fabs(sc5x5->energy()) < 5.) continue;
      superCluster5x5_E[iSC5x5]      = sc5x5->energy();
      superCluster5x5_rawE[iSC5x5]   = sc5x5->rawEnergy();
      superCluster5x5_x[iSC5x5]      = sc5x5->x();
      superCluster5x5_y[iSC5x5]      = sc5x5->y();
      superCluster5x5_z[iSC5x5]      = sc5x5->z();
      superCluster5x5_Eta[iSC5x5]    = sc5x5->eta();
      superCluster5x5_Phi[iSC5x5]    = sc5x5->phi();
      superCluster5x5_nHits[iSC5x5]  = sc5x5->size();
      iSC5x5++;
    } 
  } 
  superCluster5x5_size = iSC5x5;

  //  ************* Load CaloTowers  ********************
  edm::Handle<CaloTowerCollection> caloTowers;
  iEvent.getByLabel("towerMaker", caloTowers);
  CaloTowerCollection::const_iterator ct;
  int iCT = 0;
  if(caloTowers->size() > 0) {
    for ( ct = caloTowers->begin(); ct != caloTowers->end(); ++ct){
      if (ct->hadEt()+ct->emEt() < 2.) continue;
      caloTower_hadE [iCT]  = ct->hadEnergy();
      caloTower_emE  [iCT]  = ct->emEnergy();
      caloTower_hadEt[iCT]  = ct->hadEt();
      caloTower_emEt [iCT]  = ct->emEt();
      caloTower_Eta  [iCT]  = ct->eta();
      caloTower_Phi  [iCT]  = ct->phi();
      iCT++;
    } 
  } 
  caloTower_size = iCT;
  
  //  ************* Load Taus ********************
  typedef edm::View<pat::Tau> TauView;
  edm::Handle<TauView> taus;
  iEvent.getByLabel("selectedPatTausPFlow", taus);
  edm::View<pat::Tau>::const_iterator tau;
  int iTaus = 0;
  tau_size = (int) (taus->size());
  if(tau_size > 0) {
    for ( tau = taus->begin(); tau != taus->end(); ++tau){
      tau_px[iTaus]  = tau->px();
      tau_py[iTaus]  = tau->py();
      tau_pz[iTaus]  = tau->pz();
      tau_Pt[iTaus]  = tau->pt();
      tau_Eta[iTaus] = tau->eta();
      tau_Phi[iTaus] = tau->phi();
      iTaus++;
    } 
  } 



  // --- genParticles loop -------------------
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);
  //  const reco::Candidate *genPart=NULL;
  const reco::GenParticle* genPart=NULL;

  genPart_size = (int) (genParticles->size());  
  if(genPart_size > 0) {
    for(size_t i = 0; i < genParticles->size(); ++ i) {

      genPart = &((*genParticles)[i]);
      genPart_Charge[i]           = genPart->charge();
      genPart_Vx[i]               = genPart->vx();
      genPart_Vy[i]               = genPart->vy();
      genPart_Vz[i]               = genPart->vz();
      genPart_Y[i]                = genPart->rapidity();
      genPart_Eta[i]              = genPart->eta();    
      genPart_Phi[i]              = genPart->phi();
      genPart_E[i]                = genPart->energy();
      genPart_px[i]               = genPart->px();
      genPart_py[i]               = genPart->py();
      genPart_pz[i]               = genPart->pz();
      genPart_Pt[i]               = genPart->pt();
      genPart_Status[i]           = genPart->status();
      genPart_pdgId[i]            = genPart->pdgId();
    }// end genParticles loop
  }


  //--------------- genjet filling ----------------------------
  edm::Handle<edm::View<reco::Jet> > genjets;
  iEvent.getByLabel( "ak5GenJets", genjets ); 
  genjet_size = (int) (genjets->size());
  if( genjet_size > 0) {
    int iGenJet = 0;
    edm::View<reco::Jet>::const_iterator genjet, endpgenjets = genjets->end(); 
    for (genjet = genjets->begin();  genjet != endpgenjets;  ++genjet, ++iGenJet) {

      genjet_Y[iGenJet]               = (*genjet).rapidity();
      genjet_Eta[iGenJet]             = (*genjet).eta();
      genjet_Phi[iGenJet]             = (*genjet).phi();
      genjet_E[iGenJet]               = (*genjet).energy();
      genjet_px[iGenJet]              = (*genjet).px();
      genjet_py[iGenJet]              = (*genjet).py();
      genjet_pz[iGenJet]              = (*genjet).pz();
      genjet_Pt[iGenJet]              = (*genjet).pt();
      genjet_area[iGenJet]            = (*genjet).jetArea();
    }
  }






  tree_->Fill();
}


void ewk::SnowmassTreeProducer::endJob()
{
//   hOutputFile->SetCompressionLevel(2);
//   hOutputFile->cd();
//   tree_->Write();

//   delete tree_;
//   hOutputFile->Close();
//   delete hOutputFile;
}





////////////////// utilities, helpers ///////////////////
 
void ewk::SnowmassTreeProducer::SetBranch( float* x, std::string brName)
{
  tree_->Branch( brName.c_str(), x, ( brName+"/F").c_str() );
}


void ewk::SnowmassTreeProducer::SetBranch( int* x, std::string brName)
{
  tree_->Branch( brName.c_str(), x, ( brName+"/I").c_str() );
}


void ewk::SnowmassTreeProducer::SetBranch( bool* x, std::string brName)
{
  tree_->Branch( brName.c_str(), x, ( brName+"/O").c_str() );
}



//////////////////////////////////////////////////////////////////
/////// Helper for Effective Areas ///////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
float ewk::SnowmassTreeProducer::EAch( float x){
 float EA = 0.012;
 if(x>1.0)   EA = 0.010;
 if(x>1.479) EA = 0.014;
 if(x>2.0)   EA = 0.012;
 if(x>2.2)   EA = 0.016;
 if(x>2.3)   EA = 0.020;
 if(x>2.4)   EA = 0.012;
 return EA;
}

float ewk::SnowmassTreeProducer::EAnh( float x){
 float EA = 0.030;
 if(x>1.0)   EA = 0.057;
 if(x>1.479) EA = 0.039;
 if(x>2.0)   EA = 0.015;
 if(x>2.2)   EA = 0.024;
 if(x>2.3)   EA = 0.039;
 if(x>2.4)   EA = 0.072;
 return EA;
}

float ewk::SnowmassTreeProducer::EApho( float x){
 float EA = 0.148;
 if(x>1.0)   EA = 0.130;
 if(x>1.479) EA = 0.112;
 if(x>2.0)   EA = 0.216;
 if(x>2.2)   EA = 0.262;
 if(x>2.3)   EA = 0.260;
 if(x>2.4)   EA = 0.266;
 return EA;
}


// declare this class as a plugin
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
using ewk::SnowmassTreeProducer;
DEFINE_FWK_MODULE(SnowmassTreeProducer);
