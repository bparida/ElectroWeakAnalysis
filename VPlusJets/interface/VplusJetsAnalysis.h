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
 *   To fill W/Z + jets related quantities into a specified TTree
 *   Can work with CaloJet, GenJet, JPT jet, PF jet.
 *   Can work with jets in RECO/AOD/PAT data formats.
 * History:
 *   
 *
 * Copyright (C) 2010 FNAL 
 *****************************************************************************/

#ifndef VplusJetsAnalysis_h
#define VplusJetsAnalysis_h

// system include files
#include <memory>
#include <string>
#include <iostream>
#include <map>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h" 

#include "TFile.h"
#include "TTree.h"

#include "ElectroWeakAnalysis/VPlusJets/interface/JetTreeFiller.h"
//#include "ElectroWeakAnalysis/VPlusJets/interface/GroomedJetFiller.h"
//#include "ElectroWeakAnalysis/VPlusJets/interface/PhotonTreeFiller.h"
#include "ElectroWeakAnalysis/VPlusJets/interface/VtoElectronTreeFiller.h"
#include "ElectroWeakAnalysis/VPlusJets/interface/VtoMuonTreeFiller.h"
#include "ElectroWeakAnalysis/VPlusJets/interface/MCTreeFiller.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/GhostedAreaSpec.hh"


//
// class decleration
//

namespace ewk
{
  class VplusJetsAnalysis : public edm::EDAnalyzer {
  public:
    explicit VplusJetsAnalysis(const edm::ParameterSet&);
    ~VplusJetsAnalysis();

    virtual void beginJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup& iSetup);
    virtual void declareTreeBranches();
    virtual void endJob() ;


  private:
    // ----------member data ---------------------------
    // names of modules, producing object collections
     
    /// output ROOT file for the tree and histograms
    TTree*  myTree;
    edm::Service<TFileService> fs;
    bool runningOverMC_;
    std::string VBosonType_;
    std::string LeptonType_;
    edm::InputTag mInputBoson;
    edm::InputTag mPrimaryVertex;

    edm::InputTag mInputMet;
    edm::InputTag mInputMetMVA;
	std::string JetsFor_rho;
	edm::InputTag mInputgenMet;
	bool runoverAOD;

    /// The objects that actually computes variables and fill the tree 
    std::auto_ptr<ewk::JetTreeFiller> CorrectedPFJetFiller;
    std::auto_ptr<ewk::JetTreeFiller> CorrectedPFJetFillerVBFTag; //For VBF Tag Jets

    /*std::auto_ptr<ewk::GroomedJetFiller> AK5groomedJetFiller;
    std::auto_ptr<ewk::GroomedJetFiller> AK7groomedJetFiller;
    std::auto_ptr<ewk::GroomedJetFiller> AK8groomedJetFiller;
    std::auto_ptr<ewk::GroomedJetFiller> CA8groomedJetFiller;
    std::auto_ptr<ewk::GroomedJetFiller> CA12groomedJetFiller;

    std::auto_ptr<ewk::GroomedJetFiller> genAK5groomedJetFiller;
    std::auto_ptr<ewk::GroomedJetFiller> genAK7groomedJetFiller;
    std::auto_ptr<ewk::GroomedJetFiller> genAK8groomedJetFiller;
    std::auto_ptr<ewk::GroomedJetFiller> genCA8groomedJetFiller;
    std::auto_ptr<ewk::GroomedJetFiller> genCA12groomedJetFiller;*/


    std::auto_ptr<ewk::JetTreeFiller> GenJetFiller;
  //  std::auto_ptr<ewk::PhotonTreeFiller> PhotonFiller;
    std::auto_ptr<ewk::VtoElectronTreeFiller> recoBosonFillerE;
    std::auto_ptr<ewk::VtoMuonTreeFiller> recoBosonFillerMu;
    std::auto_ptr<ewk::MCTreeFiller> genBosonFiller;


    // private data members
    int run;
    int event; 
    int lumi; 
    int bunch; 
    int nPV; 
    int mNVB;
    float mpfMET;
    float mpfSumET;
    float mpfMETSign;
    float mpfMETPhi;
    float spfMET;
    float spfSumET;
    float spfMETSign;
    float spfMETPhi;
    float mvaMET;
    float mvaSumET;
    float mvaMETSign;
    float mvaMETPhi;
    float fastJetRho;
    float genMET;
    float genSumET;
    float genMETSign;
    float genMETPhi;

    float mcPUtrueInteractions;
    float mcPUtotnvtx;
    float mcPUbx[3];
    float mcPUnvtx[3];
  };
}
#endif
