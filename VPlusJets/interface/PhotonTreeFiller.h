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
 *   To fill jet related quantities into a specified TTree
 *   Can work with CaloJet, GenJet, JPT jet, PF jet.
 *   Can work with jets in RECO/AOD/PAT data formats.
 * History:
 *   
 *
 * Copyright (C) 2010 FNAL 
 *****************************************************************************/

#ifndef ElectroWeakAnalysis_VPlusJets_PhotonTreeFiller_h
#define ElectroWeakAnalysis_VPlusJets_PhotonTreeFiller_h

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h" 
#include <map>

#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "EGamma/EGammaAnalysisTools/interface/PFIsolationEstimator.h"

namespace ewk {

  class PhotonTreeFiller {
  public:
    /// specify the name of the TTree, and the configuration for it
    PhotonTreeFiller(const char *name, TTree* tree, 
		  const edm::ParameterSet iConfig );


    /// default constructor
    PhotonTreeFiller() {};


    /// Destructor, does nothing 
      ~PhotonTreeFiller() {};


    /// To be called once per event to fill the values for jets
    void fill(const edm::Event &iEvent);

    static const int NUM_PHO_MAX = 10;

  protected:

    /// To be called once per event, to initialize variable values to default
    void init() const ;
    /// Helper function for main constructor
    void SetBranches(); 
    void SetBranch( float* x, std::string name);
    void SetBranch( int* x, std::string name);
    void SetBranchSingle( float* x, std::string name);
    void SetBranchSingle( int* x, std::string name);

    float EAch(float x);
    float EAnh(float x);
    float EApho(float x);
    void FillBranches() const;
    void init();

    TTree* tree_;
    edm::InputTag mPrimaryVertex;
    edm::InputTag mInputPhotons;
    edm::InputTag inputTagPFCandidateMap_;

    std::vector<edm::InputTag> inputTagIsoDepPhotons_;
    std::vector<edm::InputTag> inputTagIsoValPhotonsPFId_;   

    edm::InputTag sourceByValue;
    bool runoverAOD;
    // 'mutable' because we will fill it from a 'const' method
    mutable std::vector<std::string> bnames;

  private:
    // private data members
    
    int NumPhotons; 
    float Et[NUM_PHO_MAX];
    float E[NUM_PHO_MAX];
    float Eta[NUM_PHO_MAX];
    float Phi[NUM_PHO_MAX];
    float Theta[NUM_PHO_MAX];

    float Px[NUM_PHO_MAX];
    float Py[NUM_PHO_MAX];
    float Pz[NUM_PHO_MAX];
    float Vx[NUM_PHO_MAX];
    float Vy[NUM_PHO_MAX];
    float Vz[NUM_PHO_MAX];

    float SC_Et[NUM_PHO_MAX];
    float SC_E[NUM_PHO_MAX];
    float SC_Eta[NUM_PHO_MAX];
    float SC_Phi[NUM_PHO_MAX];
    float SC_Theta[NUM_PHO_MAX];

    float SC_x[NUM_PHO_MAX];
    float SC_y[NUM_PHO_MAX];
    float SC_z[NUM_PHO_MAX];

    float PFisocharged03[NUM_PHO_MAX];
    float PFisophoton03[NUM_PHO_MAX];
    float PFisoneutral03[NUM_PHO_MAX];
    float trkSumPtHollConeDR04_Pho11[NUM_PHO_MAX];
    float ecalRecHitSumEtConeDR04_Pho11[NUM_PHO_MAX];
    float hcalTowSumEtConeDR04_Pho11[NUM_PHO_MAX];

    float HoverE[NUM_PHO_MAX];
    float HoverE2011[NUM_PHO_MAX];
    float SigmaIetaIeta[NUM_PHO_MAX];
    int hasPixelSeed[NUM_PHO_MAX];
    int passElecVeto[NUM_PHO_MAX];

    int Id2011[NUM_PHO_MAX];
    int Id2012[NUM_PHO_MAX];

    PFIsolationEstimator isolator;

  //Pfiso variables
  float  charged03;
  float photon03;
  float neutral03;
  unsigned nrecopho;

  //PFisolation
  typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;
	
  };

} //namespace

#endif


