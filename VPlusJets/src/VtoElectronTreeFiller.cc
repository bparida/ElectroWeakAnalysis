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
 *   To fill W--> enu or Z-->ee  related quantities into a specified TTree
 *   Can work with jets in RECO/AOD/PAT data formats.
 * History:
 *   
 *
 * Copyright (C) 2010 FNAL 
 *****************************************************************************/

// CMS includes
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// Header file
#include "ElectroWeakAnalysis/VPlusJets/interface/VtoElectronTreeFiller.h"
#include "ElectroWeakAnalysis/VPlusJets/interface/METzCalculator.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"



ewk::VtoElectronTreeFiller::VtoElectronTreeFiller(const char *name, TTree* tree, 
							const edm::ParameterSet iConfig)
{

  // ********** Vector boson ********** //
  if(  iConfig.existsAs<edm::InputTag>("srcVectorBoson") )
    mInputBoson = iConfig.getParameter<edm::InputTag>("srcVectorBoson"); 
  else std::cout << "***Error:" << name << 
    " Collection not specified !" << std::endl;
    if(  iConfig.existsAs<edm::InputTag>("srcMet") )
		mInputMet = iConfig.getParameter<edm::InputTag>("srcMet");
    if(  iConfig.existsAs<edm::InputTag>("srcMetMVA") )
                mInputMetMVA = iConfig.getParameter<edm::InputTag>("srcMetMVA");
    if(  iConfig.existsAs<edm::InputTag>("srcElectrons") )
		mInputElectrons  = iConfig.getParameter<edm::InputTag>("srcElectrons");
    if(  iConfig.existsAs<edm::InputTag>("srcBeamSpot") )
		mInputBeamSpot  = iConfig.getParameter<edm::InputTag>("srcBeamSpot");
	if( iConfig.existsAs<bool>("runningOverAOD"))
		runoverAOD = iConfig.getParameter<bool>("runningOverAOD");
  tree_     = tree;
  name_     = name;
  Vtype_    = iConfig.getParameter<std::string>("VBosonType"); 
  LeptonType_ = iConfig.getParameter<std::string>("LeptonType");
  if( !(tree==0) && LeptonType_=="electron") SetBranches();
}





void ewk::VtoElectronTreeFiller::SetBranches()
{
  // Declare jet branches
  std::string lept1 = "eplus";
  std::string lept2 = "eminus";
  if( !(Vtype_=="Z") ) lept1 = "electron";

  SetBranch( &V_mass,      "mass");
  SetBranch( &V_mt,        "mt");
  SetBranch( &V_mtMVA,        "mtMVA");
  SetBranch( &V_px,        "px");
  SetBranch( &V_py,        "py");
  SetBranch( &V_pz,        "pz");
  SetBranch( &V_E,         "e");
  SetBranch( &V_Pt,        "pt");
  SetBranch( &V_Et,        "et");
  SetBranch( &V_Eta,       "eta");    
  SetBranch( &V_Phi,       "phi");
  SetBranch( &V_Vx,        "vx");
  SetBranch( &V_Vy,        "vy");
  SetBranch( &V_Vz,        "vz");
  SetBranch( &V_Y,         "y");
  SetBranch( &nTightElectron, "numTightElectrons");
  SetBranch( &nLooseElectron, "numLooseElectrons");
  if(Vtype_=="W") {
    SetBranch( &V_pzNu1,     "pzNu1");
    SetBranch( &V_pzNu2,     "pzNu2");
  }
  ///////////////////////////////////////////////
  SetBranch( &e1px,             lept1+"_px" );
  SetBranch( &e1py,             lept1+"_py" );
  SetBranch( &e1pz,             lept1+"_pz" );
  SetBranch( &e1E,              lept1+"_e" );
  SetBranch( &e1Pt,             lept1+"_pt" );
  SetBranch( &e1Et,             lept1+"_et" );
  SetBranch( &e1Eta,            lept1+"_eta" ); 
  SetBranch( &e1Theta,          lept1+"_theta" ); 
  SetBranch( &e1Phi,            lept1+"_phi" );
  SetBranch( &e1Charge,         lept1+"_charge" );
  SetBranch( &e1Vx,             lept1+"_vx" );
  SetBranch( &e1Vy,             lept1+"_vy" );
  SetBranch( &e1Vz,             lept1+"_vz" );
  SetBranch( &e1Y,              lept1+"_y" );
  SetBranch( &e1_trackiso,      lept1+"_trackiso" );
  SetBranch( &e1_hcaliso,       lept1+"_hcaliso" );
  SetBranch( &e1_ecaliso,       lept1+"_ecaliso" );
  SetBranch( &e1Classification, lept1+"_classification" );
//   SetBranch( &e1_sc_x,          lept1+"_sc_x" );
//   SetBranch( &e1_sc_y,          lept1+"_sc_y" );
//   SetBranch( &e1_sc_z,          lept1+"_sc_z" );
//   SetBranch( &e1_sc_Theta,      lept1+"_sc_Theta" );
//   SetBranch( &e1_sc_Eta,        lept1+"_sc_Eta" );
//   SetBranch( &e1_sc_Phi,        lept1+"_sc_Phi" );
//   SetBranch( &e1_sc_E,          lept1+"_sc_E" );
//   SetBranch( &e1_sc_px,         lept1+"_sc_px" );
//   SetBranch( &e1_sc_py,         lept1+"_sc_py" );
//   SetBranch( &e1_sc_pz,         lept1+"_sc_pz" );
//   SetBranch( &e1_sc_Pt,         lept1+"_sc_Pt" );
//   SetBranch( &e1_sc_Et,         lept1+"_sc_Et" );	  
//   SetBranch( &e1_EoverPout,     lept1+"_eoverp_out" );
//   SetBranch( &e1_EoverPin,      lept1+"_eoverp_in" );
  //  SetBranch( &e1_numberOfBrems, lept1+"_numbrem" );
//   SetBranch( &e1_BremFraction,  lept1+"_fbrem" );
//   SetBranch( &e1_DeltaEtaIn,    lept1+"_deltaeta_in" );
//   SetBranch( &e1_DeltaPhiIn,    lept1+"_deltaphi_in" );
//   SetBranch( &e1_DeltaPhiOut,   lept1+"_deltaphi_out" );
//   SetBranch( &e1_DeltaEtaOut,   lept1+"_deltaeta_out" );
//   SetBranch( &e1_Trackmom_calo, lept1+"_trackmom_calo" );
//   SetBranch( &e1_Trackmom_vtx,  lept1+"_trackmom_vtx" );	  
//   SetBranch( &e1_HoverE,        lept1+"_hovere" );	    
//   SetBranch( &e1_E9overE25,     lept1+"_e9e25" );
//   SetBranch( &e1_SigmaEtaEta,   lept1+"_sigmaetaeta" );
//   SetBranch( &e1_SigmaIetaIeta, lept1+"_sigmaietaieta" );
  SetBranch( &e1_missingHits,   lept1+"_missingHits" );	  
  SetBranch( &e1_dist,          lept1+"_dist" );
  SetBranch( &e1_dcot,          lept1+"_dcot" );
  //  SetBranch( &e1_convradius,    lept1+"_convradius" );
//   SetBranch( &ise1WP95,         lept1+"_isWP95" );
//   SetBranch( &ise1WP80,         lept1+"_isWP80" );
//   SetBranch( &ise1WP70,         lept1+"_isWP70" );
  SetBranch( &e1_d0bsp,         lept1+"_d0bsp" );
  SetBranch( &e1_dz000,         lept1+"_dz000" );
  SetBranch( &e1_mvaTrigV0,         lept1+"_mvaTrigV0" );
  SetBranch( &e1_mvaNonTrigV0,         lept1+"_mvaNonTrigV0" );

  SetBranch( &e1_pfiso_chargedHadronIso,         lept1+"_pfiso_chargedHadronIso" );
  SetBranch( &e1_pfiso_photonIso,                lept1+"_pfiso_photonIso" );
  SetBranch( &e1_pfiso_neutralHadronIso,         lept1+"_pfiso_neutralHadronIso" );
  SetBranch( &e1_pfIsoEA,         lept1+"_pfIsoEA" );


  ////////////////////////////////////////////////////////
  if(Vtype_=="Z") {	  
    SetBranch( &e2px,             lept2+"_px" );
    SetBranch( &e2py,             lept2+"_py" );
    SetBranch( &e2pz,             lept2+"_pz" );
    SetBranch( &e2E,              lept2+"_e" );
    SetBranch( &e2Pt,             lept2+"_pt" );
    SetBranch( &e2Et,             lept2+"_et" );
    SetBranch( &e2Eta,            lept2+"_eta" ); 
    SetBranch( &e2Theta,          lept2+"_theta" );    
    SetBranch( &e2Phi,            lept2+"_phi" );
    SetBranch( &e2Charge,         lept2+"_charge" );
    SetBranch( &e2Vx,             lept2+"_vx" );
    SetBranch( &e2Vy,             lept2+"_vy" );
    SetBranch( &e2Vz,             lept2+"_vz" );
    SetBranch( &e2Y,              lept2+"_y" );
    SetBranch( &e2_trackiso,      lept2+"_trackiso" );
    SetBranch( &e2_hcaliso,       lept2+"_hcaliso" );
    SetBranch( &e2_ecaliso,       lept2+"_ecaliso");
    SetBranch( &e2Classification, lept2+"_classification" );
//     SetBranch( &e2_sc_x,          lept2+"_sc_x" );
//     SetBranch( &e2_sc_y,          lept2+"_sc_y" );
//     SetBranch( &e2_sc_z,          lept2+"_sc_z" );
//     SetBranch( &e2_sc_Theta,      lept2+"_sc_Theta" );
//     SetBranch( &e2_sc_Eta,        lept2+"_sc_Eta" );
//     SetBranch( &e2_sc_Phi,        lept2+"_sc_Phi" );
//     SetBranch( &e2_sc_E,          lept2+"_sc_E" );
//     SetBranch( &e2_sc_px,         lept2+"_sc_px" );
//     SetBranch( &e2_sc_py,         lept2+"_sc_py" );
//     SetBranch( &e2_sc_pz,         lept2+"_sc_pz" );
//     SetBranch( &e2_sc_Pt,         lept2+"_sc_Pt" );
//     SetBranch( &e2_sc_Et,         lept2+"_sc_Et" );	  
//     SetBranch( &e2_EoverPout,     lept2+"_eoverp_out" );
//     SetBranch( &e2_EoverPin,      lept2+"_eoverp_in" );
//     SetBranch( &e2_numberOfBrems, lept2+"_numbrem" );
//     SetBranch( &e2_BremFraction,  lept2+"_fbrem" );
//     SetBranch( &e2_DeltaEtaIn,    lept2+"_deltaeta_in" );
//     SetBranch( &e2_DeltaPhiIn,    lept2+"_deltaphi_in" );
//     SetBranch( &e2_DeltaPhiOut,   lept2+"_deltaphi_out" );
//     SetBranch( &e2_DeltaEtaOut,   lept2+"_deltaeta_out" );
//     SetBranch( &e2_Trackmom_calo, lept2+"_trackmom_calo" );
//     SetBranch( &e2_Trackmom_vtx,  lept2+"_trackmom_vtx" );	  
//     SetBranch( &e2_HoverE,        lept2+"_hovere" );	   	  
//     SetBranch( &e2_E9overE25,     lept2+"_e9e25" );
//     SetBranch( &e2_SigmaEtaEta,   lept2+"_sigmaetaeta" );
//     SetBranch( &e2_SigmaIetaIeta, lept2+"_sigmaietaieta" );
    SetBranch( &e2_missingHits,   lept2+"_missingHits" );
    SetBranch( &e2_dist,          lept2+"_dist" );
    SetBranch( &e2_dcot,          lept2+"_dcot" );
//     SetBranch( &e2_convradius,    lept2+"_convradius" );
//     SetBranch( &ise2WP95,         lept2+"_isWP95" );
//     SetBranch( &ise2WP80,         lept2+"_isWP80" );
//     SetBranch( &ise2WP70,         lept2+"_isWP70" );
  }	  
}
/////////////////////////////////////////////////////////////////////////






void ewk::VtoElectronTreeFiller::init()   
{
  // initialize private data members
  nLooseElectron = 0;
  nTightElectron = 0;
  V_mass                = -1.;
  V_mt                  = -1.;
  V_mtMVA                  = -1.;
  V_px                  = -99999.;
  V_py                  = -99999.;
  V_pz                  = -99999.;
  V_E                   = -1.;
  V_Pt                  = -1.;
  V_Et                  = -1.;
  V_Eta                 = -10.;
  V_Phi                 = -10.;
  V_Vx                  = -10.;
  V_Vy                  = -10.;
  V_Vz                  = -10.;
  V_Y                   = -10.;
  V_pzNu1               = -10000.0;
  V_pzNu2               = -10000.0;

  ise1WP95          = false;
  ise1WP80          = false;
  ise1WP70          = false;
  ise2WP95          = false;
  ise2WP80          = false;
  ise2WP70          = false;

  e1Classification   = -1; 
  e1Charge           = -10;
  e2Classification  = -1;
  e2Charge          = -10;

  e1px               = -99999.;
  e1py               = -99999.;
  e1pz               = -99999.;
  e1E                = -1.;
  e1Et               = -1.;
  e1Pt               = -1.;
  e1Eta              = -10.;
  e1Theta            = -99999.;
  e1Phi              = -10.;
  e1Vx               = -10.;
  e1Vy               = -10.;
  e1Vz               = -10.;
  e1Y                = -10.;
  e1_sc_x            = -10.;
  e1_sc_y            = -10.;
  e1_sc_z            = -10.;
  e1_sc_Theta        = -99999.;
  e1_sc_Eta          = -10.;    
  e1_sc_Phi          = -10.;
  e1_sc_E            = -1.;
  e1_sc_px           = -99999.;
  e1_sc_py           = -99999.;
  e1_sc_pz           = -99999.;
  e1_sc_Pt           = -1.;
  e1_sc_Et           = -1.;	  
  e1_EoverPout       = -1.;
  e1_EoverPin        = -1.;
  e1_numberOfBrems   = -10.;
  e1_BremFraction    = -1.;
  e1_DeltaEtaIn      = -10.;
  e1_DeltaPhiIn      = -10.;
  e1_DeltaPhiOut     = -10.;
  e1_DeltaEtaOut     = -10.;
  e1_Trackmom_calo   = -1.;
  e1_Trackmom_vtx    = -1.;	  
  e1_HoverE          = -1.;	   	  
  e1_E9overE25       = -10.;
  e1_SigmaEtaEta     = -1.;
  e1_SigmaIetaIeta   = -1.;	  
  e1_missingHits     = 100;
  e1_dist            = -10.;
  e1_dcot            = -10.;
  e1_convradius      = -10.;
  e1_d0bsp           = -99999.;
  e1_dz000           = -99999.;
  e1_mvaTrigV0       = -99999.;
  e1_mvaNonTrigV0    = -99999.;

  e1_pfiso_chargedHadronIso       = -99999.;
  e1_pfiso_photonIso              = -99999.;
  e1_pfiso_neutralHadronIso       = -99999.;
  e1_EffArea                      = 0.;
  e1_pfIsoEA                      = -99999.;

  e2px              = -99999.;
  e2py              = -99999.;
  e2pz              = -99999.;
  e2E               = -1.;
  e2Pt              = -1.;
  e2Et              = -1.;
  e2Eta             = -10.;
  e2Theta           = -99999.;
  e2Phi             = -10.;
  e2Vx              = -10.;
  e2Vy              = -10.;
  e2Vz              = -10.;
  e2Y               = -10.;
  e2_sc_x           = -10.;
  e2_sc_y           = -10.;
  e2_sc_z           = -10.;
  e2_sc_Theta       = -99999.;
  e2_sc_Eta         = -10.;    
  e2_sc_Phi         = -10.;
  e2_sc_E           = -1.;
  e2_sc_px          = -99999.;
  e2_sc_py          = -99999.;
  e2_sc_pz          = -99999.;
  e2_sc_Pt          = -1.;
  e2_sc_Et          = -1.;	  
  e2_EoverPout       = -1.;
  e2_EoverPin        = -1.;
  e2_numberOfBrems   = -10.;
  e2_BremFraction    = -1.;
  e2_DeltaEtaIn      = -10.;
  e2_DeltaPhiIn      = -10.;
  e2_DeltaPhiOut     = -10.;
  e2_DeltaEtaOut     = -10.;
  e2_Trackmom_calo   = -1.;
  e2_Trackmom_vtx    = -1.;	  
  e2_HoverE          = -1.;	    
  e2_E9overE25       = -10.;
  e2_SigmaEtaEta     = -1.;
  e2_SigmaIetaIeta     = -1.;	  
  e2_missingHits     = 100;
  e2_dist            = -10.;
  e2_dcot            = -10.;
  e2_convradius      = -10.;
	  
  //////////////
  e1_trackiso     = 5000.0;
  e1_hcaliso      = 5000.0;
  e1_ecaliso      = 5000.0;
  e2_trackiso     = 5000.0;
  e2_hcaliso      = 5000.0;
  e2_ecaliso      = 5000.0;

  // initialization done
}

void ewk::VtoElectronTreeFiller::fill(const edm::Event& iEvent, int vecBosonIndex)
{
  // protection
  if( (tree_==0) || !(LeptonType_=="electron") )  return;

  // first initialize to the default values
  init();

  edm::Handle<reco::CandidateView> boson;
  iEvent.getByLabel( mInputBoson, boson);
  if( boson->size()<1 ) return; // Nothing to fill
  

//  edm::Handle<reco::PFMETCollection> pfmet;
  // iEvent.getByLabel("pfMet", pfmet);

   edm::Handle<edm::View<reco::MET> > pfmet;
   iEvent.getByLabel(mInputMet, pfmet);

 /////// Pileup density "rho" in the event from fastJet pileup calculation /////
  edm::Handle<double> rho;
  const edm::InputTag eventrho("kt6PFJetsPFlow", "rho");
  iEvent.getByLabel(eventrho,rho);
  double fastJetRho = *rho;

  nTightElectron = 0;
  nLooseElectron = 0;

  // Loop over electrons for loose/tight counting
// edm::Handle<reco::GsfElectronCollection> electrons;
// iEvent.getByLabel("gsfElectrons", electrons);
 edm::Handle<edm::View<reco::GsfElectron> > electrons;
 iEvent.getByLabel(mInputElectrons, electrons);
 for(edm::View<reco::GsfElectron>::const_iterator 
       elec = electrons->begin(); elec != electrons->end();++elec) {
   if( isTightElectron( iEvent, *elec) ) nTightElectron++;
   if( isLooseElectron( iEvent, *elec) ) nLooseElectron++;
 }


  const reco::Candidate *Vboson = &((*boson)[vecBosonIndex]); 
  if( Vboson == 0) return;

  /////// MVA MET information /////
  edm::Handle<edm::View<reco::MET> > metMVA;
  iEvent.getByLabel(mInputMetMVA, metMVA);

  ////////// Vector boson quantities //////////////
  V_mass = Vboson->mass();
  V_mt = sqrt(2.0*Vboson->daughter(0)->pt()*Vboson->daughter(1)->pt()*
	      (1.0-cos(Vboson->daughter(0)->phi()-Vboson->daughter(1)->phi())));
  if (metMVA->size() > 0)
        V_mtMVA = sqrt(2.0*Vboson->daughter(0)->pt()*(*metMVA)[0].et()*
              (1.0-cos(Vboson->daughter(0)->phi()-(*metMVA)[0].phi())));
  V_Eta = Vboson->eta();   
  V_Phi = Vboson->phi();
  V_Vx = Vboson->vx();
  V_Vy = Vboson->vy();
  V_Vz = Vboson->vz();
  V_Y  = Vboson->rapidity();
  V_px = Vboson->px();
  V_py = Vboson->py();
  V_pz = Vboson->pz();
  V_E  = Vboson->energy();
  V_Pt = Vboson->pt();
  V_Et = Vboson->et();

  // now iterate over the daughters  
  if(Vboson->numberOfDaughters()<2 ) {
    throw cms::Exception( "***Error: V boson has < 2 daughters !\n");
    return;  // if no electron found, then return
  } 
  std::cout << "number of daugters " << Vboson->numberOfDaughters() <<std::endl;
  // get the two daughters
  reco::CandidateBaseRef m0 = Vboson->daughter(0)->masterClone();
  reco::CandidateBaseRef m1 = Vboson->daughter(1)->masterClone();
  std::cout << "ptr m1 is nonnull " << m1.isNonnull() <<std::endl ;

  const reco::GsfElectron* e1=NULL;
  const reco::GsfElectron* e2=NULL;
  const std::type_info & type0 = typeid(*m0);
  const std::type_info & type1 = typeid(*m1);

  if(runoverAOD){
	  if( type0 == typeid(reco::GsfElectron) ){
		  e1 = dynamic_cast<const reco::GsfElectron *>( &*m0 ); }}
  
  else{
    if( type0 == typeid(pat::Electron) ){
      //  std::cout << " boson daughter - pat el " << std::endl;      
      //	  std::cout << "ptr is nonnull " << m0.isNonnull() <<std::endl ;
      //	  std::cout << "charge from Ptr = " << m0->charge() << std::endl;
      e1 = dynamic_cast<const pat::Electron *>(&*m0) ;
      //	  std::cout << "charge from el 1 = " << e1->charge() << std::endl;
      //  std::cout << "pt from el 1 = " << e1->pt() << std::endl;
    }
  }

  if(runoverAOD){
	  if( type1 == typeid(reco::GsfElectron) ){
		  e2 = dynamic_cast<const reco::GsfElectron *>( &*m1 ); }}
  
  else{
    //std::cout << " boson daughter - pat el 2" << std::endl;
    if( type1 == typeid(pat::Electron) ){
	    //	  std::cout << " boson daughter - pat el 2" << std::endl;
      //	  std::cout << "ptr is nonnull " << m1.isNonnull() <<std::endl ;
      //  std::cout << "charge from Ptr = " << m1->charge() << std::endl;
      e2 = dynamic_cast<const pat::Electron *>(&*m1) ;
      //	  std::cout << "charge from el 2 = " << e2->charge() << std::endl;
      //	  std::cout << "pt from el 2 = " << e2->pt() << std::endl;
    }
  }

   if(0==e1 && 0==e2) {
    throw cms::Exception("***Error: couldn't" 
			 " do dynamic cast of vector boson daughters !\n");
    return;  // if no electron found, then return
  } 
  

 	  const reco::GsfElectron* ele1=NULL;
	  const reco::GsfElectron* ele2=NULL;

		 
  // if Z--> e+e- then ele1 = e+, ele2 = e-
  if(Vtype_=="Z") {
    if(e1->charge() > 0) {  ele1 = e1;   ele2 = e2; }
    else { ele1 = e2;  ele2 = e1; }
  }
  // if W--> enu then ele1 = e, ele2 = NULL 
  if(Vtype_=="W") {
    if( abs(e1->charge())==1 ) ele1  = e1;
    else if( abs(e2->charge())==1 ) ele1  = e2;

    if( !(ele1 == NULL) ) {
      // estimate Pz of neutrino
      TLorentzVector p4MET((*pfmet)[0].px(), (*pfmet)[0].py(), (*pfmet)[0].pz(), (*pfmet)[0].energy());
      TLorentzVector p4lepton(ele1->px(), ele1->py(), ele1->pz(), ele1->energy());
      METzCalculator metz;
      metz.SetMET(p4MET);
      metz.SetLepton(p4lepton);
      if (LeptonType_=="electron") metz.SetLeptonType("electron");
      if (LeptonType_=="muon")     metz.SetLeptonType("muon");
      V_pzNu1 = metz.Calculate();
      V_pzNu2 = metz.getOther();
    }
  }

	  
  ////////// electron #1 quantities //////////////
  if( !(ele1 == NULL) ) {
    e1Charge           = ele1-> charge();
    e1Vx               = ele1->vx();
    e1Vy               = ele1->vy();
    e1Vz               = ele1->vz();
    e1Y                = ele1->rapidity();
    e1Theta            = ele1->theta();
    e1Eta              = ele1->eta();    
    e1Phi              = ele1->phi();
    e1E                = ele1->energy();
    e1px               = ele1->px();
    e1py               = ele1->py();
    e1pz               = ele1->pz();
    e1Pt               = ele1->pt();
    e1Et               = ele1->et();

	  /// isolation 
    e1_trackiso       = ele1->dr03TkSumPt();
    e1_ecaliso        = ele1->dr03EcalRecHitSumEt();
    e1_hcaliso        = ele1->dr03HcalTowerSumEt();
    e1Classification  = ele1->classification();
    e1_numberOfBrems  = ele1->numberOfBrems();      
    e1_BremFraction   = ele1->fbrem();
    e1_DeltaEtaIn     = ele1->deltaEtaSuperClusterTrackAtVtx();
    e1_DeltaPhiIn     = ele1->deltaPhiSuperClusterTrackAtVtx();
    e1_DeltaPhiOut    = ele1->deltaPhiSeedClusterTrackAtCalo();
    e1_DeltaEtaOut    = ele1->deltaEtaSeedClusterTrackAtCalo();
	
    //Track Momentum information
    e1_Trackmom_calo  = sqrt(ele1->trackMomentumAtCalo().perp2());
    e1_Trackmom_vtx   = sqrt(ele1->trackMomentumAtVtx().perp2());
    //get Hcal energy over Ecal Energy
    e1_HoverE = ele1->hadronicOverEm();	  
    e1_EoverPout      = ele1->eSeedClusterOverPout();
    e1_EoverPin       = ele1->eSuperClusterOverP();
    //get SuperCluster (sc) infos
    reco::SuperClusterRef SCp = ele1->superCluster();
    e1_sc_x          = SCp->x();
    e1_sc_y          = SCp->y();
    e1_sc_z          = SCp->z();
    e1_sc_Theta      = SCp->position().Theta();
    e1_sc_Eta        = SCp->eta();    
    e1_sc_Phi        = SCp->phi();
    e1_sc_E          = SCp->energy();
    e1_sc_px         = e1E * sin(e1Theta) * cos(e1Phi);
    e1_sc_py         = e1E * sin(e1Theta) * sin(e1Phi);
    e1_sc_pz         = e1E * cos(e1Theta);
    e1_sc_Pt         = e1E / cosh(e1Eta);
    e1_sc_Et         = e1Pt;	  
    e1_SigmaEtaEta   = ele1->sigmaEtaEta();
    e1_SigmaIetaIeta = ele1->sigmaIetaIeta();
    e1_missingHits   = ele1->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
    e1_dist          = ele1->convDist();
    e1_dcot          = ele1->convDcot();
    e1_convradius    = ele1->convRadius();
    
    if( ele1->isEB() ) {
      ise1WP95      = (e1_missingHits<=1) && (e1_trackiso/e1Et<0.15) && (e1_ecaliso/e1Et<2.0) 
	&& (e1_hcaliso/e1Et<0.12) && (e1_SigmaIetaIeta<0.01) && (fabs(e1_DeltaPhiIn)<0.8) 
	&& (fabs(e1_DeltaEtaIn)<0.007) && (e1_HoverE<0.15);
      ise1WP80      = (e1_missingHits==0) && (fabs(e1_dist)>0.02 || fabs(e1_dcot)>0.02) 
	&& (e1_trackiso/e1Et<0.09) && (e1_ecaliso/e1Et<0.07) 
	&& (e1_hcaliso/e1Et<0.1) && (e1_SigmaIetaIeta<0.01) && (fabs(e1_DeltaPhiIn)<0.06) 
	&& (fabs(e1_DeltaEtaIn)<0.004) && (e1_HoverE<0.04);
      ise1WP70      = (e1_missingHits==0) && (fabs(e1_dist)>0.02 || fabs(e1_dcot)>0.02) 
	&& (e1_trackiso/e1Et<0.05) && (e1_ecaliso/e1Et<0.06) 
	&& (e1_hcaliso/e1Et<0.03) && (e1_SigmaIetaIeta<0.01) && (fabs(e1_DeltaPhiIn)<0.03) 
	&& (fabs(e1_DeltaEtaIn)<0.004) && (e1_HoverE<0.025);
    }
    if( ele1->isEE() ) {
      ise1WP95      = (e1_missingHits<=1) && (e1_trackiso/e1Et<0.08) && (e1_ecaliso/e1Et<0.06) 
	&& (e1_hcaliso/e1Et<0.05) && (e1_SigmaIetaIeta<0.03) && (fabs(e1_DeltaPhiIn)<0.7) 
	&& (fabs(e1_DeltaEtaIn)<0.01) && (e1_HoverE<0.07);
      ise1WP80      = (e1_missingHits==0)  && (fabs(e1_dist)>0.02 || fabs(e1_dcot)>0.02)  
	&& (e1_trackiso/e1Et<0.08) && (e1_ecaliso/e1Et<0.06) 
	&& (e1_hcaliso/e1Et<0.05) && (e1_SigmaIetaIeta<0.03) && (fabs(e1_DeltaPhiIn)<0.7) 
	&& (fabs(e1_DeltaEtaIn)<0.01) && (e1_HoverE<0.07);
      ise1WP70      = (e1_missingHits==0)  && (fabs(e1_dist)>0.02 || fabs(e1_dcot)>0.02)  
	&& (e1_trackiso/e1Et<0.025) && (e1_ecaliso/e1Et<0.025) 
	&& (e1_hcaliso/e1Et<0.02) && (e1_SigmaIetaIeta<0.03) && (fabs(e1_DeltaPhiIn)<0.02) 
	&& (fabs(e1_DeltaEtaIn)<0.005) && (e1_HoverE<0.025);
    }

  // IP relative to beam spot  & dz relative to vertex
   edm::Handle<reco::BeamSpot > beamSpot;
   if(runoverAOD){
     iEvent.getByLabel(mInputBeamSpot, beamSpot);
     e1_d0bsp = e1->gsfTrack()->dxy( beamSpot->position() ) ;
     e1_dz000 = e1->vertex().z(); }
   else{
     const pat::Electron* patel1 = dynamic_cast<const pat::Electron *>( &*ele1 );
	 e1_d0bsp = patel1->dB(pat::Electron::BS2D) ;
	 e1_dz000 = patel1->dB(pat::Electron::PV2D); 

         if(patel1->isElectronIDAvailable("mvaNonTrigV0"))
           e1_mvaNonTrigV0 =  patel1->electronID("mvaNonTrigV0");
         if(patel1->isElectronIDAvailable("mvaTrigV0"))
           e1_mvaTrigV0 =  patel1->electronID("mvaTrigV0");
   }
 	   
   // PF Isolation 
   e1_pfiso_chargedHadronIso = e1->pfIsolationVariables().chargedHadronIso;
   e1_pfiso_photonIso        = e1->pfIsolationVariables().photonIso;
   e1_pfiso_neutralHadronIso = e1->pfIsolationVariables().neutralHadronIso;
   e1_EffArea = ElectronEffectiveArea::GetElectronEffectiveArea( ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03 , e1_sc_Eta , ElectronEffectiveArea::kEleEAData2011);
   e1_pfIsoEA = (e1_pfiso_chargedHadronIso +
                      max(0.,e1_pfiso_neutralHadronIso +
                      e1_pfiso_photonIso  -
                      e1_EffArea*fastJetRho)) / e1Pt;

  }

  ////////// electron #2 quantities //////////////
  if( !(ele2 == NULL) ) {
    e2Charge          = ele2->charge();
    e2Vx              = ele2->vx();
    e2Vy              = ele2->vy();
    e2Vz              = ele2->vz();
    e2Y               = ele2->rapidity();
    e2Theta           = ele2->theta();
    e2Eta             = ele2->eta();    
    e2Phi             = ele2->phi();
    e2E               = ele2->energy();
    e2px              = ele2->px();
    e2py              = ele2->py();
    e2pz              = ele2->pz();
    e2Pt              = ele2->pt();
    e2Et              = ele2->et();	  

    /// isolation 
    e2_trackiso       = ele2->dr03TkSumPt();
    e2_ecaliso        = ele2->dr03EcalRecHitSumEt();
    e2_hcaliso        = ele2->dr03HcalTowerSumEt();
    e2Classification  = ele2->classification();
    e2_EoverPout = ele2->eSeedClusterOverPout();
    e2_EoverPin  = ele2->eSuperClusterOverP();
    e2_numberOfBrems  = ele2->numberOfBrems();
    e2_BremFraction   = ele2->fbrem();
    e2_DeltaEtaIn     = ele2->deltaEtaSuperClusterTrackAtVtx();
    e2_DeltaPhiIn     = ele2->deltaPhiSuperClusterTrackAtVtx();
    e2_DeltaPhiOut    = ele2->deltaPhiSeedClusterTrackAtCalo();
    e2_DeltaEtaOut    = ele2->deltaEtaSeedClusterTrackAtCalo();
    //Track Momentum information
    e2_Trackmom_calo  = sqrt(ele2->trackMomentumAtCalo().perp2());
    e2_Trackmom_vtx   = sqrt(ele2->trackMomentumAtVtx().perp2());   
    //get Hcal energy over Ecal Energy
    e2_HoverE = ele2->hadronicOverEm();	  
    //get SuperCluster (sc) infos
    reco::SuperClusterRef SCm = ele2->superCluster();
    e2_sc_x          = SCm->x();
    e2_sc_y          = SCm->y();
    e2_sc_z          = SCm->z();
    e2_sc_Theta      = SCm->position().Theta();
    e2_sc_Eta        = SCm->eta();    
    e2_sc_Phi        = SCm->phi();
    e2_sc_E          = SCm->energy();
    e2_sc_px         = e2E * sin(e2Theta) * cos(e2Phi);
    e2_sc_py         = e2E * sin(e2Theta) * sin(e2Phi);
    e2_sc_pz         = e2E * cos(e2Theta);
    e2_sc_Pt         = e2E / cosh(e2Eta);
    e2_sc_Et         = e2Pt;	 
    e2_SigmaEtaEta   = ele2->sigmaEtaEta();
    e2_SigmaIetaIeta = ele2->sigmaIetaIeta();
    e2_missingHits   = ele2->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
    e2_dist          = ele2->convDist();
    e2_dcot          = ele2->convDcot();
    e2_convradius    = ele2->convRadius();
    if( ele2->isEB() ) {
      ise2WP95      = (e2_missingHits<=1) && (e2_trackiso/e2Et<0.15) && (e2_ecaliso/e2Et<2.0) 
	&& (e2_hcaliso/e2Et<0.12) && (e2_SigmaIetaIeta<0.01) && (fabs(e2_DeltaPhiIn)<0.8) 
	&& (fabs(e2_DeltaEtaIn)<0.007) && (e2_HoverE<0.15);
      ise2WP80      = (e2_missingHits==0) && (fabs(e2_dist)>0.02 || fabs(e2_dcot)>0.02) 
	&& (e2_trackiso/e2Et<0.09) && (e2_ecaliso/e2Et<0.07) 
	&& (e2_hcaliso/e2Et<0.1) && (e2_SigmaIetaIeta<0.01) && (fabs(e2_DeltaPhiIn)<0.06) 
	&& (fabs(e2_DeltaEtaIn)<0.004) && (e2_HoverE<0.04);
      ise2WP70      = (e2_missingHits==0) && (fabs(e2_dist)>0.02 || fabs(e2_dcot)>0.02) 
	&& (e2_trackiso/e2Et<0.05) && (e2_ecaliso/e2Et<0.06) 
	&& (e2_hcaliso/e2Et<0.03) && (e2_SigmaIetaIeta<0.01) && (fabs(e2_DeltaPhiIn)<0.03) 
	&& (fabs(e2_DeltaEtaIn)<0.004) && (e2_HoverE<0.025);
    }
    if( ele2->isEE() ) {
      ise2WP95      = (e2_missingHits<=1) && (e2_trackiso/e2Et<0.08) && (e2_ecaliso/e2Et<0.06) 
	&& (e2_hcaliso/e2Et<0.05) && (e2_SigmaIetaIeta<0.03) && (fabs(e2_DeltaPhiIn)<0.7) 
	&& (fabs(e2_DeltaEtaIn)<0.01) && (e2_HoverE<0.07);
      ise2WP80      = (e2_missingHits==0)  && (fabs(e2_dist)>0.02 || fabs(e2_dcot)>0.02)  
	&& (e2_trackiso/e2Et<0.08) && (e2_ecaliso/e2Et<0.06) 
	&& (e2_hcaliso/e2Et<0.05) && (e2_SigmaIetaIeta<0.03) && (fabs(e2_DeltaPhiIn)<0.7) 
	&& (fabs(e2_DeltaEtaIn)<0.01) && (e2_HoverE<0.07);
      ise2WP70      = (e2_missingHits==0)  && (fabs(e2_dist)>0.02 || fabs(e2_dcot)>0.02)  
	&& (e2_trackiso/e2Et<0.025) && (e2_ecaliso/e2Et<0.025) 
	&& (e2_hcaliso/e2Et<0.02) && (e2_SigmaIetaIeta<0.03) && (fabs(e2_DeltaPhiIn)<0.02) 
	&& (fabs(e2_DeltaEtaIn)<0.005) && (e2_HoverE<0.025);
    }
  } 

}





// WP80: but only track iso, and cut on impact parameter |d0|
bool ewk::VtoElectronTreeFiller::isTightElectron(const edm::Event& iEvent, const reco::GsfElectron& ele) 
{
  float pt = ele.pt();
  float eta = fabs(ele.superCluster()->eta());  

  if( pt<20.0 )  return false;
  if( eta > 1.4442 && eta < 1.5660) return false;
  if ( eta > 2.5 ) return false;


  int missingHits   = ele.gsfTrack()->trackerExpectedHitsInner().numberOfHits();
  float dist          = fabs(ele.convDist());
  float dcot          = fabs(ele.convDcot());

  if( !(missingHits==0) )  return false;
  if( dist<0.02 && dcot<0.02) return false;


  float trackiso      = ele.dr03TkSumPt()/pt;
  float sigmaIetaIeta = ele.sigmaIetaIeta();
  float deltaEta      = ele.deltaEtaSuperClusterTrackAtVtx();
  float deltaPhi      = ele.deltaPhiSuperClusterTrackAtVtx();
  float HoverE        = ele.hadronicOverEm();	  
  bool isEB = ele.isEB();
  bool isEE = ele.isEE();


  bool isWP80Id    = (trackiso<0.1) && 
    ((isEB && sigmaIetaIeta<0.01 && deltaPhi<0.06 && deltaEta<0.004 && HoverE<0.04) || 
     (isEE && sigmaIetaIeta<0.03 && deltaPhi<0.03 && deltaEta<0.007 && HoverE<0.025)); 


  //////////// Beam spot //////////////
  edm::Handle<reco::BeamSpot> beamSpot;
  double dz =0.0;
	if(runoverAOD){
		iEvent.getByLabel(mInputBeamSpot, beamSpot);
		dz = fabs( ele.gsfTrack()->dxy( beamSpot->position() ) );

	}
  if( !(isWP80Id && dz<0.02) )  return false;

  return true;
}





// WP95: but only track iso
bool ewk::VtoElectronTreeFiller::isLooseElectron(const edm::Event& iEvent, const reco::GsfElectron& ele) 
{
  float pt = ele.pt();
  float eta = fabs(ele.superCluster()->eta());  

  if( pt<20.0 )  return false;
  if( eta > 1.4442 && eta < 1.5660) return false;
  if ( eta > 2.5 ) return false;


  float trackiso      = ele.dr03TkSumPt()/pt;
  float sigmaIetaIeta = ele.sigmaIetaIeta();
  float deltaEta      = ele.deltaEtaSuperClusterTrackAtVtx();
  float deltaPhi      = ele.deltaPhiSuperClusterTrackAtVtx();
  float HoverE        = ele.hadronicOverEm();	  
  bool isEB = ele.isEB();
  bool isEE = ele.isEE();


  bool isWP95Id    = (trackiso<0.2) && 
    ((isEB && sigmaIetaIeta<0.01 && deltaPhi<0.8 && deltaEta<0.007 && HoverE<0.15) || 
     (isEE && sigmaIetaIeta<0.03 && deltaPhi<0.7 && deltaEta<0.010 && HoverE<0.07)); 

  if( !(isWP95Id) )  return false;

  return true;
}





////////////////// utilities, helpers ///////////////////
 
void ewk::VtoElectronTreeFiller::SetBranch( float* x, std::string name)
{
  std::string brName = std::string(name_) + "_" + name;
  tree_->Branch( brName.c_str(), x, ( brName+"/F").c_str() );
}


void ewk::VtoElectronTreeFiller::SetBranch( int* x, std::string name)
{
  std::string brName = std::string(name_) + "_" + name;
  tree_->Branch( brName.c_str(), x, ( brName+"/I").c_str() );
}


void ewk::VtoElectronTreeFiller::SetBranch( bool* x, std::string name)
{
  std::string brName = std::string(name_) + "_" + name;
  tree_->Branch( brName.c_str(), x, ( brName+"/O").c_str() );
}

