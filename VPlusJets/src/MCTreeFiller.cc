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
 *   To fill W/Z related MC information into a specified TTree
 *   Can work with jets in RECO/AOD/PAT data formats.
 * History:
 *   
 *
 * Copyright (C) 2010 FNAL 
 *****************************************************************************/

// CMS includes
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TMath.h" 
// Header file
#include "ElectroWeakAnalysis/VPlusJets/interface/MCTreeFiller.h"


ewk::MCTreeFiller::MCTreeFiller(const char *name, TTree* tree, 
			const edm::ParameterSet iConfig)
{
	// ********** Vector boson ********** //
	tree_     = tree;
	name_     = name;
	Vtype_    = iConfig.getParameter<std::string>("VBosonType"); 
	ptype_    = iConfig.getParameter<std::string>("LeptonType");
	pdgIdDau_ = 11;
	if(  iConfig.existsAs<edm::InputTag>("srcGenParticles") )
	  mInputgenParticles  = iConfig.getParameter<edm::InputTag>("srcGenParticles");
	if(ptype_=="muon") pdgIdDau_ = 13; 

	if( !(tree==0) ) SetBranches();
}





void ewk::MCTreeFiller::SetBranches()
{
	//	std::cout << "srcgenparticles " << mInputgenParticles << std::endl;
	// Declare jet branches
	std::string lept1;
	std::string lept2;
	if( Vtype_=="Z" ) {
		if(ptype_=="muon") {
			lept1 = "muplus";
			lept2 = "muminus";
		} else {
			lept1 = "eplus";
			lept2 = "eminus";
		} 
	} else {
		if(ptype_=="muon") {
			lept1 = "muon";
			lept2 = "neutrino";
		} else {
			lept1 = "electron";
			lept2 = "neutrino";
		}
	}

	SetBranch( &Photon_pt_gen,      "Photon_pt_gen");
	//SetBranch( &H_mass,      "H_mass_gen");
	//SetBranch( &H_px,        "H_px_gen");
	//SetBranch( &H_py,        "H_py_gen");
	//SetBranch( &H_pz,        "H_pz_gen");
	//SetBranch( &H_E,         "H_e_gen");
	//SetBranch( &H_Pt,        "H_pt_gen");
	//SetBranch( &H_Et,        "H_et_gen");
	//SetBranch( &H_Eta,       "H_eta_gen");    
	//SetBranch( &H_Phi,       "H_phi_gen");
	//SetBranch( &H_Vx,        "H_vx_gen");
	//SetBranch( &H_Vy,        "H_vy_gen");
	//SetBranch( &H_Vz,        "H_vz_gen");
	//SetBranch( &H_Y,         "H_y_gen");
	//SetBranch( &H_Id,         "H_Id_gen");


	if( Vtype_=="Z" ) {
	    ///////////////////////////////////////////////
	    SetBranch( &l1px,             lept1+"_px_gen" );
	    SetBranch( &l1py,             lept1+"_py_gen" );
	    SetBranch( &l1pz,             lept1+"_pz_gen" );
	    SetBranch( &l1E,              lept1+"_e_gen" );
	    SetBranch( &l1Pt,             lept1+"_pt_gen" );
	    SetBranch( &l1Et,             lept1+"_et_gen" );
	    SetBranch( &l1Eta,            lept1+"_eta_gen" ); 
	    SetBranch( &l1Theta,          lept1+"_theta_gen" ); 
	    SetBranch( &l1Phi,            lept1+"_phi_gen" );
	    SetBranch( &l1Charge,         lept1+"_charge_gen" );
	    SetBranch( &l1Vx,             lept1+"_vx_gen" );
	    SetBranch( &l1Vy,             lept1+"_vy_gen" );
	    SetBranch( &l1Vz,             lept1+"_vz_gen" );
	    SetBranch( &l1Y,              lept1+"_y_gen" );

	    ////////////////////////////////////////////////////////
	    SetBranch( &l2px,             lept2+"_px_gen" );
	    SetBranch( &l2py,             lept2+"_py_gen" );
	    SetBranch( &l2pz,             lept2+"_pz_gen" );
	    SetBranch( &l2E,              lept2+"_e_gen" );
	    SetBranch( &l2Pt,             lept2+"_pt_gen" );
	    SetBranch( &l2Et,             lept2+"_et_gen" );
	    SetBranch( &l2Eta,            lept2+"_eta_gen" ); 
	    SetBranch( &l2Theta,          lept2+"_theta_gen" );    
	    SetBranch( &l2Phi,            lept2+"_phi_gen" );
	    SetBranch( &l2Charge,         lept2+"_charge_gen" );
	    SetBranch( &l2Vx,             lept2+"_vx_gen" );
	    SetBranch( &l2Vy,             lept2+"_vy_gen" );
	    SetBranch( &l2Vz,             lept2+"_vz_gen" );
	    SetBranch( &l2Y,              lept2+"_y_gen" );
	}

	/*////////////////////////////////////////////////////////
	SetBranch( Parton_px,             "Parton_px[2]" );
	SetBranch( Parton_py,             "Parton_py[2]" );
	SetBranch( Parton_pz,             "Parton_pz[2]" );
	SetBranch( Parton_E,              "Parton_E[2]" );
	SetBranch( Parton_Pt,             "Parton_pt[2]" );
	SetBranch( Parton_Et,             "Parton_et[2]" );
	SetBranch( Parton_Eta,            "Parton_eta[2]" ); 
	SetBranch( Parton_Theta,          "Parton_theta[2]" );    
	SetBranch( Parton_Phi,            "Parton_phi[2]" );
	SetBranch( Parton_Charge,         "Parton_charge[2]" );
	SetBranch( Parton_Vx,             "Parton_vx[2]" );
	SetBranch( Parton_Vy,             "Parton_vy[2]" );
	SetBranch( Parton_Vz,             "Parton_vz[2]" );
	SetBranch( Parton_Y,              "Parton_y[2]" );
	SetBranch( Parton_Id,              "Parton_Id[2]" );
    */

	////////////////////////////////////////////////////////
	if( Vtype_=="W" ) {
	    SetBranch( &Lepton_px,             "Lepton_px" );
	    SetBranch( &Lepton_py,             "Lepton_py" );
	    SetBranch( &Lepton_pz,             "Lepton_pz" );
	    SetBranch( &Lepton_E,              "Lepton_E" );
	    SetBranch( &Lepton_Pt,             "Lepton_pt" );
	    SetBranch( &Lepton_Et,             "Lepton_et" );
	    SetBranch( &Lepton_Eta,            "Lepton_eta" );    
	    SetBranch( &Lepton_Theta,          "Lepton_theta" );    
	    SetBranch( &Lepton_Phi,            "Lepton_phi" );
	    SetBranch( &Lepton_Charge,         "Lepton_charge" );
	    SetBranch( &Lepton_Vx,             "Lepton_vx" );
	    SetBranch( &Lepton_Vy,             "Lepton_vy" );
	    SetBranch( &Lepton_Vz,             "Lepton_vz" );
	    SetBranch( &Lepton_Y,              "Lepton_y" );
	    SetBranch( &Lepton_Id,              "Lepton_Id" );

	    ////////////////////////////////////////////////////////
	    SetBranch( &Met_px,             "Met_px" );
	    SetBranch( &Met_py,             "Met_py" );
	    SetBranch( &Met_pz,             "Met_pz" );
	    SetBranch( &Met_E,              "Met_E" );
	    SetBranch( &Met_Pt,             "Met_pt" );
	    SetBranch( &Met_Et,             "Met_et" );
	    SetBranch( &Met_Eta,            "Met_eta" );    
	    SetBranch( &Met_Theta,          "Met_theta" );    
	    SetBranch( &Met_Phi,            "Met_phi" );
	    SetBranch( &Met_Charge,         "Met_charge" );
	    SetBranch( &Met_Vx,             "Met_vx" );
	    SetBranch( &Met_Vy,             "Met_vy" );
	    SetBranch( &Met_Vz,             "Met_vz" );
	    SetBranch( &Met_Y,              "Met_y" );
	    SetBranch( &Met_Id,              "Met_Id" );
		/////////////////////////////////////////////////////////
	    SetBranch( &nParton_Winclusive,              "nParton_Winclusive" );
	}



	////////////////////////////////////////////////////////
	/*SetBranch( tParton_px,             "tParton_px[2]" );
	SetBranch( tParton_py,             "tParton_py[2]" );
	SetBranch( tParton_pz,             "tParton_pz[2]" );
	SetBranch( tParton_E,              "tParton_E[2]" );
	SetBranch( tParton_Pt,             "tParton_pt[2]" );
	SetBranch( tParton_Et,             "tParton_et[2]" );
	SetBranch( tParton_Eta,            "tParton_eta[2]" ); 
	SetBranch( tParton_Theta,          "tParton_theta[2]" );    
	SetBranch( tParton_Phi,            "tParton_phi[2]" );
	SetBranch( tParton_Charge,         "tParton_charge[2]" );
	SetBranch( tParton_Vx,             "tParton_vx[2]" );
	SetBranch( tParton_Vy,             "tParton_vy[2]" );
	SetBranch( tParton_Vz,             "tParton_vz[2]" );
	SetBranch( tParton_Y,              "tParton_y[2]" );
	SetBranch( tParton_Id,              "tParton_Id[2]" );


	////////////////////////////////////////////////////////
	SetBranch( &tLepton_px,             "tLepton_px" );
	SetBranch( &tLepton_py,             "tLepton_py" );
	SetBranch( &tLepton_pz,             "tLepton_pz" );
	SetBranch( &tLepton_E,              "tLepton_E" );
	SetBranch( &tLepton_Pt,             "tLepton_pt" );
	SetBranch( &tLepton_Et,             "tLepton_et" );
	SetBranch( &tLepton_Eta,            "tLepton_eta" );    
	SetBranch( &tLepton_Theta,          "tLepton_theta" );    
	SetBranch( &tLepton_Phi,            "tLepton_phi" );
	SetBranch( &tLepton_Charge,         "tLepton_charge" );
	SetBranch( &tLepton_Vx,             "tLepton_vx" );
	SetBranch( &tLepton_Vy,             "tLepton_vy" );
	SetBranch( &tLepton_Vz,             "tLepton_vz" );
	SetBranch( &tLepton_Y,              "tLepton_y" );
	SetBranch( &tLepton_Id,              "tLepton_Id" );

	////////////////////////////////////////////////////////
	SetBranch( &tMet_px,             "tMet_px" );
	SetBranch( &tMet_py,             "tMet_py" );
	SetBranch( &tMet_pz,             "tMet_pz" );
	SetBranch( &tMet_E,              "tMet_E" );
	SetBranch( &tMet_Pt,             "tMet_pt" );
	SetBranch( &tMet_Et,             "tMet_et" );
	SetBranch( &tMet_Eta,            "tMet_eta" );    
	SetBranch( &tMet_Theta,          "tMet_theta" );    
	SetBranch( &tMet_Phi,            "tMet_phi" );
	SetBranch( &tMet_Charge,         "tMet_charge" );
	SetBranch( &tMet_Vx,             "tMet_vx" );
	SetBranch( &tMet_Vy,             "tMet_vy" );
	SetBranch( &tMet_Vz,             "tMet_vz" );
	SetBranch( &tMet_Y,              "tMet_y" );
	SetBranch( &tMet_Id,              "tMet_Id" );

	SetBranch( &tb_px,             "tb_px" );
	SetBranch( &tb_py,             "tb_py" );
	SetBranch( &tb_pz,             "tb_pz" );
	SetBranch( &tb_E,              "tb_E" );
	SetBranch( &tb_Pt,             "tb_pt" );
	SetBranch( &tb_Et,             "tb_et" );
	SetBranch( &tb_Eta,            "tb_eta" );
	SetBranch( &tb_Theta,          "tb_theta" );
	SetBranch( &tb_Phi,            "tb_phi" );
	SetBranch( &tb_Charge,         "tb_charge" );
	SetBranch( &tb_Vx,             "tb_vx" );
	SetBranch( &tb_Vy,             "tb_vy" );
	SetBranch( &tb_Vz,             "tb_vz" );
	SetBranch( &tb_Y,              "tb_y" );
	SetBranch( &tb_Id,              "tb_Id" );

	SetBranch( &tbbar_px,             "tbbar_px" );
	SetBranch( &tbbar_py,             "tbbar_py" );
	SetBranch( &tbbar_pz,             "tbbar_pz" );
	SetBranch( &tbbar_E,              "tbbar_E" );
	SetBranch( &tbbar_Pt,             "tbbar_pt" );
	SetBranch( &tbbar_Et,             "tbbar_et" );
	SetBranch( &tbbar_Eta,            "tbbar_eta" );
	SetBranch( &tbbar_Theta,          "tbbar_theta" );
	SetBranch( &tbbar_Phi,            "tbbar_phi" );
	SetBranch( &tbbar_Charge,         "tbbar_charge" );
	SetBranch( &tbbar_Vx,             "tbbar_vx" );
	SetBranch( &tbbar_Vy,             "tbbar_vy" );
	SetBranch( &tbbar_Vz,             "tbbar_vz" );
	SetBranch( &tbbar_Y,              "tbbar_y" );
	SetBranch( &tbbar_Id,              "tbbar_Id" );

	SetBranch( &Hb_px,             "Hb_px" );
	SetBranch( &Hb_py,             "Hb_py" );
	SetBranch( &Hb_pz,             "Hb_pz" );
	SetBranch( &Hb_E,              "Hb_E" );
	SetBranch( &Hb_Pt,             "Hb_pt" );
	SetBranch( &Hb_Et,             "Hb_et" );
	SetBranch( &Hb_Eta,            "Hb_eta" );
	SetBranch( &Hb_Theta,          "Hb_theta" );
	SetBranch( &Hb_Phi,            "Hb_phi" );
	SetBranch( &Hb_Charge,         "Hb_charge" );
	SetBranch( &Hb_Vx,             "Hb_vx" );
	SetBranch( &Hb_Vy,             "Hb_vy" );
	SetBranch( &Hb_Vz,             "Hb_vz" );
	SetBranch( &Hb_Y,              "Hb_y" );
	SetBranch( &Hb_Id,              "Hb_Id" );

	SetBranch( &Hbbar_px,             "Hbbar_px" );
	SetBranch( &Hbbar_py,             "Hbbar_py" );
	SetBranch( &Hbbar_pz,             "Hbbar_pz" );
	SetBranch( &Hbbar_E,              "Hbbar_E" );
	SetBranch( &Hbbar_Pt,             "Hbbar_pt" );
	SetBranch( &Hbbar_Et,             "Hbbar_et" );
	SetBranch( &Hbbar_Eta,            "Hbbar_eta" );
	SetBranch( &Hbbar_Theta,          "Hbbar_theta" );
	SetBranch( &Hbbar_Phi,            "Hbbar_phi" );
	SetBranch( &Hbbar_Charge,         "Hbbar_charge" );
	SetBranch( &Hbbar_Vx,             "Hbbar_vx" );
	SetBranch( &Hbbar_Vy,             "Hbbar_vy" );
	SetBranch( &Hbbar_Vz,             "Hbbar_vz" );
	SetBranch( &Hbbar_Y,              "Hbbar_y" );
	SetBranch( &Hbbar_Id,              "Hbbar_Id" );

	/////////////////////////////////////////////////
	SetBranch( TagQuark_px,             "TagQuark_px[2]" );
	SetBranch( TagQuark_py,             "TagQuark_py[2]" );
	SetBranch( TagQuark_pz,             "TagQuark_pz[2]" );
	SetBranch( TagQuark_E,              "TagQuark_E[2]" );
	SetBranch( TagQuark_Pt,             "TagQuark_pt[2]" );
	SetBranch( TagQuark_Et,             "TagQuark_et[2]" );
	SetBranch( TagQuark_Eta,            "TagQuark_eta[2]" ); 
	SetBranch( TagQuark_Theta,          "TagQuark_theta[2]" );    
	SetBranch( TagQuark_Phi,            "TagQuark_phi[2]" );
	SetBranch( TagQuark_Charge,         "TagQuark_charge[2]" );
	SetBranch( TagQuark_Vx,             "TagQuark_vx[2]" );
	SetBranch( TagQuark_Vy,             "TagQuark_vy[2]" );
	SetBranch( TagQuark_Vz,             "TagQuark_vz[2]" );
	SetBranch( TagQuark_Y,              "TagQuark_y[2]" );
	SetBranch( TagQuark_Id,              "TagQuark_Id[2]" );

	////////////////////////////////////////////////

	SetBranch( &EWKW_px,             "EWKW_px" );
	SetBranch( &EWKW_py,             "EWKW_py" );
	SetBranch( &EWKW_pz,             "EWKW_pz" );
	SetBranch( &EWKW_E,              "EWKW_E" );
	SetBranch( &EWKW_Pt,             "EWKW_pt" );
	SetBranch( &EWKW_Et,             "EWKW_et" );
	SetBranch( &EWKW_Eta,            "EWKW_eta" ); 
	SetBranch( &EWKW_Theta,          "EWKW_theta" );    
	SetBranch( &EWKW_Phi,            "EWKW_phi" );
	SetBranch( &EWKW_Charge,         "EWKW_charge" );
	SetBranch( &EWKW_Vx,             "EWKW_vx" );
	SetBranch( &EWKW_Vy,             "EWKW_vy" );
	SetBranch( &EWKW_Vz,             "EWKW_vz" );
	SetBranch( &EWKW_Y,              "EWKW_y" );
	SetBranch( &EWKW_Id,              "EWKW_Id" );

	SetBranch( EWKTagQuark_px,             "EWKTagQuark_px[2]" );
	SetBranch( EWKTagQuark_py,             "EWKTagQuark_py[2]" );
	SetBranch( EWKTagQuark_pz,             "EWKTagQuark_pz[2]" );
	SetBranch( EWKTagQuark_E,              "EWKTagQuark_E[2]" );
	SetBranch( EWKTagQuark_Pt,             "EWKTagQuark_pt[2]" );
	SetBranch( EWKTagQuark_Et,             "EWKTagQuark_et[2]" );
	SetBranch( EWKTagQuark_Eta,            "EWKTagQuark_eta[2]" ); 
	SetBranch( EWKTagQuark_Theta,          "EWKTagQuark_theta[2]" );    
	SetBranch( EWKTagQuark_Phi,            "EWKTagQuark_phi[2]" );
	SetBranch( EWKTagQuark_Charge,         "EWKTagQuark_charge[2]" );
	SetBranch( EWKTagQuark_Vx,             "EWKTagQuark_vx[2]" );
	SetBranch( EWKTagQuark_Vy,             "EWKTagQuark_vy[2]" );
	SetBranch( EWKTagQuark_Vz,             "EWKTagQuark_vz[2]" );
	SetBranch( EWKTagQuark_Y,              "EWKTagQuark_y[2]" );
	SetBranch( EWKTagQuark_Id,              "EWKTagQuark_Id[2]" );
    */
}
/////////////////////////////////////////////////////////////////////////






void ewk::MCTreeFiller::init()   
{
	// initialize private data members
	Photon_pt_gen         = -1. ;
	H_mass                  = -1.;
	H_px                  = -99999.;
	H_py                  = -99999.;
	H_pz                  = -99999.;
	H_E                   = -1.;
	H_Pt                  = -1.;
	H_Et                  = -1.;
	H_Eta                 = -10.;
	H_Phi                 = -10.;
	H_Vx                  = -10.;
	H_Vy                  = -10.;
	H_Vz                  = -10.;
	H_Y                   = -10.;
	H_Id                   = 0;

	l1Charge           = -10;
	l2Charge          = -10;

	l1px               = -99999.;
	l1py               = -99999.;
	l1pz               = -99999.;
	l1E                = -1.;
	l1Et               = -1.;
	l1Pt               = -1.;
	l1Eta              = -10.;
	l1Theta            = -99999.;
	l1Phi              = -10.;
	l1Vx               = -10.;
	l1Vy               = -10.;
	l1Vz               = -10.;
	l1Y                = -10.;	  

	l2px              = -99999.;
	l2py              = -99999.;
	l2pz              = -99999.;
	l2E               = -1.;
	l2Pt              = -1.;
	l2Et              = -1.;
	l2Eta             = -10.;
	l2Theta           = -99999.;
	l2Phi             = -10.;
	l2Vx              = -10.;
	l2Vy              = -10.;
	l2Vz              = -10.;
	l2Y               = -10.;

	for ( int i =0; i<2; i++){
		Parton_px[i]              = -99999.;
		Parton_py[i]              = -99999.;
		Parton_pz[i]              = -99999.;
		Parton_E[i]               = -1.;
		Parton_Pt[i]              = -1.;
		Parton_Et[i]              = -1.;
		Parton_Eta[i]             = -10.;
		Parton_Theta[i]           = -99999.;
		Parton_Phi[i]             = -10.;
		Parton_Vx[i]              = -10.;
		Parton_Vy[i]              = -10.;
		Parton_Vz[i]              = -10.;
		Parton_Y[i]              =  -10.;
		Parton_Id[i]              = 0;
	}

	Lepton_px              = -99999.;
	Lepton_py              = -99999.;
	Lepton_pz              = -99999.;
	Lepton_E               = -1.;
	Lepton_Pt              = -1.;
	Lepton_Et              = -1.;
	Lepton_Eta             = -10.;
	Lepton_Charge          = -10;
	Lepton_Theta           = -99999.;
	Lepton_Phi             = -10.;
	Lepton_Vx              = -10.;
	Lepton_Vy              = -10.;
	Lepton_Vz              = -10.;
	Lepton_Y               = -10.;
	Lepton_Id              = 0;

	Met_px              = -99999.;
	Met_py              = -99999.;
	Met_pz              = -99999.;
	Met_E               = -1.;
	Met_Pt              = -1.;
	Met_Et              = -1.;
	Met_Eta             = -10.;
	Met_Charge          = -10;
	Met_Theta           = -99999.;
	Met_Phi             = -10.;
	Met_Vx              = -10.;
	Met_Vy              = -10.;
	Met_Vz              = -10.;
	Met_Y               = -10.;
	Met_Id              = 0;

	for ( int i =0; i<2; i++){
		tParton_px[i]              = -99999.;
		tParton_py[i]              = -99999.;
		tParton_pz[i]              = -99999.;
		tParton_E[i]               = -1.;
		tParton_Pt[i]              = -1.;
		tParton_Et[i]              = -1.;
		tParton_Eta[i]             = -10.;
		tParton_Theta[i]           = -99999.;
		tParton_Phi[i]             = -10.;
		tParton_Vx[i]              = -10.;
		tParton_Vy[i]              = -10.;
		tParton_Vz[i]              = -10.;
		tParton_Y[i]              =  -10.;
		tParton_Id[i]              = 0;
	}

	tLepton_px              = -99999.;
	tLepton_py              = -99999.;
	tLepton_pz              = -99999.;
	tLepton_E               = -1.;
	tLepton_Pt              = -1.;
	tLepton_Et              = -1.;
	tLepton_Eta             = -10.;
	tLepton_Charge          = -10;
	tLepton_Theta           = -99999.;
	tLepton_Phi             = -10.;
	tLepton_Vx              = -10.;
	tLepton_Vy              = -10.;
	tLepton_Vz              = -10.;
	tLepton_Y               = -10.;
	tLepton_Id              = 0;

	tMet_px              = -99999.;
	tMet_py              = -99999.;
	tMet_pz              = -99999.;
	tMet_E               = -1.;
	tMet_Pt              = -1.;
	tMet_Et              = -1.;
	tMet_Eta             = -10.;
	tMet_Charge          = -10;
	tMet_Theta           = -99999.;
	tMet_Phi             = -10.;
	tMet_Vx              = -10.;
	tMet_Vy              = -10.;
	tMet_Vz              = -10.;
	tMet_Y               = -10.;
	tMet_Id              = 0;


	tb_px              = -99999.;
	tb_py              = -99999.;
	tb_pz              = -99999.;
	tb_E               = -1.;
	tb_Pt              = -1.;
	tb_Et              = -1.;
	tb_Eta             = -10.;
	tb_Charge          = -10;
	tb_Theta           = -99999.;
	tb_Phi             = -10.;
	tb_Vx              = -10.;
	tb_Vy              = -10.;
	tb_Vz              = -10.;
	tb_Y               = -10.;
	tb_Id              = 0;


	tbbar_px              = -99999.;
	tbbar_py              = -99999.;
	tbbar_pz              = -99999.;
	tbbar_E               = -1.;
	tbbar_Pt              = -1.;
	tbbar_Et              = -1.;
	tbbar_Eta             = -10.;
	tbbar_Charge          = -10;
	tbbar_Theta           = -99999.;
	tbbar_Phi             = -10.;
	tbbar_Vx              = -10.;
	tbbar_Vy              = -10.;
	tbbar_Vz              = -10.;
	tbbar_Y               = -10.;
	tbbar_Id              = 0;

	Hb_px              = -99999.;
	Hb_py              = -99999.;
	Hb_pz              = -99999.;
	Hb_E               = -1.;
	Hb_Pt              = -1.;
	Hb_Et              = -1.;
	Hb_Eta             = -10.;
	Hb_Charge          = -10;
	Hb_Theta           = -99999.;
	Hb_Phi             = -10.;
	Hb_Vx              = -10.;
	Hb_Vy              = -10.;
	Hb_Vz              = -10.;
	Hb_Y               = -10.;
	Hb_Id              = 0;


	Hbbar_px              = -99999.;
	Hbbar_py              = -99999.;
	Hbbar_pz              = -99999.;
	Hbbar_E               = -1.;
	Hbbar_Pt              = -1.;
	Hbbar_Et              = -1.;
	Hbbar_Eta             = -10.;
	Hbbar_Charge          = -10;
	Hbbar_Theta           = -99999.;
	Hbbar_Phi             = -10.;
	Hbbar_Vx              = -10.;
	Hbbar_Vy              = -10.;
	Hbbar_Vz              = -10.;
	Hbbar_Y               = -10.;
	Hbbar_Id              = 0;


	for ( int i =0; i<2; i++){
		TagQuark_px[i]              = -99999.;
		TagQuark_py[i]              = -99999.;
		TagQuark_pz[i]              = -99999.;
		TagQuark_E[i]               = -1.;
		TagQuark_Pt[i]              = -1.;
		TagQuark_Et[i]              = -1.;
		TagQuark_Eta[i]             = -10.;
		TagQuark_Theta[i]           = -99999.;
		TagQuark_Phi[i]             = -10.;
		TagQuark_Vx[i]              = -10.;
		TagQuark_Vy[i]              = -10.;
		TagQuark_Vz[i]              = -10.;
		TagQuark_Y[i]              =  -10.;
		TagQuark_Id[i]              = 0;
	}


	EWKW_px              = -99999.;
	EWKW_py              = -99999.;
	EWKW_pz              = -99999.;
	EWKW_E               = -1.;
	EWKW_Pt              = -1.;
	EWKW_Et              = -1.;
	EWKW_Eta             = -10.;
	EWKW_Theta           = -99999.;
	EWKW_Phi             = -10.;
	EWKW_Vx              = -10.;
	EWKW_Vy              = -10.;
	EWKW_Vz              = -10.;
	EWKW_Y              =  -10.;
	EWKW_Id              = 0;

	for ( int i =0; i<2; i++){
		EWKTagQuark_px[i]              = -99999.;
		EWKTagQuark_py[i]              = -99999.;
		EWKTagQuark_pz[i]              = -99999.;
		EWKTagQuark_E[i]               = -1.;
		EWKTagQuark_Pt[i]              = -1.;
		EWKTagQuark_Et[i]              = -1.;
		EWKTagQuark_Eta[i]             = -10.;
		EWKTagQuark_Theta[i]           = -99999.;
		EWKTagQuark_Phi[i]             = -10.;
		EWKTagQuark_Vx[i]              = -10.;
		EWKTagQuark_Vy[i]              = -10.;
		EWKTagQuark_Vz[i]              = -10.;
		EWKTagQuark_Y[i]              =  -10.;
		EWKTagQuark_Id[i]              = 0;
	}

	nParton_Winclusive = -4; //Incoming Particles Calculated Twice

	// initialization done
}



void ewk::MCTreeFiller::fill(const edm::Event& iEvent)
{
	// first initialize to the default values
	init();


	edm::Handle<reco::GenParticleCollection> genParticles;
	iEvent.getByLabel(mInputgenParticles, genParticles);

	size_t nGen = genParticles->size();
	if( nGen < 1 ) return; // Nothing to fill


	// now iterate over the daughters  
	const reco::Candidate *V=NULL;
	const reco::Candidate* lepton1=NULL;
	const reco::Candidate* lepton2=NULL;

	const reco::Candidate *H=NULL;
	const reco::Candidate *Lepton=NULL;
	const reco::Candidate *Parton1=NULL;
	const reco::Candidate *Parton2=NULL;
	const reco::Candidate *Met=NULL;

	const reco::Candidate *t=NULL;
	const reco::Candidate *tbar=NULL;
	const reco::Candidate *aH=NULL;

	//  const reco::Candidate *ab=NULL;
	//  const reco::Candidate *abbar=NULL;

	const reco::Candidate *tb=NULL;
	const reco::Candidate *tbbar=NULL;

	const reco::Candidate *Hb=NULL;
	const reco::Candidate *Hbbar=NULL;

	const reco::Candidate *tMet=NULL;
	const reco::Candidate *tLepton=NULL;
	const reco::Candidate *tParton1=NULL;
	const reco::Candidate *tParton2=NULL;

	const reco::Candidate *TagQuark1=NULL;
	const reco::Candidate *TagQuark2=NULL;

	const reco::Candidate *EWKW = NULL;
	const reco::Candidate *EWKTagQuark1=NULL;
	const reco::Candidate *EWKTagQuark2=NULL;

	// Vector Bosons Info 
	for(size_t i = 0; i < nGen; ++ i) {

		V = &((*genParticles)[i]);

		// The vector boson must have stutus==3  
		if( !(abs(V->status())==3) ) continue;

		if(V->pdgId()==22)Photon_pt_gen = V->pt();

		size_t ndau = 0;

		if(!(V==NULL)) ndau = V->numberOfDaughters();
		// The vector boson must decay to leptons

		if(ndau<1) continue;
		if( (Vtype_=="Z") && !( V->pdgId()==22 || V->pdgId()==23) ) continue;
		if( (Vtype_=="W") && !(abs(V->pdgId())==24) ) continue;

		// Loop over daugthers
		for(size_t j = 0; j < ndau; ++ j) {
			const reco::Candidate *d = V->daughter( j );
			// first look for Z --> l+l-
			if( !(d==NULL) && (V->pdgId()==23 || V->pdgId()==22) ) {
				if (d->pdgId()==-pdgIdDau_)  lepton1  = d;
				if (d->pdgId()==pdgIdDau_) lepton2 = d;
			} // if not, then look for W-->lnu
			else if( !(d==NULL) && abs(V->pdgId())==24) {
				if ( abs(d->pdgId())==pdgIdDau_ )  lepton1  = d;
				if ( abs(d->pdgId())==(pdgIdDau_+1) )  lepton2  = d;
			} 
		} // end ndaughter loop
	}// end genParticles loop



	/*//tth gen information

	for(size_t i = 0; i < nGen; ++ i) {
		aH = &((*genParticles)[i]);
		t = &((*genParticles)[i]);
		tbar = &((*genParticles)[i]);
		// asociated higgs
		if( ((abs(aH->status())==3) && (abs(aH->pdgId())==25) )){
			//        std::cout<<"tth    "<<aH->pdgId()<<std::endl;
			size_t aHndau =0;
			if(!(aH==NULL)) aHndau = aH->numberOfDaughters();
			//        std::cout<<"bbar    "<<aHndau<<std::endl;
			if(aHndau<1) continue;
			for(size_t k =0; k< aHndau-1; k++){
				//loop over higgs daughter
				const reco::Candidate *g = aH->daughter( k );
				if( !(g==NULL) && (g->pdgId()==5) ) {Hb=g;
					//         std::cout<<"h decay b    "<<g->pdgId()<<std::endl;
				}
				else if (!(g==NULL) && (g->pdgId()==-5)) {Hbbar=g;
					//        std::cout<<" h decay bbar    "<<g->pdgId()<<std::endl;
				}
			}
		}

		//associated top 

		if( ((abs(t->status())==3) && (t->pdgId()==6) )){
			std::cout<<" associated top   "<<t->pdgId()<<std::endl;
			size_t atndau =0;
			if(!(t==NULL)) atndau = t->numberOfDaughters();
			//        std::cout<<"bbar    "<<aHndau<<std::endl;
			if(atndau<1) continue;
			for(size_t k =0; k< atndau-1; k++){
				//loop over associated top daughter
				const reco::Candidate *p = t->daughter( k );
				if (!(p==NULL) && (abs(p->pdgId())==5)) tb=p;
				else if( !(p==NULL) && (p->pdgId()==24) ) {
					//        std::cout<<"W    "<<p->pdgId()<<std::endl;
					size_t tWndau =0;
					tWndau = p->numberOfDaughters();
					if(tWndau<1) continue;
					for (size_t l =0; l<tWndau-1; l++){
						//loop over tW daughter
						const reco::Candidate *g = p->daughter( l );
						//std::cout<<" dau    "<<e->numberOfDaughters()<<std::endl;
						//  std::cout<<" dau of tW   "<<g->pdgId()<<std::endl;
						if( !(g==NULL)&& (abs(g->pdgId())<=4)){
							tParton1 =g->daughter(0);
							tParton2 =g->daughter(1);}
						else if (!(g==NULL) && ((abs(g->pdgId())==12) ||(abs(g->pdgId())==14))) tMet =g;
						else if (!(g==NULL) && ((abs(g->pdgId())==11) ||(abs(g->pdgId())==13))) tLepton =g;

					}// loop over dau of tW

				} // W condidate
			} // loop over top daughter
		}// associated top

		// associated anti top

		if( ((abs(tbar->status())==3) && (tbar->pdgId()==-6) )){
			std::cout<<" associated anti top    "<<tbar->pdgId()<<std::endl;
			size_t atbarndau =0;
			if(!(tbar==NULL)) atbarndau = tbar->numberOfDaughters();
			//std::cout<<"bbar    "<<aHndau<<std::endl;
			if(atbarndau<1) continue;
			for(size_t k =0; k< atbarndau-1; k++){
				//loop over associated anti top daughter
				const reco::Candidate *q = tbar->daughter( k );
				if (!(q==NULL) && (q->pdgId()==-5)){ tbbar=q;
					//       std::cout<<" anti top bbar    "<<q->pdgId()<<std::endl;
				}
				else if( !(q==NULL) && (q->pdgId()==-24) ) {
					//        std::cout<<" anti top W    "<<q->pdgId()<<std::endl;
					size_t tbarWndau =0;
					tbarWndau = q->numberOfDaughters();
					if(tbarWndau<1) continue;
					for (size_t l =0; l<tbarWndau-1; l++){
						//loop over tbarW daughter
						const reco::Candidate *h = q->daughter( l );
						//std::cout<<" dau    "<<e->numberOfDaughters()<<std::endl;
						//       std::cout<<" dau of anti top W    "<<h->pdgId()<<std::endl;
						if( !(h==NULL)&& (abs(h->pdgId())<=4)){
							tParton1 =h->daughter(0);
							tParton2 =h->daughter(1);}
						else if (!(h==NULL) && ((abs(h->pdgId())==12) ||(abs(h->pdgId())==14))) tMet =h;
						else if (!(h==NULL) && ((abs(h->pdgId())==11) ||(abs(h->pdgId())==13))) tLepton =h;

					} //loop over dau of tbar W 	
				}	// W conditate loop	

			} // loop over tbar daughter
		}// anti top


	}// gen Particles loop
     
	// tth info ends


	//;;;;;;;;;;;;;;;;;;;;;;;;;;;my stuff;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	for(size_t i = 0; i < nGen; ++ i) {
		H = &((*genParticles)[i]);
		if( !((abs(H->status())==3) && (abs(H->pdgId())==25) ))continue;
		//Higgs must decay to W
		//	std::cout<<"hhhhhhhhhhhh    "<<H->pdgId()<<std::endl;
		size_t Hndau =0;
		size_t Wndau =0;
		if(!(H==NULL)) Hndau = H->numberOfDaughters();
		//        std::cout<<"WWWW    "<<Hndau<<std::endl;
		if(Hndau<1) continue;
		for(size_t k =0; k< Hndau-1; k++){
			//loop over higgs daughter
			const reco::Candidate *e = H->daughter( k );
			if( !(e==NULL) && (abs(e->pdgId())==24) ) {
				//if higgs dau is W
				//      std::cout<<"WWWW    "<<e->pdgId()<<std::endl;
				Wndau = e->numberOfDaughters();
				if(Wndau<1) continue;
				for (size_t l =0; l<Wndau-1; l++){
					//loop over W daughter
					const reco::Candidate *f = e->daughter( l );
					//std::cout<<" dau    "<<e->numberOfDaughters()<<std::endl;
					//    std::cout<<" dau    "<<f->pdgId()<<std::endl;
					if( !(f==NULL)&& (abs(f->pdgId())<=4)){ 
						Parton1 =e->daughter(0);
						Parton2 =e->daughter(1);}
					else if (!(f==NULL) && ((abs(f->pdgId())==12) ||(abs(f->pdgId())==14))) Met =f;
					else if (!(f==NULL) && ((abs(f->pdgId())==11) ||(abs(f->pdgId())==13))) Lepton =f;
				}
			}	
		} 

		// save tag quark information for vbf case in HWW topology 
		int counter = 0;
		if(H!=NULL){
			for(reco::GenParticleCollection::const_iterator p = genParticles -> begin();p != genParticles -> end(); ++p,counter++){

				const reco::Candidate* pCurrent = &(*p);

				if(counter == 7) TagQuark1 = pCurrent;
				if(counter == 8) TagQuark2 = pCurrent;

			}
		}



		////////// Higgs boson quantities //////////////
		if(H==NULL) return;

		H_mass = H->mass();
		H_Eta = H->eta();   
		H_Phi = H->phi();
		H_Vx = H->vx();
		H_Vy = H->vy();
		H_Vz = H->vz();
		H_Y  = H->rapidity();
		H_px = H->px();
		H_py = H->py();
		H_pz = H->pz();
		H_E  = H->energy();
		H_Pt = H->pt();
		H_Et = H->et();
		H_Id = H->pdgId();

	} //nGen loop end

	//#################################EWKW2Jets###################################
	for(size_t i = 0; i < nGen; ++ i) {
		EWKW = &((*genParticles)[i]);
		//Generated W from Matrix Element
		if( !((abs(EWKW->status())==3) && (abs(EWKW->pdgId())==24) ))continue;

		if( EWKW!=NULL && (abs(EWKW->status())==3) && (abs(EWKW->pdgId())==24) ) {
			EWKW_Charge          = EWKW->charge();
			EWKW_Vx              = EWKW->vx();
			EWKW_Vy              = EWKW->vy();
			EWKW_Vz              = EWKW->vz();
			EWKW_Y               = EWKW->rapidity();
			EWKW_Theta           = EWKW->theta();
			EWKW_Eta             = EWKW->eta();
			EWKW_Phi             = EWKW->phi();
			EWKW_E               = EWKW->energy();
			EWKW_px              = EWKW->px();
			EWKW_py              = EWKW->py();
			EWKW_pz              = EWKW->pz();
			EWKW_Pt              = EWKW->pt();
			EWKW_Et              = EWKW->et(); 
			EWKW_Id              = EWKW->pdgId();

			int tmp1 = i + 2;
			int tmp2 = i + 1;
			EWKTagQuark1 = &((*genParticles)[tmp1]);
			EWKTagQuark2 = &((*genParticles)[tmp2]);

		}

	} //nGen loop end
    */

	for(size_t i = 0; i < nGen; ++ i) {
		const reco::Candidate* tmp = &((*genParticles)[i]);
		//Generated W from Matrix Element
		if( !((abs(tmp->status())==3) && (abs(tmp->pdgId())==24) )) continue;

		if((abs(tmp->status())==3) && (abs(tmp->pdgId())==24) )
		{
			for(reco::GenParticleCollection::const_iterator p = genParticles -> begin();p != genParticles -> end(); ++p){

				const reco::Candidate* pCurrent = &(*p);

				if((pCurrent->status())==3 && (abs(pCurrent->pdgId())<=5 || abs(pCurrent->pdgId()) == 21) )
				{
					nParton_Winclusive++;
				}
			}
		}
	} //nGen loop end

	//std::cout << "nPartons: " << nParton_Winclusive << std::endl;

	////////// lepton #1 quantities //////////////
	if( !(lepton1 == NULL) ) {
		l1Charge           = lepton1-> charge();
		l1Vx               = lepton1->vx();
		l1Vy               = lepton1->vy();
		l1Vz               = lepton1->vz();
		l1Y                = lepton1->rapidity();
		l1Theta            = lepton1->theta();
		l1Eta              = lepton1->eta();    
		l1Phi              = lepton1->phi();
		l1E                = lepton1->energy();
		l1px               = lepton1->px();
		l1py               = lepton1->py();
		l1pz               = lepton1->pz();
		l1Pt               = lepton1->pt();
		l1Et               = lepton1->et();	  
	}

	////////// lepton #2 quantities: in case of Z ///////
	if( !(lepton2 == NULL) ) {
		l2Charge          = lepton2->charge();
		l2Vx              = lepton2->vx();
		l2Vy              = lepton2->vy();
		l2Vz              = lepton2->vz();
		l2Y               = lepton2->rapidity();
		l2Theta           = lepton2->theta();
		l2Eta             = lepton2->eta();    
		l2Phi             = lepton2->phi();
		l2E               = lepton2->energy();
		l2px              = lepton2->px();
		l2py              = lepton2->py();
		l2pz              = lepton2->pz();
		l2Pt              = lepton2->pt();
		l2Et              = lepton2->et();	 
	} 

	/*//Parton filling
	if( !(Parton1 == NULL) ) {

		Parton_Charge[0]          = Parton1->charge();
		Parton_Vx[0]              = Parton1->vx();
		Parton_Vy[0]              = Parton1->vy();
		Parton_Vz[0]              = Parton1->vz();
		Parton_Y[0]               = Parton1->rapidity();
		Parton_Theta[0]           = Parton1->theta();
		Parton_Eta[0]             = Parton1->eta();
		Parton_Phi[0]             = Parton1->phi();
		Parton_E[0]               = Parton1->energy();
		Parton_px[0]              = Parton1->px();
		Parton_py[0]              = Parton1->py();
		Parton_pz[0]              = Parton1->pz();
		Parton_Pt[0]              = Parton1->pt();
		Parton_Et[0]              = Parton1->et(); 
		Parton_Id[0]              = Parton1->pdgId();
	}

	//Parton filling
	if( !(Parton2 == NULL) ) {

		Parton_Charge[1]          = Parton2->charge();
		Parton_Vx[1]              = Parton2->vx();
		Parton_Vy[1]              = Parton2->vy();
		Parton_Vz[1]              = Parton2->vz();
		Parton_Y[1]               = Parton2->rapidity();
		Parton_Theta[1]           = Parton2->theta();
		Parton_Eta[1]             = Parton2->eta();
		Parton_Phi[1]             = Parton2->phi();
		Parton_E[1]               = Parton2->energy();
		Parton_px[1]              = Parton2->px();
		Parton_py[1]              = Parton2->py();
		Parton_pz[1]              = Parton2->pz();
		Parton_Pt[1]              = Parton2->pt();
		Parton_Et[1]              = Parton2->et();
		Parton_Id[1]              = Parton2->pdgId();
	}*/

	//Lepton filling

	if( !(Lepton == NULL) ) {
		Lepton_Charge          = Lepton->charge();
		Lepton_Vx              = Lepton->vx();
		Lepton_Vy              = Lepton->vy();
		Lepton_Vz              = Lepton->vz();
		Lepton_Y               = Lepton->rapidity();
		Lepton_Theta           = Lepton->theta();
		Lepton_Eta             = Lepton->eta();
		Lepton_Phi             = Lepton->phi();
		Lepton_E               = Lepton->energy();
		Lepton_px              = Lepton->px();
		Lepton_py              = Lepton->py();
		Lepton_pz              = Lepton->pz();
		Lepton_Pt              = Lepton->pt();
		Lepton_Et              = Lepton->et();
		Lepton_Id              = Lepton->pdgId();
	}


	//Met filling

	if( !(Met == NULL) ) {
		Met_Charge          = Met->charge();
		Met_Vx              = Met->vx();
		Met_Vy              = Met->vy();
		Met_Vz              = Met->vz();
		Met_Y               = Met->rapidity();
		Met_Theta           = Met->theta();
		Met_Eta             = Met->eta();
		Met_Phi             = Met->phi();
		Met_E               = Met->energy();
		Met_px              = Met->px();
		Met_py              = Met->py();
		Met_pz              = Met->pz();
		Met_Pt              = Met->pt();
		Met_Et              = Met->et();
		Met_Id              = Met->pdgId();
	}


/*	//tParton filling
	if( !(tParton1 == NULL) ) {

		tParton_Charge[0]          = tParton1->charge();
		tParton_Vx[0]              = tParton1->vx();
		tParton_Vy[0]              = tParton1->vy();
		tParton_Vz[0]              = tParton1->vz();
		tParton_Y[0]               = tParton1->rapidity();
		tParton_Theta[0]           = tParton1->theta();
		tParton_Eta[0]             = tParton1->eta();
		tParton_Phi[0]             = tParton1->phi();
		tParton_E[0]               = tParton1->energy();
		tParton_px[0]              = tParton1->px();
		tParton_py[0]              = tParton1->py();
		tParton_pz[0]              = tParton1->pz();
		tParton_Pt[0]              = tParton1->pt();
		tParton_Et[0]              = tParton1->et(); 
		tParton_Id[0]              = tParton1->pdgId();
	}

	//tParton filling
	if( !(tParton2 == NULL) ) {

		tParton_Charge[1]          = tParton2->charge();
		tParton_Vx[1]              = tParton2->vx();
		tParton_Vy[1]              = tParton2->vy();
		tParton_Vz[1]              = tParton2->vz();
		tParton_Y[1]               = tParton2->rapidity();
		tParton_Theta[1]           = tParton2->theta();
		tParton_Eta[1]             = tParton2->eta();
		tParton_Phi[1]             = tParton2->phi();
		tParton_E[1]               = tParton2->energy();
		tParton_px[1]              = tParton2->px();
		tParton_py[1]              = tParton2->py();
		tParton_pz[1]              = tParton2->pz();
		tParton_Pt[1]              = tParton2->pt();
		tParton_Et[1]              = tParton2->et();
		tParton_Id[1]              = tParton2->pdgId();
	}

	//Lepton filling

	if( !(tLepton == NULL) ) {
		tLepton_Charge          = tLepton->charge();
		tLepton_Vx              = tLepton->vx();
		tLepton_Vy              = tLepton->vy();
		tLepton_Vz              = tLepton->vz();
		tLepton_Y               = tLepton->rapidity();
		tLepton_Theta           = tLepton->theta();
		tLepton_Eta             = tLepton->eta();
		tLepton_Phi             = tLepton->phi();
		tLepton_E               = tLepton->energy();
		tLepton_px              = tLepton->px();
		tLepton_py              = tLepton->py();
		tLepton_pz              = tLepton->pz();
		tLepton_Pt              = tLepton->pt();
		tLepton_Et              = tLepton->et();
		tLepton_Id              = tLepton->pdgId();
	}


	//Met filling

	if( !(tMet == NULL) ) {
		tMet_Charge          = tMet->charge();
		tMet_Vx              = tMet->vx();
		tMet_Vy              = tMet->vy();
		tMet_Vz              = tMet->vz();
		tMet_Y               = tMet->rapidity();
		tMet_Theta           = tMet->theta();
		tMet_Eta             = tMet->eta();
		tMet_Phi             = tMet->phi();
		tMet_E               = tMet->energy();
		tMet_px              = tMet->px();
		tMet_py              = tMet->py();
		tMet_pz              = tMet->pz();
		tMet_Pt              = tMet->pt();
		tMet_Et              = tMet->et();
		tMet_Id              = tMet->pdgId();
	}



	if( !(tb == NULL) ) {
		tb_Charge          = tb->charge();
		tb_Vx              = tb->vx();
		tb_Vy              = tb->vy();
		tb_Vz              = tb->vz();
		tb_Y               = tb->rapidity();
		tb_Theta           = tb->theta();
		tb_Eta             = tb->eta();
		tb_Phi             = tb->phi();
		tb_E               = tb->energy();
		tb_px              = tb->px();
		tb_py              = tb->py();
		tb_pz              = tb->pz();
		tb_Pt              = tb->pt();
		tb_Et              = tb->et();
		tb_Id              = tb->pdgId();
	}


	if( !(tbbar == NULL) ) {
		tbbar_Charge          = tbbar->charge();
		tbbar_Vx              = tbbar->vx();
		tbbar_Vy              = tbbar->vy();
		tbbar_Vz              = tbbar->vz();
		tbbar_Y               = tbbar->rapidity();
		tbbar_Theta           = tbbar->theta();
		tbbar_Eta             = tbbar->eta();
		tbbar_Phi             = tbbar->phi();
		tbbar_E               = tbbar->energy();
		tbbar_px              = tbbar->px();
		tbbar_py              = tbbar->py();
		tbbar_pz              = tbbar->pz();
		tbbar_Pt              = tbbar->pt();
		tbbar_Et              = tbbar->et();
		tbbar_Id              = tbbar->pdgId();
	}


	if( !(Hb == NULL) ) {
		Hb_Charge          = Hb->charge();
		Hb_Vx              = Hb->vx();
		Hb_Vy              = Hb->vy();
		Hb_Vz              = Hb->vz();
		Hb_Y               = Hb->rapidity();
		Hb_Theta           = Hb->theta();
		Hb_Eta             = Hb->eta();
		Hb_Phi             = Hb->phi();
		Hb_E               = Hb->energy();
		Hb_px              = Hb->px();
		Hb_py              = Hb->py();
		Hb_pz              = Hb->pz();
		Hb_Pt              = Hb->pt();
		Hb_Et              = Hb->et();
		Hb_Id              = Hb->pdgId();
	}


	if( !(Hbbar == NULL) ) {
		Hbbar_Charge          = Hbbar->charge();
		Hbbar_Vx              = Hbbar->vx();
		Hbbar_Vy              = Hbbar->vy();
		Hbbar_Vz              = Hbbar->vz();
		Hbbar_Y               = Hbbar->rapidity();
		Hbbar_Theta           = Hbbar->theta();
		Hbbar_Eta             = Hbbar->eta();
		Hbbar_Phi             = Hbbar->phi();
		Hbbar_E               = Hbbar->energy();
		Hbbar_px              = Hbbar->px();
		Hbbar_py              = Hbbar->py();
		Hbbar_pz              = Hbbar->pz();
		Hbbar_Pt              = Hbbar->pt();
		Hbbar_Et              = Hbbar->et();
		Hbbar_Id              = Hbbar->pdgId();
	}

	///////////////////////////////// Fill Tag Quark info for VBF case

	//Parton filling
	if( H!=NULL && TagQuark1 != NULL && TagQuark2!=NULL) {

		TagQuark_Charge[0]          = TagQuark1->charge();
		TagQuark_Vx[0]              = TagQuark1->vx();
		TagQuark_Vy[0]              = TagQuark1->vy();
		TagQuark_Vz[0]              = TagQuark1->vz();
		TagQuark_Y[0]               = TagQuark1->rapidity();
		TagQuark_Theta[0]           = TagQuark1->theta();
		TagQuark_Eta[0]             = TagQuark1->eta();
		TagQuark_Phi[0]             = TagQuark1->phi();
		TagQuark_E[0]               = TagQuark1->energy();
		TagQuark_px[0]              = TagQuark1->px();
		TagQuark_py[0]              = TagQuark1->py();
		TagQuark_pz[0]              = TagQuark1->pz();
		TagQuark_Pt[0]              = TagQuark1->pt();
		TagQuark_Et[0]              = TagQuark1->et(); 
		TagQuark_Id[0]              = TagQuark1->pdgId();

		TagQuark_Charge[1]          = TagQuark2->charge();
		TagQuark_Vx[1]              = TagQuark2->vx();
		TagQuark_Vy[1]              = TagQuark2->vy();
		TagQuark_Vz[1]              = TagQuark2->vz();
		TagQuark_Y[1]               = TagQuark2->rapidity();
		TagQuark_Theta[1]           = TagQuark2->theta();
		TagQuark_Eta[1]             = TagQuark2->eta();
		TagQuark_Phi[1]             = TagQuark2->phi();
		TagQuark_E[1]               = TagQuark2->energy();
		TagQuark_px[1]              = TagQuark2->px();
		TagQuark_py[1]              = TagQuark2->py();
		TagQuark_pz[1]              = TagQuark2->pz();
		TagQuark_Pt[1]              = TagQuark2->pt();
		TagQuark_Et[1]              = TagQuark2->et();
		TagQuark_Id[1]              = TagQuark2->pdgId();

	}


	if( EWKTagQuark1 != NULL && EWKTagQuark2!=NULL) {

		EWKTagQuark_Charge[0]          = EWKTagQuark1->charge();
		EWKTagQuark_Vx[0]              = EWKTagQuark1->vx();
		EWKTagQuark_Vy[0]              = EWKTagQuark1->vy();
		EWKTagQuark_Vz[0]              = EWKTagQuark1->vz();
		EWKTagQuark_Y[0]               = EWKTagQuark1->rapidity();
		EWKTagQuark_Theta[0]           = EWKTagQuark1->theta();
		EWKTagQuark_Eta[0]             = EWKTagQuark1->eta();
		EWKTagQuark_Phi[0]             = EWKTagQuark1->phi();
		EWKTagQuark_E[0]               = EWKTagQuark1->energy();
		EWKTagQuark_px[0]              = EWKTagQuark1->px();
		EWKTagQuark_py[0]              = EWKTagQuark1->py();
		EWKTagQuark_pz[0]              = EWKTagQuark1->pz();
		EWKTagQuark_Pt[0]              = EWKTagQuark1->pt();
		EWKTagQuark_Et[0]              = EWKTagQuark1->et(); 
		EWKTagQuark_Id[0]              = EWKTagQuark1->pdgId();

		EWKTagQuark_Charge[1]          = EWKTagQuark2->charge();
		EWKTagQuark_Vx[1]              = EWKTagQuark2->vx();
		EWKTagQuark_Vy[1]              = EWKTagQuark2->vy();
		EWKTagQuark_Vz[1]              = EWKTagQuark2->vz();
		EWKTagQuark_Y[1]               = EWKTagQuark2->rapidity();
		EWKTagQuark_Theta[1]           = EWKTagQuark2->theta();
		EWKTagQuark_Eta[1]             = EWKTagQuark2->eta();
		EWKTagQuark_Phi[1]             = EWKTagQuark2->phi();
		EWKTagQuark_E[1]               = EWKTagQuark2->energy();
		EWKTagQuark_px[1]              = EWKTagQuark2->px();
		EWKTagQuark_py[1]              = EWKTagQuark2->py();
		EWKTagQuark_pz[1]              = EWKTagQuark2->pz();
		EWKTagQuark_Pt[1]              = EWKTagQuark2->pt();
		EWKTagQuark_Et[1]              = EWKTagQuark2->et();
		EWKTagQuark_Id[1]              = EWKTagQuark2->pdgId();

	}
*/
}



////////////////// utilities, helpers ///////////////////

void ewk::MCTreeFiller::SetBranch( float* x, std::string name)
{
	std::string brName = std::string(name_) + "_" + name;
	tree_->Branch( brName.c_str(), x, ( brName+"/F").c_str() );
}


void ewk::MCTreeFiller::SetBranch( int* x, std::string name)
{
	std::string brName = std::string(name_) + "_" + name;
	tree_->Branch( brName.c_str(), x, ( brName+"/I").c_str() );
}


void ewk::MCTreeFiller::SetBranch( bool* x, std::string name)
{
	std::string brName = std::string(name_) + "_" + name;
	tree_->Branch( brName.c_str(), x, ( brName+"/O").c_str() );
}

