
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
 *   Can work with GenJet, PF jet.
 *   Can work with jets in RECO/AOD/PAT data formats.
 * History:
 *   
 *
 * Copyright (C) 2010 FNAL 
 *****************************************************************************/

// CMS includes
#include "TMath.h" 
#include "JetMETCorrections/MCJet/plugins/JetUtilMC.h" // needed for dPhi,dR


// Header file
#include "ElectroWeakAnalysis/VPlusJets/interface/JetTreeFiller.h"
#include "ElectroWeakAnalysis/VPlusJets/interface/METzCalculator.h"
#include "ElectroWeakAnalysis/VPlusJets/interface/AngularVars.h"
#include "ElectroWeakAnalysis/VPlusJets/interface/ColorCorrel.h"

const float BTAG_DISCRIM_DEFAULT=-999.;

ewk::JetTreeFiller::JetTreeFiller(const char *name, TTree* tree, const std::string jetType, const edm::ParameterSet iConfig) {
	// ********** GenJets ********** //
	if( jetType=="Gen" && iConfig.existsAs<edm::InputTag>("srcGen") )
	  mInputJets = iConfig.getParameter<edm::InputTag>("srcGen"); 
	// ********** PFJets ********** //
	if( jetType=="PF" && iConfig.existsAs<edm::InputTag>("srcPFJets") )
	  mInputJets = iConfig.getParameter<edm::InputTag>("srcPFJets"); 
	// ********** Corrected PFJets ********** //
	if( jetType=="PFCor" && iConfig.existsAs<edm::InputTag>("srcPFCor") )
	  mInputJets = iConfig.getParameter<edm::InputTag>("srcPFCor"); 
	// ********** Corrected PFJets for VBF Tag ********** //
	if( jetType=="PFCorVBFTag" && iConfig.existsAs<edm::InputTag>("srcPFCorVBFTag") )
	  mInputJets = iConfig.getParameter<edm::InputTag>("srcPFCorVBFTag"); 


	if(  iConfig.existsAs<edm::InputTag>("srcMet") )
	  mInputMet = iConfig.getParameter<edm::InputTag>("srcMet");
	if(  iConfig.existsAs<edm::InputTag>("srcMetMVA") )
	  mInputMetMVA = iConfig.getParameter<edm::InputTag>("srcMetMVA");

	// ********** Vector boson ********** //
	if(  iConfig.existsAs<edm::InputTag>("srcVectorBoson") )
	  mInputBoson = iConfig.getParameter<edm::InputTag>("srcVectorBoson"); 

	//*********************  Run Over AOD or PAT  ***********//
	if( iConfig.existsAs<bool>("runningOverAOD"))
	  runoverAOD = iConfig.getParameter<bool>("runningOverAOD");


	tree_     = tree;
	jetType_ = jetType;

	Vtype_    = iConfig.getParameter<std::string>("VBosonType"); 
	LeptonType_ = iConfig.getParameter<std::string>("LeptonType");

	// ---- Quark Gluon Likelihood
	qglikeli = new QGLikelihoodCalculator();  

	if( !(tree==0) ) SetBranches();
}


void ewk::JetTreeFiller::SetBranches()
{
	// Declare jet branches
	SetBranchSingle( &NumJets, "num" + jetType_ + "Jets");
	SetBranchSingle( &numBTags, "num" + jetType_ + "JetBTags");
	SetBranch( Et, "Jet" + jetType_ + "_Et");
	SetBranch( Pt, "Jet" + jetType_ + "_Pt");
	SetBranch( Eta, "Jet" + jetType_ + "_Eta");
	SetBranch( Phi, "Jet" + jetType_ + "_Phi");
	SetBranch( Theta, "Jet" + jetType_ + "_Theta");
	SetBranch( Px, "Jet" + jetType_ + "_Px");
	SetBranch( Py, "Jet" + jetType_ + "_Py");
	SetBranch( Pz, "Jet" + jetType_ + "_Pz");
	SetBranch( E, "Jet" + jetType_ + "_E");
	SetBranch( Y, "Jet" + jetType_ + "_Y");
	SetBranch( Mass, "Jet" + jetType_ + "_Mass");

	/// eta-eta second moment, ET weighted
	SetBranch( etaetaMoment, "Jet" + jetType_ + "_etaetaMoment");
	/// phi-phi second moment, ET weighted
	SetBranch( phiphiMoment, "Jet" + jetType_ + "_phiphiMoment");
	/// eta-phi second moment, ET weighted
	SetBranch( etaphiMoment, "Jet" + jetType_ + "_etaphiMoment");
	/// maximum distance from jet to constituent
	SetBranch( maxDistance, "Jet" + jetType_ + "_maxDistance");
	/// # of constituents
	SetBranch( nConstituents, "Jet" + jetType_ + "_nConstituents");
	/// Area of the jet
	SetBranch( Area, "Jet" + jetType_ + "_Area");

	SetBranch( VjetMass, "Vplus" + jetType_ + "Jet_Mass");
	SetBranch( Dphi, "Jet" + jetType_ + "_dphiBoson");
	SetBranch( Deta, "Jet" + jetType_ + "_detaBoson");
	SetBranch( DR, "Jet" + jetType_ + "_dRBoson");
	SetBranch( DphiMET, "Jet" + jetType_ + "_dphiMET");
	// SetBranch( Response, "Jet" + jetType_ + "_Response");
	SetBranch( bDiscriminator, "Jet" + jetType_ + "_bDiscriminator");
	SetBranch( bDiscriminatorSSVHE, "Jet" + jetType_ + "_bDiscriminatorSSVHE");
	SetBranch( bDiscriminatorTCHE, "Jet" + jetType_ + "_bDiscriminatorTCHE");
	SetBranch( bDiscriminatorCSV, "Jet" + jetType_ + "_bDiscriminatorCSV");
	SetBranch( bDiscriminatorJP, "Jet" + jetType_ + "_bDiscriminatorJP");
	SetBranch( bDiscriminatorSSVHP, "Jet" + jetType_ + "_bDiscriminatorSSVHP");
	SetBranch( bDiscriminatorTCHP, "Jet" + jetType_ + "_bDiscriminatorTCHP");
	SetBranch( secVertexMass, "Jet" + jetType_ + "_secVertexMass");
	//   SetBranch( Dphi2, "Jet" + jetType_ + "_dphiBoson2");
	//   SetBranch( Deta2, "Jet" + jetType_ + "_detaBoson2");
	//   SetBranch( DR2, "Jet" + jetType_ + "_dRBoson2");
	//   SetBranch( VjetMass2, "Vplus" + jetType_ + "Jet_Mass2");
	//   SetBranch( Response2, "Jet" + jetType_ + "_Response2");

	/////////////////////////////////////////////////////////////////////////

	if( jetType_ == "PF" || jetType_ == "PFCor" || jetType_ == "PFCorVBFTag") {
		/// chargedHadronEnergy 
		SetBranch( PFChargedHadronEnergy, "Jet" + jetType_ + "_ChargedHadronEnergy");
		///  chargedHadronEnergyFraction
		SetBranch( PFChargedHadronEnergyFraction, "Jet" + jetType_ + "_ChargedHadronEnergyFrac");
		/// neutralHadronEnergy
		SetBranch( PFNeutralHadronEnergy, "Jet" + jetType_ + "_NeutralHadronEnergy");
		/// neutralHadronEnergyFraction
		SetBranch( PFNeutralHadronEnergyFraction, "Jet" + jetType_ + "_NeutralHadronEnergyFrac");
		/// chargedEmEnergy
		SetBranch( PFChargedEmEnergy, "Jet" + jetType_ + "_ChargedEmEnergy");
		/// chargedEmEnergyFraction
		SetBranch( PFChargedEmEnergyFraction, "Jet" + jetType_ + "_ChargedEmEnergyFrac");
		/// chargedMuEnergy
		SetBranch( PFChargedMuEnergy, "Jet" + jetType_ + "_ChargedMuEnergy");
		/// chargedMuEnergyFraction
		SetBranch( PFChargedMuEnergyFraction, "Jet" + jetType_ + "_ChargedMuEnergyFrac");
		/// neutralEmEnergy
		SetBranch( PFNeutralEmEnergy, "Jet" + jetType_ + "_NeutralEmEnergy");
		/// neutralEmEnergyFraction
		SetBranch( PFNeutralEmEnergyFraction, "Jet" + jetType_ + "_NeutralEmEnergyFrac");
		/// chargedMultiplicity
		SetBranch( PFChargedMultiplicity, "Jet" + jetType_ + "_ChargedMultiplicity");
		/// neutralMultiplicity
		SetBranch( PFNeutralMultiplicity, "Jet" + jetType_ + "_NeutralMultiplicity");
		/// muonMultiplicity
		SetBranch( PFMuonMultiplicity, "Jet" + jetType_ + "_MuonMultiplicity");
		/// photonEnergy 
		SetBranch( PFPhotonEnergy, "Jet" + jetType_ + "_PhotonEnergy");
		/// photonEnergyFraction 
		SetBranch( PFPhotonEnergyFraction, "Jet" + jetType_ + "_PhotonEnergyFraction");
		/// electronEnergy 
		SetBranch( PFElectronEnergy, "Jet" + jetType_ + "_ElectronEnergy"); 
		/// electronEnergyFraction 
		SetBranch( PFElectronEnergyFraction, "Jet" + jetType_ + "_ElectronEnergyFraction");
		/// muonEnergy 
		SetBranch( PFMuonEnergy, "Jet" + jetType_ + "_MuonEnergy"); 
		/// muonEnergyFraction 
		SetBranch( PFMuonEnergyFraction, "Jet" + jetType_ + "_MuonEnergyFraction");
		/// HFHadronEnergy  
		SetBranch( PFHFHadronEnergy, "Jet" + jetType_ + "_HFHadronEnergy");
		/// HFHadronEnergyFraction 
		SetBranch( PFHFHadronEnergyFraction, "Jet" + jetType_ + "_HFHadronEnergyFraction");
		/// HFEMEnergy  
		SetBranch( PFHFEMEnergy, "Jet" + jetType_ + "_HFEMEnergy");
		/// HFEMEnergyFraction 
		SetBranch( PFHFEMEnergyFraction, "Jet" + jetType_ + "_HFEMEnergyFraction");
		/// chargedHadronMultiplicity 
		SetBranch( PFChargedHadronMultiplicity, "Jet" + jetType_ + "_ChargedHadronMultiplicity");
		/// neutralHadronMultiplicity 
		SetBranch( PFNeutralHadronMultiplicity, "Jet" + jetType_ + "_NeutralHadronMultiplicity");
		/// photonMultiplicity 
		SetBranch( PFPhotonMultiplicity, "Jet" + jetType_ + "_PhotonMultiplicity");
		/// electronMultiplicity 
		SetBranch( PFElectronMultiplicity, "Jet" + jetType_ + "_ElectronMultiplicity");
		/// HFHadronMultiplicity 
		SetBranch( PFHFHadronMultiplicity, "Jet" + jetType_ + "_HFHadronMultiplicity");
		/// HFEMMultiplicity 
		SetBranch( PFHFEMMultiplicity, "Jet" + jetType_ + "_HFEMMultiplicity");
		/// Sum pt of all pf constituents
		SetBranch( PFsumPtCands, "Jet" + jetType_ + "_SumPtCands");
		/// Sum pt^2 of all pf constituents
		SetBranch( PFsumPt2Cands, "Jet" + jetType_ + "_SumPt2Cands");
		/// [Sum pt^2*deltaR(cand, jet)^2] / Sum pt^2 of all pf constituents
		SetBranch( PFrmsCands, "Jet" + jetType_ + "_rmsCands");
		/// pt_D variable for QG likelihood:  pt_D = sqrt(Sum_i{pt^2})/ Sum_i(pt), over all PF cands
		// ------- See for details: CMS AN-2011/215
		// -- https://indico.cern.ch/getFile.py/access?contribId=3&resId=0&materialId=slides&confId=129897
		// -- https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=135378
		// -- https://indico.cern.ch/getFile.py/access?contribId=0&resId=0&materialId=slides&confId=144396
		// ------- JetMETCorrections/GammaJet/src/GammaJetAnalyzer.cc
		// 
		SetBranch( PFptD, "Jet" + jetType_ + "_PtD");
		// ---- Quark Gluon Likelihood
		// ---- first check out: cvs co   UserCode/pandolf/QGLikelihood
		// ---- then instantiate:  QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator();
		// ---- finally, compute:   float qgLH = qglikeli->computeQGLikelihoodPU(jet_pt, 
		// -----                         fastjet_rho, chargedMultiplicity, neutralMultiplicity, ptD);
		SetBranch( PFqgLikelihood, "Jet" + jetType_ + "_QGLikelihood");
	}

	SetBranchSingle( &V2jMassMVAMET,  "MassV2j_" + jetType_ + "_MVAMET");
	SetBranchSingle( &V2jMass,  "MassV2j_" + jetType_);
	SetBranchSingle( &V3jMass, "MassV3j_" + jetType_);
	SetBranchSingle( &V4jMass, "MassV4j_" + jetType_);
	SetBranchSingle( &V5jMass, "MassV5j_" + jetType_);
	SetBranchSingle( &V6jMass, "MassV6j_" + jetType_);
	SetBranchSingle( &c2jMass, "Mass2j_" + jetType_);
	SetBranchSingle( &c3jMass, "Mass3j_" + jetType_);
	SetBranchSingle( &c4jMass, "Mass4j_" + jetType_);
	SetBranchSingle( &c5jMass, "Mass5j_" + jetType_);
	SetBranchSingle( &c6jMass, "Mass6j_" + jetType_);


	SetBranchSingle( &V2jCosJacksonAngle, "cosJacksonAngleV2j_" + jetType_);
	SetBranchSingle( &c2jCosJacksonAngle, "cosJacksonAngle2j_" + jetType_);
	SetBranchSingle( &V3jCosJacksonAngle, "cosJacksonAngleV3j_" + jetType_);
	SetBranchSingle( &c3jCosJacksonAngle12, "cosJacksonAngle3j12_" + jetType_);
	SetBranchSingle( &c3jCosJacksonAngle23, "cosJacksonAngle3j23_" + jetType_);
	SetBranchSingle( &c3jCosJacksonAngle31, "cosJacksonAngle3j31_" + jetType_);
	SetBranchSingle( &cosphiDecayPlane, "cosphiDecayPlane_" + jetType_); 
	SetBranchSingle( &cosThetaLnu, "cosThetaLnu_" + jetType_); 
	SetBranchSingle( &cosThetaJJ, "cosThetaJJ_" + jetType_);


	if( jetType_ == "PF" || jetType_ == "PFCor" || jetType_ == "PFCorVBFTag") {
		/// Color Correlation between W jets ( jets pull )
		SetBranchSingle( &colorCorr01, "colorCorrPull01" + jetType_);
		SetBranchSingle( &colorCorr02, "colorCorrPull02" + jetType_);
		SetBranchSingle( &colorCorr12, "colorCorrPull12" + jetType_);
		SetBranchSingle( &colorCorr03, "colorCorrPull03" + jetType_);
		SetBranchSingle( &colorCorr13, "colorCorrPull13" + jetType_);
		SetBranchSingle( &colorCorr23, "colorCorrPull23" + jetType_);
		SetBranchSingle( &colorCorr04, "colorCorrPull04" + jetType_);
		SetBranchSingle( &colorCorr14, "colorCorrPull14" + jetType_);
		SetBranchSingle( &colorCorr24, "colorCorrPull24" + jetType_);
		SetBranchSingle( &colorCorr34, "colorCorrPull34" + jetType_);
		SetBranchSingle( &colorCorr05, "colorCorrPull05" + jetType_);
		SetBranchSingle( &colorCorr15, "colorCorrPull15" + jetType_);
		SetBranchSingle( &colorCorr25, "colorCorrPull25" + jetType_);
		SetBranchSingle( &colorCorr35, "colorCorrPull35" + jetType_);
		SetBranchSingle( &colorCorr45, "colorCorrPull45" + jetType_);

		/// Helicity angles in the Higgs rest frame
		SetBranchSingle( &j1Hel_HiggsCM, "cosThetaJ1HiggsCM_" + jetType_);
		SetBranchSingle( &j2Hel_HiggsCM, "cosThetaJ2HiggsCM_" + jetType_);
		SetBranchSingle( &l1Hel_HiggsCM, "cosThetaL1HiggsCM_" + jetType_);
		SetBranchSingle( &l2Hel_HiggsCM, "cosThetaL2HiggsCM_" + jetType_);
		SetBranchSingle( &b1Hel_HiggsCM, "cosThetaV1HiggsCM_" + jetType_);
		SetBranchSingle( &b2Hel_HiggsCM, "cosThetaV2HiggsCM_" + jetType_);

	}

	if( jetType_ == "PF" || jetType_ == "PFCor" || jetType_ == "PFCorVBFTag") {
		SetBranch( isPileUpJetLoose, "Jet" + jetType_ + "_isPileUpJetLoose");
		SetBranch( isPileUpJetMedium, "Jet" + jetType_ + "_isPileUpJetMedium");
		SetBranch( isPileUpJetTight, "Jet" + jetType_ + "_isPileUpJetTight");
	}
}


//////////////////////////////////////////////////////////////////
/////// Helper for above function ////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

void ewk::JetTreeFiller::SetBranchSingle( float* x, std::string name)
{
	tree_->Branch( name.c_str(), x, ( name+"/F").c_str() );
	bnames.push_back( name );

}

void ewk::JetTreeFiller::SetBranchSingle( int* x, std::string name)
{
	tree_->Branch( name.c_str(), x, ( name+"/I").c_str() );
	bnames.push_back( name );
}

void ewk::JetTreeFiller::SetBranchSingle( bool* x, std::string name){
	tree_->Branch( name.c_str(), x, ( name+"/O").c_str() );
	bnames.push_back( name );
}


void ewk::JetTreeFiller::SetBranch( float* x, std::string name)
{
	tree_->Branch( name.c_str(), x, ( name+"[8]/F").c_str() );
	bnames.push_back( name );
}


void ewk::JetTreeFiller::SetBranch( int* x, std::string name)
{
	tree_->Branch( name.c_str(), x, ( name+"[8]/I").c_str() );
	bnames.push_back( name );
}


void ewk::JetTreeFiller::SetBranch( bool* x, std::string name){
	tree_->Branch( name.c_str(), x, ( name+"[8]/O").c_str() );
	bnames.push_back( name );
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////




void ewk::JetTreeFiller::FillBranches() const 
{
	for(std::vector<std::string>::iterator it = bnames.begin(); it != bnames.end(); ++it) {
		if(TBranch *br = tree_->GetBranch( (*it).c_str() ) ) br->Fill();
	}
}



void ewk::JetTreeFiller::init()   {
	// initialize private data members
	NumJets = 0; 
	numBTags = 0;

	for (int j =0; j< NUM_JET_MAX; ++j) {
		Et[j] = -1.0;
		Pt[j] = -1.0;
		Eta[j] = -10.0;
		Phi[j] = -10.0;
		Theta[j] = -10.0;
		E[j] = -1.0;
		Y[j] = -10.0;
		etaetaMoment[j]  = -10.0;  
		phiphiMoment[j]  = -10.0;      
		etaphiMoment[j]  = -10.0;      
		maxDistance[j]  = -10.0;
		nConstituents[j]  = -1;
		Px[j] = -999999.9;
		Py[j] = -999999.9;
		Pz[j] = -999999.9;
		Mass[j] = -1.0;
		Area[j] = -10.;

		VjetMass[j] = -1.0;
		Dphi[j] = -10.0;
		Deta[j] = -10.0;
		DR[j] = -10.0;
		DphiMET[j] = -10.0;
		Response[j] = -1.0;
		bDiscriminator[j] = BTAG_DISCRIM_DEFAULT;
		bDiscriminatorSSVHE[j] = BTAG_DISCRIM_DEFAULT;
		bDiscriminatorTCHE[j] = BTAG_DISCRIM_DEFAULT;
		bDiscriminatorCSV[j] = BTAG_DISCRIM_DEFAULT;
		bDiscriminatorJP[j] = BTAG_DISCRIM_DEFAULT;
		bDiscriminatorSSVHP[j] = BTAG_DISCRIM_DEFAULT;
		bDiscriminatorTCHP[j] = BTAG_DISCRIM_DEFAULT;
		secVertexMass[j] = -1.0;

		VjetMass2[j] = -1.0;
		DR2[j] = -10.0;
		Dphi2[j] = -10.0;
		Deta2[j] = -10.0;
		Response2[j] = -1.0;

		PFChargedHadronEnergy[j] = -1.0;
		PFChargedHadronEnergyFraction[j] = -1.0;
		PFNeutralHadronEnergy[j] = -1.0;
		PFNeutralHadronEnergyFraction[j] = -1.0;
		PFChargedEmEnergy[j] = -1.0;
		PFChargedEmEnergyFraction[j] = -1.0;
		PFChargedMuEnergy[j] = -1.0;
		PFChargedMuEnergyFraction[j] = -1.0;
		PFNeutralEmEnergy[j] = -1.0;
		PFNeutralEmEnergyFraction[j] = -1.0;
		PFChargedMultiplicity[j] = -1;
		PFNeutralMultiplicity[j] = -1;
		PFMuonMultiplicity[j] = -1;
		PFPhotonEnergy[j] = -1.0;
		PFPhotonEnergyFraction[j] = -1.0;
		PFElectronEnergy[j] = -1.0;
		PFElectronEnergyFraction[j] = -1.0;
		PFMuonEnergy[j] = -1.0;
		PFMuonEnergyFraction[j] = -1.0;
		PFHFHadronEnergy[j] = -1.0;
		PFHFHadronEnergyFraction[j] = -1.0;
		PFHFEMEnergy[j] = -1.0;
		PFHFEMEnergyFraction[j] = -1.0;	 
		PFChargedHadronMultiplicity[j] = -1;
		PFNeutralHadronMultiplicity[j] = -1;
		PFPhotonMultiplicity[j] = -1;
		PFElectronMultiplicity[j] = -1;
		PFHFHadronMultiplicity[j] = -1;
		PFHFEMMultiplicity[j] = -1;
		PFsumPtCands[j]=0.;
		PFsumPt2Cands[j]=0.;
		PFrmsCands[j]=0.;
		PFptD[j] = -1.0;
		PFqgLikelihood[j] = -1.0;

		isPileUpJetLoose[j]=false;
		isPileUpJetMedium[j]=false;
		isPileUpJetTight[j]=false;

	}
	// initialization done

	V2jMassMVAMET = -1.0;
	V2jMass = -1.0;
	V3jMass = -1.0;
	V4jMass = -1.0;
	V5jMass = -1.0;
	V6jMass = -1.0;
	c2jMass = -1.0;
	c3jMass = -1.0;
	c4jMass = -1.0;
	c5jMass = -1.0;
	c6jMass = -1.0;

	V2jCosJacksonAngle = -10.0;
	c2jCosJacksonAngle = -10.0;
	V3jCosJacksonAngle = -10.0;
	c3jCosJacksonAngle12 = -10.0;
	c3jCosJacksonAngle23 = -10.0;
	c3jCosJacksonAngle31 = -10.0;

	cosphiDecayPlane = 10.0; 
	cosThetaLnu = 10.0; 
	cosThetaJJ = 10.0;

	colorCorr01 = -10.0;
	colorCorr02 = -10.0;
	colorCorr12 = -10.0;
	colorCorr03 = -10.0;
	colorCorr13 = -10.0;
	colorCorr23 = -10.0;
	colorCorr04 = -10.0;
	colorCorr14 = -10.0;
	colorCorr24 = -10.0;
	colorCorr34 = -10.0;
	colorCorr05 = -10.0;
	colorCorr15 = -10.0;
	colorCorr25 = -10.0;
	colorCorr35 = -10.0;
	colorCorr45 = -10.0;

	j1Hel_HiggsCM = -10.0;;
	j2Hel_HiggsCM = -10.0;;
	l1Hel_HiggsCM = -10.0;;
	l2Hel_HiggsCM = -10.0;;
	b1Hel_HiggsCM = -10.0;;
	b2Hel_HiggsCM = -10.0;;

}


void ewk::JetTreeFiller::fill(const edm::Event& iEvent){
	// first initialize to the default values
	init();

	edm::Handle<reco::CandidateView> boson;
	iEvent.getByLabel( mInputBoson, boson);
	const int nBoson = boson->size();
	if( nBoson<1 ) return; // Nothing to fill

	const reco::Candidate *Vboson = &((*boson)[0]); 

	// in case when we have two candidates for the W boson in the event
	const reco::Candidate *Vboson2(0);
	if( nBoson==2) Vboson2  = &((*boson)[1]);

	edm::Handle<edm::View<reco::Jet> > jets;
	iEvent.getByLabel( mInputJets, jets ); 

	if(jets->size() < 1) return;

	/////// PfMET information /////
	edm::Handle<edm::View<reco::MET> > pfmet;
	iEvent.getByLabel(mInputMet, pfmet);

	/////// MVA MET information /////
	edm::Handle<edm::View<reco::MET> > metMVA;
	iEvent.getByLabel(mInputMetMVA, metMVA);


	// ----- b-tags ------------
	edm::Handle<reco::JetTagCollection> bTagHandle1;
	edm::Handle<reco::JetTagCollection> bTagHandle2;
	edm::Handle<reco::JetTagCollection> bTagHandle3;
	edm::Handle<reco::JetTagCollection> bTagHandle4;
	edm::Handle<reco::JetTagCollection> bTagHandle5;
	edm::Handle<reco::JetTagCollection> bTagHandle6;

	if(runoverAOD){
		// Get b tag information
		// ------------- SSV-HE ------------------------
		iEvent.getByLabel("simpleSecondaryVertexHighEffBJetTags", bTagHandle1);
		// ------------- TC-HE ------------------------
		iEvent.getByLabel("trackCountingHighEffBJetTags", bTagHandle2);
		// ------------- CSV ------------------------
		iEvent.getByLabel("combinedSecondaryVertexBJetTags", bTagHandle3);
		// ------------- JP ------------------------
		iEvent.getByLabel("jetProbabilityBJetTags", bTagHandle4);
		// ------------- SSV-HP ------------------------
		iEvent.getByLabel("simpleSecondaryVertexHighPurBJetTags", bTagHandle5);
		// ------------- TC-HP ------------------------
		iEvent.getByLabel("trackCountingHighPurBJetTags", bTagHandle6);
	}
	const reco::JetTagCollection  &  bTagsSSVHE = *(bTagHandle1.product()) ;
	const reco::JetTagCollection  &  bTagsTCHE = *(bTagHandle2.product()) ;
	const reco::JetTagCollection  &  bTagsCSV = *(bTagHandle3.product()) ;
	const reco::JetTagCollection  &  bTagsJP = *(bTagHandle4.product());
	const reco::JetTagCollection  &  bTagsSSVHP = *(bTagHandle5.product());
	const reco::JetTagCollection  &  bTagsTCHP = *(bTagHandle6.product()) ;

	edm::Handle<reco::SecondaryVertexTagInfoCollection> svTagInfos;
	if(runoverAOD){
		iEvent.getByLabel("secondaryVertexTagInfos", svTagInfos);
	}

	size_t iJet = 0;
	NumJets = 0;
	numBTags = 0;


	/////// Pileup density "rho" in the event from fastJet pileup calculation /////
	float fastjet_rho = -999999.9;
	edm::Handle<double> rho;
	const edm::InputTag eventrho("kt6PFJetsPFlow", "rho");
	iEvent.getByLabel(eventrho,rho);
	if( *rho == *rho) fastjet_rho = *rho;

	//   // get PFCandidates
	//   edm::Handle<reco::PFCandidateCollection>  PFCandidates;
	//   if(runoverAOD) iEvent.getByLabel("particleFlow", PFCandidates);

	// Fill pile Up jet id info
	if(mInputJets.label()!="Gen") fillPileUpJetID (jets);     


	// Loop over reco jets 
	edm::View<reco::Jet>::const_iterator jet, endpjets = jets->end(); 
	for (jet = jets->begin();  jet != endpjets;  ++jet, ++iJet) {
		if( !(iJet< (unsigned int) NUM_JET_MAX) ) break;

		//--------- store jet 4-vectors and related stuff ---
		fillBasicJetQuantities(iJet, *jet, Vboson, Vboson2, (*pfmet)[0]);
		// --------- Fill b-tag information -----------
		if(runoverAOD){
		  fillBtagInfoRECO( iJet, svTagInfos, bTagsSSVHE, bTagsTCHE, bTagsCSV, bTagsJP, bTagsSSVHP, bTagsTCHP);
		}else { 
			edm::Ptr<reco::Jet> ptrJet = jets->ptrAt( jet - jets->begin() );		  
			if ( ptrJet.isNonnull() && ptrJet.isAvailable() ) {
				const pat::Jet* pjet = dynamic_cast<const pat::Jet *>(ptrJet.get()) ;
				fillBtagInfoPAT( iJet, pjet);
			}
		}

		// ------- fill energy fractions --------------------	
		const std::type_info & type = typeid(*jet);
		if ( type == typeid(reco::PFJet) || type == typeid(pat::Jet)) {

			// PFJet specific quantities
			std::vector<reco::PFCandidatePtr> pfCandidates;
			if(type == typeid(reco::PFJet)) {
				reco::PFJet pfjet  = static_cast<const reco::PFJet &>(*jet);
				fillEnergyFractionsPFjets(pfjet, iJet);
				pfCandidates = pfjet.getPFConstituents();
			}
			if(type == typeid(pat::Jet)) {
				pat::Jet pfjet  = static_cast<const pat::Jet &>(*jet);
				if(pfjet.isPFJet()) {
					fillEnergyFractionsPFjets(pfjet, iJet);
					pfCandidates = pfjet.getPFConstituents();
				}
			}

			// ------- Compute pt_D and Quark Gluon Likelihood		 
			fillQGLH(iJet, fastjet_rho, pfCandidates);

		}// close PF jets loop

	}// close jets iteration loop

	NumJets = (int) iJet;


	// get 4-vectors for the two daughters of vector boson
	TLorentzVector p4lepton1;
	TLorentzVector p4lepton2;


	// 4body mass with MVA MET
	if(metMVA->size() > 0 ) {
		computeLeptonAndNu4Vectors(Vboson, (*metMVA)[0], p4lepton1, p4lepton2);
		fillInvariantMasses(p4lepton1, p4lepton2);
	} // second pass will overwrite everithing, but 4-body mass with MVA MET. 

	computeLeptonAndNu4Vectors(Vboson, (*pfmet)[0], p4lepton1, p4lepton2);



	// Now compute all the invariant masses
	fillInvariantMasses(p4lepton1, p4lepton2);


	// Angle between the decay planes of two W
	cosphiDecayPlane = 10.0; 
	cosThetaLnu = 10.0; 
	cosThetaJJ = 10.0;
	if( NumJets>1 ) { 
		TLorentzVector p4j1;
		TLorentzVector p4j2;
		p4j1.SetPxPyPzE( Px[0], Py[0], Pz[0], E[0] );
		p4j2.SetPxPyPzE( Px[1], Py[1], Pz[1], E[1] );
		dg_kin_Wuv_Wjj( p4lepton1, p4lepton2, 
					p4j1, p4j2, cosphiDecayPlane, 
					cosThetaLnu, cosThetaJJ);
	}


	// Color correlation between two leading jets ( jets pull ) 
	if( jetType_!="Gen") {
		if(NumJets>1) colorCorr01 = TMath::Abs( getDeltaTheta(  &(*jets)[0], &(*jets)[1]) );
		if(NumJets>2) {
			colorCorr02 = TMath::Abs( getDeltaTheta(  &(*jets)[0], &(*jets)[2]) );
			colorCorr12 = TMath::Abs( getDeltaTheta(  &(*jets)[1], &(*jets)[2]) );
		}
		if(NumJets>3) {
			colorCorr03 = TMath::Abs( getDeltaTheta(  &(*jets)[0], &(*jets)[3]) );
			colorCorr13 = TMath::Abs( getDeltaTheta(  &(*jets)[1], &(*jets)[3]) );
			colorCorr23 = TMath::Abs( getDeltaTheta(  &(*jets)[2], &(*jets)[3]) );
		}
		if(NumJets>4) {
			colorCorr04 = TMath::Abs( getDeltaTheta(  &(*jets)[0], &(*jets)[4]) );
			colorCorr14 = TMath::Abs( getDeltaTheta(  &(*jets)[1], &(*jets)[4]) );
			colorCorr24 = TMath::Abs( getDeltaTheta(  &(*jets)[2], &(*jets)[4]) );
			colorCorr34 = TMath::Abs( getDeltaTheta(  &(*jets)[3], &(*jets)[4]) );
		}
		if(NumJets>5) {
			colorCorr05 = TMath::Abs( getDeltaTheta(  &(*jets)[0], &(*jets)[5]) );
			colorCorr15 = TMath::Abs( getDeltaTheta(  &(*jets)[1], &(*jets)[5]) );
			colorCorr25 = TMath::Abs( getDeltaTheta(  &(*jets)[2], &(*jets)[5]) );
			colorCorr35 = TMath::Abs( getDeltaTheta(  &(*jets)[3], &(*jets)[5]) );
			colorCorr45 = TMath::Abs( getDeltaTheta(  &(*jets)[4], &(*jets)[5]) );
		}
	} // not Gen Jet


	// Cos(theta*) or Helicity Angles in Higgs rest frame
	fillHelicityIn4bodyFrame( p4lepton1, p4lepton2);


	//FillBranches();
}




//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////--------- store jet 4-vectors ---
template <typename T1> 
void ewk::JetTreeFiller::fillBasicJetQuantities(int iJet, 
			const T1& pfjet,
			const reco::Candidate* Vboson, 
			const reco::Candidate* Vboson2, 
			const reco::MET met) 
{
	Et[iJet] = pfjet.et();
	Pt[iJet] = pfjet.pt();
	Eta[iJet] = pfjet.eta();
	Phi[iJet] = pfjet.phi();
	Theta[iJet] = pfjet.theta();
	Px[iJet] = pfjet.px();
	Py[iJet] = pfjet.py();
	Pz[iJet] = pfjet.pz();
	E[iJet]  = pfjet.energy();
	Y[iJet]  = pfjet.rapidity();
	Mass[iJet] = pfjet.mass();
	Area[iJet] = pfjet.jetArea();



	Dphi[iJet] = dPhi( pfjet.phi(), Vboson->phi() );
	Deta[iJet] = fabs( pfjet.eta() - Vboson->eta() );
	DR[iJet] = radius( pfjet.eta(), pfjet.phi(), 
				Vboson->eta(), Vboson->phi());
	DphiMET[iJet] = dPhi( pfjet.phi(), met.phi() );

	Response[iJet] = 10.0;
	float vpt = Vboson->pt();
	if( vpt>0.0 ) Response[iJet] = pfjet.pt() / vpt;

	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > qstar 
		= pfjet.p4() + Vboson->p4();
	VjetMass[iJet] = qstar.M();



	// ------- Compute deltaR (jet, vector boson) -----
	if( !(Vboson2==0) ) {
		DR2[iJet] = radius( pfjet.eta(), pfjet.phi(), 
					Vboson2->eta(), Vboson2->phi());
		Dphi2[iJet] = dPhi( pfjet.phi(), Vboson2->phi() );
		Deta2[iJet] = fabs( pfjet.eta() - Vboson2->eta() );
		Response2[iJet] = 10.0;
		vpt = Vboson2->pt();
		if( vpt>0.0 ) Response2[iJet] = pfjet.pt() / vpt;	
		qstar = pfjet.p4() + Vboson2->p4();
		VjetMass2[iJet] = qstar.M();	
	}



	if(runoverAOD){
		etaetaMoment[iJet]  = pfjet.etaetaMoment();
		phiphiMoment[iJet]  = pfjet.phiphiMoment();      
		etaphiMoment[iJet]  = pfjet.etaphiMoment();
		maxDistance[iJet]  = pfjet.maxDistance();
		nConstituents[iJet]  = pfjet.nConstituents();
	}
}
template 
void ewk::JetTreeFiller::fillBasicJetQuantities(int, const reco::PFJet&,
			const reco::Candidate* Vboson, 
			const reco::Candidate* Vboson2, 
			const reco::MET met);
template 
void ewk::JetTreeFiller::fillBasicJetQuantities(int, const pat::Jet&,
			const reco::Candidate* Vboson, 
			const reco::Candidate* Vboson2, 
			const reco::MET met); 


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////--------- Energy fractions for generic PF jets (RECO or PAT) ---

template <typename T1>
void ewk::JetTreeFiller::fillEnergyFractionsPFjets(const T1& pfjet, 
			int iJet) 
{

	PFChargedHadronEnergy[iJet] = pfjet.chargedHadronEnergy();
	PFChargedHadronEnergyFraction[iJet] = pfjet.chargedHadronEnergyFraction ();
	PFNeutralHadronEnergy[iJet] = pfjet.neutralHadronEnergy();
	PFNeutralHadronEnergyFraction[iJet] = pfjet.neutralHadronEnergyFraction ();
	PFChargedEmEnergy[iJet] = pfjet.chargedEmEnergy ();
	PFChargedEmEnergyFraction[iJet] = pfjet.chargedEmEnergyFraction ();
	PFChargedMuEnergy[iJet] = pfjet.chargedMuEnergy ();
	PFChargedMuEnergyFraction[iJet] = pfjet.chargedMuEnergyFraction ();
	PFNeutralEmEnergy[iJet] = pfjet.neutralEmEnergy ();
	PFNeutralEmEnergyFraction[iJet] = pfjet.neutralEmEnergyFraction ();
	PFChargedMultiplicity[iJet] = pfjet.chargedMultiplicity();
	PFNeutralMultiplicity[iJet] = pfjet.neutralMultiplicity();
	PFMuonMultiplicity[iJet] = pfjet.muonMultiplicity();
	PFPhotonEnergy[iJet] = pfjet.photonEnergy();
	PFPhotonEnergyFraction[iJet] = pfjet.photonEnergyFraction();
	PFElectronEnergy[iJet] = pfjet.electronEnergy();
	PFElectronEnergyFraction[iJet] = pfjet.electronEnergyFraction();
	PFMuonEnergy[iJet] = pfjet.muonEnergy();
	PFMuonEnergyFraction[iJet] = pfjet.muonEnergyFraction();
	PFHFHadronEnergy[iJet] = pfjet.HFHadronEnergy();
	PFHFHadronEnergyFraction[iJet] = pfjet.HFHadronEnergyFraction();
	PFHFEMEnergy[iJet] = pfjet.HFEMEnergy();
	PFHFEMEnergyFraction[iJet] = pfjet.HFEMEnergyFraction();	 
	PFChargedHadronMultiplicity[iJet] = pfjet.chargedHadronMultiplicity();
	PFNeutralHadronMultiplicity[iJet] = pfjet.neutralHadronMultiplicity();
	PFPhotonMultiplicity[iJet] = pfjet.photonMultiplicity();
	PFElectronMultiplicity[iJet] = pfjet.electronMultiplicity();
	PFHFHadronMultiplicity[iJet] = pfjet.HFHadronMultiplicity();
	PFHFEMMultiplicity[iJet] = pfjet.HFEMMultiplicity();
	//  pfCandidates = pfjet.getPFConstituents();
}
// ----- explicit function template instantiation(s) ---
template 
void ewk::JetTreeFiller::fillEnergyFractionsPFjets(const reco::PFJet&, int);
template 
void ewk::JetTreeFiller::fillEnergyFractionsPFjets(const pat::Jet&, int);






//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////--------- Compute Quark-Gluon likelihood ---

void ewk::JetTreeFiller::fillQGLH(int iJet, float fastjet_rho, 
			std::vector<reco::PFCandidatePtr> pfCandidates) 
{
	// ------- Compute pt_D and Quark Gluon Likelihood		  
	float sumPt_cands=0.;
	float sumPt2_cands=0.;
	float rms_cands=0.;

	typedef std::vector<reco::PFCandidatePtr>::const_iterator IC;

	for (IC jt = pfCandidates.begin(); // If pfCandidates has no constituents then the loop simply won't execute
				jt != pfCandidates.end(); jt++) { // and so no segmentation fault should occur
		const reco::PFCandidatePtr pfCandPtr = *jt;
		if ( !(pfCandPtr.isNonnull() && pfCandPtr.isAvailable()) ) continue;

		// reco::PFCandidate::ParticleType id = (*jt)->particleId();
		// Convert particle momentum to normal TLorentzVector, wrong type :(
		math::XYZTLorentzVectorD const& p4t = pfCandPtr->p4();

		TLorentzVector p4(p4t.px(), p4t.py(), p4t.pz(), p4t.energy());
		TLorentzVector jetp4;
		jetp4.SetPtEtaPhiE(Pt[iJet], Eta[iJet], Phi[iJet], E[iJet]);
		if(p4.Pt()!=0){
			sumPt_cands += p4.Pt();
			sumPt2_cands += (p4.Pt()*p4.Pt());
			float deltaR = jetp4.DeltaR(p4);
			rms_cands += (p4.Pt()*p4.Pt()*deltaR*deltaR);
		}
	}			  
	PFsumPtCands[iJet]  = sumPt_cands;
	PFsumPt2Cands[iJet] = sumPt2_cands;
	if(sumPt_cands != 0)  PFptD[iJet] = sqrt( sumPt2_cands )/sumPt_cands;
	if(rms_cands  != 0)   PFrmsCands[iJet] = rms_cands/sumPt2_cands;
	PFqgLikelihood[iJet]= qglikeli->computeQGLikelihoodPU( Pt[iJet], fastjet_rho, 
				PFChargedMultiplicity[iJet], 
				PFNeutralMultiplicity[iJet], 
				PFptD[iJet]);	 
}





//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////--------- Fill b-tag information from RECO ---

void ewk::JetTreeFiller::fillBtagInfoRECO(int iJet, 
			edm::Handle<reco::SecondaryVertexTagInfoCollection> svTagInfos,
			const reco::JetTagCollection  &  bTagsSSVHE,
			const reco::JetTagCollection  &  bTagsTCHE,
			const reco::JetTagCollection  &  bTagsCSV,
			const reco::JetTagCollection  &  bTagsJP,
			const reco::JetTagCollection  &  bTagsSSVHP,
			const reco::JetTagCollection  &  bTagsTCHP) 
{
	// study b tag info.

	double dist = 100000.0;
	double closestDistance = 100000.0;
	unsigned int closestIndex = 10000;

	// ------------- SSV-HE ------------------------
	// compute B-tag discriminator
	for (unsigned int i = 0; i != bTagsSSVHE.size(); ++i) {
		edm::RefToBase<reco::Jet> aJet  = bTagsSSVHE[i].first;   
		dist = radius(aJet->eta(), aJet->phi(), Eta[iJet], Phi[iJet]);
		if( dist < closestDistance ) { 
			closestDistance = dist;
			closestIndex = i;
		}
	}

	// std::cout << " ++++++++++++++++++++  jbTagsSSVHE.size() " <<  bTagsSSVHE.size() <<std::endl;
	if( closestDistance<0.3 && closestIndex<bTagsSSVHE.size() ) {
		bDiscriminator[iJet] = bTagsSSVHE[closestIndex].second;
		bDiscriminatorSSVHE[iJet] = bTagsSSVHE[closestIndex].second;
		const reco::SecondaryVertexTagInfo svTagInfo = (*svTagInfos)[closestIndex];
		if (  svTagInfo.nVertices() > 0  && 
					bDiscriminator[iJet]>-1.0) {
			if(bDiscriminator[iJet]>1.74) numBTags ++;

			///////////////////////////
			// Calculate SecVtx Mass //
			///////////////////////////

			ROOT::Math::LorentzVector< ROOT::Math::PxPyPzM4D<double> > sumVec;
			reco::CompositeCandidate vertexCand;
			reco::Vertex::trackRef_iterator 
				kEndTracks = svTagInfo.secondaryVertex(0).tracks_end();
			for (reco::Vertex::trackRef_iterator trackIter = 
						svTagInfo.secondaryVertex(0).tracks_begin(); 
						trackIter != kEndTracks; 
						++trackIter ) 
			{
				const double kPionMass = 0.13957018;
				ROOT::Math::LorentzVector< ROOT::Math::PxPyPzM4D<double> > vec;
				vec.SetPx( (*trackIter)->px() );
				vec.SetPy( (*trackIter)->py() );
				vec.SetPz( (*trackIter)->pz() );
				vec.SetM (kPionMass);
				sumVec += vec;
			} // for trackIter
			secVertexMass[iJet] = sumVec.M();
		} // endif svTagInfo.nVertices condition
	}// endif closestDistance condition				      


	// ------------- TC-HE ------------------------
	bDiscriminatorTCHE[iJet] = BTAG_DISCRIM_DEFAULT;
	closestDistance = 100000.0;
	closestIndex = 10000;

	for (unsigned int i = 0; i != bTagsTCHE.size(); ++i) {
		edm::RefToBase<reco::Jet> aJet  = bTagsTCHE[i].first;   
		dist = radius(aJet->eta(), aJet->phi(), Eta[iJet], Phi[iJet]);
		if( dist < closestDistance ) { 
			closestDistance = dist;
			closestIndex = i;
		}
	}
	if( closestDistance<0.3 && closestIndex<bTagsTCHE.size() )
	  bDiscriminatorTCHE[iJet] = bTagsTCHE[closestIndex].second;

	// ------------- CSV ------------------------
	bDiscriminatorCSV[iJet] = BTAG_DISCRIM_DEFAULT;
	closestDistance = 100000.0;
	closestIndex = 10000;

	for (unsigned int i = 0; i != bTagsCSV.size(); ++i) {
		edm::RefToBase<reco::Jet> aJet  = bTagsCSV[i].first;   
		dist = radius(aJet->eta(), aJet->phi(), Eta[iJet], Phi[iJet]);
		if( dist < closestDistance ) { 
			closestDistance = dist;
			closestIndex = i;
		}
	}
	if( closestDistance<0.3 && closestIndex<bTagsCSV.size() )
	  bDiscriminatorCSV[iJet] = bTagsCSV[closestIndex].second;

	// ------------- JP ------------------------
	bDiscriminatorJP[iJet] = BTAG_DISCRIM_DEFAULT;
	closestDistance = 100000.0;
	closestIndex = 10000;

	for (unsigned int i = 0; i != bTagsJP.size(); ++i) {
		edm::RefToBase<reco::Jet> aJet  = bTagsJP[i].first;   
		dist = radius(aJet->eta(), aJet->phi(), Eta[iJet], Phi[iJet]);
		if( dist < closestDistance ) { 
			closestDistance = dist;
			closestIndex = i;
		}
	}
	if( closestDistance<0.3 && closestIndex<bTagsJP.size() )
	  bDiscriminatorJP[iJet] = bTagsJP[closestIndex].second;

	// ------------- SSV-HP ------------------------
	bDiscriminatorSSVHP[iJet] = BTAG_DISCRIM_DEFAULT;
	closestDistance = 100000.0;
	closestIndex = 10000;

	for (unsigned int i = 0; i != bTagsSSVHP.size(); ++i) {
		edm::RefToBase<reco::Jet> aJet  = bTagsSSVHP[i].first;   
		dist = radius(aJet->eta(), aJet->phi(), Eta[iJet], Phi[iJet]);
		if( dist < closestDistance ) { 
			closestDistance = dist;
			closestIndex = i;
		}
	}
	if( closestDistance<0.3 && closestIndex<bTagsSSVHP.size() )
	  bDiscriminatorSSVHP[iJet] = bTagsSSVHP[closestIndex].second;

	// ------------- TC-HP ------------------------
	bDiscriminatorTCHP[iJet] = BTAG_DISCRIM_DEFAULT;
	closestDistance = 100000.0;
	closestIndex = 10000;

	for (unsigned int i = 0; i != bTagsTCHP.size(); ++i) {
		edm::RefToBase<reco::Jet> aJet  = bTagsTCHP[i].first;   
		dist = radius(aJet->eta(), aJet->phi(), Eta[iJet], Phi[iJet]);
		if( dist < closestDistance ) { 
			closestDistance = dist;
			closestIndex = i;
		}
	}
	if( closestDistance<0.3 && closestIndex<bTagsTCHP.size() )
	  bDiscriminatorTCHP[iJet] = bTagsTCHP[closestIndex].second;

}





//////--------- Fill b-tag information from PAT ---

void ewk::JetTreeFiller::fillBtagInfoPAT(int iJet, const pat::Jet* pjet) 
{
	if(pjet !=0)
	{
		bDiscriminatorSSVHE[iJet] = (*pjet).bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
		bDiscriminatorTCHE[iJet] = (*pjet).bDiscriminator("trackCountingHighEffBJetTags");
		bDiscriminatorCSV[iJet] = (*pjet).bDiscriminator("combinedSecondaryVertexBJetTags");
		bDiscriminatorJP[iJet] = (*pjet).bDiscriminator("jetProbabilityBJetTags");
		bDiscriminatorSSVHP[iJet] = (*pjet).bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
		bDiscriminatorTCHP[iJet] = (*pjet).bDiscriminator("trackCountingHighPurBJetTags");
		bDiscriminator[iJet] = bDiscriminatorSSVHE[iJet];
		if(bDiscriminator[iJet]>1.74) numBTags ++;
	}
}






//////--------- Helicity in 4-body frame ---

void ewk::JetTreeFiller::fillHelicityIn4bodyFrame(TLorentzVector& p4lepton1, 
			TLorentzVector& p4lepton2) 
{
	if( NumJets<2 ) return;

	TLorentzVector p4j1;
	TLorentzVector p4j2;
	p4j1.SetPxPyPzE( Px[0], Py[0], Pz[0], E[0] );
	p4j2.SetPxPyPzE( Px[1], Py[1], Pz[1], E[1] );

	TLorentzVector p4boson1 = p4lepton1 + p4lepton2;
	TLorentzVector p4boson2 = p4j1 + p4j2;
	TLorentzVector p4higgs  = p4boson1 + p4boson2;
	TVector3 higgsBoost = p4higgs.BoostVector();

	j1Hel_HiggsCM = getHelicity( p4j1 , higgsBoost );
	j2Hel_HiggsCM = getHelicity( p4j2 , higgsBoost );
	l1Hel_HiggsCM = getHelicity( p4lepton1 , higgsBoost );
	l2Hel_HiggsCM = getHelicity( p4lepton2 , higgsBoost );
	b1Hel_HiggsCM = getHelicity( p4boson1 , higgsBoost );
	b2Hel_HiggsCM = getHelicity( p4boson2 , higgsBoost );
}






//////--------- Invariant masses of Njets (N>1) and V+Njets (N>0) ---

void ewk::JetTreeFiller::fillInvariantMasses(TLorentzVector& p4lepton1, 
			TLorentzVector& p4lepton2) 
{

	// Now compute all the invariant masses
	TLorentzVector Vj;
	TLorentzVector V2j;
	TLorentzVector V3j;
	TLorentzVector V4j;
	TLorentzVector V5j;
	TLorentzVector V6j;

	TLorentzVector c2j;
	TLorentzVector c3j;
	TLorentzVector c4j;
	TLorentzVector c5j;
	TLorentzVector c6j;

	// 4-vectors for the first four jets
	TLorentzVector p4j1;
	TLorentzVector p4j2;
	TLorentzVector p4j3;
	TLorentzVector p4j4;
	TLorentzVector p4j5;
	TLorentzVector p4j6;
	TLorentzVector p4V = p4lepton1 + p4lepton2;

	if( NumJets>1 ) { 
		p4j1.SetPxPyPzE( Px[0], Py[0], Pz[0], E[0] );
		p4j2.SetPxPyPzE( Px[1], Py[1], Pz[1], E[1] );
		c2j =  p4j1 + p4j2;
		c2jMass =  c2j.M();
		V2j =  p4V + c2j;
		V2jMass =  V2j.M();
		if(V2jMassMVAMET<0) V2jMassMVAMET = V2jMass;  // fill only at first pass 
		V2jCosJacksonAngle = JacksonAngle( p4V, V2j);
		c2jCosJacksonAngle = JacksonAngle( p4j1, p4j2);
	}

	if( NumJets>2 ) {
		p4j3.SetPxPyPzE( Px[2], Py[2], Pz[2], E[2] );
		c3j =  p4j1 + p4j2 + p4j3;
		c3jMass =  c3j.M();
		V3j =  p4V + c3j;
		V3jMass =  V3j.M();
		V3jCosJacksonAngle = JacksonAngle( p4V, V3j);
		c3jCosJacksonAngle12 = JacksonAngle( p4j1,  p4j2 );
		c3jCosJacksonAngle23 = JacksonAngle( p4j2,  p4j3 );
		c3jCosJacksonAngle31 = JacksonAngle( p4j3,  p4j1 );
	}

	if( NumJets>3 ) { 
		p4j4.SetPxPyPzE( Px[3], Py[3], Pz[3], E[3] );
		c4j =  p4j1 + p4j2 + p4j3 + p4j4;
		c4jMass =  c4j.M();
		V4j =  p4V + c4j;
		V4jMass =  V4j.M();

	}

	if( NumJets>4 ) { 
		p4j5.SetPxPyPzE( Px[4], Py[4], Pz[4], E[4] );
		c5j =  p4j1 + p4j2 + p4j3 + p4j4 + p4j5;
		c5jMass =  c5j.M();
		V5j =  p4V + c5j;
		V5jMass =  V5j.M();

	}

	if( NumJets>5 ) { 
		p4j6.SetPxPyPzE( Px[5], Py[5], Pz[5], E[5] );

		c6j =  p4j1 + p4j2 + p4j3 + p4j4 + p4j5 + p4j6;
		c6jMass =  c6j.M();
		V6j =  p4V + c6j;
		V6jMass =  V6j.M();
	}
}





//////--------- compute lepton and neutrino 4-vectors --------

void ewk::JetTreeFiller::computeLeptonAndNu4Vectors(const reco::Candidate* Vboson, 
			const reco::MET met, 
			TLorentzVector& p4lepton1, 
			TLorentzVector& p4lepton2) 
{
	const reco::Candidate* m0 = Vboson->daughter(0);
	const reco::Candidate* m1 = Vboson->daughter(1);
	TLorentzVector p4MET;
	TLorentzVector p4lepton;
	METzCalculator* metz = new METzCalculator();
	if (LeptonType_=="electron") metz->SetLeptonType("electron");
	double nupz;


	// Compute pz if one of the lepton is neutrino
	if( m0->isElectron() || m0->isMuon() ) 
	  p4lepton1.SetPxPyPzE(m0->px(), m0->py(), m0->pz(), m0->energy());
	else {
		p4MET.SetPxPyPzE(met.px(), met.py(), met.pz(), met.energy());
		p4lepton.SetPxPyPzE(m1->px(), m1->py(), m1->pz(), m1->energy());
		metz->SetMET(p4MET);
		metz->SetLepton(p4lepton);
		nupz = metz->Calculate();
		p4lepton1.SetPxPyPzE( m0->px(), m0->py(), nupz, sqrt(m0->px()*m0->px()+m0->py()*m0->py()+nupz*nupz) );
	}

	if( m1->isElectron() || m1->isMuon() ) 
	  p4lepton2.SetPxPyPzE(m1->px(), m1->py(), m1->pz(), m1->energy());
	else {
		p4MET.SetPxPyPzE(met.px(), met.py(), met.pz(), met.energy());
		p4lepton.SetPxPyPzE(m0->px(), m0->py(), m0->pz(), m0->energy());
		metz->SetMET(p4MET);
		metz->SetLepton(p4lepton);
		nupz = metz->Calculate();
		p4lepton2.SetPxPyPzE( m1->px(), m1->py(), nupz, sqrt(m1->px()*m1->px()+m1->py()*m1->py()+nupz*nupz) );
	}

	delete metz;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////--------- Pile Up ID information for generic PF jets (RECO or PAT) ---

template <typename T1>
void ewk::JetTreeFiller::fillPileUpJetID (const edm::Handle<edm::View<T1> > & jets) {

	unsigned int  iJet = 0;

	typename edm::View<T1>::const_iterator itJet, endpjets = jets->end();

	for (itJet = jets->begin();  itJet != endpjets;  ++iJet, ++itJet) {

		if( !(iJet< (unsigned int) NUM_JET_MAX) ) break;

		const std::type_info & type = typeid(*itJet);

		if ( type == typeid(pat::Jet)){

			edm::Ptr<T1> ptrJet = jets->ptrAt(itJet - jets->begin());

			if (ptrJet.isNonnull() && ptrJet.isAvailable() ) {


				const pat::Jet* pjet = dynamic_cast<const pat::Jet *>(ptrJet.get()) ;


				if((*pjet).userInt("PUChargedWorkingPoint") == 3 )                                     isPileUpJetTight[iJet] = true ;
				if((*pjet).userInt("PUChargedWorkingPoint") == 2 || isPileUpJetTight[iJet]  == true )  isPileUpJetMedium[iJet] = true ;
				if((*pjet).userInt("PUChargedWorkingPoint") == 1 || isPileUpJetMedium[iJet] == true )  isPileUpJetLoose[iJet] = true ;

			}
		}  
	}

}

// ----- explicit function template instantiation(s) ---
template
void ewk::JetTreeFiller::fillPileUpJetID (const edm::Handle<edm::View<pat::Jet> > &);



// ////////--------- Get constituents of a jet---------
// void 
// ewk::JetTreeFiller::getConstituents(const reco::Jet* jet, 
// 				    std::vector<reco::PFCandidatePtr> &jetpfc) {
//   const reco::PFJet* pfjet  = static_cast<const reco::PFJet *>(jet);
//   jetpfc = pfjet->getPFConstituents();
// }


// void  
// ewk::JetTreeFiller::getConstituents(const pat::Jet* jet, 
// 				    std::vector<reco::PFCandidatePtr> &jetpfc) {
//   const pat::Jet* pfjet  = static_cast<const pat::Jet *>(jet);
//   if(pfjet->isPFJet()) jetpfc = pfjet->getPFConstituents();
// }





// template <typename T> 
// std::vector<reco::PFCandidatePtr> ewk::JetTreeFiller::getConstituents(const T&  jet, bool isReco) {

//   std::vector<reco::PFCandidatePtr> jetpfc ;

//   if(isReco)  {
//     const reco::PFJet* pfjet  = static_cast<const reco::PFJet *>(jet);
//     jetpfc = pfjet->getPFConstituents(); 
//   }
//   else {
//     const pat::Jet* pfjet  = static_cast<const pat::Jet *>(jet);
//     if(pfjet->isPFJet()) jetpfc = pfjet->getPFConstituents();
//   } 
//   return jetpfc;
// }

// template 
// std::vector<reco::PFCandidatePtr> ewk::JetTreeFiller::getConstituents(const reco::Jet*, bool);
// template 
// std::vector<reco::PFCandidatePtr> ewk::JetTreeFiller::getConstituents(const pat::Jet*, bool);
