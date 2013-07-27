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


// user include files
#include "ElectroWeakAnalysis/VPlusJets/interface/VplusJetsAnalysis.h" 
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "JetSubstructure/SubstructureTools/interface/PseudoJetUserInfo.h"
#include "JetSubstructure/SubstructureTools/interface/JetSubstructureTools.h" 


ewk::VplusJetsAnalysis::VplusJetsAnalysis(const edm::ParameterSet& iConfig) :
	myTree ( fs -> mkdir("../").make<TTree>(iConfig.getParameter<std::string>("TreeName").c_str(),"V+jets Tree") ), 
	//CorrectedPFJetFiller ( iConfig.existsAs<edm::InputTag>("srcPFCor") ?  new JetTreeFiller("CorrectedPFJetFiller", myTree, "PFCor", iConfig) : 0),
	//PhotonFiller (  iConfig.existsAs<edm::InputTag>("srcPhoton") ?  new PhotonTreeFiller("PhotonFiller", myTree,  iConfig) : 0),
	recoBosonFillerE( new VtoElectronTreeFiller( iConfig.getParameter<std::string>("VBosonType").c_str(), myTree, iConfig) ),
	recoBosonFillerMu( new VtoMuonTreeFiller( iConfig.getParameter<std::string>("VBosonType").c_str(), myTree, iConfig) ),
	genBosonFiller( (iConfig.existsAs<bool>("runningOverMC") && iConfig.getParameter<bool>("runningOverMC")) ?  new MCTreeFiller(iConfig.getParameter<std::string>("VBosonType").c_str(), myTree, iConfig) : 0)
{
	// Are we running over Monte Carlo ?
	if( iConfig.existsAs<bool>("runningOverMC") ) runningOverMC_=iConfig.getParameter< bool >("runningOverMC");
	else runningOverMC_= false;

	if(  iConfig.existsAs<edm::InputTag>("srcVectorBoson") ) mInputBoson = iConfig.getParameter<edm::InputTag>("srcVectorBoson"); 
	LeptonType_ = iConfig.getParameter<std::string>("LeptonType");
	VBosonType_ = iConfig.getParameter<std::string>("VBosonType");

	if(  iConfig.existsAs<edm::InputTag>("srcPrimaryVertex") ) mPrimaryVertex = iConfig.getParameter<edm::InputTag>("srcPrimaryVertex"); 
	else mPrimaryVertex =  edm::InputTag("offlinePrimaryVertices");

	if(  iConfig.existsAs<edm::InputTag>("srcMet") ) mInputMet = iConfig.getParameter<edm::InputTag>("srcMet");

	//if(  iConfig.existsAs<edm::InputTag>("srcMetMVA") ) mInputMetMVA = iConfig.getParameter<edm::InputTag>("srcMetMVA");

	//*********************  Run Over AOD or PAT  ***********//
	if( iConfig.existsAs<bool>("runningOverAOD")) runoverAOD = iConfig.getParameter<bool>("runningOverAOD");

	//JetsFor_rho =  iConfig.getParameter<std::string>("srcJetsforRho") ; 

	//if(  iConfig.existsAs<edm::InputTag>("srcgenMet") ) mInputgenMet =  iConfig.getParameter<edm::InputTag>("srcgenMet") ; 

}



ewk::VplusJetsAnalysis::~VplusJetsAnalysis() {}


void ewk::VplusJetsAnalysis::beginJob() {

	// Declare all the branches of the tree
	declareTreeBranches();
	std::cout<<"VplusJetsAnalysis begin work!"<<std::endl;
}



// ------------ method called to produce the data  ------------
void ewk::VplusJetsAnalysis::analyze(const edm::Event& iEvent, 
			const edm::EventSetup& iSetup) {

	// write event information: run, event, bunch crossing, ....
	run   = iEvent.id().run();
	event = iEvent.id().event();
	lumi  = iEvent.luminosityBlock();
	bunch = iEvent.bunchCrossing();


	// primary/secondary vertices
	// edm::Handle<reco::VertexCollection > recVtxs;
	edm::Handle <edm::View<reco::Vertex> > recVtxs;
	iEvent.getByLabel( mPrimaryVertex, recVtxs);
	nPV = recVtxs->size();


	/////// PfMET information /////
	edm::Handle<edm::View<reco::MET> > pfmet;
	iEvent.getByLabel(mInputMet, pfmet);
	if (pfmet->size() == 0) {
		mpfMET   = -1;
		mpfSumET = -1;
		mpfMETSign = -1;
		mpfMETPhi   = -10.0;
	} else {
		mpfMET   = (*pfmet)[0].et();
		mpfSumET = (*pfmet)[0].sumEt();
		mpfMETSign = (*pfmet)[0].significance();
		mpfMETPhi   = (*pfmet)[0].phi();
	}


	/*/////// MVA MET information /////
	edm::Handle<edm::View<reco::MET> > metMVA;
	iEvent.getByLabel(mInputMetMVA, metMVA);
	if (metMVA->size() == 0) {
		mvaMET   = -1;
		mvaSumET = -1;
		mvaMETSign = -1;
		mvaMETPhi   = -10.0;
	} else {
		mvaMET   = (*metMVA)[0].et();
		mvaSumET = (*metMVA)[0].sumEt();
		mvaMETSign = (*metMVA)[0].significance();
		mvaMETPhi   = (*metMVA)[0].phi();
	}*/

	
	/*/////// Pileup density "rho" in the event from fastJet pileup calculation /////
	edm::Handle<double> rho;
	const edm::InputTag eventrho(JetsFor_rho, "rho");
	iEvent.getByLabel(eventrho,rho);
	if( *rho == *rho) fastJetRho = *rho;
	else  fastJetRho =  -999999.9;*/
	



	/////////// GenMET information & MC Pileup Summary Info  //////////
	mcPUtrueInteractions = 0;
	mcPUtotnvtx = 0;
	mcPUbx[0]   = -999; mcPUbx[1]   = -999; mcPUbx[2]   = -999;
	mcPUnvtx[0] = -999; mcPUnvtx[1] = -999; mcPUnvtx[2] = -999;
	if ( runningOverMC_ ){
		/*edm::Handle<reco::GenMETCollection> genMETs;
		if(!runoverAOD){
			iEvent.getByLabel(mInputgenMet,genMETs);
			if ( genMETs->size() == 0) {
				genMET   = -1.0;
				genSumET = -1.0;
				genMETSign  = -1.0;
				genMETPhi   = -10.0;
			} else {
				genMET = (*genMETs)[0].et();
				genSumET = (*genMETs)[0].sumEt();  
				genMETSign = (*genMETs)[0].significance();  
				genMETPhi = (*genMETs)[0].phi();
			}
		}*/
		// MC Pileup Summary Info
		const edm::InputTag PileupSrc("addPileupInfo");
		edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
		iEvent.getByLabel(PileupSrc, PupInfo);
		std::vector<PileupSummaryInfo>::const_iterator PVI;
		int ctid = 0;
		for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
			if (ctid>2) break;
			mcPUbx[ctid]   =  PVI->getBunchCrossing();
			mcPUnvtx[ctid] =  PVI->getPU_NumInteractions();
			mcPUtotnvtx   +=  PVI->getPU_NumInteractions();
			if(PVI->getBunchCrossing() == 0)mcPUtrueInteractions = PVI->getTrueNumInteractions();
			ctid++;
			cout<<"bx="<<PVI->getBunchCrossing()<<" getTrueNumInteractions="<<PVI->getTrueNumInteractions()<<" getPU_NumInteractions()="<<PVI->getPU_NumInteractions()<<endl;
		}
	}
	// fill jet branches
	edm::Handle<edm::View< reco::Candidate> > boson;
	iEvent.getByLabel( mInputBoson, boson);
	mNVB = boson->size();
	if( mNVB<1 ) return; // Nothing to fill


	if(GenJetFiller.get()) GenJetFiller->fill(iEvent);
	//if(PhotonFiller.get()) PhotonFiller->fill(iEvent);

	if(CorrectedPFJetFiller.get()) CorrectedPFJetFiller->fill(iEvent);
	if(CorrectedPFJetFillerVBFTag.get()) CorrectedPFJetFillerVBFTag->fill(iEvent);//For VBF Tag Jets


	/**  Store groomed jet information */
	//if(AK5groomedJetFiller.get()) AK5groomedJetFiller->fill(iEvent);
	//if(AK7groomedJetFiller.get()) AK7groomedJetFiller->fill(iEvent);
	//if(AK8groomedJetFiller.get()) AK8groomedJetFiller->fill(iEvent);
	//if(CA8groomedJetFiller.get()) CA8groomedJetFiller->fill(iEvent);
	//if(CA12groomedJetFiller.get()) CA12groomedJetFiller->fill(iEvent);
	//if(genAK5groomedJetFiller.get()) genAK5groomedJetFiller->fill(iEvent);
	//if(genAK7groomedJetFiller.get()) genAK7groomedJetFiller->fill(iEvent);
	//if(genAK8groomedJetFiller.get()) genAK8groomedJetFiller->fill(iEvent);
	//if(genCA8groomedJetFiller.get()) genCA8groomedJetFiller->fill(iEvent);
	//if(genCA12groomedJetFiller.get()) genCA12groomedJetFiller->fill(iEvent);



	/**  Store reconstructed vector boson information */
	recoBosonFillerE->fill(iEvent, 0);
	// if(mNVB==2) recoBosonFillerE->fill(iEvent, 1);

	recoBosonFillerMu->fill(iEvent,0);
	// if(mNVB==2) recoBosonFillerMu->fill(iEvent, 1);


	/**  Store generated vector boson information */
	if(genBosonFiller.get()) genBosonFiller->fill(iEvent);


	//Begin Jet substructure work
    // ---------------------------------------------------------------------
    // reading in different collections
    // ---------------------------------------------------------------------
    
    edm::Handle< std::vector< reco::PFCandidate > >  pfs;
    iEvent.getByLabel("pfNoElectronPFlow", pfs);
    std::cout << "pfs->size() = " << pfs->size() << std::endl;
    
    edm::Handle< std::vector< reco::GenParticle > >  gens;
    iEvent.getByLabel("genParticles", gens);
    std::cout << "gens->size() = " << gens->size() << std::endl;

    edm::Handle < reco::GenParticleRefVector > gens_nonu;
    iEvent.getByLabel("genParticlesForJetsNoNu", gens_nonu);
    std::cout << "gens_nonu->size() = " << gens_nonu->size() << std::endl;
    
    // ---------------------------------------------------------------------
    // filling in pseudojets to be used as input for clustering
    // ---------------------------------------------------------------------
    
    std::vector<fastjet::PseudoJet> fjinputs;
    fjinputs.clear();

    for (unsigned int i = 0; i < gens_nonu->size(); i++){
        
        const reco::GenParticle& P = *(gens_nonu->at(i));
        fastjet::PseudoJet tmp_psjet( P.px(), P.py(), P.pz(), P.energy() );
        // add user info about charge and pdgId
        tmp_psjet.set_user_info( new PseudoJetUserInfo(P.pdgId(), P.charge()) );
        fjinputs.push_back( tmp_psjet );
        
    }
    
    // ---------------------------------------------------------------------
    // recluster on the fly
    // ---------------------------------------------------------------------
    double mJetRadius = 0.5;
    //fastjet::JetDefinition jetDef(fastjet::cambridge_algorithm, mJetRadius);
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, mJetRadius);
    
    int activeAreaRepeats = 1;
    double ghostArea = 0.01;
    double ghostEtaMax = 7.0;
    
    fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
    fastjet::AreaDefinition fjAreaDefinition( fastjet::active_area, fjActiveArea );
    fastjet::AreaDefinition fjAreaDefinition_wGhosts( fastjet::active_area_explicit_ghosts, fjActiveArea );
    
    fastjet::ClusterSequenceArea* thisClustering_wGhosts_ = new fastjet::ClusterSequenceArea(fjinputs, jetDef, fjAreaDefinition_wGhosts);
    fastjet::ClusterSequence* thisClustering_basic_ = new fastjet::ClusterSequence(fjinputs, jetDef);
    
    std::vector<fastjet::PseudoJet> out_jets_wGhosts_ = sorted_by_pt(thisClustering_wGhosts_->inclusive_jets(25.0));    
    std::vector<fastjet::PseudoJet> out_jets_basic_ = sorted_by_pt(thisClustering_basic_->inclusive_jets(25.0));
    
    std::cout << "out_jets_wGhosts_->size() = " << out_jets_wGhosts_.size() << std::endl;
    std::cout << "out_jets_basic_->size() = " << out_jets_basic_.size() << std::endl;    
    
    // ---------------------------------------------------------------------
    // now pass on to a class which stores all the substructure information you would want...
    // inputs are jet definition, vector< PseudoJet >, area 
    //    - area, because you need to compute the area using all the constituents in an event
    //    - jet definition, because need to recluster the single jet to recover the cluster sequence
    // ---------------------------------------------------------------------
    
    //JetSubstructureTools *tmp = new JetSubstructureTools( jetDef, out_jets_basic_[0].constituents(), out_jets_wGhosts_[0].area() );

    for ( unsigned int i = 0; i < out_jets_basic_.size(); i++ ){
     
        std::cout << "out_jets_wGhosts_[i].constituents().size() = " << out_jets_wGhosts_[i].constituents().size() << std::endl;
        std::cout << "out_jets_basic_[i].constituents().size() = " << out_jets_basic_[i].constituents().size() << std::endl;        
        
        std::cout << "pt = " << out_jets_basic_[i].pt() << ", " << out_jets_wGhosts_[i].pt() << std::endl;
        
        JetSubstructureTools *tmp = new JetSubstructureTools( jetDef, out_jets_basic_[i].constituents(), out_jets_wGhosts_[i].area() );

        std::cout << "obj pt = " << tmp->getJet().pt() << ", jet charge = " << tmp->getJetCharge( 0.3 ) << std::endl;
        
    }
    
   
    
    



	myTree->Fill();

} // analyze method






void ewk::VplusJetsAnalysis::endJob()
{
	std::cout<<"VplusJetsAnalysis end!"<<std::endl;
}




//  **** Utility: declare TTree branches for ntuple variables ***
void ewk::VplusJetsAnalysis::declareTreeBranches() {

	myTree->Branch("event_runNo",  &run,   "event_runNo/I");
	myTree->Branch("event_evtNo",  &event, "event_evtNo/I");
	myTree->Branch("event_lumi",   &lumi,  "event_lumi/I"); 
	myTree->Branch("event_bunch",  &bunch, "event_bunch/I"); 
	myTree->Branch("event_nPV",    &nPV,   "event_nPV/I"); 
	myTree->Branch("event_met_pfmet",             &mpfMET,      "event_met_pfmet/F"); 
	myTree->Branch("event_met_pfsumet",           &mpfSumET,    "event_met_pfsumet/F"); 
	myTree->Branch("event_met_pfmetsignificance", &mpfMETSign,  "event_met_pfmetsignificance/F"); 
	myTree->Branch("event_met_pfmetPhi",          &mpfMETPhi,   "event_met_pfmetPhi/F"); 
	myTree->Branch("event_metMVA_met",             &mvaMET,      "event_metMVA_met/F"); 
	myTree->Branch("event_metMVA_sumet",           &mvaSumET,    "event_metMVA_sumet/F"); 
	myTree->Branch("event_metMVA_metsignificance", &mvaMETSign,  "event_metMVA_metsignificance/F"); 
	myTree->Branch("event_metMVA_metPhi",          &mvaMETPhi,   "event_metMVA_metPhi/F"); 

	myTree->Branch("event_fastJetRho",                &fastJetRho,    "event_fastJetRho/F"); 

	if ( runningOverMC_ ){
		myTree->Branch("event_met_genmet",    &genMET,  "event_met_genmet/F"); 
		myTree->Branch("event_met_gensumet",  &genSumET,"event_met_gensumet/F"); 
		myTree->Branch("event_met_genmetsignificance", &genMETSign,  "event_met_genmetsignificance/F"); 
		myTree->Branch("event_met_genmetPhi",    &genMETPhi,  "event_met_genmetPhi/F"); 	  
		myTree->Branch("event_mcPU_totnvtx",    &mcPUtotnvtx,  "event_mcPU_totnvtx/F"); 
		myTree->Branch("event_mcPU_trueInteractions",    &mcPUtrueInteractions,  "event_mcPU_trueInteractions/F"); 
		myTree->Branch("event_mcPU_bx",         mcPUbx ,       "event_mcPU_bx[3]/F"); 
		myTree->Branch("event_mcPU_nvtx",       mcPUnvtx,      "event_mcPU_nvtx[3]/F"); 
	}

}  





// declare this class as a plugin
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
using ewk::VplusJetsAnalysis;
DEFINE_FWK_MODULE(VplusJetsAnalysis);
