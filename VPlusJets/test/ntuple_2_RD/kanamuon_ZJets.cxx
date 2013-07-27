#define kanamuon_ZJets_cxx
#include "kanamuon_ZJets.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>
#include <algorithm>
#include <string>
#include <sstream>
#include <map>
#include <TMap.h>
//#include "LOTable.h"

#include "Resolution.h"
#include "PhysicsTools/KinFitter/interface/TFitConstraintMGaus.h"
#include "PhysicsTools/KinFitter/interface/TFitConstraintM.h"
#include "PhysicsTools/KinFitter/interface/TFitConstraintEp.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEtEtaPhi.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleCart.h"
#include "PhysicsTools/KinFitter/interface/TKinFitter.h"

#include "ElectroWeakAnalysis/VPlusJets/interface/AngularVars.h"

#include "ElectroWeakAnalysis/VPlusJets/interface/METzCalculator.h"

#include "ClassifierOut/TMVAClassification_noqg_nJ2_mu_BDT.class.C"
#include "ClassifierOut/TMVAClassification_noqg_nJ3_mu_BDT.class.C"
#include "ClassifierOut/TMVAClassification_withqg_nJ2_mu_BDT.class.C"
#include "ClassifierOut/TMVAClassification_withqg_nJ3_mu_BDT.class.C"

#include "EffTableReader.h"
#include "EffTableLoader.h"

#include "ElectroWeakAnalysis/VPlusJets/interface/QGLikelihoodCalculator.h"
#include "MMozer/powhegweight/interface/pwhg_wrapper.h"

//const TString inQCDDir   = "/eos/uscms/store/user/lnujj/HCP2012METfix/MergedNtuples_v1/";
const TString inDataDir  = "/uscms_data/d3/zixu/BoostJet/ZPlusJets_jetsubstructure/Data/muon/ntuple/";
//const TString outDataDir = "/uscms_data/d3/zixu/BoostJet/ZPlusJets_jetsubstructure/Data/muon/RDTree/";
const TString outDataDir = "./";
const std::string fDir   = "EffTable2012/";
const std::string fInterferenceDir   = "InterferenceTable2012/";

bool large(const double &a, const double &b) { return a > b; }

struct sortPt {
	bool operator()(TLorentzVector* s1, TLorentzVector* s2) const { return s1->Pt() >= s2->Pt(); }
} mysortPt;

void kanamuon::myana(double myflag, bool isQCD, int runflag) {
	//Prepare the histogram for the cut-flow control : 8 presel + 7 sel
	const int n_step = 15;
	TH1F* h_events          = new TH1F("h_events", "h_events", n_step, 0, n_step);
	TH1F* h_events_weighted = new TH1F("h_events_weighted", "h_events_weighted", n_step, 0, n_step);

	string step_names[n_step] = {
		"all",
		"no scraping",
		"HBHE noise filter",
		"good PV",
		"tight lepton",
		"loose ele veto",
		"loose mu veto",
		"loose jets",
		"P_{T}(WJ1) > 30", 
		"P_{T}(WJ2) > 30",
		"M_{T}(W^{lep}) > 30",
		"tighter lepton",
		"#Delta#eta(W^{had}) < 1.5",
		"#Delta#phi(WJ1,MET) > 0.4",
		"P_{T}(W^{had}) > 40"};

	for ( int istep = 0; istep < n_step; istep++ ) {
		h_events -> GetXaxis() -> SetBinLabel( istep + 1, step_names[istep].c_str() );
		h_events_weighted -> GetXaxis() -> SetBinLabel( istep + 1, step_names[istep].c_str() );
	}

	TChain * myChain;

	if (myflag == 20121001 || myflag == -1){
		InitCounters( inDataDir + "DYJets_ntuple.root", h_events, h_events_weighted);
		myChain = new TChain("ZJet");
		myChain->Add( inDataDir + "DYJets_ntuple.root");
		Init(myChain);Loop( h_events, h_events_weighted, 20121001,runflag, outDataDir + "DYJets_RD");
	}

	if (myflag == 20121002 || myflag == -1){
		InitCounters( inDataDir + "ZZ_ntuple.root", h_events, h_events_weighted);
		myChain = new TChain("ZJet");
		myChain->Add( inDataDir + "ZZ_ntuple.root");
		Init(myChain);Loop( h_events, h_events_weighted, 20121002,runflag, outDataDir + "ZZ_RD");
	}

	if (myflag == 20121003 || myflag == -1){
		InitCounters( inDataDir + "WZ_ntuple.root", h_events, h_events_weighted);
		myChain = new TChain("ZJet");
		myChain->Add( inDataDir + "WZ_ntuple.root");
		Init(myChain);Loop( h_events, h_events_weighted, 20121003,runflag, outDataDir + "WZ_RD");
	}

	if (myflag == 20121004 || myflag == -1){
		InitCounters( inDataDir + "WW_ntuple.root", h_events, h_events_weighted);
		myChain = new TChain("ZJet");
		myChain->Add( inDataDir + "WW_ntuple.root");
		Init(myChain);Loop( h_events, h_events_weighted, 20121004,runflag, outDataDir + "WW_RD");
	}

	if (myflag == 20121005 || myflag == -1){
		InitCounters( inDataDir + "TtW_ntuple.root", h_events, h_events_weighted);
		myChain = new TChain("ZJet");
		myChain->Add( inDataDir + "TtW_ntuple.root");
		Init(myChain);Loop( h_events, h_events_weighted, 20121005,runflag, outDataDir + "TtW_RD");
	}

	if (myflag == 20121006 || myflag == -1){
		InitCounters( inDataDir + "TbartW_ntuple.root", h_events, h_events_weighted);
		myChain = new TChain("ZJet");
		myChain->Add( inDataDir + "TbartW_ntuple.root");
		Init(myChain);Loop( h_events, h_events_weighted, 20121006,runflag, outDataDir + "TbartW_RD");
	}
}



void kanamuon::Loop(TH1F* h_events, TH1F* h_events_weighted, int wda, int runflag, const char *outfilename, bool isQCD)
{
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntries();
	// Out Put File Here
	char rootfn[200]; 
	if (runflag ==0 ) {sprintf(rootfn, "%s.root",outfilename);}
	else {             sprintf(rootfn, "%s-VS-%i.root",outfilename,runflag);}
	TFile fresults= TFile(rootfn,"RECREATE");

	// Disable some variables never used to reduce the size of file
	fChain->SetBranchStatus("JetPFCor_etaetaMoment",    0);
	fChain->SetBranchStatus("JetPFCor_phiphiMoment",    0);
	fChain->SetBranchStatus("JetPFCor_etaphiMoment",    0);
	fChain->SetBranchStatus("JetPFCor_maxDistance",    0);
	fChain->SetBranchStatus("JetPFCor_nConstituents",    0);
	fChain->SetBranchStatus("JetPFCor_ChargedHadronEnergy",    0);
	fChain->SetBranchStatus("JetPFCor_ChargedHadronEnergyFrac",    0);
	fChain->SetBranchStatus("JetPFCor_NeutralHadronEnergy",    0);
	fChain->SetBranchStatus("JetPFCor_NeutralHadronEnergyFrac",    0);
	fChain->SetBranchStatus("JetPFCor_ChargedEmEnergy",    0);
	fChain->SetBranchStatus("JetPFCor_ChargedEmEnergyFrac",    0);
	fChain->SetBranchStatus("JetPFCor_ChargedMuEnergy",    0);
	fChain->SetBranchStatus("JetPFCor_ChargedMuEnergyFrac",    0);
	fChain->SetBranchStatus("JetPFCor_NeutralEmEnergy",    0);
	fChain->SetBranchStatus("JetPFCor_NeutralEmEnergyFrac",    0);
	fChain->SetBranchStatus("JetPFCor_MuonMultiplicity",    0);
	fChain->SetBranchStatus("JetPFCor_PhotonEnergy",    0);
	fChain->SetBranchStatus("JetPFCor_PhotonEnergyFraction",    0);
	fChain->SetBranchStatus("JetPFCor_ElectronEnergy",    0);
	fChain->SetBranchStatus("JetPFCor_ElectronEnergyFraction",    0);
	fChain->SetBranchStatus("JetPFCor_MuonEnergy",    0);
	fChain->SetBranchStatus("JetPFCor_MuonEnergyFraction",    0);
	fChain->SetBranchStatus("JetPFCor_HFHadronEnergy",    0);
	fChain->SetBranchStatus("JetPFCor_HFHadronEnergyFraction",    0);
	fChain->SetBranchStatus("JetPFCor_HFEMEnergy",    0);
	fChain->SetBranchStatus("JetPFCor_HFEMEnergyFraction",    0);
	fChain->SetBranchStatus("JetPFCor_ChargedHadronMultiplicity",    0);
	fChain->SetBranchStatus("JetPFCor_NeutralHadronMultiplicity",    0);
	fChain->SetBranchStatus("JetPFCor_PhotonMultiplicity",    0);
	fChain->SetBranchStatus("JetPFCor_ElectronMultiplicity",    0);
	fChain->SetBranchStatus("JetPFCor_HFHadronMultiplicity",    0);
	fChain->SetBranchStatus("JetPFCor_HFEMMultiplicity",    0);
	fChain->SetBranchStatus("JetPFCor_SumPtCands",    0);
	fChain->SetBranchStatus("JetPFCor_SumPt2Cands",    0);
	fChain->SetBranchStatus("JetPFCor_rmsCands",    0);
	// Drop gen jet information
	//Drop Some Groomed information
	//fChain->SetBranchStatus("GroomedJet_*_pt_uncorr" , 0);
	//fChain->SetBranchStatus("GroomedJet_*_tau1" , 0);
	//fChain->SetBranchStatus("GroomedJet_*_tau2" , 0);
	//fChain->SetBranchStatus("GroomedJet_*_tau3" , 0);
	//fChain->SetBranchStatus("GroomedJet_*_tau4" , 0);
	//   fChain->SetBranchStatus("GroomedJet_*_pt_*_uncorr" , 0);
	//   fChain->SetBranchStatus("GroomedJet_*_area", 0);
	//   fChain->SetBranchStatus("GroomedJet_*_area_tr", 0);
	//   fChain->SetBranchStatus("GroomedJet_*_area_ft", 0);
	//   fChain->SetBranchStatus("GroomedJet_*_area_pr", 0);
	//   fChain->SetBranchStatus("GroomedJet_*_jetcharge", 0); //This may be useful later
	fChain->SetBranchStatus("GroomedJet_*_constituents0_eta", 0);
	fChain->SetBranchStatus("GroomedJet_*_constituents0_phi", 0);
	fChain->SetBranchStatus("GroomedJet_*_constituents0_e", 0);
	fChain->SetBranchStatus("GroomedJet_*_nconstituents0", 0);
	fChain->SetBranchStatus("GroomedJet_*_constituents0pr_eta", 0);
	fChain->SetBranchStatus("GroomedJet_*_constituents0pr_phi", 0);
	fChain->SetBranchStatus("GroomedJet_*_constituents0pr_e", 0);
	fChain->SetBranchStatus("GroomedJet_*_nconstituents0pr", 0);
	fChain->SetBranchStatus("Mass*_PFCor*",0);
	fChain->SetBranchStatus("cosJacksonAngle*",0);
	fChain->SetBranchStatus("colorCorr*",0);
	fChain->SetBranchStatus("cosTheta*",0);
	fChain->SetBranchStatus("cosJackson*",0);
	fChain->SetBranchStatus("cosphiDecayPlane*",0);
	fChain->SetBranchStatus("colorCorrPull*",0);


	TTree *newtree = fChain->CloneTree();
	//char textfn[100]; 
	//sprintf(textfn,"%s.txt", rootfn);
	//FILE *textfile = fopen(textfn,"w");

	Int_t   ggdevt   =0,   evtNJ     =0;
	Int_t   ggdevtinclusive   =0; //For inclusive Jet Bin

	TBranch *branch_ggdevt= newtree->Branch("ggdevt",    &ggdevt,     "ggdevt/I");
	TBranch *branch_evtNJ = newtree->Branch("evtNJ",     &evtNJ,      "evtNJ/I");
	TBranch *branch_ggdevtinclusive = newtree->Branch("ggdevtinclusive", &ggdevtinclusive, "ggdevtinclusive/I");

	Float_t fit_muplus_px=0,   fit_muplus_py =0,   fit_muplus_pz=0,   fit_muplus_e=0;
	Float_t fit_muminus_px=0,  fit_muminus_py =0,  fit_muminus_pz=0,  fit_muminus_e=0;
	Float_t fit_aj_px=0,   fit_aj_py =0,   fit_aj_pz=0,   fit_aj_e=0;
	Float_t fit_bj_px=0,   fit_bj_py =0,   fit_bj_pz=0,   fit_bj_e=0;
	Float_t fit_mlljj=0,   fit_chi2  =999;  
	Int_t   fit_NDF  =999, fit_status=999;
	Float_t fit_mll  =0,   fit_mjj   =0;

	TBranch *branch_muplus_px = newtree->Branch("fit_muplus_px", &fit_muplus_px,  "fit_muplus_px/F");
	TBranch *branch_muplus_py = newtree->Branch("fit_muplus_py", &fit_muplus_py,  "fit_muplus_py/F");
	TBranch *branch_muplus_pz = newtree->Branch("fit_muplus_pz", &fit_muplus_pz,  "fit_muplus_pz/F");
	TBranch *branch_muplus_e  = newtree->Branch("fit_muplus_e",  &fit_muplus_e,   "fit_muplus_e/F");

	TBranch *branch_muminus_px = newtree->Branch("fit_muminus_px", &fit_muminus_px,  "fit_muminus_px/F");
	TBranch *branch_muminus_py = newtree->Branch("fit_muminus_py", &fit_muminus_py,  "fit_muminus_py/F");
	TBranch *branch_muminus_pz = newtree->Branch("fit_muminus_pz", &fit_muminus_pz,  "fit_muminus_pz/F");
	TBranch *branch_muminus_e  = newtree->Branch("fit_muminus_e",  &fit_muminus_e,   "fit_muminus_e/F");

	TBranch *branch_aj_px = newtree->Branch("fit_aj_px", &fit_aj_px,  "fit_aj_px/F");
	TBranch *branch_aj_py = newtree->Branch("fit_aj_py", &fit_aj_py,  "fit_aj_py/F");
	TBranch *branch_aj_pz = newtree->Branch("fit_aj_pz", &fit_aj_pz,  "fit_aj_pz/F");
	TBranch *branch_aj_e  = newtree->Branch("fit_aj_e",  &fit_aj_e,   "fit_aj_e/F");

	TBranch *branch_bj_px = newtree->Branch("fit_bj_px", &fit_bj_px,  "fit_bj_px/F");
	TBranch *branch_bj_py = newtree->Branch("fit_bj_py", &fit_bj_py,  "fit_bj_py/F");
	TBranch *branch_bj_pz = newtree->Branch("fit_bj_pz", &fit_bj_pz,  "fit_bj_pz/F");
	TBranch *branch_bj_e  = newtree->Branch("fit_bj_e",  &fit_bj_e,   "fit_bj_e/F");

	TBranch *branch_mlljj = newtree->Branch("fit_mlljj", &fit_mlljj,  "fit_mlljj/F");
	TBranch *branch_mll   = newtree->Branch("fit_mll",   &fit_mll,    "fit_mll/F");
	TBranch *branch_mjj   = newtree->Branch("fit_mjj",   &fit_mjj,    "fit_mjj/F");
	TBranch *branch_chi2  = newtree->Branch("fit_chi2",  &fit_chi2,   "fit_chi2/F");
	TBranch *branch_NDF   = newtree->Branch("fit_NDF",   &fit_NDF,    "fit_NDF/I");
	TBranch *branch_status= newtree->Branch("fit_status",&fit_status, "fit_status/I");

	Float_t TopWm=0,   TopWm5j=0;
	Float_t Tchi2=999, Tchi25j=999;
	TBranch *branch_TopWm   = newtree->Branch("TopWm",       &TopWm,      "TopWm/F");
	TBranch *branch_TopWm5j = newtree->Branch("TopWm5j",     &TopWm5j,    "TopWm5j/F");
	TBranch *branch_Tchi2   = newtree->Branch("Tchi2",       &Tchi2,      "Tchi2/F");
	TBranch *branch_Tchi25j = newtree->Branch("Tchi25j",     &Tchi25j,    "Tchi25j/F");

	Float_t ang_ha   = 999, ang_hb = 999, ang_hs = 999, ang_phi = 999, ang_phia = 999, ang_phib = 999;
	Float_t masslljj =-999, ptlljj =-999,  ylljj = -999,philljj = -999;
	TBranch * branch_ha   =  newtree->Branch("ang_ha",   &ang_ha,    "ang_ha/F");
	TBranch * branch_hb   =  newtree->Branch("ang_hb",   &ang_hb,    "ang_hb/F");
	TBranch * branch_hs   =  newtree->Branch("ang_hs",   &ang_hs,    "ang_hs/F");
	TBranch * branch_phi  =  newtree->Branch("ang_phi",  &ang_phi,   "ang_phi/F");
	TBranch * branch_phia =  newtree->Branch("ang_phia", &ang_phia,  "ang_phia/F");
	TBranch * branch_phib =  newtree->Branch("ang_phib", &ang_phib,  "ang_phib/F");
	TBranch * branch_orgm =  newtree->Branch("masslljj", &masslljj,  "masslljj/F");
	TBranch * branch_orgpt=  newtree->Branch("ptlljj",   &ptlljj,    "ptlljj/F");
	TBranch * branch_orgy =  newtree->Branch("ylljj",    &ylljj,     "ylljj/F");
	TBranch * branch_orgph=  newtree->Branch("philljj",  &philljj,   "philljj/F");



	Float_t effwt = 1.0, puwt = 1.0, puwt_up = 1.0, puwt_down = 1.0;
	TBranch * branch_effwt          =  newtree->Branch("effwt",       &effwt,        "effwt/F");
	TBranch * branch_puwt           =  newtree->Branch("puwt",        &puwt,         "puwt/F");
	TBranch * branch_puwt_up        =  newtree->Branch("puwt_up",     &puwt_up,      "puwt_up/F");
	TBranch * branch_puwt_down      =  newtree->Branch("puwt_down",   &puwt_down,    "puwt_down/F");


	Float_t qgld_Spring11[6]={-1,-1,-1,-1,-1,-1}; 
	Float_t qgld_Summer11[6]={-1,-1,-1,-1,-1,-1};
	Float_t qgld_Summer11CHS[6]={-1,-1,-1,-1,-1,-1};

	TBranch * branch_qgld_Spring11     =  newtree->Branch("qgld_Spring11",     qgld_Spring11,        "qgld_Spring11[6]/F");
	TBranch * branch_qgld_Summer11     =  newtree->Branch("qgld_Summer11",     qgld_Summer11,        "qgld_Summer11[6]/F");
	TBranch * branch_qgld_Summer11CHS  =  newtree->Branch("qgld_Summer11CHS",  qgld_Summer11CHS,     "qgld_Summer11CHS[6]/F");

	//Event Flag for Boosted W Analysis
	Int_t isgengdboostedZevt = 0;
	TBranch *branch_isgengdboostedZevt = newtree->Branch("isgengdboostedZevt", &isgengdboostedZevt, "isgengdboostedZevt/I");

	Int_t   ggdboostedZevt =0;
	TBranch *branch_ggdboostedZevt = newtree->Branch("ggdboostedZevt", &ggdboostedZevt, "ggdboostedZevt/I"); 

	Int_t   GroomedJet_numberbjets_csvl = 0;
	TBranch *branch_GroomedJet_numberbjets_csvl = newtree->Branch("GroomedJet_numberbjets_csvl", &GroomedJet_numberbjets_csvl,"GroomedJet_numberbjets_csvl/I");

	Int_t   GroomedJet_numberbjets_csvm = 0;
	TBranch *branch_GroomedJet_numberbjets_csvm = newtree->Branch("GroomedJet_numberbjets_csvm", &GroomedJet_numberbjets_csvm,"GroomedJet_numberbjets_csvm/I");

	Int_t   GroomedJet_numberbjets_csvl_veto = 0;
	TBranch *branch_GroomedJet_numberbjets_csvl_veto = newtree->Branch("GroomedJet_numberbjets_csvl_veto", &GroomedJet_numberbjets_csvl_veto,"GroomedJet_numberbjets_csvl_veto/I");

	Int_t   GroomedJet_numberbjets_csvm_veto = 0;
	TBranch *branch_GroomedJet_numberbjets_csvm_veto = newtree->Branch("GroomedJet_numberbjets_csvm_veto", &GroomedJet_numberbjets_csvm_veto,"GroomedJet_numberbjets_csvm_veto/I");

	Int_t   GroomedJet_numberjets = 0;
	TBranch *branch_GroomedJet_numberjets = newtree->Branch("GroomedJet_numberjets", &GroomedJet_numberjets,"GroomedJet_numberjets/I");

	//lepton, MET angular information
	Float_t GroomedJet_CA8_deltaR_lplusca8jet = -999, GroomedJet_CA8_deltaR_lminusca8jet = -999, GroomedJet_CA8_deltaphi_Vca8jet = -999;
	TBranch *branch_GroomedJet_CA8_deltaR_lplusca8jet = newtree->Branch("GroomedJet_CA8_deltaR_lplusca8jet", &GroomedJet_CA8_deltaR_lplusca8jet, "GroomedJet_CA8_deltaR_lplusca8jet/F");
	TBranch *branch_GroomedJet_CA8_deltaR_lminusca8jet = newtree->Branch("GroomedJet_CA8_deltaR_lminusca8jet", &GroomedJet_CA8_deltaR_lminusca8jet, "GroomedJet_CA8_deltaR_lminusca8jet/F");
	TBranch *branch_GroomedJet_CA8_deltaphi_Vca8jet = newtree->Branch("GroomedJet_CA8_deltaphi_Vca8jet", &GroomedJet_CA8_deltaphi_Vca8jet, "GroomedJet_CA8_deltaphi_Vca8jet/F");

	//Some More Variables To be Added in the Reduced Tree Or used in the TMVA Training
	Float_t GroomedJet_CA8_rcores01 = -1, GroomedJet_CA8_rcores02 = -1, GroomedJet_CA8_rcores03 = -1, GroomedJet_CA8_rcores04 = -1;
	Float_t GroomedJet_CA8_rcores05 = -1, GroomedJet_CA8_rcores06 = -1, GroomedJet_CA8_rcores07 = -1, GroomedJet_CA8_rcores08 = -1;
	Float_t GroomedJet_CA8_rcores09 = -1, GroomedJet_CA8_rcores10 = -1, GroomedJet_CA8_rcores11 = -1;

	TBranch *branch_GroomedJet_CA8_rcores01 = newtree->Branch("GroomedJet_CA8_rcores01", &GroomedJet_CA8_rcores01, "GroomedJet_CA8_rcores01/F");
	TBranch *branch_GroomedJet_CA8_rcores02 = newtree->Branch("GroomedJet_CA8_rcores02", &GroomedJet_CA8_rcores02, "GroomedJet_CA8_rcores02/F");
	TBranch *branch_GroomedJet_CA8_rcores03 = newtree->Branch("GroomedJet_CA8_rcores03", &GroomedJet_CA8_rcores03, "GroomedJet_CA8_rcores03/F");
	TBranch *branch_GroomedJet_CA8_rcores04 = newtree->Branch("GroomedJet_CA8_rcores04", &GroomedJet_CA8_rcores04, "GroomedJet_CA8_rcores04/F");
	TBranch *branch_GroomedJet_CA8_rcores05 = newtree->Branch("GroomedJet_CA8_rcores05", &GroomedJet_CA8_rcores05, "GroomedJet_CA8_rcores05/F");
	TBranch *branch_GroomedJet_CA8_rcores06 = newtree->Branch("GroomedJet_CA8_rcores06", &GroomedJet_CA8_rcores06, "GroomedJet_CA8_rcores06/F");
	TBranch *branch_GroomedJet_CA8_rcores07 = newtree->Branch("GroomedJet_CA8_rcores07", &GroomedJet_CA8_rcores07, "GroomedJet_CA8_rcores07/F");
	TBranch *branch_GroomedJet_CA8_rcores08 = newtree->Branch("GroomedJet_CA8_rcores08", &GroomedJet_CA8_rcores08, "GroomedJet_CA8_rcores08/F");
	TBranch *branch_GroomedJet_CA8_rcores09 = newtree->Branch("GroomedJet_CA8_rcores09", &GroomedJet_CA8_rcores09, "GroomedJet_CA8_rcores09/F");
	TBranch *branch_GroomedJet_CA8_rcores10 = newtree->Branch("GroomedJet_CA8_rcores10", &GroomedJet_CA8_rcores10, "GroomedJet_CA8_rcores10/F");
	TBranch *branch_GroomedJet_CA8_rcores11 = newtree->Branch("GroomedJet_CA8_rcores11", &GroomedJet_CA8_rcores11, "GroomedJet_CA8_rcores11/F");

	Float_t GroomedJet_CA8_ptcores01 = -1, GroomedJet_CA8_ptcores02 = -1, GroomedJet_CA8_ptcores03 = -1, GroomedJet_CA8_ptcores04 = -1;
	Float_t GroomedJet_CA8_ptcores05 = -1, GroomedJet_CA8_ptcores06 = -1, GroomedJet_CA8_ptcores07 = -1, GroomedJet_CA8_ptcores08 = -1;
	Float_t GroomedJet_CA8_ptcores09 = -1, GroomedJet_CA8_ptcores10 = -1, GroomedJet_CA8_ptcores11 = -1;

	TBranch *branch_GroomedJet_CA8_ptcores01 = newtree->Branch("GroomedJet_CA8_ptcores01", &GroomedJet_CA8_ptcores01, "GroomedJet_CA8_ptcores01/F");
	TBranch *branch_GroomedJet_CA8_ptcores02 = newtree->Branch("GroomedJet_CA8_ptcores02", &GroomedJet_CA8_ptcores02, "GroomedJet_CA8_ptcores02/F");
	TBranch *branch_GroomedJet_CA8_ptcores03 = newtree->Branch("GroomedJet_CA8_ptcores03", &GroomedJet_CA8_ptcores03, "GroomedJet_CA8_ptcores03/F");
	TBranch *branch_GroomedJet_CA8_ptcores04 = newtree->Branch("GroomedJet_CA8_ptcores04", &GroomedJet_CA8_ptcores04, "GroomedJet_CA8_ptcores04/F");
	TBranch *branch_GroomedJet_CA8_ptcores05 = newtree->Branch("GroomedJet_CA8_ptcores05", &GroomedJet_CA8_ptcores05, "GroomedJet_CA8_ptcores05/F");
	TBranch *branch_GroomedJet_CA8_ptcores06 = newtree->Branch("GroomedJet_CA8_ptcores06", &GroomedJet_CA8_ptcores06, "GroomedJet_CA8_ptcores06/F");
	TBranch *branch_GroomedJet_CA8_ptcores07 = newtree->Branch("GroomedJet_CA8_ptcores07", &GroomedJet_CA8_ptcores07, "GroomedJet_CA8_ptcores07/F");
	TBranch *branch_GroomedJet_CA8_ptcores08 = newtree->Branch("GroomedJet_CA8_ptcores08", &GroomedJet_CA8_ptcores08, "GroomedJet_CA8_ptcores08/F");
	TBranch *branch_GroomedJet_CA8_ptcores09 = newtree->Branch("GroomedJet_CA8_ptcores09", &GroomedJet_CA8_ptcores09, "GroomedJet_CA8_ptcores09/F");
	TBranch *branch_GroomedJet_CA8_ptcores10 = newtree->Branch("GroomedJet_CA8_ptcores10", &GroomedJet_CA8_ptcores10, "GroomedJet_CA8_ptcores10/F");
	TBranch *branch_GroomedJet_CA8_ptcores11 = newtree->Branch("GroomedJet_CA8_ptcores11", &GroomedJet_CA8_ptcores11, "GroomedJet_CA8_ptcores11/F");

	Float_t GroomedJet_CA8_planarflow01 = -1, GroomedJet_CA8_planarflow02 = -1, GroomedJet_CA8_planarflow03 = -1, GroomedJet_CA8_planarflow04 = -1;
	Float_t GroomedJet_CA8_planarflow05 = -1, GroomedJet_CA8_planarflow06 = -1, GroomedJet_CA8_planarflow07 = -1, GroomedJet_CA8_planarflow08 = -1;
	Float_t GroomedJet_CA8_planarflow09 = -1, GroomedJet_CA8_planarflow10 = -1, GroomedJet_CA8_planarflow11 = -1;

	TBranch *branch_GroomedJet_CA8_planarflow01 = newtree->Branch("GroomedJet_CA8_planarflow01", &GroomedJet_CA8_planarflow01, "GroomedJet_CA8_planarflow01/F");
	TBranch *branch_GroomedJet_CA8_planarflow02 = newtree->Branch("GroomedJet_CA8_planarflow02", &GroomedJet_CA8_planarflow02, "GroomedJet_CA8_planarflow02/F");
	TBranch *branch_GroomedJet_CA8_planarflow03 = newtree->Branch("GroomedJet_CA8_planarflow03", &GroomedJet_CA8_planarflow03, "GroomedJet_CA8_planarflow03/F");
	TBranch *branch_GroomedJet_CA8_planarflow04 = newtree->Branch("GroomedJet_CA8_planarflow04", &GroomedJet_CA8_planarflow04, "GroomedJet_CA8_planarflow04/F");
	TBranch *branch_GroomedJet_CA8_planarflow05 = newtree->Branch("GroomedJet_CA8_planarflow05", &GroomedJet_CA8_planarflow05, "GroomedJet_CA8_planarflow05/F");
	TBranch *branch_GroomedJet_CA8_planarflow06 = newtree->Branch("GroomedJet_CA8_planarflow06", &GroomedJet_CA8_planarflow06, "GroomedJet_CA8_planarflow06/F");
	TBranch *branch_GroomedJet_CA8_planarflow07 = newtree->Branch("GroomedJet_CA8_planarflow07", &GroomedJet_CA8_planarflow07, "GroomedJet_CA8_planarflow07/F");
	TBranch *branch_GroomedJet_CA8_planarflow08 = newtree->Branch("GroomedJet_CA8_planarflow08", &GroomedJet_CA8_planarflow08, "GroomedJet_CA8_planarflow08/F");
	TBranch *branch_GroomedJet_CA8_planarflow09 = newtree->Branch("GroomedJet_CA8_planarflow09", &GroomedJet_CA8_planarflow09, "GroomedJet_CA8_planarflow09/F");
	TBranch *branch_GroomedJet_CA8_planarflow10 = newtree->Branch("GroomedJet_CA8_planarflow10", &GroomedJet_CA8_planarflow10, "GroomedJet_CA8_planarflow10/F");
	TBranch *branch_GroomedJet_CA8_planarflow11 = newtree->Branch("GroomedJet_CA8_planarflow11", &GroomedJet_CA8_planarflow11, "GroomedJet_CA8_planarflow11/F");

	Float_t GroomedJet_CA8_mass_sensi_tr = -1, GroomedJet_CA8_mass_sensi_ft = -1, GroomedJet_CA8_mass_sensi_pr = -1;
	TBranch *branch_GroomedJet_CA8_mass_sensi_tr = newtree->Branch("GroomedJet_CA8_mass_sensi_tr", &GroomedJet_CA8_mass_sensi_tr, "GroomedJet_CA8_mass_sensi_tr/F");
	TBranch *branch_GroomedJet_CA8_mass_sensi_ft = newtree->Branch("GroomedJet_CA8_mass_sensi_ft", &GroomedJet_CA8_mass_sensi_ft, "GroomedJet_CA8_mass_sensi_ft/F");
	TBranch *branch_GroomedJet_CA8_mass_sensi_pr = newtree->Branch("GroomedJet_CA8_mass_sensi_pr", &GroomedJet_CA8_mass_sensi_pr, "GroomedJet_CA8_mass_sensi_pr/F");

	Float_t GroomedJet_CA8_qjetmassvolatility = -1;
	TBranch *branch_GroomedJet_CA8_qjetmassvolatility = newtree->Branch("GroomedJet_CA8_qjetmassvolatility", &GroomedJet_CA8_qjetmassvolatility, "GroomedJet_CA8_qjetmassvolatility/F");

	Float_t GroomedJet_CA8_prsubjet1ptoverjetpt = -1, GroomedJet_CA8_prsubjet2ptoverjetpt = -1;
	Float_t GroomedJet_CA8_prsubjet1subjet2_deltaR = -1;
	Float_t GroomedJet_CA8_JetResp=-1;
	Float_t GroomedJet_CA8_JetResp1=-1;
	Float_t GroomedJet_CA8_JetResp2=-1;

	TBranch *branch_GroomedJet_CA8_prsubjet1ptoverjetpt = newtree->Branch("GroomedJet_CA8_prsubjet1ptoverjetpt", &GroomedJet_CA8_prsubjet1ptoverjetpt, "GroomedJet_CA8_prsubjet1ptoverjetpt/F");
	TBranch *branch_GroomedJet_CA8_prsubjet2ptoverjetpt = newtree->Branch("GroomedJet_CA8_prsubjet2ptoverjetpt", &GroomedJet_CA8_prsubjet2ptoverjetpt, "GroomedJet_CA8_prsubjet2ptoverjetpt/F");
	TBranch *branch_GroomedJet_CA8_prsubjet1subjet2_deltaR = newtree->Branch("GroomedJet_CA8_prsubjet1subjet2_deltaR", &GroomedJet_CA8_prsubjet1subjet2_deltaR, "GroomedJet_CA8_prsubjet1subjet2_deltaR/F");
	TBranch *branch_GroomedJet_CA8_JetResp = newtree->Branch("GroomedJet_CA8_JetResp",&GroomedJet_CA8_JetResp,"GroomedJet_CA8_JetResp/F");
	TBranch *branch_GroomedJet_CA8_JetResp1 = newtree->Branch("GroomedJet_CA8_JetResp1",&GroomedJet_CA8_JetResp1,"GroomedJet_CA8_JetResp1/F");
	TBranch *branch_GroomedJet_CA8_JetResp2 = newtree->Branch("GroomedJet_CA8_JetResp2",&GroomedJet_CA8_JetResp2,"GroomedJet_CA8_JetResp2/F");


	Float_t boostedZ_llj_e=-999,   boostedZ_llj_pt=-999,   boostedZ_llj_eta=-999,   boostedZ_llj_phi=-999,   boostedZ_llj_m=-999,   boostedZ_llj_y=-999;
	TBranch *branch_boostedZ_llj_e    = newtree->Branch("boostedZ_llj_e",    &boostedZ_llj_e,     "boostedZ_llj_e/F");
	TBranch *branch_boostedZ_llj_pt   = newtree->Branch("boostedZ_llj_pt",   &boostedZ_llj_pt,    "boostedZ_llj_pt/F");
	TBranch *branch_boostedZ_llj_eta  = newtree->Branch("boostedZ_llj_eta",  &boostedZ_llj_eta,   "boostedZ_llj_eta/F");
	TBranch *branch_boostedZ_llj_phi  = newtree->Branch("boostedZ_llj_phi",  &boostedZ_llj_phi,   "boostedZ_llj_phi/F");
	TBranch *branch_boostedZ_llj_m    = newtree->Branch("boostedZ_llj_m",    &boostedZ_llj_m,     "boostedZ_llj_m/F");
	TBranch *branch_boostedZ_llj_y    = newtree->Branch("boostedZ_llj_y",    &boostedZ_llj_y,     "boostedZ_llj_y/F");

	Float_t boostedZ_zjj_ang_ha   = 999, boostedZ_zjj_ang_hb = 999, boostedZ_zjj_ang_hs = 999, boostedZ_zjj_ang_phi = 999, boostedZ_zjj_ang_phia = 999, boostedZ_zjj_ang_phib = 999;

	TBranch * branch_boostedZ_zjj_ang_ha   = newtree->Branch("boostedZ_zjj_ang_ha",   &boostedZ_zjj_ang_ha,    "boostedZ_zjj_ang_ha/F");
	TBranch * branch_boostedZ_zjj_ang_hb   = newtree->Branch("boostedZ_zjj_ang_hb",   &boostedZ_zjj_ang_hb,    "boostedZ_zjj_ang_hb/F");
	TBranch * branch_boostedZ_zjj_ang_hs   = newtree->Branch("boostedZ_zjj_ang_hs",   &boostedZ_zjj_ang_hs,    "boostedZ_zjj_ang_hs/F");
	TBranch * branch_boostedZ_zjj_ang_phi  = newtree->Branch("boostedZ_zjj_ang_phi",  &boostedZ_zjj_ang_phi,   "boostedZ_zjj_ang_phi/F");
	TBranch * branch_boostedZ_zjj_ang_phia = newtree->Branch("boostedZ_zjj_ang_phia", &boostedZ_zjj_ang_phia,  "boostedZ_zjj_ang_phia/F");
	TBranch * branch_boostedZ_zjj_ang_phib = newtree->Branch("boostedZ_zjj_ang_phib", &boostedZ_zjj_ang_phib,  "boostedZ_zjj_ang_phib/F");

	/*//Event Flag for Boosted W Analysis With AK7
	Int_t   ggdboostedZevt_ak7 =0;
	TBranch *branch_ggdboostedZevt_ak7 = newtree->Branch("ggdboostedZevt_ak7", &ggdboostedZevt_ak7, "ggdboostedZevt_ak7/I"); 

	Int_t   GroomedJet_numberbjets_csvl_ak7 = 0;
	TBranch *branch_GroomedJet_numberbjets_csvl_ak7 = newtree->Branch("GroomedJet_numberbjets_csvl_ak7", &GroomedJet_numberbjets_csvl_ak7,"GroomedJet_numberbjets_csvl_ak7/I");

	Int_t   GroomedJet_numberbjets_csvm_ak7 = 0;
	TBranch *branch_GroomedJet_numberbjets_csvm_ak7 = newtree->Branch("GroomedJet_numberbjets_csvm_ak7", &GroomedJet_numberbjets_csvm_ak7,"GroomedJet_numberbjets_csvm_ak7/I");

	Int_t   GroomedJet_numberbjets_csvl_veto_ak7 = 0;
	TBranch *branch_GroomedJet_numberbjets_csvl_veto_ak7 = newtree->Branch("GroomedJet_numberbjets_csvl_veto_ak7", &GroomedJet_numberbjets_csvl_veto_ak7,"GroomedJet_numberbjets_csvl_veto_ak7/I");

	Int_t   GroomedJet_numberbjets_csvm_veto_ak7 = 0;
	TBranch *branch_GroomedJet_numberbjets_csvm_veto_ak7 = newtree->Branch("GroomedJet_numberbjets_csvm_veto_ak7", &GroomedJet_numberbjets_csvm_veto_ak7,"GroomedJet_numberbjets_csvm_veto_ak7/I");

	Int_t   GroomedJet_numberjets_ak7 = 0;
	TBranch *branch_GroomedJet_numberjets_ak7 = newtree->Branch("GroomedJet_numberjets_ak7", &GroomedJet_numberjets_ak7,"GroomedJet_numberjets_ak7/I");

	//lepton, MET angular information
	Float_t GroomedJet_AK7_deltaR_lplusak7jet = -999, GroomedJet_AK7_deltaR_lminusak7jet = -999, GroomedJet_AK7_deltaphi_Vak7jet = -999;
	TBranch *branch_GroomedJet_AK7_deltaR_lplusak7jet = newtree->Branch("GroomedJet_AK7_deltaR_lplusak7jet", &GroomedJet_AK7_deltaR_lplusak7jet, "GroomedJet_AK7_deltaR_lplusak7jet/F");
	TBranch *branch_GroomedJet_AK7_deltaR_lminusak7jet = newtree->Branch("GroomedJet_AK7_deltaR_lminusak7jet", &GroomedJet_AK7_deltaR_lminusak7jet, "GroomedJet_AK7_deltaR_lminusak7jet/F");
	TBranch *branch_GroomedJet_AK7_deltaphi_Vak7jet = newtree->Branch("GroomedJet_AK7_deltaphi_Vak7jet", &GroomedJet_AK7_deltaphi_Vak7jet, "GroomedJet_AK7_deltaphi_Vak7jet/F");

	//Some More Variables To be Added in the Reduced Tree Or used in the TMVA Training
	Float_t GroomedJet_AK7_rcores01 = -1, GroomedJet_AK7_rcores02 = -1, GroomedJet_AK7_rcores03 = -1, GroomedJet_AK7_rcores04 = -1;
	Float_t GroomedJet_AK7_rcores05 = -1, GroomedJet_AK7_rcores06 = -1, GroomedJet_AK7_rcores07 = -1, GroomedJet_AK7_rcores08 = -1;
	Float_t GroomedJet_AK7_rcores09 = -1, GroomedJet_AK7_rcores10 = -1, GroomedJet_AK7_rcores11 = -1;

	TBranch *branch_GroomedJet_AK7_rcores01 = newtree->Branch("GroomedJet_AK7_rcores01", &GroomedJet_AK7_rcores01, "GroomedJet_AK7_rcores01/F");
	TBranch *branch_GroomedJet_AK7_rcores02 = newtree->Branch("GroomedJet_AK7_rcores02", &GroomedJet_AK7_rcores02, "GroomedJet_AK7_rcores02/F");
	TBranch *branch_GroomedJet_AK7_rcores03 = newtree->Branch("GroomedJet_AK7_rcores03", &GroomedJet_AK7_rcores03, "GroomedJet_AK7_rcores03/F");
	TBranch *branch_GroomedJet_AK7_rcores04 = newtree->Branch("GroomedJet_AK7_rcores04", &GroomedJet_AK7_rcores04, "GroomedJet_AK7_rcores04/F");
	TBranch *branch_GroomedJet_AK7_rcores05 = newtree->Branch("GroomedJet_AK7_rcores05", &GroomedJet_AK7_rcores05, "GroomedJet_AK7_rcores05/F");
	TBranch *branch_GroomedJet_AK7_rcores06 = newtree->Branch("GroomedJet_AK7_rcores06", &GroomedJet_AK7_rcores06, "GroomedJet_AK7_rcores06/F");
	TBranch *branch_GroomedJet_AK7_rcores07 = newtree->Branch("GroomedJet_AK7_rcores07", &GroomedJet_AK7_rcores07, "GroomedJet_AK7_rcores07/F");
	TBranch *branch_GroomedJet_AK7_rcores08 = newtree->Branch("GroomedJet_AK7_rcores08", &GroomedJet_AK7_rcores08, "GroomedJet_AK7_rcores08/F");
	TBranch *branch_GroomedJet_AK7_rcores09 = newtree->Branch("GroomedJet_AK7_rcores09", &GroomedJet_AK7_rcores09, "GroomedJet_AK7_rcores09/F");
	TBranch *branch_GroomedJet_AK7_rcores10 = newtree->Branch("GroomedJet_AK7_rcores10", &GroomedJet_AK7_rcores10, "GroomedJet_AK7_rcores10/F");
	TBranch *branch_GroomedJet_AK7_rcores11 = newtree->Branch("GroomedJet_AK7_rcores11", &GroomedJet_AK7_rcores11, "GroomedJet_AK7_rcores11/F");

	Float_t GroomedJet_AK7_ptcores01 = -1, GroomedJet_AK7_ptcores02 = -1, GroomedJet_AK7_ptcores03 = -1, GroomedJet_AK7_ptcores04 = -1;
	Float_t GroomedJet_AK7_ptcores05 = -1, GroomedJet_AK7_ptcores06 = -1, GroomedJet_AK7_ptcores07 = -1, GroomedJet_AK7_ptcores08 = -1;
	Float_t GroomedJet_AK7_ptcores09 = -1, GroomedJet_AK7_ptcores10 = -1, GroomedJet_AK7_ptcores11 = -1;

	TBranch *branch_GroomedJet_AK7_ptcores01 = newtree->Branch("GroomedJet_AK7_ptcores01", &GroomedJet_AK7_ptcores01, "GroomedJet_AK7_ptcores01/F");
	TBranch *branch_GroomedJet_AK7_ptcores02 = newtree->Branch("GroomedJet_AK7_ptcores02", &GroomedJet_AK7_ptcores02, "GroomedJet_AK7_ptcores02/F");
	TBranch *branch_GroomedJet_AK7_ptcores03 = newtree->Branch("GroomedJet_AK7_ptcores03", &GroomedJet_AK7_ptcores03, "GroomedJet_AK7_ptcores03/F");
	TBranch *branch_GroomedJet_AK7_ptcores04 = newtree->Branch("GroomedJet_AK7_ptcores04", &GroomedJet_AK7_ptcores04, "GroomedJet_AK7_ptcores04/F");
	TBranch *branch_GroomedJet_AK7_ptcores05 = newtree->Branch("GroomedJet_AK7_ptcores05", &GroomedJet_AK7_ptcores05, "GroomedJet_AK7_ptcores05/F");
	TBranch *branch_GroomedJet_AK7_ptcores06 = newtree->Branch("GroomedJet_AK7_ptcores06", &GroomedJet_AK7_ptcores06, "GroomedJet_AK7_ptcores06/F");
	TBranch *branch_GroomedJet_AK7_ptcores07 = newtree->Branch("GroomedJet_AK7_ptcores07", &GroomedJet_AK7_ptcores07, "GroomedJet_AK7_ptcores07/F");
	TBranch *branch_GroomedJet_AK7_ptcores08 = newtree->Branch("GroomedJet_AK7_ptcores08", &GroomedJet_AK7_ptcores08, "GroomedJet_AK7_ptcores08/F");
	TBranch *branch_GroomedJet_AK7_ptcores09 = newtree->Branch("GroomedJet_AK7_ptcores09", &GroomedJet_AK7_ptcores09, "GroomedJet_AK7_ptcores09/F");
	TBranch *branch_GroomedJet_AK7_ptcores10 = newtree->Branch("GroomedJet_AK7_ptcores10", &GroomedJet_AK7_ptcores10, "GroomedJet_AK7_ptcores10/F");
	TBranch *branch_GroomedJet_AK7_ptcores11 = newtree->Branch("GroomedJet_AK7_ptcores11", &GroomedJet_AK7_ptcores11, "GroomedJet_AK7_ptcores11/F");

	Float_t GroomedJet_AK7_planarflow01 = -1, GroomedJet_AK7_planarflow02 = -1, GroomedJet_AK7_planarflow03 = -1, GroomedJet_AK7_planarflow04 = -1;
	Float_t GroomedJet_AK7_planarflow05 = -1, GroomedJet_AK7_planarflow06 = -1, GroomedJet_AK7_planarflow07 = -1, GroomedJet_AK7_planarflow08 = -1;
	Float_t GroomedJet_AK7_planarflow09 = -1, GroomedJet_AK7_planarflow10 = -1, GroomedJet_AK7_planarflow11 = -1;

	TBranch *branch_GroomedJet_AK7_planarflow01 = newtree->Branch("GroomedJet_AK7_planarflow01", &GroomedJet_AK7_planarflow01, "GroomedJet_AK7_planarflow01/F");
	TBranch *branch_GroomedJet_AK7_planarflow02 = newtree->Branch("GroomedJet_AK7_planarflow02", &GroomedJet_AK7_planarflow02, "GroomedJet_AK7_planarflow02/F");
	TBranch *branch_GroomedJet_AK7_planarflow03 = newtree->Branch("GroomedJet_AK7_planarflow03", &GroomedJet_AK7_planarflow03, "GroomedJet_AK7_planarflow03/F");
	TBranch *branch_GroomedJet_AK7_planarflow04 = newtree->Branch("GroomedJet_AK7_planarflow04", &GroomedJet_AK7_planarflow04, "GroomedJet_AK7_planarflow04/F");
	TBranch *branch_GroomedJet_AK7_planarflow05 = newtree->Branch("GroomedJet_AK7_planarflow05", &GroomedJet_AK7_planarflow05, "GroomedJet_AK7_planarflow05/F");
	TBranch *branch_GroomedJet_AK7_planarflow06 = newtree->Branch("GroomedJet_AK7_planarflow06", &GroomedJet_AK7_planarflow06, "GroomedJet_AK7_planarflow06/F");
	TBranch *branch_GroomedJet_AK7_planarflow07 = newtree->Branch("GroomedJet_AK7_planarflow07", &GroomedJet_AK7_planarflow07, "GroomedJet_AK7_planarflow07/F");
	TBranch *branch_GroomedJet_AK7_planarflow08 = newtree->Branch("GroomedJet_AK7_planarflow08", &GroomedJet_AK7_planarflow08, "GroomedJet_AK7_planarflow08/F");
	TBranch *branch_GroomedJet_AK7_planarflow09 = newtree->Branch("GroomedJet_AK7_planarflow09", &GroomedJet_AK7_planarflow09, "GroomedJet_AK7_planarflow09/F");
	TBranch *branch_GroomedJet_AK7_planarflow10 = newtree->Branch("GroomedJet_AK7_planarflow10", &GroomedJet_AK7_planarflow10, "GroomedJet_AK7_planarflow10/F");
	TBranch *branch_GroomedJet_AK7_planarflow11 = newtree->Branch("GroomedJet_AK7_planarflow11", &GroomedJet_AK7_planarflow11, "GroomedJet_AK7_planarflow11/F");

	Float_t GroomedJet_AK7_mass_sensi_tr = -1, GroomedJet_AK7_mass_sensi_ft = -1, GroomedJet_AK7_mass_sensi_pr = -1;
	TBranch *branch_GroomedJet_AK7_mass_sensi_tr = newtree->Branch("GroomedJet_AK7_mass_sensi_tr", &GroomedJet_AK7_mass_sensi_tr, "GroomedJet_AK7_mass_sensi_tr/F");
	TBranch *branch_GroomedJet_AK7_mass_sensi_ft = newtree->Branch("GroomedJet_AK7_mass_sensi_ft", &GroomedJet_AK7_mass_sensi_ft, "GroomedJet_AK7_mass_sensi_ft/F");
	TBranch *branch_GroomedJet_AK7_mass_sensi_pr = newtree->Branch("GroomedJet_AK7_mass_sensi_pr", &GroomedJet_AK7_mass_sensi_pr, "GroomedJet_AK7_mass_sensi_pr/F");

	Float_t GroomedJet_AK7_qjetmassvolatility = -1;
	TBranch *branch_GroomedJet_AK7_qjetmassvolatility = newtree->Branch("GroomedJet_AK7_qjetmassvolatility", &GroomedJet_AK7_qjetmassvolatility, "GroomedJet_AK7_qjetmassvolatility/F");

	Float_t GroomedJet_AK7_prsubjet1ptoverjetpt = -1, GroomedJet_AK7_prsubjet2ptoverjetpt = -1;
	Float_t GroomedJet_AK7_prsubjet1subjet2_deltaR = -1;

	TBranch *branch_GroomedJet_AK7_prsubjet1ptoverjetpt = newtree->Branch("GroomedJet_AK7_prsubjet1ptoverjetpt", &GroomedJet_AK7_prsubjet1ptoverjetpt, "GroomedJet_AK7_prsubjet1ptoverjetpt/F");
	TBranch *branch_GroomedJet_AK7_prsubjet2ptoverjetpt = newtree->Branch("GroomedJet_AK7_prsubjet2ptoverjetpt", &GroomedJet_AK7_prsubjet2ptoverjetpt, "GroomedJet_AK7_prsubjet2ptoverjetpt/F");
	TBranch *branch_GroomedJet_AK7_prsubjet1subjet2_deltaR = newtree->Branch("GroomedJet_AK7_prsubjet1subjet2_deltaR", &GroomedJet_AK7_prsubjet1subjet2_deltaR, "GroomedJet_AK7_prsubjet1subjet2_deltaR/F");

	Float_t boostedZ_llj_e_ak7 =-999,   boostedZ_llj_pt_ak7=-999,   boostedZ_llj_eta_ak7=-999,   boostedZ_llj_phi_ak7=-999,   boostedZ_llj_m_ak7=-999,   boostedZ_llj_y_ak7=-999;
	TBranch *branch_boostedZ_llj_e_ak7    = newtree->Branch("boostedZ_llj_e_ak7",    &boostedZ_llj_e_ak7,     "boostedZ_llj_e_ak7/F");
	TBranch *branch_boostedZ_llj_pt_ak7   = newtree->Branch("boostedZ_llj_pt_ak7",   &boostedZ_llj_pt_ak7,    "boostedZ_llj_pt_ak7/F");
	TBranch *branch_boostedZ_llj_eta_ak7  = newtree->Branch("boostedZ_llj_eta_ak7",  &boostedZ_llj_eta_ak7,   "boostedZ_llj_eta_ak7/F");
	TBranch *branch_boostedZ_llj_phi_ak7  = newtree->Branch("boostedZ_llj_phi_ak7",  &boostedZ_llj_phi_ak7,   "boostedZ_llj_phi_ak7/F");
	TBranch *branch_boostedZ_llj_m_ak7    = newtree->Branch("boostedZ_llj_m_ak7",    &boostedZ_llj_m_ak7,     "boostedZ_llj_m_ak7/F");
	TBranch *branch_boostedZ_llj_y_ak7    = newtree->Branch("boostedZ_llj_y_ak7",    &boostedZ_llj_y_ak7,     "boostedZ_llj_y_ak7/F");

	Float_t boostedZ_zjj_ang_ha_ak7   = 999, boostedZ_zjj_ang_hb_ak7 = 999, boostedZ_zjj_ang_hs_ak7 = 999, boostedZ_zjj_ang_phi_ak7 = 999, boostedZ_zjj_ang_phia_ak7 = 999, boostedZ_zjj_ang_phib_ak7 = 999;

	TBranch * branch_boostedZ_zjj_ang_ha_ak7   = newtree->Branch("boostedZ_zjj_ang_ha_ak7",   &boostedZ_zjj_ang_ha_ak7,    "boostedZ_zjj_ang_ha_ak7/F");
	TBranch * branch_boostedZ_zjj_ang_hb_ak7   = newtree->Branch("boostedZ_zjj_ang_hb_ak7",   &boostedZ_zjj_ang_hb_ak7,    "boostedZ_zjj_ang_hb_ak7/F");
	TBranch * branch_boostedZ_zjj_ang_hs_ak7   = newtree->Branch("boostedZ_zjj_ang_hs_ak7",   &boostedZ_zjj_ang_hs_ak7,    "boostedZ_zjj_ang_hs_ak7/F");
	TBranch * branch_boostedZ_zjj_ang_phi_ak7  = newtree->Branch("boostedZ_zjj_ang_phi_ak7",  &boostedZ_zjj_ang_phi_ak7,   "boostedZ_zjj_ang_phi_ak7/F");
	TBranch * branch_boostedZ_zjj_ang_phia_ak7 = newtree->Branch("boostedZ_zjj_ang_phia_ak7", &boostedZ_zjj_ang_phia_ak7,  "boostedZ_zjj_ang_phia_ak7/F");
	TBranch * branch_boostedZ_zjj_ang_phib_ak7 = newtree->Branch("boostedZ_zjj_ang_phib_ak7", &boostedZ_zjj_ang_phib_ak7,  "boostedZ_zjj_ang_phib_ak7/F");
	*/
	//End Some More Variables To be Added in the Reduced Tree Or used in the TMVA Training
	


	EffTableLoader muIDEff(            fDir + "scaleFactor-Run2012ABCD-RecoToIso.txt");
	EffTableLoader muHLTEff(           fDir + "efficiency-Run2012ABCD-IsoToIsoMuHLT.txt");

	// Pile up Re-weighting
	/*
	   edm::Lumi3DReWeighting LumiWeights_ = edm::Lumi3DReWeighting("PUMC_dist.root", "PUData_dist.root", "pileup", "pileup", "Weight_3D.root");
	   LumiWeights_.weight3D_init( 1.08 );

	   edm::Lumi3DReWeighting up_LumiWeights_ = edm::Lumi3DReWeighting("PUMC_dist.root", "PUData_dist.root", "pileup", "pileup", "Weight_3D_up.root");
	   up_LumiWeights_.weight3D_init( 1.16 );

	   edm::Lumi3DReWeighting dn_LumiWeights_ = edm::Lumi3DReWeighting("PUMC_dist.root", "PUData_dist.root", "pileup", "pileup", "Weight_3D_down.root");
	   dn_LumiWeights_.weight3D_init( 1.00 );
	   */  
	// S7 MC PU True profile - hardcoded, wow
	// https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupMCReweightingUtilities
	// TFile *dataFile_      = new TFile( "PileupHistogramGold_190456-196531_8TeV_PromptReco_Collisions12_true.root" );
	//TFile *dataFile_      = new TFile( "Data190389-200041_PileupHistogram.root" );
	TFile *dataFile_      = new TFile( "Data190456-208686_PileupHistogram.root" );
	TH1F* PU_intended = new TH1F(  *(static_cast<TH1F*>(dataFile_->Get( "pileup" )->Clone() )) );
	TH1F* PU_generated = new TH1F("PU_generated","Generated pileup distribution (i.e., MC)",60,0.,60);
	Double_t Summer2012[60] = { //Pile Up For S10
		2.560E-06,
		5.239E-06,
		1.420E-05,
		5.005E-05,
		1.001E-04,
		2.705E-04,
		1.999E-03,
		6.097E-03,
		1.046E-02,
		1.383E-02,
		1.685E-02,
		2.055E-02,
		2.572E-02,
		3.262E-02,
		4.121E-02,
		4.977E-02,
		5.539E-02,
		5.725E-02,
		5.607E-02,
		5.312E-02,
		5.008E-02,
		4.763E-02,
		4.558E-02,
		4.363E-02,
		4.159E-02,
		3.933E-02,
		3.681E-02,
		3.406E-02,
		3.116E-02,
		2.818E-02,
		2.519E-02,
		2.226E-02,
		1.946E-02,
		1.682E-02,
		1.437E-02,
		1.215E-02,
		1.016E-02,
		8.400E-03,
		6.873E-03,
		5.564E-03,
		4.457E-03,
		3.533E-03,
		2.772E-03,
		2.154E-03,
		1.656E-03,
		1.261E-03,
		9.513E-04,
		7.107E-04,
		5.259E-04,
		3.856E-04,
		2.801E-04,
		2.017E-04,
		1.439E-04,
		1.017E-04,
		7.126E-05,
		4.948E-05,
		3.405E-05,
		2.322E-05,
		1.570E-05,
		5.005E-06
			//Pile Up for S7
			/*2.344E-05,
			  2.344E-05,
			  2.344E-05,
			  2.344E-05,
			  4.687E-04,
			  4.687E-04,
			  7.032E-04,
			  9.414E-04,
			  1.234E-03,
			  1.603E-03,
			  2.464E-03,
			  3.250E-03,
			  5.021E-03,
			  6.644E-03,
			  8.502E-03,
			  1.121E-02,
			  1.518E-02,
			  2.033E-02,
			  2.608E-02,
			  3.171E-02,
			  3.667E-02,
			  4.060E-02,
			  4.338E-02,
			  4.520E-02,
			  4.641E-02,
			  4.735E-02,
			  4.816E-02,
			  4.881E-02,
			  4.917E-02,
			  4.909E-02,
			  4.842E-02,
			  4.707E-02,
			  4.501E-02,
			  4.228E-02,
			  3.896E-02,
			  3.521E-02,
			  3.118E-02,
			  2.702E-02,
			  2.287E-02,
			  1.885E-02,
			  1.508E-02,
			  1.166E-02,
			  8.673E-03,
			  6.190E-03,
			  4.222E-03,
			  2.746E-03,
			  1.698E-03,
			  9.971E-04,
			  5.549E-04,
			  2.924E-04,
			  1.457E-04,
			  6.864E-05,
			  3.054E-05,
			  1.282E-05,
			  5.081E-06,
			  1.898E-06,
			  6.688E-07,
			  2.221E-07,
			  6.947E-08,
			  2.047E-08
			  */
	};   
	for (int i=1;i<=60;i++)  {
		PU_generated->SetBinContent(i,Summer2012[i-1]);
	}
	PU_intended->Scale( 1.0/ PU_intended->Integral() );
	PU_generated->Scale( 1.0/ PU_generated->Integral() );

	TH1F *weights_ = new TH1F( *(PU_intended)) ;

	weights_->Divide(PU_generated);


	//Re-calculate Q/G Likelihood
	//QGLikelihoodCalculator *qglikeli_Spring11    = new QGLikelihoodCalculator("./QG_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1.root");  
	//QGLikelihoodCalculator *qglikeli_Summer11    = new QGLikelihoodCalculator("./QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root");  
	//QGLikelihoodCalculator *qglikeli_Summer11CHS = new QGLikelihoodCalculator("./QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2_CHS.root");  

	// Parameter Setup
	const unsigned int jetsize         = 6;
	const double Jpt                   = 30;    // Jet pt threshold
	const double boostedZJpt           = 80;   // Conservative boosted Jet cut
	const double boostedZtranpt        = 150;
	const double btssv                 = 1.74;  // BTagging SSVHE Medium
	const double btcsvm                = 0.679; //CSVM
	const double btcsvl                = 0.244; //CSVL
	// loop over all events
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		//Long64_t ientry = LoadTree(jentry);
		//if (ientry < 0) break;
		if(jentry%100000==0){cout<< "jentry: " << jentry << endl;}
		nb = newtree->GetEntry(jentry);   nbytes += nb;
		// Cut variable definitions
		double jess    = 1.00; // control the jet energy scale
		//    double muoniso = (Z_muplus_pfiso_sumChargedHadronPt+Z_muplus_pfiso_sumNeutralHadronEt+Z_muplus_pfiso_sumPhotonEt-event_RhoForLeptonIsolation*3.141592653589*0.09)/Z_muplus_pt;
		//double muoniso = (Z_muplus_pfiso_sumChargedHadronPt+max(0.,Z_muplus_pfiso_sumNeutralHadronEt+Z_muplus_pfiso_sumPhotonEt-0.5*Z_muplus_pfiso_sumPUPt))/Z_muplus_pt;
		double dijetpt = sqrt(JetPFCor_Pt[0]*JetPFCor_Pt[0]+
					JetPFCor_Pt[1]*JetPFCor_Pt[1]+
					2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]));

		// Save variable initialization
		ggdevt    = 0;
		evtNJ     = 0;
		ggdevtinclusive    = 0;

		fit_muplus_px = 0; fit_muplus_py = 0; fit_muplus_pz = 0;  fit_muplus_e  = 0; 
		fit_muminus_px = 0; fit_muminus_py = 0; fit_muminus_pz = 0;  fit_muminus_e  = 0; 
		fit_aj_px = 0; fit_aj_py = 0; fit_aj_pz = 0;  fit_aj_e  = 0; 
		fit_bj_px = 0; fit_bj_py = 0; fit_bj_pz = 0;  fit_bj_e  = 0; 
		fit_mlljj = 0; fit_chi2  =999;fit_NDF   =999; fit_status=999;
		fit_mll   = 0; fit_mjj   = 0;

		TopWm     = 0; TopWm5j   = 0; Tchi2     =999; Tchi25j   =999;

		ang_ha  = 999; ang_hb    =999;ang_hs    =999; ang_phi   =999; 
		ang_phia= 999; ang_phib  =999;
		masslljj=-999; ptlljj    =-999; ylljj   =-999;philljj   =-999;


		effwt = 1.0; puwt = 1.0; puwt_up = 1.0; puwt_down = 1.0;
		qgld_Spring11[0]= -1;       qgld_Spring11[1]= -1;       qgld_Spring11[2]= -1;       qgld_Spring11[3]= -1;       qgld_Spring11[4]= -1;       qgld_Spring11[5]= -1;
		qgld_Summer11[0]= -1;       qgld_Summer11[1]= -1;       qgld_Summer11[2]= -1;       qgld_Summer11[3]= -1;       qgld_Summer11[4]= -1;       qgld_Summer11[5]= -1;
		qgld_Summer11CHS[0]= -1;    qgld_Summer11CHS[1]= -1;    qgld_Summer11CHS[2]= -1;    qgld_Summer11CHS[3]= -1;    qgld_Summer11CHS[4]= -1;    qgld_Summer11CHS[5]= -1;



		//CA8
		isgengdboostedZevt = 0; ggdboostedZevt = 0; GroomedJet_numberbjets_csvm = 0; GroomedJet_numberbjets_csvl = 0; GroomedJet_numberbjets_csvl_veto = 0; GroomedJet_numberbjets_csvm_veto = 0; GroomedJet_numberjets = 0;

		GroomedJet_CA8_deltaR_lplusca8jet = -999;
		GroomedJet_CA8_deltaR_lminusca8jet = -999;
		GroomedJet_CA8_deltaphi_Vca8jet = -999;

		GroomedJet_CA8_rcores01 = -1; GroomedJet_CA8_rcores02 = -1; GroomedJet_CA8_rcores03 = -1; GroomedJet_CA8_rcores04 = -1;
		GroomedJet_CA8_rcores05 = -1; GroomedJet_CA8_rcores06 = -1; GroomedJet_CA8_rcores07 = -1; GroomedJet_CA8_rcores08 = -1;
		GroomedJet_CA8_rcores09 = -1; GroomedJet_CA8_rcores10 = -1; GroomedJet_CA8_rcores11 = -1;

		GroomedJet_CA8_ptcores01 = -1; GroomedJet_CA8_ptcores02 = -1; GroomedJet_CA8_ptcores03 = -1; GroomedJet_CA8_ptcores04 = -1;
		GroomedJet_CA8_ptcores05 = -1; GroomedJet_CA8_ptcores06 = -1; GroomedJet_CA8_ptcores07 = -1; GroomedJet_CA8_ptcores08 = -1;
		GroomedJet_CA8_ptcores09 = -1; GroomedJet_CA8_ptcores10 = -1; GroomedJet_CA8_ptcores11 = -1;

		GroomedJet_CA8_planarflow01 = -1; GroomedJet_CA8_planarflow02 = -1; GroomedJet_CA8_planarflow03 = -1; GroomedJet_CA8_planarflow04 = -1;
		GroomedJet_CA8_planarflow05 = -1; GroomedJet_CA8_planarflow06 = -1; GroomedJet_CA8_planarflow07 = -1; GroomedJet_CA8_planarflow08 = -1;
		GroomedJet_CA8_planarflow09 = -1; GroomedJet_CA8_planarflow10 = -1; GroomedJet_CA8_planarflow11 = -1;

		GroomedJet_CA8_mass_sensi_tr = -1; GroomedJet_CA8_mass_sensi_ft = -1; GroomedJet_CA8_mass_sensi_pr = -1;

		GroomedJet_CA8_qjetmassvolatility = -1;

		GroomedJet_CA8_prsubjet1ptoverjetpt = -1; GroomedJet_CA8_prsubjet2ptoverjetpt = -1;
		GroomedJet_CA8_prsubjet1subjet2_deltaR = -1;
		GroomedJet_CA8_JetResp = -1;
		GroomedJet_CA8_JetResp1 = -1;
		GroomedJet_CA8_JetResp2 = -1;



		boostedZ_llj_e=-999;   boostedZ_llj_pt=-999;   boostedZ_llj_eta=-999;   boostedZ_llj_phi=-999;   boostedZ_llj_m=-999;   boostedZ_llj_y=-999;

		boostedZ_zjj_ang_ha = 999; boostedZ_zjj_ang_hb = 999; boostedZ_zjj_ang_hs = 999; boostedZ_zjj_ang_phi = 999; boostedZ_zjj_ang_phia = 999; boostedZ_zjj_ang_phib = 999;

		/*//AK7
		ggdboostedZevt_ak7 = 0; GroomedJet_numberbjets_csvm_ak7 = 0; GroomedJet_numberbjets_csvl_ak7 = 0; GroomedJet_numberbjets_csvl_veto_ak7 = 0; GroomedJet_numberbjets_csvm_veto_ak7 = 0; GroomedJet_numberjets_ak7 = 0;

		GroomedJet_AK7_deltaR_lplusak7jet = -999;
		GroomedJet_AK7_deltaR_lminusak7jet = -999;
		GroomedJet_AK7_deltaphi_Vak7jet = -999;

		GroomedJet_AK7_rcores01 = -1; GroomedJet_AK7_rcores02 = -1; GroomedJet_AK7_rcores03 = -1; GroomedJet_AK7_rcores04 = -1;
		GroomedJet_AK7_rcores05 = -1; GroomedJet_AK7_rcores06 = -1; GroomedJet_AK7_rcores07 = -1; GroomedJet_AK7_rcores08 = -1;
		GroomedJet_AK7_rcores09 = -1; GroomedJet_AK7_rcores10 = -1; GroomedJet_AK7_rcores11 = -1;

		GroomedJet_AK7_ptcores01 = -1; GroomedJet_AK7_ptcores02 = -1; GroomedJet_AK7_ptcores03 = -1; GroomedJet_AK7_ptcores04 = -1;
		GroomedJet_AK7_ptcores05 = -1; GroomedJet_AK7_ptcores06 = -1; GroomedJet_AK7_ptcores07 = -1; GroomedJet_AK7_ptcores08 = -1;
		GroomedJet_AK7_ptcores09 = -1; GroomedJet_AK7_ptcores10 = -1; GroomedJet_AK7_ptcores11 = -1;

		GroomedJet_AK7_planarflow01 = -1; GroomedJet_AK7_planarflow02 = -1; GroomedJet_AK7_planarflow03 = -1; GroomedJet_AK7_planarflow04 = -1;
		GroomedJet_AK7_planarflow05 = -1; GroomedJet_AK7_planarflow06 = -1; GroomedJet_AK7_planarflow07 = -1; GroomedJet_AK7_planarflow08 = -1;
		GroomedJet_AK7_planarflow09 = -1; GroomedJet_AK7_planarflow10 = -1; GroomedJet_AK7_planarflow11 = -1;

		GroomedJet_AK7_mass_sensi_tr = -1; GroomedJet_AK7_mass_sensi_ft = -1; GroomedJet_AK7_mass_sensi_pr = -1;

		GroomedJet_AK7_qjetmassvolatility = -1;

		GroomedJet_AK7_prsubjet1ptoverjetpt = -1; GroomedJet_AK7_prsubjet2ptoverjetpt = -1;
		GroomedJet_AK7_prsubjet1subjet2_deltaR = -1;

		boostedZ_llj_e_ak7=-999;   boostedZ_llj_pt_ak7=-999;   boostedZ_llj_eta_ak7=-999;   boostedZ_llj_phi_ak7=-999;   boostedZ_llj_m_ak7=-999;   boostedZ_llj_y_ak7=-999;

		boostedZ_zjj_ang_ha_ak7 = 999; boostedZ_zjj_ang_hb_ak7 = 999; boostedZ_zjj_ang_hs_ak7 = 999; boostedZ_zjj_ang_phi_ak7 = 999; boostedZ_zjj_ang_phia_ak7 = 999; boostedZ_zjj_ang_phib_ak7 = 999;
        */


		// Calculate efficiency
		effwt = 
			muIDEff.GetEfficiency(Z_muplus_pt, Z_muplus_eta) * 
			muHLTEff.GetEfficiency(Z_muplus_pt, Z_muplus_eta);

		// Pile up Re-weighting
		if (wda>20120999) {
			//      puwt      =    LumiWeights_.weight3D(event_mcPU_nvtx[0], event_mcPU_nvtx[1], event_mcPU_nvtx[2]);   
			puwt      =    weights_->GetBinContent(int(event_mcPU_trueInteractions+0.01)+1);   
			//      puwt_up   = up_LumiWeights_.weight3D(event_mcPU_nvtx[0], event_mcPU_nvtx[1], event_mcPU_nvtx[2]);   
			puwt_up   = puwt;
			//      puwt_down = dn_LumiWeights_.weight3D(event_mcPU_nvtx[0], event_mcPU_nvtx[1], event_mcPU_nvtx[2]);   
			puwt_down = puwt;   
		} else {effwt=1.0;puwt=1.0;puwt_up=1.0;puwt_down=1.0;} // if data, always put 1 as the weighting factor


		// Good Event Selection Requirement for all events
		int istep = 8; //starting selection step after preselection
		bool  isgengdevt = 0;
		if (JetPFCor_Pt[0]>Jpt 
					&& JetPFCor_Pt[1]>Jpt 
					&& Z_mt>30. //Move to MVA MET Later
					//&& Z_mtMVA>30.
					&& Z_muplus_pt>25.
					&& fabs(Z_muplus_dz000)<0.02
					&& fabs(Z_muplus_dzPV)<0.5
					&& fabs(Z_muplus_eta)<2.1 //Fix the Muon Eta Range to 2.1
		   ) {isgengdevt = 1;}


		if (JetPFCor_Pt[0]>Jpt ) {
			h_events          -> Fill ( istep ); 
			h_events_weighted -> Fill ( istep, effwt*puwt ); 
			istep++;
			if ( JetPFCor_Pt[1]>Jpt ) {
				h_events          -> Fill ( istep ); 
				h_events_weighted -> Fill ( istep, effwt*puwt ); 
				istep++;
				if ( Z_mt>30. ) {
					h_events          -> Fill ( istep ); 
					h_events_weighted -> Fill ( istep, effwt*puwt ); 
					istep++;
					if ( Z_muplus_pt>30. 
								&& fabs(Z_muplus_dz000)<0.02
								&& fabs(Z_muplus_dzPV)<0.5
								&& fabs(Z_muplus_eta)<2.1
					   ) { 
						h_events          -> Fill ( istep ); 
						h_events_weighted -> Fill ( istep, effwt*puwt ); 
						istep++;
					}
				}
			}
		}

		// Event Selection Requirement for Standard vs QCD events
		if ( !isQCD ) {
			//keep muons with iso<0.20(loose=0.20;tight=0.12) && event_met_pfmet>25.
			//if ( !(muoniso<0.12)         ) {isgengdevt=0;}
			if ( !(event_met_pfmet>25.0) ) {isgengdevt=0;}
			// if ( !(event_metMVA_met >25.0) ) isgengdevt=0; //Move to MVA MET Later
		} else {
			//keep muons with iso>0.20 (loose=0.20;tight=0.12)
			//if ( !(muoniso>0.12)         ) {isgengdevt=0;}
		}


		// Fill lepton information
		TLorentzVector  mup, mum;
		mup.SetPtEtaPhiE(Z_muplus_pt,              Z_muplus_eta,       Z_muplus_phi,       Z_muplus_e               );
		mum.SetPtEtaPhiE(Z_muminus_pt,              Z_muminus_eta,       Z_muminus_phi,       Z_muminus_e               );
		/*TLorentzVector b_metpt; b_metpt.SetPxPyPzE(event_met_pfmet * cos(event_met_pfmetPhi), event_met_pfmet * sin(event_met_pfmetPhi), 0, sqrt(event_met_pfmet*event_met_pfmet) );
		//TLorentzVector b_metpt; b_metpt.SetPxPyPzE(event_metMVA_met * cos(event_metMVA_metPhi),  event_metMVA_met * sin(event_metMVA_metPhi), 0, sqrt(event_metMVA_met * event_metMVA_met) );//Move to MVA MET Later
		METzCalculator b_metpz;
		b_metpz.SetMET(b_metpt);
		b_metpz.SetLepton(mup);
		b_metpz.SetLeptonType("muon");
		double b_mumz = b_metpz.Calculate(); // Default one
		TLorentzVector b_mum; b_mum.SetPxPyPzE(b_metpt.Px(), b_metpt.Py(), b_mumz, sqrt(b_metpt.Px()*b_metpt.Px() + b_metpt.Py()*b_metpt.Py() + b_mumz*b_mumz) );
		if (b_metpz.IsComplex()) {// if this is a complix, change MET
			double nu_pt1 = b_metpz.getPtneutrino(1);
			double nu_pt2 = b_metpz.getPtneutrino(2);
			TLorentzVector tmpp1; 
			tmpp1.SetPxPyPzE(nu_pt1 * cos(event_met_pfmetPhi), nu_pt1 * sin(event_met_pfmetPhi), b_mumz, sqrt(nu_pt1*nu_pt1 + b_mumz*b_mumz) );
			//tmpp1.SetPxPyPzE(nu_pt1 * cos(event_metMVA_metPhi), nu_pt1 * sin(event_metMVA_metPhi), b_mumz, sqrt(nu_pt1*nu_pt1 + b_mumz*b_mumz));//Move to MVA MET Later
			TLorentzVector tmpp2; 
			tmpp2.SetPxPyPzE(nu_pt2 * cos(event_met_pfmetPhi), nu_pt2 * sin(event_met_pfmetPhi), b_mumz, sqrt(nu_pt2*nu_pt2 + b_mumz*b_mumz) );
			//tmpp2.SetPxPyPzE(nu_pt2 * cos(event_metMVA_metPhi), nu_pt2 * sin(event_metMVA_metPhi), b_mumz, sqrt(nu_pt2*nu_pt2 + b_mumz*b_mumz) ); //Move to MVA MET Later
			b_mum = tmpp1;	if ( fabs((mup+tmpp1).M()-80.4) > fabs((mup+tmpp2).M()-80.4) ) 	b_mum = tmpp2;
		}
		*/

		//###########Begin Boosted W analysis Flag ###################################
		//Event Flag  Passed The Boosted W analysis Leptonic W Pt > 150GeV
		TLorentzVector lepwtransversep = mup + mum;
		float lepwtransversept = lepwtransversep.Pt();

		//isgengdboostedZevt = 0; //Initialized Already

		if(   GroomedJet_CA8_pt[0] > boostedZJpt 
					&& Z_mt>30. //Move to MVA MET Later
					// && Z_mtMVA>30.
					&& Z_muplus_pt>25.
					&& fabs(Z_muplus_dz000)<0.02
					&& fabs(Z_muplus_dzPV)<0.5
					&& fabs(Z_muplus_eta)<2.1
					&& lepwtransversept > boostedZtranpt 
		  ) {isgengdboostedZevt = 1;}

		// Event Selection Requirement for Standard vs QCD events
		if ( !isQCD ) {
			//keep muons with iso<0.20(loose=0.20;tight=0.12) && event_met_pfmet>25.
			//if ( !(muoniso<0.12)         ) {isgengdboostedZevt = 0;}
			if ( !(event_met_pfmet>25.0) ) {isgengdboostedZevt = 0;}
			//if ( !(event_metMVA_met >25.0) ) isgengdboostedZevt = 0; //Move to MVA MET Later
		} else {
			//keep muons with iso>0.20 (loose=0.20;tight=0.12)
			//if ( !(muoniso>0.12)         ) {isgengdboostedZevt = 0;}
		}
		//###########End Boosted W analysis Flag ###################################


		// 2 and 3 jet event for Mjj
		if (isgengdevt 
					&& fabs(JetPFCor_Eta[0]-JetPFCor_Eta[1])<1.5 ) { 
			h_events          -> Fill ( istep );  
			h_events_weighted -> Fill ( istep, effwt*puwt );  
			istep++; 
			if ( fabs(JetPFCor_dphiMET[0])>0.4 ) { 
				h_events          -> Fill ( istep );  
				h_events_weighted -> Fill ( istep, effwt*puwt );  
				istep++; 
				if ( dijetpt>40.){ 
					h_events          -> Fill ( istep );  
					h_events_weighted -> Fill ( istep, effwt*puwt );  
					istep++; 
					if ( JetPFCor_Pt[1] > Jpt && JetPFCor_Pt[2] < Jpt ) {evtNJ = 2;} 
					if ( JetPFCor_Pt[2] > Jpt && JetPFCor_Pt[3] < Jpt ) {evtNJ = 3;} 
				} 
			} 
		} 
		// 2 and 3 jet event for Hww
		if (isgengdevt) { ggdevt = 4;// Do the kinematic fit for all event!!!
			if ( JetPFCor_Pt[1] > Jpt && JetPFCor_Pt[2] < Jpt ) {ggdevt = 2;}
			if ( JetPFCor_Pt[2] > Jpt && JetPFCor_Pt[3] < Jpt ) {ggdevt = 3;}
			//Add the inclusive Jet bin
			for(size_t i = 0; i < jetsize; i++)
			{
				if(JetPFCor_Pt[i] > Jpt)
				{
					ggdevtinclusive++;
				}
			}

			int Aj = 0, Bj = 1;    TLorentzVector ajp, bjp; 
			ajp.SetPtEtaPhiE(jess * JetPFCor_Pt[Aj], JetPFCor_Eta[Aj], JetPFCor_Phi[Aj], jess * JetPFCor_E[Aj]  );
			bjp.SetPtEtaPhiE(jess * JetPFCor_Pt[Bj], JetPFCor_Eta[Bj], JetPFCor_Phi[Bj], jess * JetPFCor_E[Bj]  );

			// Do kinematic fit
			TLorentzVector fit_mup(0,0,0,0), fit_mum(0,0,0,0), fit_ajp(0,0,0,0), fit_bjp(0,0,0,0) ;
			doKinematicFit( 1, mup, mum, ajp, bjp,  fit_mup, fit_mum, fit_ajp, fit_bjp, fit_chi2, fit_NDF, fit_status);
			fit_muplus_px = fit_mup.Px(); fit_muplus_py = fit_mup.Py(); fit_muplus_pz = fit_mup.Pz(); fit_muplus_e = fit_mup.E(); 
			fit_muminus_px = fit_mum.Px(); fit_muminus_py = fit_mum.Py(); fit_muminus_pz = fit_mum.Pz(); fit_muminus_e = fit_mum.E(); 
			fit_aj_px = fit_ajp.Px(); fit_aj_py = fit_ajp.Py(); fit_aj_pz = fit_ajp.Pz(); fit_aj_e = fit_ajp.E(); 
			fit_bj_px = fit_bjp.Px(); fit_bj_py = fit_bjp.Py(); fit_bj_pz = fit_bjp.Pz(); fit_bj_e = fit_bjp.E(); 
			fit_mlljj = (fit_mup+fit_mum+fit_ajp+fit_bjp).M();
			fit_mll   = (fit_mup+fit_mum).M();
			fit_mjj   = (fit_ajp+fit_bjp).M(); 

			// Calculate angular distribution
			masslljj = (mup+mum+ajp+bjp).M();
			ptlljj   = (mup+mum+ajp+bjp).Pt();
			ylljj    = (mup+mum+ajp+bjp).Rapidity();
			philljj  = (mup+mum+ajp+bjp).Phi();
			double a_costheta1, a_costheta2, a_phi, a_costhetastar, a_phistar1, a_phistar2;
			if (Z_muplus_charge < 0){
				calculateAngles(mup, mum, ajp, bjp, a_costheta1, a_costheta2, a_phi, a_costhetastar, a_phistar1, a_phistar2);
			}
			else{
				calculateAngles(mum, mup, ajp, bjp, a_costheta1, a_costheta2, a_phi, a_costhetastar, a_phistar1, a_phistar2);
			}
			ang_ha = a_costheta1; ang_hb = fabs(a_costheta2); ang_hs = a_costhetastar;  ang_phi = a_phi; ang_phia = a_phistar1; ang_phib = a_phistar2;

		}
		// For Hadronic W in Top sample
		if (isgengdevt)
		{
			if (JetPFCor_Pt[3] > Jpt && JetPFCor_Pt[4] < Jpt){
				int nbjet = 0;
				int nbnot = 0;
				int Aj    = -999;
				int Bj    = -999;
				/*if (JetPFCor_bDiscriminator[0]>btssv) { nbjet++; } else { nbnot++; if (nbnot==1) Aj=0; if (nbnot==2) Bj=0;}
				  if (JetPFCor_bDiscriminator[1]>btssv) { nbjet++; } else { nbnot++; if (nbnot==1) Aj=1; if (nbnot==2) Bj=1;}
				  if (JetPFCor_bDiscriminator[2]>btssv) { nbjet++; } else { nbnot++; if (nbnot==1) Aj=2; if (nbnot==2) Bj=2;}
				  if (JetPFCor_bDiscriminator[3]>btssv) { nbjet++; } else { nbnot++; if (nbnot==1) Aj=3; if (nbnot==2) Bj=3;}
				  */
				if (JetPFCor_bDiscriminatorCSV[0]>btcsvm) { nbjet++; } else { nbnot++; if (nbnot==1) Aj=0; if (nbnot==2) Bj=0;}
				if (JetPFCor_bDiscriminatorCSV[1]>btcsvm) { nbjet++; } else { nbnot++; if (nbnot==1) Aj=1; if (nbnot==2) Bj=1;}
				if (JetPFCor_bDiscriminatorCSV[2]>btcsvm) { nbjet++; } else { nbnot++; if (nbnot==1) Aj=2; if (nbnot==2) Bj=2;}
				if (JetPFCor_bDiscriminatorCSV[3]>btcsvm) { nbjet++; } else { nbnot++; if (nbnot==1) Aj=3; if (nbnot==2) Bj=3;}

				if (nbjet==2 && nbnot==2 && Aj!=-999 && Bj!=-999){
					TLorentzVector  ajp, bjp; 
					ajp.SetPtEtaPhiE(jess * JetPFCor_Pt[Aj], JetPFCor_Eta[Aj], JetPFCor_Phi[Aj], jess * JetPFCor_E[Aj]  );
					bjp.SetPtEtaPhiE(jess * JetPFCor_Pt[Bj], JetPFCor_Eta[Bj], JetPFCor_Phi[Bj], jess * JetPFCor_E[Bj]  );
					TopWm   = (ajp+bjp).M(); 

					TLorentzVector fit_mup(0,0,0,0), fit_mum(0,0,0,0), fit_ajp(0,0,0,0), fit_bjp(0,0,0,0) ; Int_t tmpa =0, tmpb=0;
					doKinematicFit( 1, mup, mum, ajp, bjp,  fit_mup, fit_mum, fit_ajp, fit_bjp, Tchi2, tmpa, tmpb);
				}
			}
		}
		if (isgengdevt)
		{
			if (JetPFCor_Pt[4] > Jpt && JetPFCor_Pt[5] < Jpt){
				int nbjet = 0;
				int nbnot = 0;
				int Aj    = -999;
				int Bj    = -999;
				/*if (JetPFCor_bDiscriminator[0]>btssv) { nbjet++; } else { nbnot++; if (nbnot==1) Aj=0; if (nbnot==2) Bj=0;}
				  if (JetPFCor_bDiscriminator[1]>btssv) { nbjet++; } else { nbnot++; if (nbnot==1) Aj=1; if (nbnot==2) Bj=1;}
				  if (JetPFCor_bDiscriminator[2]>btssv) { nbjet++; } else { nbnot++; if (nbnot==1) Aj=2; if (nbnot==2) Bj=2;}
				  if (JetPFCor_bDiscriminator[3]>btssv) { nbjet++; } else { nbnot++; if (nbnot==1) Aj=3; if (nbnot==2) Bj=3;}
				  if (JetPFCor_bDiscriminator[4]>btssv) { nbjet++; } else { nbnot++; if (nbnot==1) Aj=4; if (nbnot==2) Bj=4;}
				  */
				if (JetPFCor_bDiscriminatorCSV[0]>btcsvm) { nbjet++; } else { nbnot++; if (nbnot==1) Aj=0; if (nbnot==2) Bj=0;}
				if (JetPFCor_bDiscriminatorCSV[1]>btcsvm) { nbjet++; } else { nbnot++; if (nbnot==1) Aj=1; if (nbnot==2) Bj=1;}
				if (JetPFCor_bDiscriminatorCSV[2]>btcsvm) { nbjet++; } else { nbnot++; if (nbnot==1) Aj=2; if (nbnot==2) Bj=2;}
				if (JetPFCor_bDiscriminatorCSV[3]>btcsvm) { nbjet++; } else { nbnot++; if (nbnot==1) Aj=3; if (nbnot==2) Bj=3;}
				if (JetPFCor_bDiscriminatorCSV[4]>btcsvm) { nbjet++; } else { nbnot++; if (nbnot==1) Aj=4; if (nbnot==2) Bj=4;}

				if (nbjet==2 && nbnot==3 && Aj!=-999 && Bj!=-999){
					TLorentzVector  ajp, bjp; 
					ajp.SetPtEtaPhiE(jess * JetPFCor_Pt[Aj], JetPFCor_Eta[Aj], JetPFCor_Phi[Aj], jess * JetPFCor_E[Aj]  );
					bjp.SetPtEtaPhiE(jess * JetPFCor_Pt[Bj], JetPFCor_Eta[Bj], JetPFCor_Phi[Bj], jess * JetPFCor_E[Bj]  );
					TopWm5j = (ajp+bjp).M(); 

					TLorentzVector fit_mup(0,0,0,0), fit_mum(0,0,0,0), fit_ajp(0,0,0,0), fit_bjp(0,0,0,0) ; Int_t tmpa =0, tmpb=0;
					doKinematicFit( 1, mup, mum, ajp, bjp,  fit_mup, fit_mum, fit_ajp, fit_bjp, Tchi25j, tmpa, tmpb);
				}
			}
		}

		//################Begin Boosted W Analysis########################################
		if(isgengdboostedZevt){

			//CA8 Jet And we also need to check the AK7 Jet
			TLorentzVector ca8jetp4;
			ca8jetp4.SetPtEtaPhiE(GroomedJet_CA8_pt[0], GroomedJet_CA8_eta[0], GroomedJet_CA8_phi[0], GroomedJet_CA8_e[0]);
			double deltaR_lplusca8jet = mup.DeltaR(ca8jetp4);
			double deltaR_lminusca8jet = mum.DeltaR(ca8jetp4);
			TLorentzVector wbosonp = mup + mum;
			//double deltaphi_Vca8jet = wbosonp.DeltaPhi(ca8jetp4);
			double deltaphi_Vca8jet = getDeltaPhi(wbosonp.Phi(),ca8jetp4.Phi());

			GroomedJet_CA8_deltaR_lplusca8jet = deltaR_lplusca8jet;
			GroomedJet_CA8_deltaR_lminusca8jet = deltaR_lminusca8jet;
			GroomedJet_CA8_deltaphi_Vca8jet = deltaphi_Vca8jet;

			//Count the number of B tag jet, for ttbar and contral plots
			//Synchrinize the Cuts with the Z'->tt group using the semi-leptonic channel selection
			for(int i = 0; i < numPFCorJets; i++)
			{
				if(JetPFCor_Pt[i] > Jpt)
				{
					TLorentzVector  ajp;
					ajp.SetPtEtaPhiE(jess * JetPFCor_Pt[i], JetPFCor_Eta[i], JetPFCor_Phi[i], jess * JetPFCor_E[i]  );

					double tmpdelatR = ca8jetp4.DeltaR(ajp);

					double tmpdeltaRlj = ajp.DeltaR(mup);

					if(tmpdelatR > 0.8)//Veto the AK5 jet in the CA8 jet cone
					{
						GroomedJet_numberjets = GroomedJet_numberjets + 1;
					}

					//if(JetPFCor_bDiscriminator[i] > btssv && tmpdelatR > 0.8)//Veto the AK5 jet in the CA8 jet cone
					//if(JetPFCor_bDiscriminatorCSV[i] > btcsvm && tmpdelatR > 0.8)//Veto the AK5 jet in the CA8 jet cone and Move to CSVM tagger
					if(JetPFCor_bDiscriminatorCSV[i] > btcsvm && tmpdelatR > 0.8 && tmpdeltaRlj > TMath::Pi()/2.0)//Veto the AK5 jet in the CA8 jet cone and Move to CSVM tagger and leptonic and hadronic top should be in the different hemisphere for ttbar control
					{
						GroomedJet_numberbjets_csvm = GroomedJet_numberbjets_csvm + 1;
					}

					if(JetPFCor_bDiscriminatorCSV[i] > btcsvl && tmpdelatR > 0.8 && tmpdeltaRlj > TMath::Pi()/2.0)//Veto the AK5 jet in the CA8 jet cone and Move to CSVM tagger and leptonic and hadronic top should be in the different hemisphere for ttbar control
					{
						GroomedJet_numberbjets_csvl = GroomedJet_numberbjets_csvl + 1;
					}

					/*if(JetPFCor_bDiscriminator[i] > btssv && tmpdelatR > 0.8 && tmpdeltaRlj > TMath::Pi()/2.0)//Veto the AK5 jet in the CA8 jet cone and Move to CSVM tagger and leptonic and hadronic top should be in the different hemisphere for ttbar control
					{
						GroomedJet_numberbjets_ssvhem = GroomedJet_numberbjets_ssvhem + 1;
					}*/

					if(JetPFCor_bDiscriminatorCSV[i] > btcsvl)//ttbar veto
					{
						GroomedJet_numberbjets_csvl_veto = GroomedJet_numberbjets_csvl_veto + 1;
					}

					if(JetPFCor_bDiscriminatorCSV[i] > btcsvm)//ttbar veto
					{
						GroomedJet_numberbjets_csvm_veto = GroomedJet_numberbjets_csvm_veto + 1;
					}

					/*if(JetPFCor_bDiscriminator[i] > btssv)//ttbar veto
					{
						GroomedJet_numberbjets_ssvhem_veto = GroomedJet_numberbjets_ssvhem_veto + 1;
					}*/
				}
			}

			//if(GroomedJet_CA8_pt[0] > boostedZJpt && GroomedJet_CA8_pt[1] < boostedZJpt && deltaR_lca8jet > 1.0 && deltaphi_METca8jet > 0.4 && deltaphi_Vca8jet > 2.0) 
			//if(deltaR_lca8jet > 1.0 && deltaphi_METca8jet > 0.4 && deltaphi_Vca8jet > 2.0) 
			if(deltaR_lplusca8jet > TMath::Pi()/ 2.0 && deltaR_lminusca8jet > TMath::Pi()/ 2.0 && deltaphi_Vca8jet > 2.0)//Tighter Angular Cuts 
			{
				ggdboostedZevt = 1;
			}

			GroomedJet_CA8_rcores01 = GroomedJet_CA8_rcores[0][0];
			GroomedJet_CA8_ptcores01 = GroomedJet_CA8_ptcores[0][0];
			GroomedJet_CA8_planarflow01 = GroomedJet_CA8_planarflow[0][0];

			GroomedJet_CA8_rcores02 = GroomedJet_CA8_rcores[1][0];
			GroomedJet_CA8_ptcores02 = GroomedJet_CA8_ptcores[1][0];
			GroomedJet_CA8_planarflow02 = GroomedJet_CA8_planarflow[1][0];

			GroomedJet_CA8_rcores03 = GroomedJet_CA8_rcores[2][0];
			GroomedJet_CA8_ptcores03 = GroomedJet_CA8_ptcores[2][0];
			GroomedJet_CA8_planarflow03 = GroomedJet_CA8_planarflow[2][0];

			GroomedJet_CA8_rcores04 = GroomedJet_CA8_rcores[3][0];
			GroomedJet_CA8_ptcores04 = GroomedJet_CA8_ptcores[3][0];
			GroomedJet_CA8_planarflow04 = GroomedJet_CA8_planarflow[3][0];

			GroomedJet_CA8_rcores05 = GroomedJet_CA8_rcores[4][0];
			GroomedJet_CA8_ptcores05 = GroomedJet_CA8_ptcores[4][0];
			GroomedJet_CA8_planarflow05 = GroomedJet_CA8_planarflow[4][0];

			GroomedJet_CA8_rcores06 = GroomedJet_CA8_rcores[5][0];
			GroomedJet_CA8_ptcores06 = GroomedJet_CA8_ptcores[5][0];
			GroomedJet_CA8_planarflow06 = GroomedJet_CA8_planarflow[5][0];

			GroomedJet_CA8_rcores07 = GroomedJet_CA8_rcores[6][0];
			GroomedJet_CA8_ptcores07 = GroomedJet_CA8_ptcores[6][0];
			GroomedJet_CA8_planarflow07 = GroomedJet_CA8_planarflow[6][0];

			GroomedJet_CA8_rcores08 = GroomedJet_CA8_rcores[7][0];
			GroomedJet_CA8_ptcores08 = GroomedJet_CA8_ptcores[7][0];
			GroomedJet_CA8_planarflow08 = GroomedJet_CA8_planarflow[7][0];

			GroomedJet_CA8_rcores09 = GroomedJet_CA8_rcores[8][0];
			GroomedJet_CA8_ptcores09 = GroomedJet_CA8_ptcores[8][0];
			GroomedJet_CA8_planarflow09 = GroomedJet_CA8_planarflow[8][0];

			GroomedJet_CA8_rcores10 = GroomedJet_CA8_rcores[9][0];
			GroomedJet_CA8_ptcores10 = GroomedJet_CA8_ptcores[9][0];
			GroomedJet_CA8_planarflow10 = GroomedJet_CA8_planarflow[9][0];

			GroomedJet_CA8_rcores11 = GroomedJet_CA8_rcores[10][0];
			GroomedJet_CA8_ptcores11 = GroomedJet_CA8_ptcores[10][0];
			GroomedJet_CA8_planarflow11 = GroomedJet_CA8_planarflow[10][0];

			GroomedJet_CA8_mass_sensi_tr = GroomedJet_CA8_mass_tr[0]/GroomedJet_CA8_mass[0];
			GroomedJet_CA8_mass_sensi_ft = GroomedJet_CA8_mass_ft[0]/GroomedJet_CA8_mass[0];
			GroomedJet_CA8_mass_sensi_pr = GroomedJet_CA8_mass_pr[0]/GroomedJet_CA8_mass[0];

			//QJet mass Volatility
			Int_t qjetsize = 50;
			double averagemsquare = 0;
			double averagem = 0;
			for(Int_t i = 0; i < qjetsize; i++)
			{
				averagemsquare = averagemsquare + GroomedJet_CA8_qjetmass[i] * GroomedJet_CA8_qjetmass[i];
				averagem = averagem + GroomedJet_CA8_qjetmass[i];
			} 

			averagemsquare = averagemsquare / qjetsize;
			averagem = averagem / qjetsize;

			GroomedJet_CA8_qjetmassvolatility = TMath::Sqrt(averagemsquare - TMath::Power(averagem,2))/averagem;

			TLorentzVector ca8subjet1p4;
			TLorentzVector ca8subjet2p4;
			TLorentzVector ca8prjetp4;

			ca8subjet1p4.SetPxPyPzE(GroomedJet_CA8_prsubjet1_px[0],GroomedJet_CA8_prsubjet1_py[0],GroomedJet_CA8_prsubjet1_pz[0],GroomedJet_CA8_prsubjet1_e[0]);
			ca8subjet2p4.SetPxPyPzE(GroomedJet_CA8_prsubjet2_px[0],GroomedJet_CA8_prsubjet2_py[0],GroomedJet_CA8_prsubjet2_pz[0],GroomedJet_CA8_prsubjet2_e[0]);
			ca8prjetp4.SetPtEtaPhiE(GroomedJet_CA8_pt_pr[0],GroomedJet_CA8_eta_pr[0],GroomedJet_CA8_phi_pr[0],GroomedJet_CA8_e_pr[0]);

			GroomedJet_CA8_prsubjet1ptoverjetpt = ca8subjet1p4.Pt()/ca8prjetp4.Pt();
			GroomedJet_CA8_prsubjet2ptoverjetpt = ca8subjet2p4.Pt()/ca8prjetp4.Pt();

			if((ca8subjet1p4.Pt() > 0.001) && (ca8subjet2p4.Pt() > 0.001)) //Avoid Too Samll Pt
			{GroomedJet_CA8_prsubjet1subjet2_deltaR = ca8subjet1p4.DeltaR(ca8subjet2p4);}

			GroomedJet_CA8_JetResp = GroomedJet_CA8_pt[0]/Z_pt;
			double subptxx=ca8subjet1p4.Pt();
			if(ca8subjet1p4.Pt()<ca8subjet2p4.Pt())  subptxx=ca8subjet2p4.Pt();
			GroomedJet_CA8_JetResp1 = subptxx/Z_pt;
			if( TMath::Min( ca8subjet1p4.Pt(), ca8subjet2p4.Pt() )/Z_pt < 0.1) 
			{GroomedJet_CA8_JetResp2 = subptxx/Z_pt;}


			//Angular Correlation For the Boosted W Analysis
			boostedZ_llj_e      = (mup+mum+ca8jetp4).E();
			boostedZ_llj_pt     = (mup+mum+ca8jetp4).Pt();
			boostedZ_llj_eta    = (mup+mum+ca8jetp4).Eta();
			boostedZ_llj_phi    = (mup+mum+ca8jetp4).Phi();
			boostedZ_llj_m      = (mup+mum+ca8jetp4).M();
			boostedZ_llj_y      = (mup+mum+ca8jetp4).Rapidity();

			double a_costheta1, a_costheta2, a_phi, a_costhetastar, a_phistar1, a_phistar2;
			//Use the Subjet in the Boosted W Analyisis
			if (Z_muplus_charge < 0){
				calculateAngles(mup, mum, ca8subjet1p4, ca8subjet2p4, a_costheta1, a_costheta2, a_phi, a_costhetastar, a_phistar1, a_phistar2);
			}
			else{
				calculateAngles(mum, mup, ca8subjet1p4, ca8subjet2p4, a_costheta1, a_costheta2, a_phi, a_costhetastar, a_phistar1, a_phistar2);
			}
			boostedZ_zjj_ang_ha = a_costheta1; boostedZ_zjj_ang_hb = fabs(a_costheta2); boostedZ_zjj_ang_hs = a_costhetastar;  boostedZ_zjj_ang_phi = a_phi; boostedZ_zjj_ang_phia = a_phistar1; boostedZ_zjj_ang_phib = a_phistar2;

			//Input For the TMVA Training And Classification TODO
/*
			//check the AK7 Jet Used in our Reduced Tree
			TLorentzVector ak7jetp4;
			ak7jetp4.SetPtEtaPhiE(GroomedJet_AK7_pt[0], GroomedJet_AK7_eta[0], GroomedJet_AK7_phi[0], GroomedJet_AK7_e[0]);
			double deltaR_lplusak7jet = mup.DeltaR(ak7jetp4);
			double deltaR_lminusak7jet = mum.DeltaR(ak7jetp4);
			//TLorentzVector wbosonp = mup + mum;
			//double deltaphi_Vak7jet = wbosonp.DeltaPhi(ak7jetp4);
			double deltaphi_Vak7jet = getDeltaPhi(wbosonp.Phi(),ak7jetp4.Phi());

			GroomedJet_AK7_deltaR_lplusak7jet = deltaR_lplusak7jet;
			GroomedJet_AK7_deltaR_lminusak7jet = deltaR_lminusak7jet;
			GroomedJet_AK7_deltaphi_Vak7jet = deltaphi_Vak7jet;

			//Count the number of B tag jet, for ttbar and contral plots
			for(int i = 0; i < numPFCorJets; i++)
			{
				if(JetPFCor_Pt[i] > Jpt)
				{
					TLorentzVector  ajp;
					ajp.SetPtEtaPhiE(jess * JetPFCor_Pt[i], JetPFCor_Eta[i], JetPFCor_Phi[i], jess * JetPFCor_E[i]  );

					double tmpdelatR = ak7jetp4.DeltaR(ajp);
					double tmpdeltaRlj = ajp.DeltaR(mup);

					if(tmpdelatR > 0.7)//Veto the AK5 jet in the AK7 jet cone
					{
						GroomedJet_numberjets_ak7 = GroomedJet_numberjets_ak7 + 1;
					}

					//if(JetPFCor_bDiscriminator[i] > btssv && tmpdelatR > 0.8)//Veto the AK5 jet in the AK7 jet cone
					if(JetPFCor_bDiscriminatorCSV[i] > btcsvm && tmpdelatR > 0.7 && tmpdeltaRlj > TMath::Pi()/2.0)//Veto the AK5 jet in the AK7 jet cone and Move to CSVM tagger
					{
						GroomedJet_numberbjets_csvm_ak7 = GroomedJet_numberbjets_csvm_ak7 + 1;
					}
					if(JetPFCor_bDiscriminatorCSV[i] > btcsvl && tmpdelatR > 0.7 && tmpdeltaRlj > TMath::Pi()/2.0)//Veto the AK5 jet in the AK7 jet cone and Move to CSVM tagger
					{
						GroomedJet_numberbjets_csvl_ak7 = GroomedJet_numberbjets_csvl_ak7 + 1;
					}
					if(JetPFCor_bDiscriminatorCSV[i] > btcsvl)//ttbar veto
					{
						GroomedJet_numberbjets_csvl_veto_ak7 = GroomedJet_numberbjets_csvl_veto_ak7 + 1;
					}
					if(JetPFCor_bDiscriminatorCSV[i] > btcsvm)//ttbar veto
					{
						GroomedJet_numberbjets_csvm_veto_ak7 = GroomedJet_numberbjets_csvm_veto_ak7 + 1;
					}
				}
			}

			//if(deltaR_lak7jet > 1.0 && deltaphi_METak7jet > 0.4 && deltaphi_Vak7jet > 2.0) 
			if(deltaR_lplusak7jet > TMath::Pi()/2.0 && deltaR_lminusak7jet > 2.0 && deltaphi_Vak7jet > 2.0)//Tighter Angular Cuts 
			{
				ggdboostedZevt_ak7 = 1;
			}

			GroomedJet_AK7_rcores01 = GroomedJet_AK7_rcores[0][0];
			GroomedJet_AK7_ptcores01 = GroomedJet_AK7_ptcores[0][0];
			GroomedJet_AK7_planarflow01 = GroomedJet_AK7_planarflow[0][0];

			GroomedJet_AK7_rcores02 = GroomedJet_AK7_rcores[1][0];
			GroomedJet_AK7_ptcores02 = GroomedJet_AK7_ptcores[1][0];
			GroomedJet_AK7_planarflow02 = GroomedJet_AK7_planarflow[1][0];

			GroomedJet_AK7_rcores03 = GroomedJet_AK7_rcores[2][0];
			GroomedJet_AK7_ptcores03 = GroomedJet_AK7_ptcores[2][0];
			GroomedJet_AK7_planarflow03 = GroomedJet_AK7_planarflow[2][0];

			GroomedJet_AK7_rcores04 = GroomedJet_AK7_rcores[3][0];
			GroomedJet_AK7_ptcores04 = GroomedJet_AK7_ptcores[3][0];
			GroomedJet_AK7_planarflow04 = GroomedJet_AK7_planarflow[3][0];

			GroomedJet_AK7_rcores05 = GroomedJet_AK7_rcores[4][0];
			GroomedJet_AK7_ptcores05 = GroomedJet_AK7_ptcores[4][0];
			GroomedJet_AK7_planarflow05 = GroomedJet_AK7_planarflow[4][0];

			GroomedJet_AK7_rcores06 = GroomedJet_AK7_rcores[5][0];
			GroomedJet_AK7_ptcores06 = GroomedJet_AK7_ptcores[5][0];
			GroomedJet_AK7_planarflow06 = GroomedJet_AK7_planarflow[5][0];

			GroomedJet_AK7_rcores07 = GroomedJet_AK7_rcores[6][0];
			GroomedJet_AK7_ptcores07 = GroomedJet_AK7_ptcores[6][0];
			GroomedJet_AK7_planarflow07 = GroomedJet_AK7_planarflow[6][0];

			GroomedJet_AK7_rcores08 = GroomedJet_AK7_rcores[7][0];
			GroomedJet_AK7_ptcores08 = GroomedJet_AK7_ptcores[7][0];
			GroomedJet_AK7_planarflow08 = GroomedJet_AK7_planarflow[7][0];

			GroomedJet_AK7_rcores09 = GroomedJet_AK7_rcores[8][0];
			GroomedJet_AK7_ptcores09 = GroomedJet_AK7_ptcores[8][0];
			GroomedJet_AK7_planarflow09 = GroomedJet_AK7_planarflow[8][0];

			GroomedJet_AK7_rcores10 = GroomedJet_AK7_rcores[9][0];
			GroomedJet_AK7_ptcores10 = GroomedJet_AK7_ptcores[9][0];
			GroomedJet_AK7_planarflow10 = GroomedJet_AK7_planarflow[9][0];

			GroomedJet_AK7_rcores11 = GroomedJet_AK7_rcores[10][0];
			GroomedJet_AK7_ptcores11 = GroomedJet_AK7_ptcores[10][0];
			GroomedJet_AK7_planarflow11 = GroomedJet_AK7_planarflow[10][0];

			GroomedJet_AK7_mass_sensi_tr = GroomedJet_AK7_mass_tr[0]/GroomedJet_AK7_mass[0];
			GroomedJet_AK7_mass_sensi_ft = GroomedJet_AK7_mass_ft[0]/GroomedJet_AK7_mass[0];
			GroomedJet_AK7_mass_sensi_pr = GroomedJet_AK7_mass_pr[0]/GroomedJet_AK7_mass[0];

			//QJet mass Volatility
			Int_t qjetsizeak7 = 50;
			double averagemsquareak7 = 0;
			double averagemak7 = 0;
			for(Int_t i = 0; i < qjetsizeak7; i++)
			{
				averagemsquareak7 = averagemsquareak7 + GroomedJet_AK7_qjetmass[i] * GroomedJet_AK7_qjetmass[i];
				averagemak7 = averagemak7 + GroomedJet_AK7_qjetmass[i];
			} 

			averagemsquareak7 = averagemsquareak7 / qjetsizeak7;
			averagemak7 = averagemak7 / qjetsizeak7;

			GroomedJet_AK7_qjetmassvolatility = TMath::Sqrt(averagemsquareak7 - TMath::Power(averagemak7,2))/averagemak7;

			TLorentzVector ak7subjet1p4;
			TLorentzVector ak7subjet2p4;
			TLorentzVector ak7prjetp4;

			ak7subjet1p4.SetPxPyPzE(GroomedJet_AK7_prsubjet1_px[0],GroomedJet_AK7_prsubjet1_py[0],GroomedJet_AK7_prsubjet1_pz[0],GroomedJet_AK7_prsubjet1_e[0]);
			ak7subjet2p4.SetPxPyPzE(GroomedJet_AK7_prsubjet2_px[0],GroomedJet_AK7_prsubjet2_py[0],GroomedJet_AK7_prsubjet2_pz[0],GroomedJet_AK7_prsubjet2_e[0]);
			ak7prjetp4.SetPtEtaPhiE(GroomedJet_AK7_pt_pr[0],GroomedJet_AK7_eta_pr[0],GroomedJet_AK7_phi_pr[0],GroomedJet_AK7_e_pr[0]);

			GroomedJet_AK7_prsubjet1ptoverjetpt = ak7subjet1p4.Pt()/ak7prjetp4.Pt();
			GroomedJet_AK7_prsubjet2ptoverjetpt = ak7subjet2p4.Pt()/ak7prjetp4.Pt();

			if((ak7subjet1p4.Pt() > 0.001) && (ak7subjet2p4.Pt() > 0.001)) //Avoid Too Samll Pt
			{GroomedJet_AK7_prsubjet1subjet2_deltaR = ak7subjet1p4.DeltaR(ak7subjet2p4);}

			//Angular Correlation For the Boosted W Analysis
			boostedZ_llj_e_ak7      = (mup+mum+ak7jetp4).E();
			boostedZ_llj_pt_ak7     = (mup+mum+ak7jetp4).Pt();
			boostedZ_llj_eta_ak7    = (mup+mum+ak7jetp4).Eta();
			boostedZ_llj_phi_ak7    = (mup+mum+ak7jetp4).Phi();
			boostedZ_llj_m_ak7      = (mup+mum+ak7jetp4).M();
			boostedZ_llj_y_ak7      = (mup+mum+ak7jetp4).Rapidity();

			double tmpa_costheta1, tmpa_costheta2, tmpa_phi, tmpa_costhetastar, tmpa_phistar1, tmpa_phistar2;
			//Use the Subjet in the Boosted W Analyisis
			if (Z_muplus_charge < 0){
				calculateAngles(mup, mum, ak7subjet1p4, ak7subjet2p4, tmpa_costheta1, tmpa_costheta2, tmpa_phi, tmpa_costhetastar, tmpa_phistar1, tmpa_phistar2);
			}
			else{
				calculateAngles(mum, mup, ak7subjet1p4, ak7subjet2p4, tmpa_costheta1, tmpa_costheta2, tmpa_phi, tmpa_costhetastar, tmpa_phistar1, tmpa_phistar2);
			}
			boostedZ_zjj_ang_ha_ak7 = tmpa_costheta1; boostedZ_zjj_ang_hb_ak7 = fabs(tmpa_costheta2); boostedZ_zjj_ang_hs_ak7 = tmpa_costhetastar;  boostedZ_zjj_ang_phi_ak7 = tmpa_phi; boostedZ_zjj_ang_phia_ak7 = tmpa_phistar1; boostedZ_zjj_ang_phib_ak7 = tmpa_phistar2;
*/
		}
		//###############End Boosted W Analysis########################################

		branch_ggdevt->Fill();
		branch_evtNJ ->Fill();
		branch_ggdevtinclusive->Fill();

		branch_muplus_px->Fill();
		branch_muplus_py->Fill();
		branch_muplus_pz->Fill();
		branch_muplus_e ->Fill();

		branch_muminus_px->Fill();
		branch_muminus_py->Fill();
		branch_muminus_pz->Fill();
		branch_muminus_e ->Fill();

		branch_aj_px->Fill();
		branch_aj_py->Fill();
		branch_aj_pz->Fill();
		branch_aj_e ->Fill();

		branch_bj_px->Fill();
		branch_bj_py->Fill();
		branch_bj_pz->Fill();
		branch_bj_e ->Fill();

		branch_mlljj->Fill();
		branch_mll  ->Fill();
		branch_mjj  ->Fill();
		branch_chi2 ->Fill();
		branch_NDF  ->Fill();
		branch_status->Fill();

		branch_TopWm->Fill();
		branch_TopWm5j->Fill();
		branch_Tchi2->Fill();
		branch_Tchi25j->Fill();

		branch_ha->Fill();   
		branch_hb->Fill();   
		branch_hs->Fill();  
		branch_phi->Fill(); 
		branch_phia->Fill();
		branch_phib->Fill();
		branch_orgm->Fill();
		branch_orgpt->Fill();
		branch_orgy->Fill();
		branch_orgph->Fill();


		branch_effwt->Fill();
		branch_puwt->Fill();
		branch_puwt_up->Fill();
		branch_puwt_down->Fill();


		branch_qgld_Spring11->Fill();
		branch_qgld_Summer11->Fill();
		branch_qgld_Summer11CHS->Fill();


		//Boosted W Fill
		branch_isgengdboostedZevt->Fill();
		branch_ggdboostedZevt->Fill();
		branch_GroomedJet_CA8_deltaR_lplusca8jet->Fill();
		branch_GroomedJet_CA8_deltaR_lminusca8jet->Fill();
		//branch_GroomedJet_CA8_deltaphi_METca8jet->Fill();
		branch_GroomedJet_CA8_deltaphi_Vca8jet->Fill();
		branch_GroomedJet_numberbjets_csvl->Fill();
		branch_GroomedJet_numberbjets_csvm->Fill();
		//branch_GroomedJet_numberbjets_ssvhem->Fill();
		branch_GroomedJet_numberbjets_csvl_veto->Fill();
		branch_GroomedJet_numberbjets_csvm_veto->Fill();
		//branch_GroomedJet_numberbjets_ssvhem_veto->Fill();
		branch_GroomedJet_numberjets->Fill();
		branch_GroomedJet_CA8_rcores01->Fill();
		branch_GroomedJet_CA8_rcores02->Fill();
		branch_GroomedJet_CA8_rcores03->Fill();
		branch_GroomedJet_CA8_rcores04->Fill();
		branch_GroomedJet_CA8_rcores05->Fill();
		branch_GroomedJet_CA8_rcores06->Fill();
		branch_GroomedJet_CA8_rcores07->Fill();
		branch_GroomedJet_CA8_rcores08->Fill();
		branch_GroomedJet_CA8_rcores09->Fill();
		branch_GroomedJet_CA8_rcores10->Fill();
		branch_GroomedJet_CA8_rcores11->Fill();

		branch_GroomedJet_CA8_ptcores01->Fill();
		branch_GroomedJet_CA8_ptcores02->Fill();
		branch_GroomedJet_CA8_ptcores03->Fill();
		branch_GroomedJet_CA8_ptcores04->Fill();
		branch_GroomedJet_CA8_ptcores05->Fill();
		branch_GroomedJet_CA8_ptcores06->Fill();
		branch_GroomedJet_CA8_ptcores07->Fill();
		branch_GroomedJet_CA8_ptcores08->Fill();
		branch_GroomedJet_CA8_ptcores09->Fill();
		branch_GroomedJet_CA8_ptcores10->Fill();
		branch_GroomedJet_CA8_ptcores11->Fill();

		branch_GroomedJet_CA8_planarflow01->Fill();
		branch_GroomedJet_CA8_planarflow02->Fill();
		branch_GroomedJet_CA8_planarflow03->Fill();
		branch_GroomedJet_CA8_planarflow04->Fill();
		branch_GroomedJet_CA8_planarflow05->Fill();
		branch_GroomedJet_CA8_planarflow06->Fill();
		branch_GroomedJet_CA8_planarflow07->Fill();
		branch_GroomedJet_CA8_planarflow08->Fill();
		branch_GroomedJet_CA8_planarflow09->Fill();
		branch_GroomedJet_CA8_planarflow10->Fill();
		branch_GroomedJet_CA8_planarflow11->Fill();

		branch_GroomedJet_CA8_mass_sensi_tr->Fill();
		branch_GroomedJet_CA8_mass_sensi_ft->Fill();
		branch_GroomedJet_CA8_mass_sensi_pr->Fill();

		branch_GroomedJet_CA8_qjetmassvolatility->Fill();

		branch_GroomedJet_CA8_prsubjet1ptoverjetpt->Fill();
		branch_GroomedJet_CA8_prsubjet2ptoverjetpt->Fill();
		branch_GroomedJet_CA8_prsubjet1subjet2_deltaR->Fill();

		branch_GroomedJet_CA8_JetResp->Fill();
		branch_GroomedJet_CA8_JetResp1->Fill();
		branch_GroomedJet_CA8_JetResp2->Fill();


		branch_boostedZ_llj_e->Fill();
		branch_boostedZ_llj_pt->Fill();
		branch_boostedZ_llj_eta->Fill();
		branch_boostedZ_llj_phi->Fill();
		branch_boostedZ_llj_m->Fill();
		branch_boostedZ_llj_y->Fill();

		branch_boostedZ_zjj_ang_ha->Fill();
		branch_boostedZ_zjj_ang_hb->Fill();
		branch_boostedZ_zjj_ang_hs->Fill();
		branch_boostedZ_zjj_ang_phi->Fill();
		branch_boostedZ_zjj_ang_phia->Fill();
		branch_boostedZ_zjj_ang_phib->Fill();

/*		//AK7 Jet Algorithm
		branch_ggdboostedZevt_ak7->Fill();
		branch_GroomedJet_AK7_deltaR_lplusak7jet->Fill();
		branch_GroomedJet_AK7_deltaR_lminusak7jet->Fill();
		branch_GroomedJet_AK7_deltaphi_Vak7jet->Fill();
		branch_GroomedJet_numberbjets_csvl_ak7->Fill();
		branch_GroomedJet_numberbjets_csvm_ak7->Fill();
		//branch_GroomedJet_numberbjets_ssvhem_ak7->Fill();
		branch_GroomedJet_numberbjets_csvl_veto_ak7->Fill();
		branch_GroomedJet_numberbjets_csvm_veto_ak7->Fill();
		//branch_GroomedJet_numberbjets_ssvhem_veto_ak7->Fill();
		branch_GroomedJet_numberjets_ak7->Fill();
		branch_GroomedJet_AK7_rcores01->Fill();
		branch_GroomedJet_AK7_rcores02->Fill();
		branch_GroomedJet_AK7_rcores03->Fill();
		branch_GroomedJet_AK7_rcores04->Fill();
		branch_GroomedJet_AK7_rcores05->Fill();
		branch_GroomedJet_AK7_rcores06->Fill();
		branch_GroomedJet_AK7_rcores07->Fill();
		branch_GroomedJet_AK7_rcores08->Fill();
		branch_GroomedJet_AK7_rcores09->Fill();
		branch_GroomedJet_AK7_rcores10->Fill();
		branch_GroomedJet_AK7_rcores11->Fill();

		branch_GroomedJet_AK7_ptcores01->Fill();
		branch_GroomedJet_AK7_ptcores02->Fill();
		branch_GroomedJet_AK7_ptcores03->Fill();
		branch_GroomedJet_AK7_ptcores04->Fill();
		branch_GroomedJet_AK7_ptcores05->Fill();
		branch_GroomedJet_AK7_ptcores06->Fill();
		branch_GroomedJet_AK7_ptcores07->Fill();
		branch_GroomedJet_AK7_ptcores08->Fill();
		branch_GroomedJet_AK7_ptcores09->Fill();
		branch_GroomedJet_AK7_ptcores10->Fill();
		branch_GroomedJet_AK7_ptcores11->Fill();

		branch_GroomedJet_AK7_planarflow01->Fill();
		branch_GroomedJet_AK7_planarflow02->Fill();
		branch_GroomedJet_AK7_planarflow03->Fill();
		branch_GroomedJet_AK7_planarflow04->Fill();
		branch_GroomedJet_AK7_planarflow05->Fill();
		branch_GroomedJet_AK7_planarflow06->Fill();
		branch_GroomedJet_AK7_planarflow07->Fill();
		branch_GroomedJet_AK7_planarflow08->Fill();
		branch_GroomedJet_AK7_planarflow09->Fill();
		branch_GroomedJet_AK7_planarflow10->Fill();
		branch_GroomedJet_AK7_planarflow11->Fill();

		branch_GroomedJet_AK7_mass_sensi_tr->Fill();
		branch_GroomedJet_AK7_mass_sensi_ft->Fill();
		branch_GroomedJet_AK7_mass_sensi_pr->Fill();

		branch_GroomedJet_AK7_qjetmassvolatility->Fill();

		branch_GroomedJet_AK7_prsubjet1ptoverjetpt->Fill();
		branch_GroomedJet_AK7_prsubjet2ptoverjetpt->Fill();
		branch_GroomedJet_AK7_prsubjet1subjet2_deltaR->Fill();

		branch_boostedZ_llj_e_ak7->Fill();
		branch_boostedZ_llj_pt_ak7->Fill();
		branch_boostedZ_llj_eta_ak7->Fill();
		branch_boostedZ_llj_phi_ak7->Fill();
		branch_boostedZ_llj_m_ak7->Fill();
		branch_boostedZ_llj_y_ak7->Fill();

		branch_boostedZ_zjj_ang_ha_ak7->Fill();
		branch_boostedZ_zjj_ang_hb_ak7->Fill();
		branch_boostedZ_zjj_ang_hs_ak7->Fill();
		branch_boostedZ_zjj_ang_phi_ak7->Fill();
		branch_boostedZ_zjj_ang_phia_ak7->Fill();
		branch_boostedZ_zjj_ang_phib_ak7->Fill();
*/
	} // end event loop
	fresults.cd();
	newtree->Write("ZJet",TObject::kOverwrite);
	h_events->Write();
	h_events_weighted->Write();
	delete newtree;
	fresults.Close();
	//fclose(textfile);
	std::cout <<  wda << " Finish :: " << outfilename << "    "<< nentries  << std::endl;
}

double kanamuon::getDeltaPhi(double phi1, double phi2  )
{
	const double PI = 3.14159265;
	double result = phi1 - phi2;

	if(result > PI) {result = result - 2 * PI;}
	if(result <= (-1 * PI)) {result = result + 2 * PI;}
	result = TMath::Abs(result);
	return result;
}

bool kanamuon::doKinematicFit(Int_t                 fflage,
			const TLorentzVector     mup, 
			const TLorentzVector     mum, 
			const TLorentzVector     ajp, 
			const TLorentzVector     bjp, 
			TLorentzVector     & fit_mup, 
			TLorentzVector     & fit_mum,
			TLorentzVector     & fit_ajp, 
			TLorentzVector     & fit_bjp, 
			Float_t            & fit_chi2,
			Int_t              & fit_NDF, 
			Int_t              & fit_status)
{

	bool OK                     = false;
	Resolution* resolution      = new Resolution();

	TMatrixD m1(3,3);
	TMatrixD m2(3,3);
	TMatrixD m3(3,3);
	TMatrixD m4(3,3);
	m1.Zero();
	m2.Zero();
	m3.Zero();
	m4.Zero();

	double etRes, etaRes, phiRes;
	// lepton resolution
	const std::string& leptonName = "muon";  const TLorentzVector lepton   = mup;
	if(leptonName == "electron") {
		OK = resolution->electronResolution(lepton.Et(), lepton.Eta(), etRes, etaRes, phiRes);
		if(!OK) return OK;
	} else {
		OK = resolution->muonResolution(    lepton.Et(), lepton.Eta(), etRes, etaRes, phiRes);
		if(!OK) return OK;
	}
	m1(0,0) = resolution->square(etRes);
	m1(1,1) = resolution->square(etaRes);
	m1(2,2) = resolution->square(phiRes);
	// MET resolution
	OK = resolution->PFMETResolution(     mum.Et(),            etRes, etaRes, phiRes);
	if(!OK) return OK;
	m2(0,0) = resolution->square(etRes);
	m2(1,1) = 0.01; // resolution->square(etaRes)
	m2(2,2) = resolution->square(phiRes);
	// Leading Jet resolution
	OK = resolution->udscPFJetResolution( ajp.Et(), ajp.Eta(), etRes, etaRes, phiRes);
	if(!OK) return OK;
	m3(0,0) = resolution->square(etRes);
	m3(1,1) = resolution->square(etaRes);
	m3(2,2) = resolution->square(phiRes);
	// Leading Jet resolution
	OK = resolution->udscPFJetResolution( bjp.Et(), bjp.Eta(), etRes, etaRes, phiRes);
	if(!OK) return OK;
	m4(0,0) = resolution->square(etRes);
	m4(1,1) = resolution->square(etaRes);
	m4(2,2) = resolution->square(phiRes);

	TLorentzVector tmp_mup = mup;
	TLorentzVector tmp_mum = mum;
	TLorentzVector tmp_ajp = ajp;
	TLorentzVector tmp_bjp = bjp;

	// Fit Particle
	TFitParticleEtEtaPhi* particle1 = new TFitParticleEtEtaPhi( "Lepton",   "Lepton",   &tmp_mup,    &m1 );
	TFitParticleEtEtaPhi* particle2 = new TFitParticleEtEtaPhi( "Neutrino", "Neutrino", &tmp_mum,    &m2 );
	TFitParticleEtEtaPhi* particle3 = new TFitParticleEtEtaPhi( "Jeta",     "Jeta",     &tmp_ajp,    &m3 );
	TFitParticleEtEtaPhi* particle4 = new TFitParticleEtEtaPhi( "Jetb",     "Jetb",     &tmp_bjp,    &m4 );

	// Constraint
	TFitConstraintMGaus* mCons1 = new TFitConstraintMGaus( "W1MassConstraint", "W1Mass-Constraint", 0, 0 , 80.399, 2.085);
	//TFitConstraintM *mCons1 = new TFitConstraintM( "WMassConstrainta", "WMass-Constrainta", 0, 0 , 80.4);
	mCons1->addParticles1( particle1, particle2 );

	TFitConstraintMGaus* mCons2 = new TFitConstraintMGaus( "W2MassConstraint", "W2Mass-Constraint", 0, 0 , 80.399, 2.085);
	//TFitConstraintM *mCons2 = new TFitConstraintM( "WMassConstraintb", "WMass-Constraintb", 0, 0 , 80.4);
	mCons2->addParticles1( particle3, particle4 );

	TFitConstraintEp *pxCons = new TFitConstraintEp( "PxConstraint", "Px-Constraint", 0, TFitConstraintEp::pX , (mup+mum+ajp+bjp).Px() );
	pxCons->addParticles( particle1, particle2, particle3, particle4 );

	TFitConstraintEp *pyCons = new TFitConstraintEp( "PyConstraint", "Py-Constraint", 0, TFitConstraintEp::pY , (mup+mum+ajp+bjp).Py() );
	pyCons->addParticles( particle1, particle2, particle3, particle4 );

	//Definition of the fitter
	TKinFitter* fitter = new TKinFitter("fitter", "fitter");
	if        (fflage == 1 ){
		fitter->addMeasParticle( particle1 );
		fitter->addMeasParticle( particle2 );
		fitter->addMeasParticle( particle3 );
		fitter->addMeasParticle( particle4 );
		fitter->addConstraint( mCons1 );
		fitter->addConstraint( mCons2 );
	}else   if(fflage == 2 ){
		fitter->addMeasParticle( particle1 );
		fitter->addMeasParticle( particle2 );
		fitter->addMeasParticle( particle3 );
		fitter->addMeasParticle( particle4 );
		fitter->addConstraint( pxCons );
		fitter->addConstraint( pyCons );
		fitter->addConstraint( mCons1 );
		fitter->addConstraint( mCons2 );
	}else   if(fflage == 3 ){
		fitter->addMeasParticle( particle3 );
		fitter->addMeasParticle( particle4 );
		fitter->addConstraint( mCons2 );
	}else {return false;}

	//Set convergence criteria
	fitter->setMaxNbIter( 50 );
	fitter->setMaxDeltaS( 1e-2 );
	fitter->setMaxF( 1e-1 );
	fitter->setVerbosity(1);
	fitter->fit();

	//Return the kinematic fit results
	fit_status   = fitter->getStatus();
	fit_chi2     = fitter->getS();
	fit_NDF      = fitter->getNDF();
	fit_mup      = *(particle1->getCurr4Vec()); 
	fit_mum      = *(particle2->getCurr4Vec()); 
	fit_ajp      = *(particle3->getCurr4Vec()); 
	fit_bjp      = *(particle4->getCurr4Vec()); 

	if(fitter->getStatus() == 0) { OK = true;  } else { OK = false;  }
	delete resolution;
	delete particle1;
	delete particle2;
	delete particle3;
	delete particle4;
	delete mCons1;
	delete mCons2;
	delete pxCons;
	delete pyCons;
	delete fitter;

	return OK;
}

void kanamuon::calculateAngles(TLorentzVector& thep4M11, TLorentzVector& thep4M12, TLorentzVector& thep4M21, TLorentzVector& thep4M22, double& costheta1, double& costheta2, double& phi, double& costhetastar, double& phistar1, double& phistar2){


	TLorentzVector thep4H = thep4M11 + thep4M12 + thep4M21 + thep4M22;
	TLorentzVector thep4Z1 = thep4M11 + thep4M12;
	TLorentzVector thep4Z2 = thep4M21 + thep4M22;

	double norm;

	TVector3 boostX = -(thep4H.BoostVector());
	TLorentzVector thep4Z1inXFrame( thep4Z1 );
	TLorentzVector thep4Z2inXFrame( thep4Z2 );      
	thep4Z1inXFrame.Boost( boostX );
	thep4Z2inXFrame.Boost( boostX );
	TVector3 theZ1X_p3 = TVector3( thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z() );
	TVector3 theZ2X_p3 = TVector3( thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z() );

	// calculate phi1, phi2, costhetastar
	///phi1 = theZ1X_p3.Phi();
	///phi2 = theZ2X_p3.Phi();

	///////////////////////////////////////////////
	// check for z1/z2 convention, redefine all 4 vectors with convention
	/////////////////////////////////////////////// 
	TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2, p4M21, p4M22;
	p4Z1 = thep4Z1; p4M11 = thep4M11; p4M12 = thep4M12;
	p4Z2 = thep4Z2; p4M21 = thep4M21; p4M22 = thep4M22;
	costhetastar = theZ1X_p3.CosTheta();

	// now helicity angles................................
	// ...................................................
	TVector3 boostZ1 = -(p4Z1.BoostVector());
	TLorentzVector p4Z2Z1(p4Z2);
	p4Z2Z1.Boost(boostZ1);
	//find the decay axis
	/////TVector3 unitx_1 = -Hep3Vector(p4Z2Z1);
	TVector3 unitx_1( -p4Z2Z1.X(), -p4Z2Z1.Y(), -p4Z2Z1.Z() );
	norm = 1/(unitx_1.Mag());
	unitx_1*=norm;
	//boost daughters of z2
	TLorentzVector p4M21Z1(p4M21);
	TLorentzVector p4M22Z1(p4M22);
	p4M21Z1.Boost(boostZ1);
	p4M22Z1.Boost(boostZ1);
	//create z and y axes
	/////TVector3 unitz_1 = Hep3Vector(p4M21Z1).cross(Hep3Vector(p4M22Z1));
	TVector3 p4M21Z1_p3( p4M21Z1.X(), p4M21Z1.Y(), p4M21Z1.Z() );
	TVector3 p4M22Z1_p3( p4M22Z1.X(), p4M22Z1.Y(), p4M22Z1.Z() );
	TVector3 unitz_1 = p4M21Z1_p3.Cross( p4M22Z1_p3 );
	norm = 1/(unitz_1.Mag());
	unitz_1 *= norm;
	TVector3 unity_1 = unitz_1.Cross(unitx_1);

	//caculate theta1
	TLorentzVector p4M11Z1(p4M11);
	p4M11Z1.Boost(boostZ1);
	TVector3 p3M11( p4M11Z1.X(), p4M11Z1.Y(), p4M11Z1.Z() );
	TVector3 unitM11 = p3M11.Unit();
	double x_m11 = unitM11.Dot(unitx_1); double y_m11 = unitM11.Dot(unity_1); double z_m11 = unitM11.Dot(unitz_1);
	TVector3 M11_Z1frame(y_m11, z_m11, x_m11);
	costheta1 = M11_Z1frame.CosTheta();
	//std::cout << "theta1: " << M11_Z1frame.Theta() << std::endl;
	//////-----------------------old way of calculating phi---------------/////////
	phi = M11_Z1frame.Phi();

	//set axes for other system
	TVector3 boostZ2 = -(p4Z2.BoostVector());
	TLorentzVector p4Z1Z2(p4Z1);
	p4Z1Z2.Boost(boostZ2);
	TVector3 unitx_2( -p4Z1Z2.X(), -p4Z1Z2.Y(), -p4Z1Z2.Z() );
	norm = 1/(unitx_2.Mag());
	unitx_2*=norm;
	//boost daughters of z2
	TLorentzVector p4M11Z2(p4M11);
	TLorentzVector p4M12Z2(p4M12);
	p4M11Z2.Boost(boostZ2);
	p4M12Z2.Boost(boostZ2);
	TVector3 p4M11Z2_p3( p4M11Z2.X(), p4M11Z2.Y(), p4M11Z2.Z() );
	TVector3 p4M12Z2_p3( p4M12Z2.X(), p4M12Z2.Y(), p4M12Z2.Z() );
	TVector3 unitz_2 = p4M11Z2_p3.Cross( p4M12Z2_p3 );
	norm = 1/(unitz_2.Mag());
	unitz_2*=norm;
	TVector3 unity_2 = unitz_2.Cross(unitx_2);
	//calcuate theta2
	TLorentzVector p4M21Z2(p4M21);
	p4M21Z2.Boost(boostZ2);
	TVector3 p3M21( p4M21Z2.X(), p4M21Z2.Y(), p4M21Z2.Z() );
	TVector3 unitM21 = p3M21.Unit();
	double x_m21 = unitM21.Dot(unitx_2); double y_m21 = unitM21.Dot(unity_2); double z_m21 = unitM21.Dot(unitz_2);
	TVector3 M21_Z2frame(y_m21, z_m21, x_m21);
	costheta2 = M21_Z2frame.CosTheta();

	// calculate phi
	//calculating phi_n
	TLorentzVector n_p4Z1inXFrame( p4Z1 );
	TLorentzVector n_p4M11inXFrame( p4M11 );
	n_p4Z1inXFrame.Boost( boostX );
	n_p4M11inXFrame.Boost( boostX );        
	TVector3 n_p4Z1inXFrame_unit = n_p4Z1inXFrame.Vect().Unit();
	TVector3 n_p4M11inXFrame_unit = n_p4M11inXFrame.Vect().Unit();  
	TVector3 n_unitz_1( n_p4Z1inXFrame_unit );
	//// y-axis is defined by neg lepton cross z-axis
	//// the subtle part is here...
	//////////TVector3 n_unity_1 = n_p4M11inXFrame_unit.Cross( n_unitz_1 );
	TVector3 n_unity_1 = n_unitz_1.Cross( n_p4M11inXFrame_unit );
	TVector3 n_unitx_1 = n_unity_1.Cross( n_unitz_1 );

	TLorentzVector n_p4M21inXFrame( p4M21 );
	n_p4M21inXFrame.Boost( boostX );
	TVector3 n_p4M21inXFrame_unit = n_p4M21inXFrame.Vect().Unit();
	//rotate into other plane
	TVector3 n_p4M21inXFrame_unitprime( n_p4M21inXFrame_unit.Dot(n_unitx_1), n_p4M21inXFrame_unit.Dot(n_unity_1), n_p4M21inXFrame_unit.Dot(n_unitz_1) );

	///////-----------------new way of calculating phi-----------------///////
	//double phi_n =  n_p4M21inXFrame_unitprime.Phi();
	/*
	   std::cout << "---------------------------" << std::endl;
	   std::cout << "phi: " << phi << std::endl;
	   std::cout << "phi_n: " << phi_n << std::endl;
	   std::cout << "phi + phi_n: " << (phi+phi_n) << std::endl;
	   */
	/// and then calculate phistar1
	TVector3 n_p4PartoninXFrame_unit( 0.0, 0.0, 1.0 );
	TVector3 n_p4PartoninXFrame_unitprime( n_p4PartoninXFrame_unit.Dot(n_unitx_1), n_p4PartoninXFrame_unit.Dot(n_unity_1), n_p4PartoninXFrame_unit.Dot(n_unitz_1) );
	// negative sign is for arrow convention in paper
	phistar1 = (n_p4PartoninXFrame_unitprime.Phi());

	// and the calculate phistar2
	TLorentzVector n_p4Z2inXFrame( p4Z2 );
	n_p4Z2inXFrame.Boost( boostX );
	TVector3 n_p4Z2inXFrame_unit = n_p4Z2inXFrame.Vect().Unit();
	///////TLorentzVector n_p4M21inXFrame( p4M21 );
	//////n_p4M21inXFrame.Boost( boostX );        
	////TVector3 n_p4M21inXFrame_unit = n_p4M21inXFrame.Vect().Unit();  
	TVector3 n_unitz_2( n_p4Z2inXFrame_unit );
	//// y-axis is defined by neg lepton cross z-axis
	//// the subtle part is here...
	//////TVector3 n_unity_2 = n_p4M21inXFrame_unit.Cross( n_unitz_2 );
	TVector3 n_unity_2 = n_unitz_2.Cross( n_p4M21inXFrame_unit );
	TVector3 n_unitx_2 = n_unity_2.Cross( n_unitz_2 );
	TVector3 n_p4PartoninZ2PlaneFrame_unitprime( n_p4PartoninXFrame_unit.Dot(n_unitx_2), n_p4PartoninXFrame_unit.Dot(n_unity_2), n_p4PartoninXFrame_unit.Dot(n_unitz_2) );
	phistar2 = (n_p4PartoninZ2PlaneFrame_unitprime.Phi());

	/*
	   double phistar12_0 = phistar1 + phistar2;
	   if (phistar12_0 > TMath::Pi()) phistar12 = phistar12_0 - 2*TMath::Pi();
	   else if (phistar12_0 < (-1.)*TMath::Pi()) phistar12 = phistar12_0 + 2*TMath::Pi();
	   else phistar12 = phistar12_0;
	   */

}


// function used to fill the counters with preselction level cuts
void kanamuon::InitCounters( const char* input_file_name, TH1F* h_events, TH1F* h_events_weighted)
{
	TFile* f = new TFile(input_file_name, "READ");
	std::vector<float> events;

	//get the counters from the FNAL NT
	events.push_back(((TH1F*) f->Get("AllEventsStep/totalEvents"))->GetEntries());
	events.push_back(((TH1F*) f->Get("noscrapingStep/totalEvents"))->GetEntries());
	events.push_back(((TH1F*) f->Get("HBHENoiseStep/totalEvents"))->GetEntries());
	events.push_back(((TH1F*) f->Get("primaryVertexStep/totalEvents"))->GetEntries());
	events.push_back(((TH1F*) f->Get("tightLeptonStep/totalEvents"))->GetEntries());
	events.push_back(((TH1F*) f->Get("looseElectronStep/totalEvents"))->GetEntries());
	events.push_back(((TH1F*) f->Get("looseMuonStep/totalEvents"))->GetEntries());
	events.push_back(((TH1F*) f->Get("RequireTwoJetsORboostedVStep/totalEvents"))->GetEntries());


	//put the counters in the counter histos
	for ( unsigned int istep = 0; istep < events.size(); istep++ ) {
		h_events -> SetBinContent( istep + 1, events[istep] );
		h_events_weighted -> SetBinContent( istep + 1, events[istep] );
	}
	f -> Close();
}
