#ifndef kanamuon_ZJets_h
#define kanamuon_ZJets_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include "iostream"
using namespace std;

#include "TLorentzVector.h"

class kanamuon {
	public :
		TTree          *fChain;   //!pointer to the analyzed TTree or TChain
		Int_t           fCurrent; //!current Tree number in a TChain

		// Declaration of leaf types
		Int_t           numPFCorJets;
		Int_t           numPFCorJetBTags;
		Float_t         JetPFCor_Et[8];
		Float_t         JetPFCor_Pt[8];
		Float_t         JetPFCor_Eta[8];
		Float_t         JetPFCor_Phi[8];
		Float_t         JetPFCor_Theta[8];
		Float_t         JetPFCor_Px[8];
		Float_t         JetPFCor_Py[8];
		Float_t         JetPFCor_Pz[8];
		Float_t         JetPFCor_E[8];
		Float_t         JetPFCor_Y[8];
		Float_t         JetPFCor_Mass[8];
		Float_t         JetPFCor_etaetaMoment[8];
		Float_t         JetPFCor_phiphiMoment[8];
		Float_t         JetPFCor_etaphiMoment[8];
		Float_t         JetPFCor_maxDistance[8];
		Int_t           JetPFCor_nConstituents[8];
		Float_t         JetPFCor_Area[8];
		Float_t         VplusPFCorJet_Mass[8];
		Float_t         JetPFCor_dphiBoson[8];
		Float_t         JetPFCor_detaBoson[8];
		Float_t         JetPFCor_dRBoson[8];
		Float_t         JetPFCor_dphiMET[8];
		Float_t         JetPFCor_bDiscriminator[8];
		Float_t         JetPFCor_bDiscriminatorSSVHE[8];
		Float_t         JetPFCor_bDiscriminatorTCHE[8];
		Float_t         JetPFCor_bDiscriminatorCSV[8];
		Float_t         JetPFCor_bDiscriminatorJP[8];
		Float_t         JetPFCor_bDiscriminatorSSVHP[8];
		Float_t         JetPFCor_bDiscriminatorTCHP[8];
		Float_t         JetPFCor_secVertexMass[8];
		Float_t         JetPFCor_ChargedHadronEnergy[8];
		Float_t         JetPFCor_ChargedHadronEnergyFrac[8];
		Float_t         JetPFCor_NeutralHadronEnergy[8];
		Float_t         JetPFCor_NeutralHadronEnergyFrac[8];
		Float_t         JetPFCor_ChargedEmEnergy[8];
		Float_t         JetPFCor_ChargedEmEnergyFrac[8];
		Float_t         JetPFCor_ChargedMuEnergy[8];
		Float_t         JetPFCor_ChargedMuEnergyFrac[8];
		Float_t         JetPFCor_NeutralEmEnergy[8];
		Float_t         JetPFCor_NeutralEmEnergyFrac[8];
		Int_t           JetPFCor_ChargedMultiplicity[8];
		Int_t           JetPFCor_NeutralMultiplicity[8];
		Int_t           JetPFCor_MuonMultiplicity[8];
		Float_t         JetPFCor_PhotonEnergy[8];
		Float_t         JetPFCor_PhotonEnergyFraction[8];
		Float_t         JetPFCor_ElectronEnergy[8];
		Float_t         JetPFCor_ElectronEnergyFraction[8];
		Float_t         JetPFCor_MuonEnergy[8];
		Float_t         JetPFCor_MuonEnergyFraction[8];
		Float_t         JetPFCor_HFHadronEnergy[8];
		Float_t         JetPFCor_HFHadronEnergyFraction[8];
		Float_t         JetPFCor_HFEMEnergy[8];
		Float_t         JetPFCor_HFEMEnergyFraction[8];
		Int_t           JetPFCor_ChargedHadronMultiplicity[8];
		Int_t           JetPFCor_NeutralHadronMultiplicity[8];
		Int_t           JetPFCor_PhotonMultiplicity[8];
		Int_t           JetPFCor_ElectronMultiplicity[8];
		Int_t           JetPFCor_HFHadronMultiplicity[8];
		Int_t           JetPFCor_HFEMMultiplicity[8];
		Float_t         JetPFCor_SumPtCands[8];
		Float_t         JetPFCor_SumPt2Cands[8];
		Float_t         JetPFCor_rmsCands[8];
		Float_t         JetPFCor_PtD[8];
		Float_t         JetPFCor_QGLikelihood[8];
		Float_t         MassV2j_PFCor_MVAMET;
		Float_t         MassV2j_PFCor;
		Float_t         MassV3j_PFCor;
		Float_t         MassV4j_PFCor;
		Float_t         MassV5j_PFCor;
		Float_t         MassV6j_PFCor;
		Float_t         Mass2j_PFCor;
		Float_t         Mass3j_PFCor;
		Float_t         Mass4j_PFCor;
		Float_t         Mass5j_PFCor;
		Float_t         Mass6j_PFCor;
		Float_t         cosJacksonAngleV2j_PFCor;
		Float_t         cosJacksonAngle2j_PFCor;
		Float_t         cosJacksonAngleV3j_PFCor;
		Float_t         cosJacksonAngle3j12_PFCor;
		Float_t         cosJacksonAngle3j23_PFCor;
		Float_t         cosJacksonAngle3j31_PFCor;
		Float_t         cosphiDecayPlane_PFCor;
		Float_t         cosThetaLnu_PFCor;
		Float_t         cosThetaJJ_PFCor;
		Float_t         colorCorrPull01PFCor;
		Float_t         colorCorrPull02PFCor;
		Float_t         colorCorrPull12PFCor;
		Float_t         colorCorrPull03PFCor;
		Float_t         colorCorrPull13PFCor;
		Float_t         colorCorrPull23PFCor;
		Float_t         colorCorrPull04PFCor;
		Float_t         colorCorrPull14PFCor;
		Float_t         colorCorrPull24PFCor;
		Float_t         colorCorrPull34PFCor;
		Float_t         colorCorrPull05PFCor;
		Float_t         colorCorrPull15PFCor;
		Float_t         colorCorrPull25PFCor;
		Float_t         colorCorrPull35PFCor;
		Float_t         colorCorrPull45PFCor;
		Float_t         cosThetaJ1HiggsCM_PFCor;
		Float_t         cosThetaJ2HiggsCM_PFCor;
		Float_t         cosThetaL1HiggsCM_PFCor;
		Float_t         cosThetaL2HiggsCM_PFCor;
		Float_t         cosThetaV1HiggsCM_PFCor;
		Float_t         cosThetaV2HiggsCM_PFCor;
		Bool_t          JetPFCor_isPileUpJetLoose[8];
		Bool_t          JetPFCor_isPileUpJetMedium[8];
		Bool_t          JetPFCor_isPileUpJetTight[8];
		Float_t         GroomedJet_CA8_pt_uncorr[6];
		Float_t         GroomedJet_CA8_mass_uncorr[6];
		Float_t         GroomedJet_CA8_mass_tr_uncorr[6];
		Float_t         GroomedJet_CA8_mass_ft_uncorr[6];
		Float_t         GroomedJet_CA8_mass_pr_uncorr[6];
		Float_t         GroomedJet_CA8_tau2tau1[6];
		Float_t         GroomedJet_CA8_tau1[6];
		Float_t         GroomedJet_CA8_tau2[6];
		Float_t         GroomedJet_CA8_tau3[6];
		Float_t         GroomedJet_CA8_tau4[6];
		Float_t         GroomedJet_CA8_massdrop_pr_uncorr[6];
		Float_t         GroomedJet_CA8_pt[6];
		Float_t         GroomedJet_CA8_eta[6];
		Float_t         GroomedJet_CA8_phi[6];
		Float_t         GroomedJet_CA8_e[6];
		Float_t         GroomedJet_CA8_pt_tr_uncorr[6];
		Float_t         GroomedJet_CA8_pt_tr[6];
		Float_t         GroomedJet_CA8_eta_tr[6];
		Float_t         GroomedJet_CA8_phi_tr[6];
		Float_t         GroomedJet_CA8_e_tr[6];
		Float_t         GroomedJet_CA8_pt_ft_uncorr[6];
		Float_t         GroomedJet_CA8_pt_ft[6];
		Float_t         GroomedJet_CA8_eta_ft[6];
		Float_t         GroomedJet_CA8_phi_ft[6];
		Float_t         GroomedJet_CA8_e_ft[6];
		Float_t         GroomedJet_CA8_pt_pr_uncorr[6];
		Float_t         GroomedJet_CA8_pt_pr[6];
		Float_t         GroomedJet_CA8_eta_pr[6];
		Float_t         GroomedJet_CA8_phi_pr[6];
		Float_t         GroomedJet_CA8_e_pr[6];
		Float_t         GroomedJet_CA8_prsubjet1_px[6];
		Float_t         GroomedJet_CA8_prsubjet1_py[6];
		Float_t         GroomedJet_CA8_prsubjet1_pz[6];
		Float_t         GroomedJet_CA8_prsubjet1_e[6];
		Float_t         GroomedJet_CA8_prsubjet2_px[6];
		Float_t         GroomedJet_CA8_prsubjet2_py[6];
		Float_t         GroomedJet_CA8_prsubjet2_pz[6];
		Float_t         GroomedJet_CA8_prsubjet2_e[6];
		Float_t         GroomedJet_CA8_mass[6];
		Float_t         GroomedJet_CA8_mass_tr[6];
		Float_t         GroomedJet_CA8_mass_ft[6];
		Float_t         GroomedJet_CA8_mass_pr[6];
		Float_t         GroomedJet_CA8_massdrop_pr[6];
		Float_t         GroomedJet_CA8_area[6];
		Float_t         GroomedJet_CA8_area_tr[6];
		Float_t         GroomedJet_CA8_area_ft[6];
		Float_t         GroomedJet_CA8_area_pr[6];
		Float_t         GroomedJet_CA8_jetconstituents[6];
		Float_t         GroomedJet_CA8_jetcharge[6];
		Float_t         GroomedJet_CA8_rcores[11][6];
		Float_t         GroomedJet_CA8_ptcores[11][6];
		Float_t         GroomedJet_CA8_planarflow[11][6];
		Float_t         GroomedJet_CA8_qjetmass[50];
		Float_t         GroomedJet_CA8_qjetmassdrop[50];
		Float_t         GroomedJet_CA8_constituents0_eta[100];
		Float_t         GroomedJet_CA8_constituents0_phi[100];
		Float_t         GroomedJet_CA8_constituents0_e[100];
		Int_t           GroomedJet_CA8_nconstituents0;
		Float_t         GroomedJet_CA8_constituents0pr_eta[100];
		Float_t         GroomedJet_CA8_constituents0pr_phi[100];
		Float_t         GroomedJet_CA8_constituents0pr_e[100];
		Int_t           GroomedJet_CA8_nconstituents0pr;
		Float_t         GenGroomedJet_CA8_pt_uncorr[6];
		Float_t         GenGroomedJet_CA8_mass_uncorr[6];
		Float_t         GenGroomedJet_CA8_mass_tr_uncorr[6];
		Float_t         GenGroomedJet_CA8_mass_ft_uncorr[6];
		Float_t         GenGroomedJet_CA8_mass_pr_uncorr[6];
		Float_t         GenGroomedJet_CA8_tau2tau1[6];
		Float_t         GenGroomedJet_CA8_tau1[6];
		Float_t         GenGroomedJet_CA8_tau2[6];
		Float_t         GenGroomedJet_CA8_tau3[6];
		Float_t         GenGroomedJet_CA8_tau4[6];
		Float_t         GenGroomedJet_CA8_massdrop_pr_uncorr[6];
		Float_t         GenGroomedJet_CA8_pt[6];
		Float_t         GenGroomedJet_CA8_eta[6];
		Float_t         GenGroomedJet_CA8_phi[6];
		Float_t         GenGroomedJet_CA8_e[6];
		Float_t         GenGroomedJet_CA8_pt_tr_uncorr[6];
		Float_t         GenGroomedJet_CA8_pt_tr[6];
		Float_t         GenGroomedJet_CA8_eta_tr[6];
		Float_t         GenGroomedJet_CA8_phi_tr[6];
		Float_t         GenGroomedJet_CA8_e_tr[6];
		Float_t         GenGroomedJet_CA8_pt_ft_uncorr[6];
		Float_t         GenGroomedJet_CA8_pt_ft[6];
		Float_t         GenGroomedJet_CA8_eta_ft[6];
		Float_t         GenGroomedJet_CA8_phi_ft[6];
		Float_t         GenGroomedJet_CA8_e_ft[6];
		Float_t         GenGroomedJet_CA8_pt_pr_uncorr[6];
		Float_t         GenGroomedJet_CA8_pt_pr[6];
		Float_t         GenGroomedJet_CA8_eta_pr[6];
		Float_t         GenGroomedJet_CA8_phi_pr[6];
		Float_t         GenGroomedJet_CA8_e_pr[6];
		Float_t         GenGroomedJet_CA8_prsubjet1_px[6];
		Float_t         GenGroomedJet_CA8_prsubjet1_py[6];
		Float_t         GenGroomedJet_CA8_prsubjet1_pz[6];
		Float_t         GenGroomedJet_CA8_prsubjet1_e[6];
		Float_t         GenGroomedJet_CA8_prsubjet2_px[6];
		Float_t         GenGroomedJet_CA8_prsubjet2_py[6];
		Float_t         GenGroomedJet_CA8_prsubjet2_pz[6];
		Float_t         GenGroomedJet_CA8_prsubjet2_e[6];
		Float_t         GenGroomedJet_CA8_mass[6];
		Float_t         GenGroomedJet_CA8_mass_tr[6];
		Float_t         GenGroomedJet_CA8_mass_ft[6];
		Float_t         GenGroomedJet_CA8_mass_pr[6];
		Float_t         GenGroomedJet_CA8_massdrop_pr[6];
		Float_t         GenGroomedJet_CA8_area[6];
		Float_t         GenGroomedJet_CA8_area_tr[6];
		Float_t         GenGroomedJet_CA8_area_ft[6];
		Float_t         GenGroomedJet_CA8_area_pr[6];
		Float_t         GenGroomedJet_CA8_jetconstituents[6];
		Float_t         GenGroomedJet_CA8_jetcharge[6];
		Float_t         GenGroomedJet_CA8_rcores[11][6];
		Float_t         GenGroomedJet_CA8_ptcores[11][6];
		Float_t         GenGroomedJet_CA8_planarflow[11][6];
		Float_t         GenGroomedJet_CA8_qjetmass[50];
		Float_t         GenGroomedJet_CA8_qjetmassdrop[50];
		Float_t         GenGroomedJet_CA8_constituents0_eta[100];
		Float_t         GenGroomedJet_CA8_constituents0_phi[100];
		Float_t         GenGroomedJet_CA8_constituents0_e[100];
		Int_t           GenGroomedJet_CA8_nconstituents0;
		Float_t         GenGroomedJet_CA8_constituents0pr_eta[100];
		Float_t         GenGroomedJet_CA8_constituents0pr_phi[100];
		Float_t         GenGroomedJet_CA8_constituents0pr_e[100];
		Int_t           GenGroomedJet_CA8_nconstituents0pr;
		Int_t           numGenJets;
		Int_t           numGenJetBTags;
		Float_t         JetGen_Et[8];
		Float_t         JetGen_Pt[8];
		Float_t         JetGen_Eta[8];
		Float_t         JetGen_Phi[8];
		Float_t         JetGen_Theta[8];
		Float_t         JetGen_Px[8];
		Float_t         JetGen_Py[8];
		Float_t         JetGen_Pz[8];
		Float_t         JetGen_E[8];
		Float_t         JetGen_Y[8];
		Float_t         JetGen_Mass[8];
		Float_t         JetGen_etaetaMoment[8];
		Float_t         JetGen_phiphiMoment[8];
		Float_t         JetGen_etaphiMoment[8];
		Float_t         JetGen_maxDistance[8];
		Int_t           JetGen_nConstituents[8];
		Float_t         JetGen_Area[8];
		Float_t         VplusGenJet_Mass[8];
		Float_t         JetGen_dphiBoson[8];
		Float_t         JetGen_detaBoson[8];
		Float_t         JetGen_dRBoson[8];
		Float_t         JetGen_dphiMET[8];
		Float_t         JetGen_bDiscriminator[8];
		Float_t         JetGen_bDiscriminatorSSVHE[8];
		Float_t         JetGen_bDiscriminatorTCHE[8];
		Float_t         JetGen_bDiscriminatorCSV[8];
		Float_t         JetGen_bDiscriminatorJP[8];
		Float_t         JetGen_bDiscriminatorSSVHP[8];
		Float_t         JetGen_bDiscriminatorTCHP[8];
		Float_t         JetGen_secVertexMass[8];
		Float_t         MassV2j_Gen_MVAMET;
		Float_t         MassV2j_Gen;
		Float_t         MassV3j_Gen;
		Float_t         MassV4j_Gen;
		Float_t         MassV5j_Gen;
		Float_t         MassV6j_Gen;
		Float_t         Mass2j_Gen;
		Float_t         Mass3j_Gen;
		Float_t         Mass4j_Gen;
		Float_t         Mass5j_Gen;
		Float_t         Mass6j_Gen;
		Float_t         cosJacksonAngleV2j_Gen;
		Float_t         cosJacksonAngle2j_Gen;
		Float_t         cosJacksonAngleV3j_Gen;
		Float_t         cosJacksonAngle3j12_Gen;
		Float_t         cosJacksonAngle3j23_Gen;
		Float_t         cosJacksonAngle3j31_Gen;
		Float_t         cosphiDecayPlane_Gen;
		Float_t         cosThetaLnu_Gen;
		Float_t         cosThetaJJ_Gen;
		Int_t           NumPhotons;
		Float_t         Photon_Et[5];   //[NumPhotons]
		Float_t         Photon_E[5];   //[NumPhotons]
		Float_t         Photon_Eta[5];   //[NumPhotons]
		Float_t         Photon_Phi[5];   //[NumPhotons]
		Float_t         Photon_Theta[5];   //[NumPhotons]
		Float_t         Photon_Px[5];   //[NumPhotons]
		Float_t         Photon_Py[5];   //[NumPhotons]
		Float_t         Photon_Pz[5];   //[NumPhotons]
		Float_t         Photon_Vx[5];   //[NumPhotons]
		Float_t         Photon_Vy[5];   //[NumPhotons]
		Float_t         Photon_Vz[5];   //[NumPhotons]
		Float_t         Photon_SC_Et[5];   //[NumPhotons]
		Float_t         Photon_SC_E[5];   //[NumPhotons]
		Float_t         Photon_SC_Eta[5];   //[NumPhotons]
		Float_t         Photon_SC_Phi[5];   //[NumPhotons]
		Float_t         Photon_SC_Theta[5];   //[NumPhotons]
		Float_t         Photon_SC_x[5];   //[NumPhotons]
		Float_t         Photon_SC_y[5];   //[NumPhotons]
		Float_t         Photon_SC_z[5];   //[NumPhotons]
		Float_t         PFisocharged03[5];   //[NumPhotons]
		Float_t         PFisophoton03[5];   //[NumPhotons]
		Float_t         PFisoneutral03[5];   //[NumPhotons]
		Float_t         trkSumPtHollowConeDR04_Photon11[5];   //[NumPhotons]
		Float_t         ecalRecHitSumEtConeDR04_Photon11[5];   //[NumPhotons]
		Float_t         hcalTowerSumEtConeDR04_Photon11[5];   //[NumPhotons]
		Float_t         Photon_HoverE[5];   //[NumPhotons]
		Float_t         Photon_HoverE2011[5];   //[NumPhotons]
		Float_t         Photon_SigmaIetaIeta[5];   //[NumPhotons]
		Int_t           Photon_hasPixelSeed[5];   //[NumPhotons]
		Int_t           Photon_passElecVeto[5];   //[NumPhotons]
		Int_t           Photon_Id2011[5];   //[NumPhotons]
		Int_t           Photon_Id2012[5];   //[NumPhotons]
		Float_t         Z_mass;
		Float_t         Z_mt;
		Float_t         Z_mtMVA;
		Float_t         Z_px;
		Float_t         Z_py;
		Float_t         Z_pz;
		Float_t         Z_e;
		Float_t         Z_pt;
		Float_t         Z_et;
		Float_t         Z_eta;
		Float_t         Z_phi;
		Float_t         Z_vx;
		Float_t         Z_vy;
		Float_t         Z_vz;
		Float_t         Z_y;
		Float_t         Z_muplus_px;
		Float_t         Z_muplus_py;
		Float_t         Z_muplus_pz;
		Float_t         Z_muplus_e;
		Float_t         Z_muplus_pt;
		Float_t         Z_muplus_et;
		Float_t         Z_muplus_eta;
		Float_t         Z_muplus_theta;
		Float_t         Z_muplus_phi;
		Int_t           Z_muplus_charge;
		Float_t         Z_muplus_vx;
		Float_t         Z_muplus_vy;
		Float_t         Z_muplus_vz;
		Float_t         Z_muplus_y;
		Float_t         Z_muplus_trackiso;
		Float_t         Z_muplus_hcaliso;
		Float_t         Z_muplus_ecaliso;
		Int_t           Z_muplus_type;
		Int_t           Z_muplus_numberOfChambers;
		Int_t           Z_muplus_numberOfMatches;
		Float_t         Z_muplus_d0bsp;
		Float_t         Z_muplus_dz000;
		Float_t         Z_muplus_dzPV;
		Float_t         Z_muplus_pfiso_sumChargedHadronPt;
		Float_t         Z_muplus_pfiso_sumChargedParticlePt;
		Float_t         Z_muplus_pfiso_sumNeutralHadronEt;
		Float_t         Z_muplus_pfiso_sumPhotonEt;
		Float_t         Z_muplus_pfiso_sumPUPt;
		Float_t         Z_muminus_px;
		Float_t         Z_muminus_py;
		Float_t         Z_muminus_pz;
		Float_t         Z_muminus_e;
		Float_t         Z_muminus_pt;
		Float_t         Z_muminus_et;
		Float_t         Z_muminus_eta;
		Float_t         Z_muminus_theta;
		Float_t         Z_muminus_phi;
		Int_t           Z_muminus_charge;
		Float_t         Z_muminus_vx;
		Float_t         Z_muminus_vy;
		Float_t         Z_muminus_vz;
		Float_t         Z_muminus_y;
		Float_t         Z_muminus_trackiso;
		Float_t         Z_muminus_hcaliso;
		Float_t         Z_muminus_ecaliso;
		Int_t           Z_muminus_type;
		Int_t           Z_muminus_numberOfChambers;
		Int_t           Z_muminus_numberOfMatches;
		Float_t         Z_muminus_d0bsp;
		Float_t         Z_muminus_dz000;
		Float_t         Z_muminus_dzPV;
		Float_t         Z_muminus_pfiso_sumChargedHadronPt;
		Float_t         Z_muminus_pfiso_sumChargedParticlePt;
		Float_t         Z_muminus_pfiso_sumNeutralHadronEt;
		Float_t         Z_muminus_pfiso_sumPhotonEt;
		Float_t         Z_muminus_pfiso_sumPUPt;
		Float_t         Z_Photon_pt_gen;
		Float_t         Z_muplus_px_gen;
		Float_t         Z_muplus_py_gen;
		Float_t         Z_muplus_pz_gen;
		Float_t         Z_muplus_e_gen;
		Float_t         Z_muplus_pt_gen;
		Float_t         Z_muplus_et_gen;
		Float_t         Z_muplus_eta_gen;
		Float_t         Z_muplus_theta_gen;
		Float_t         Z_muplus_phi_gen;
		Int_t           Z_muplus_charge_gen;
		Float_t         Z_muplus_vx_gen;
		Float_t         Z_muplus_vy_gen;
		Float_t         Z_muplus_vz_gen;
		Float_t         Z_muplus_y_gen;
		Float_t         Z_muminus_px_gen;
		Float_t         Z_muminus_py_gen;
		Float_t         Z_muminus_pz_gen;
		Float_t         Z_muminus_e_gen;
		Float_t         Z_muminus_pt_gen;
		Float_t         Z_muminus_et_gen;
		Float_t         Z_muminus_eta_gen;
		Float_t         Z_muminus_theta_gen;
		Float_t         Z_muminus_phi_gen;
		Int_t           Z_muminus_charge_gen;
		Float_t         Z_muminus_vx_gen;
		Float_t         Z_muminus_vy_gen;
		Float_t         Z_muminus_vz_gen;
		Float_t         Z_muminus_y_gen;
		Int_t           event_runNo;
		Int_t           event_evtNo;
		Int_t           event_lumi;
		Int_t           event_bunch;
		Int_t           event_nPV;
		Float_t         event_met_pfmet;
		Float_t         event_met_pfsumet;
		Float_t         event_met_pfmetsignificance;
		Float_t         event_met_pfmetPhi;
		Float_t         event_metMVA_met;
		Float_t         event_metMVA_sumet;
		Float_t         event_metMVA_metsignificance;
		Float_t         event_metMVA_metPhi;
		Float_t         event_fastJetRho;
		Float_t         event_met_genmet;
		Float_t         event_met_gensumet;
		Float_t         event_met_genmetsignificance;
		Float_t         event_met_genmetPhi;
		Float_t         event_mcPU_totnvtx;
		Float_t         event_mcPU_trueInteractions;
		Float_t         event_mcPU_bx[3];
		Float_t         event_mcPU_nvtx[3];
		// List of branches
		TBranch        *b_numPFCorJets;   //!
		TBranch        *b_numPFCorJetBTags;   //!
		TBranch        *b_JetPFCor_Et;   //!
		TBranch        *b_JetPFCor_Pt;   //!
		TBranch        *b_JetPFCor_Eta;   //!
		TBranch        *b_JetPFCor_Phi;   //!
		TBranch        *b_JetPFCor_Theta;   //!
		TBranch        *b_JetPFCor_Px;   //!
		TBranch        *b_JetPFCor_Py;   //!
		TBranch        *b_JetPFCor_Pz;   //!
		TBranch        *b_JetPFCor_E;   //!
		TBranch        *b_JetPFCor_Y;   //!
		TBranch        *b_JetPFCor_Mass;   //!
		TBranch        *b_JetPFCor_etaetaMoment;   //!
		TBranch        *b_JetPFCor_phiphiMoment;   //!
		TBranch        *b_JetPFCor_etaphiMoment;   //!
		TBranch        *b_JetPFCor_maxDistance;   //!
		TBranch        *b_JetPFCor_nConstituents;   //!
		TBranch        *b_JetPFCor_Area;   //!
		TBranch        *b_VplusPFCorJet_Mass;   //!
		TBranch        *b_JetPFCor_dphiBoson;   //!
		TBranch        *b_JetPFCor_detaBoson;   //!
		TBranch        *b_JetPFCor_dRBoson;   //!
		TBranch        *b_JetPFCor_dphiMET;   //!
		TBranch        *b_JetPFCor_bDiscriminator;   //!
		TBranch        *b_JetPFCor_bDiscriminatorSSVHE;   //!
		TBranch        *b_JetPFCor_bDiscriminatorTCHE;   //!
		TBranch        *b_JetPFCor_bDiscriminatorCSV;   //!
		TBranch        *b_JetPFCor_bDiscriminatorJP;   //!
		TBranch        *b_JetPFCor_bDiscriminatorSSVHP;   //!
		TBranch        *b_JetPFCor_bDiscriminatorTCHP;   //!
		TBranch        *b_JetPFCor_secVertexMass;   //!
		TBranch        *b_JetPFCor_ChargedHadronEnergy;   //!
		TBranch        *b_JetPFCor_ChargedHadronEnergyFrac;   //!
		TBranch        *b_JetPFCor_NeutralHadronEnergy;   //!
		TBranch        *b_JetPFCor_NeutralHadronEnergyFrac;   //!
		TBranch        *b_JetPFCor_ChargedEmEnergy;   //!
		TBranch        *b_JetPFCor_ChargedEmEnergyFrac;   //!
		TBranch        *b_JetPFCor_ChargedMuEnergy;   //!
		TBranch        *b_JetPFCor_ChargedMuEnergyFrac;   //!
		TBranch        *b_JetPFCor_NeutralEmEnergy;   //!
		TBranch        *b_JetPFCor_NeutralEmEnergyFrac;   //!
		TBranch        *b_JetPFCor_ChargedMultiplicity;   //!
		TBranch        *b_JetPFCor_NeutralMultiplicity;   //!
		TBranch        *b_JetPFCor_MuonMultiplicity;   //!
		TBranch        *b_JetPFCor_PhotonEnergy;   //!
		TBranch        *b_JetPFCor_PhotonEnergyFraction;   //!
		TBranch        *b_JetPFCor_ElectronEnergy;   //!
		TBranch        *b_JetPFCor_ElectronEnergyFraction;   //!
		TBranch        *b_JetPFCor_MuonEnergy;   //!
		TBranch        *b_JetPFCor_MuonEnergyFraction;   //!
		TBranch        *b_JetPFCor_HFHadronEnergy;   //!
		TBranch        *b_JetPFCor_HFHadronEnergyFraction;   //!
		TBranch        *b_JetPFCor_HFEMEnergy;   //!
		TBranch        *b_JetPFCor_HFEMEnergyFraction;   //!
		TBranch        *b_JetPFCor_ChargedHadronMultiplicity;   //!
		TBranch        *b_JetPFCor_NeutralHadronMultiplicity;   //!
		TBranch        *b_JetPFCor_PhotonMultiplicity;   //!
		TBranch        *b_JetPFCor_ElectronMultiplicity;   //!
		TBranch        *b_JetPFCor_HFHadronMultiplicity;   //!
		TBranch        *b_JetPFCor_HFEMMultiplicity;   //!
		TBranch        *b_JetPFCor_SumPtCands;   //!
		TBranch        *b_JetPFCor_SumPt2Cands;   //!
		TBranch        *b_JetPFCor_rmsCands;   //!
		TBranch        *b_JetPFCor_PtD;   //!
		TBranch        *b_JetPFCor_QGLikelihood;   //!
		TBranch        *b_MassV2j_PFCor_MVAMET;   //!
		TBranch        *b_MassV2j_PFCor;   //!
		TBranch        *b_MassV3j_PFCor;   //!
		TBranch        *b_MassV4j_PFCor;   //!
		TBranch        *b_MassV5j_PFCor;   //!
		TBranch        *b_MassV6j_PFCor;   //!
		TBranch        *b_Mass2j_PFCor;   //!
		TBranch        *b_Mass3j_PFCor;   //!
		TBranch        *b_Mass4j_PFCor;   //!
		TBranch        *b_Mass5j_PFCor;   //!
		TBranch        *b_Mass6j_PFCor;   //!
		TBranch        *b_cosJacksonAngleV2j_PFCor;   //!
		TBranch        *b_cosJacksonAngle2j_PFCor;   //!
		TBranch        *b_cosJacksonAngleV3j_PFCor;   //!
		TBranch        *b_cosJacksonAngle3j12_PFCor;   //!
		TBranch        *b_cosJacksonAngle3j23_PFCor;   //!
		TBranch        *b_cosJacksonAngle3j31_PFCor;   //!
		TBranch        *b_cosphiDecayPlane_PFCor;   //!
		TBranch        *b_cosThetaLnu_PFCor;   //!
		TBranch        *b_cosThetaJJ_PFCor;   //!
		TBranch        *b_colorCorrPull01PFCor;   //!
		TBranch        *b_colorCorrPull02PFCor;   //!
		TBranch        *b_colorCorrPull12PFCor;   //!
		TBranch        *b_colorCorrPull03PFCor;   //!
		TBranch        *b_colorCorrPull13PFCor;   //!
		TBranch        *b_colorCorrPull23PFCor;   //!
		TBranch        *b_colorCorrPull04PFCor;   //!
		TBranch        *b_colorCorrPull14PFCor;   //!
		TBranch        *b_colorCorrPull24PFCor;   //!
		TBranch        *b_colorCorrPull34PFCor;   //!
		TBranch        *b_colorCorrPull05PFCor;   //!
		TBranch        *b_colorCorrPull15PFCor;   //!
		TBranch        *b_colorCorrPull25PFCor;   //!
		TBranch        *b_colorCorrPull35PFCor;   //!
		TBranch        *b_colorCorrPull45PFCor;   //!
		TBranch        *b_cosThetaJ1HiggsCM_PFCor;   //!
		TBranch        *b_cosThetaJ2HiggsCM_PFCor;   //!
		TBranch        *b_cosThetaL1HiggsCM_PFCor;   //!
		TBranch        *b_cosThetaL2HiggsCM_PFCor;   //!
		TBranch        *b_cosThetaV1HiggsCM_PFCor;   //!
		TBranch        *b_cosThetaV2HiggsCM_PFCor;   //!
		TBranch        *b_JetPFCor_isPileUpJetLoose;   //!
		TBranch        *b_JetPFCor_isPileUpJetMedium;   //!
		TBranch        *b_JetPFCor_isPileUpJetTight;   //!
		TBranch        *b_GroomedJet_CA8_pt_uncorr;   //!
		TBranch        *b_GroomedJet_CA8_mass_uncorr;   //!
		TBranch        *b_GroomedJet_CA8_mass_tr_uncorr;   //!
		TBranch        *b_GroomedJet_CA8_mass_ft_uncorr;   //!
		TBranch        *b_GroomedJet_CA8_mass_pr_uncorr;   //!
		TBranch        *b_GroomedJet_CA8_tau2tau1;   //!
		TBranch        *b_GroomedJet_CA8_tau1;   //!
		TBranch        *b_GroomedJet_CA8_tau2;   //!
		TBranch        *b_GroomedJet_CA8_tau3;   //!
		TBranch        *b_GroomedJet_CA8_tau4;   //!
		TBranch        *b_GroomedJet_CA8_massdrop_pr_uncorr;   //!
		TBranch        *b_GroomedJet_CA8_pt;   //!
		TBranch        *b_GroomedJet_CA8_eta;   //!
		TBranch        *b_GroomedJet_CA8_phi;   //!
		TBranch        *b_GroomedJet_CA8_e;   //!
		TBranch        *b_GroomedJet_CA8_pt_tr_uncorr;   //!
		TBranch        *b_GroomedJet_CA8_pt_tr;   //!
		TBranch        *b_GroomedJet_CA8_eta_tr;   //!
		TBranch        *b_GroomedJet_CA8_phi_tr;   //!
		TBranch        *b_GroomedJet_CA8_e_tr;   //!
		TBranch        *b_GroomedJet_CA8_pt_ft_uncorr;   //!
		TBranch        *b_GroomedJet_CA8_pt_ft;   //!
		TBranch        *b_GroomedJet_CA8_eta_ft;   //!
		TBranch        *b_GroomedJet_CA8_phi_ft;   //!
		TBranch        *b_GroomedJet_CA8_e_ft;   //!
		TBranch        *b_GroomedJet_CA8_pt_pr_uncorr;   //!
		TBranch        *b_GroomedJet_CA8_pt_pr;   //!
		TBranch        *b_GroomedJet_CA8_eta_pr;   //!
		TBranch        *b_GroomedJet_CA8_phi_pr;   //!
		TBranch        *b_GroomedJet_CA8_e_pr;   //!
		TBranch        *b_GroomedJet_CA8_prsubjet1_px;   //!
		TBranch        *b_GroomedJet_CA8_prsubjet1_py;   //!
		TBranch        *b_GroomedJet_CA8_prsubjet1_pz;   //!
		TBranch        *b_GroomedJet_CA8_prsubjet1_e;   //!
		TBranch        *b_GroomedJet_CA8_prsubjet2_px;   //!
		TBranch        *b_GroomedJet_CA8_prsubjet2_py;   //!
		TBranch        *b_GroomedJet_CA8_prsubjet2_pz;   //!
		TBranch        *b_GroomedJet_CA8_prsubjet2_e;   //!
		TBranch        *b_GroomedJet_CA8_mass;   //!
		TBranch        *b_GroomedJet_CA8_mass_tr;   //!
		TBranch        *b_GroomedJet_CA8_mass_ft;   //!
		TBranch        *b_GroomedJet_CA8_mass_pr;   //!
		TBranch        *b_GroomedJet_CA8_massdrop_pr;   //!
		TBranch        *b_GroomedJet_CA8_area;   //!
		TBranch        *b_GroomedJet_CA8_area_tr;   //!
		TBranch        *b_GroomedJet_CA8_area_ft;   //!
		TBranch        *b_GroomedJet_CA8_area_pr;   //!
		TBranch        *b_GroomedJet_CA8_jetconstituents;   //!
		TBranch        *b_GroomedJet_CA8_jetcharge;   //!
		TBranch        *b_GroomedJet_CA8_rcores;   //!
		TBranch        *b_GroomedJet_CA8_ptcores;   //!
		TBranch        *b_GroomedJet_CA8_planarflow;   //!
		TBranch        *b_GroomedJet_CA8_qjetmass;   //!
		TBranch        *b_GroomedJet_CA8_qjetmassdrop;   //!
		TBranch        *b_GroomedJet_CA8_constituents0_eta;   //!
		TBranch        *b_GroomedJet_CA8_constituents0_phi;   //!
		TBranch        *b_GroomedJet_CA8_constituents0_e;   //!
		TBranch        *b_GroomedJet_CA8_nconstituents0;   //!
		TBranch        *b_GroomedJet_CA8_constituents0pr_eta;   //!
		TBranch        *b_GroomedJet_CA8_constituents0pr_phi;   //!
		TBranch        *b_GroomedJet_CA8_constituents0pr_e;   //!
		TBranch        *b_GroomedJet_CA8_nconstituents0pr;   //!
		TBranch        *b_GenGroomedJet_CA8_pt_uncorr;   //!
		TBranch        *b_GenGroomedJet_CA8_mass_uncorr;   //!
		TBranch        *b_GenGroomedJet_CA8_mass_tr_uncorr;   //!
		TBranch        *b_GenGroomedJet_CA8_mass_ft_uncorr;   //!
		TBranch        *b_GenGroomedJet_CA8_mass_pr_uncorr;   //!
		TBranch        *b_GenGroomedJet_CA8_tau2tau1;   //!
		TBranch        *b_GenGroomedJet_CA8_tau1;   //!
		TBranch        *b_GenGroomedJet_CA8_tau2;   //!
		TBranch        *b_GenGroomedJet_CA8_tau3;   //!
		TBranch        *b_GenGroomedJet_CA8_tau4;   //!
		TBranch        *b_GenGroomedJet_CA8_massdrop_pr_uncorr;   //!
		TBranch        *b_GenGroomedJet_CA8_pt;   //!
		TBranch        *b_GenGroomedJet_CA8_eta;   //!
		TBranch        *b_GenGroomedJet_CA8_phi;   //!
		TBranch        *b_GenGroomedJet_CA8_e;   //!
		TBranch        *b_GenGroomedJet_CA8_pt_tr_uncorr;   //!
		TBranch        *b_GenGroomedJet_CA8_pt_tr;   //!
		TBranch        *b_GenGroomedJet_CA8_eta_tr;   //!
		TBranch        *b_GenGroomedJet_CA8_phi_tr;   //!
		TBranch        *b_GenGroomedJet_CA8_e_tr;   //!
		TBranch        *b_GenGroomedJet_CA8_pt_ft_uncorr;   //!
		TBranch        *b_GenGroomedJet_CA8_pt_ft;   //!
		TBranch        *b_GenGroomedJet_CA8_eta_ft;   //!
		TBranch        *b_GenGroomedJet_CA8_phi_ft;   //!
		TBranch        *b_GenGroomedJet_CA8_e_ft;   //!
		TBranch        *b_GenGroomedJet_CA8_pt_pr_uncorr;   //!
		TBranch        *b_GenGroomedJet_CA8_pt_pr;   //!
		TBranch        *b_GenGroomedJet_CA8_eta_pr;   //!
		TBranch        *b_GenGroomedJet_CA8_phi_pr;   //!
		TBranch        *b_GenGroomedJet_CA8_e_pr;   //!
		TBranch        *b_GenGroomedJet_CA8_prsubjet1_px;   //!
		TBranch        *b_GenGroomedJet_CA8_prsubjet1_py;   //!
		TBranch        *b_GenGroomedJet_CA8_prsubjet1_pz;   //!
		TBranch        *b_GenGroomedJet_CA8_prsubjet1_e;   //!
		TBranch        *b_GenGroomedJet_CA8_prsubjet2_px;   //!
		TBranch        *b_GenGroomedJet_CA8_prsubjet2_py;   //!
		TBranch        *b_GenGroomedJet_CA8_prsubjet2_pz;   //!
		TBranch        *b_GenGroomedJet_CA8_prsubjet2_e;   //!
		TBranch        *b_GenGroomedJet_CA8_mass;   //!
		TBranch        *b_GenGroomedJet_CA8_mass_tr;   //!
		TBranch        *b_GenGroomedJet_CA8_mass_ft;   //!
		TBranch        *b_GenGroomedJet_CA8_mass_pr;   //!
		TBranch        *b_GenGroomedJet_CA8_massdrop_pr;   //!
		TBranch        *b_GenGroomedJet_CA8_area;   //!
		TBranch        *b_GenGroomedJet_CA8_area_tr;   //!
		TBranch        *b_GenGroomedJet_CA8_area_ft;   //!
		TBranch        *b_GenGroomedJet_CA8_area_pr;   //!
		TBranch        *b_GenGroomedJet_CA8_jetconstituents;   //!
		TBranch        *b_GenGroomedJet_CA8_jetcharge;   //!
		TBranch        *b_GenGroomedJet_CA8_rcores;   //!
		TBranch        *b_GenGroomedJet_CA8_ptcores;   //!
		TBranch        *b_GenGroomedJet_CA8_planarflow;   //!
		TBranch        *b_GenGroomedJet_CA8_qjetmass;   //!
		TBranch        *b_GenGroomedJet_CA8_qjetmassdrop;   //!
		TBranch        *b_GenGroomedJet_CA8_constituents0_eta;   //!
		TBranch        *b_GenGroomedJet_CA8_constituents0_phi;   //!
		TBranch        *b_GenGroomedJet_CA8_constituents0_e;   //!
		TBranch        *b_GenGroomedJet_CA8_nconstituents0;   //!
		TBranch        *b_GenGroomedJet_CA8_constituents0pr_eta;   //!
		TBranch        *b_GenGroomedJet_CA8_constituents0pr_phi;   //!
		TBranch        *b_GenGroomedJet_CA8_constituents0pr_e;   //!
		TBranch        *b_GenGroomedJet_CA8_nconstituents0pr;   //!
		TBranch        *b_numGenJets;   //!
		TBranch        *b_numGenJetBTags;   //!
		TBranch        *b_JetGen_Et;   //!
		TBranch        *b_JetGen_Pt;   //!
		TBranch        *b_JetGen_Eta;   //!
		TBranch        *b_JetGen_Phi;   //!
		TBranch        *b_JetGen_Theta;   //!
		TBranch        *b_JetGen_Px;   //!
		TBranch        *b_JetGen_Py;   //!
		TBranch        *b_JetGen_Pz;   //!
		TBranch        *b_JetGen_E;   //!
		TBranch        *b_JetGen_Y;   //!
		TBranch        *b_JetGen_Mass;   //!
		TBranch        *b_JetGen_etaetaMoment;   //!
		TBranch        *b_JetGen_phiphiMoment;   //!
		TBranch        *b_JetGen_etaphiMoment;   //!
		TBranch        *b_JetGen_maxDistance;   //!
		TBranch        *b_JetGen_nConstituents;   //!
		TBranch        *b_JetGen_Area;   //!
		TBranch        *b_VplusGenJet_Mass;   //!
		TBranch        *b_JetGen_dphiBoson;   //!
		TBranch        *b_JetGen_detaBoson;   //!
		TBranch        *b_JetGen_dRBoson;   //!
		TBranch        *b_JetGen_dphiMET;   //!
		TBranch        *b_JetGen_bDiscriminator;   //!
		TBranch        *b_JetGen_bDiscriminatorSSVHE;   //!
		TBranch        *b_JetGen_bDiscriminatorTCHE;   //!
		TBranch        *b_JetGen_bDiscriminatorCSV;   //!
		TBranch        *b_JetGen_bDiscriminatorJP;   //!
		TBranch        *b_JetGen_bDiscriminatorSSVHP;   //!
		TBranch        *b_JetGen_bDiscriminatorTCHP;   //!
		TBranch        *b_JetGen_secVertexMass;   //!
		TBranch        *b_MassV2j_Gen_MVAMET;   //!
		TBranch        *b_MassV2j_Gen;   //!
		TBranch        *b_MassV3j_Gen;   //!
		TBranch        *b_MassV4j_Gen;   //!
		TBranch        *b_MassV5j_Gen;   //!
		TBranch        *b_MassV6j_Gen;   //!
		TBranch        *b_Mass2j_Gen;   //!
		TBranch        *b_Mass3j_Gen;   //!
		TBranch        *b_Mass4j_Gen;   //!
		TBranch        *b_Mass5j_Gen;   //!
		TBranch        *b_Mass6j_Gen;   //!
		TBranch        *b_cosJacksonAngleV2j_Gen;   //!
		TBranch        *b_cosJacksonAngle2j_Gen;   //!
		TBranch        *b_cosJacksonAngleV3j_Gen;   //!
		TBranch        *b_cosJacksonAngle3j12_Gen;   //!
		TBranch        *b_cosJacksonAngle3j23_Gen;   //!
		TBranch        *b_cosJacksonAngle3j31_Gen;   //!
		TBranch        *b_cosphiDecayPlane_Gen;   //!
		TBranch        *b_cosThetaLnu_Gen;   //!
		TBranch        *b_cosThetaJJ_Gen;   //!
		TBranch        *b_NumPhotons;   //!
		TBranch        *b_Photon_Et;   //!
		TBranch        *b_Photon_E;   //!
		TBranch        *b_Photon_Eta;   //!
		TBranch        *b_Photon_Phi;   //!
		TBranch        *b_Photon_Theta;   //!
		TBranch        *b_Photon_Px;   //!
		TBranch        *b_Photon_Py;   //!
		TBranch        *b_Photon_Pz;   //!
		TBranch        *b_Photon_Vx;   //!
		TBranch        *b_Photon_Vy;   //!
		TBranch        *b_Photon_Vz;   //!
		TBranch        *b_Photon_SC_Et;   //!
		TBranch        *b_Photon_SC_E;   //!
		TBranch        *b_Photon_SC_Eta;   //!
		TBranch        *b_Photon_SC_Phi;   //!
		TBranch        *b_Photon_SC_Theta;   //!
		TBranch        *b_Photon_SC_x;   //!
		TBranch        *b_Photon_SC_y;   //!
		TBranch        *b_Photon_SC_z;   //!
		TBranch        *b_PFisocharged03;   //!
		TBranch        *b_PFisophoton03;   //!
		TBranch        *b_PFisoneutral03;   //!
		TBranch        *b_trkSumPtHollowConeDR04_Photon11;   //!
		TBranch        *b_ecalRecHitSumEtConeDR04_Photon11;   //!
		TBranch        *b_hcalTowerSumEtConeDR04_Photon11;   //!
		TBranch        *b_Photon_HoverE;   //!
		TBranch        *b_Photon_HoverE2011;   //!
		TBranch        *b_Photon_SigmaIetaIeta;   //!
		TBranch        *b_Photon_hasPixelSeed;   //!
		TBranch        *b_Photon_passElecVeto;   //!
		TBranch        *b_Photon_Id2011;   //!
		TBranch        *b_Photon_Id2012;   //!
		TBranch        *b_Z_mass;   //!
		TBranch        *b_Z_mt;   //!
		TBranch        *b_Z_mtMVA;   //!
		TBranch        *b_Z_px;   //!
		TBranch        *b_Z_py;   //!
		TBranch        *b_Z_pz;   //!
		TBranch        *b_Z_e;   //!
		TBranch        *b_Z_pt;   //!
		TBranch        *b_Z_et;   //!
		TBranch        *b_Z_eta;   //!
		TBranch        *b_Z_phi;   //!
		TBranch        *b_Z_vx;   //!
		TBranch        *b_Z_vy;   //!
		TBranch        *b_Z_vz;   //!
		TBranch        *b_Z_y;   //!
		TBranch        *b_Z_muplus_px;   //!
		TBranch        *b_Z_muplus_py;   //!
		TBranch        *b_Z_muplus_pz;   //!
		TBranch        *b_Z_muplus_e;   //!
		TBranch        *b_Z_muplus_pt;   //!
		TBranch        *b_Z_muplus_et;   //!
		TBranch        *b_Z_muplus_eta;   //!
		TBranch        *b_Z_muplus_theta;   //!
		TBranch        *b_Z_muplus_phi;   //!
		TBranch        *b_Z_muplus_charge;   //!
		TBranch        *b_Z_muplus_vx;   //!
		TBranch        *b_Z_muplus_vy;   //!
		TBranch        *b_Z_muplus_vz;   //!
		TBranch        *b_Z_muplus_y;   //!
		TBranch        *b_Z_muplus_trackiso;   //!
		TBranch        *b_Z_muplus_hcaliso;   //!
		TBranch        *b_Z_muplus_ecaliso;   //!
		TBranch        *b_Z_muplus_type;   //!
		TBranch        *b_Z_muplus_numberOfChambers;   //!
		TBranch        *b_Z_muplus_numberOfMatches;   //!
		TBranch        *b_Z_muplus_d0bsp;   //!
		TBranch        *b_Z_muplus_dz000;   //!
		TBranch        *b_Z_muplus_dzPV;   //!
		TBranch        *b_Z_muplus_pfiso_sumChargedHadronPt;   //!
		TBranch        *b_Z_muplus_pfiso_sumChargedParticlePt;   //!
		TBranch        *b_Z_muplus_pfiso_sumNeutralHadronEt;   //!
		TBranch        *b_Z_muplus_pfiso_sumPhotonEt;   //!
		TBranch        *b_Z_muplus_pfiso_sumPUPt;   //!
		TBranch        *b_Z_muminus_px;   //!
		TBranch        *b_Z_muminus_py;   //!
		TBranch        *b_Z_muminus_pz;   //!
		TBranch        *b_Z_muminus_e;   //!
		TBranch        *b_Z_muminus_pt;   //!
		TBranch        *b_Z_muminus_et;   //!
		TBranch        *b_Z_muminus_eta;   //!
		TBranch        *b_Z_muminus_theta;   //!
		TBranch        *b_Z_muminus_phi;   //!
		TBranch        *b_Z_muminus_charge;   //!
		TBranch        *b_Z_muminus_vx;   //!
		TBranch        *b_Z_muminus_vy;   //!
		TBranch        *b_Z_muminus_vz;   //!
		TBranch        *b_Z_muminus_y;   //!
		TBranch        *b_Z_muminus_trackiso;   //!
		TBranch        *b_Z_muminus_hcaliso;   //!
		TBranch        *b_Z_muminus_ecaliso;   //!
		TBranch        *b_Z_muminus_type;   //!
		TBranch        *b_Z_muminus_numberOfChambers;   //!
		TBranch        *b_Z_muminus_numberOfMatches;   //!
		TBranch        *b_Z_muminus_d0bsp;   //!
		TBranch        *b_Z_muminus_dz000;   //!
		TBranch        *b_Z_muminus_dzPV;   //!
		TBranch        *b_Z_muminus_pfiso_sumChargedHadronPt;   //!
		TBranch        *b_Z_muminus_pfiso_sumChargedParticlePt;   //!
		TBranch        *b_Z_muminus_pfiso_sumNeutralHadronEt;   //!
		TBranch        *b_Z_muminus_pfiso_sumPhotonEt;   //!
		TBranch        *b_Z_muminus_pfiso_sumPUPt;   //!
		TBranch        *b_Z_Photon_pt_gen;   //!
		TBranch        *b_Z_muplus_px_gen;   //!
		TBranch        *b_Z_muplus_py_gen;   //!
		TBranch        *b_Z_muplus_pz_gen;   //!
		TBranch        *b_Z_muplus_e_gen;   //!
		TBranch        *b_Z_muplus_pt_gen;   //!
		TBranch        *b_Z_muplus_et_gen;   //!
		TBranch        *b_Z_muplus_eta_gen;   //!
		TBranch        *b_Z_muplus_theta_gen;   //!
		TBranch        *b_Z_muplus_phi_gen;   //!
		TBranch        *b_Z_muplus_charge_gen;   //!
		TBranch        *b_Z_muplus_vx_gen;   //!
		TBranch        *b_Z_muplus_vy_gen;   //!
		TBranch        *b_Z_muplus_vz_gen;   //!
		TBranch        *b_Z_muplus_y_gen;   //!
		TBranch        *b_Z_muminus_px_gen;   //!
		TBranch        *b_Z_muminus_py_gen;   //!
		TBranch        *b_Z_muminus_pz_gen;   //!
		TBranch        *b_Z_muminus_e_gen;   //!
		TBranch        *b_Z_muminus_pt_gen;   //!
		TBranch        *b_Z_muminus_et_gen;   //!
		TBranch        *b_Z_muminus_eta_gen;   //!
		TBranch        *b_Z_muminus_theta_gen;   //!
		TBranch        *b_Z_muminus_phi_gen;   //!
		TBranch        *b_Z_muminus_charge_gen;   //!
		TBranch        *b_Z_muminus_vx_gen;   //!
		TBranch        *b_Z_muminus_vy_gen;   //!
		TBranch        *b_Z_muminus_vz_gen;   //!
		TBranch        *b_Z_muminus_y_gen;   //!
		TBranch        *b_event_runNo;   //!
		TBranch        *b_event_evtNo;   //!
		TBranch        *b_event_lumi;   //!
		TBranch        *b_event_bunch;   //!
		TBranch        *b_event_nPV;   //!
		TBranch        *b_event_met_pfmet;   //!
		TBranch        *b_event_met_pfsumet;   //!
		TBranch        *b_event_met_pfmetsignificance;   //!
		TBranch        *b_event_met_pfmetPhi;   //!
		TBranch        *b_event_metMVA_met;   //!
		TBranch        *b_event_metMVA_sumet;   //!
		TBranch        *b_event_metMVA_metsignificance;   //!
		TBranch        *b_event_metMVA_metPhi;   //!
		TBranch        *b_event_fastJetRho;   //!
		TBranch        *b_event_met_genmet;   //!
		TBranch        *b_event_met_gensumet;   //!
		TBranch        *b_event_met_genmetsignificance;   //!
		TBranch        *b_event_met_genmetPhi;   //!
		TBranch        *b_event_mcPU_totnvtx;   //!
		TBranch        *b_event_mcPU_trueInteractions;   //!
		TBranch        *b_event_mcPU_bx;   //!
		TBranch        *b_event_mcPU_nvtx;   //!

		kanamuon(TTree *tree=0);
		virtual ~kanamuon();
		virtual Int_t    Cut(Long64_t entry);
		virtual Int_t    GetEntry(Long64_t entry);
		virtual Long64_t LoadTree(Long64_t entry);
		virtual void     Init(TTree *tree);
		virtual Bool_t   Notify();
		virtual void     Show(Long64_t entry = -1);

		virtual void     myana(double myflag = -999, bool isQCD = false, int runflag=0);
		virtual void     Loop(TH1F* h_events,
					TH1F* h_events_weighted,
					int wda,
					int runflag,
					const char * outfilename,
					bool isQCD = false);
		virtual double   getDeltaPhi(double phi1, double phi2);
		//virtual bool     large(const double &a, const double &b);
		virtual bool     doKinematicFit(Int_t                 fflage,
					const TLorentzVector     mup, 
					const TLorentzVector     nvp,
					const TLorentzVector     ajp, 
					const TLorentzVector     bjp, 
					TLorentzVector     & fit_mup, 
					TLorentzVector     & fit_nvp,
					TLorentzVector     & fit_ajp, 
					TLorentzVector     & fit_bjp, 
					Float_t            & fit_chi2,
					Int_t              & fit_NDF, 
					Int_t              & fit_status);
		/*      virtual bool    dottHKinematicFit(const TLorentzVector     mup, 
				const TLorentzVector     nvp, 
				const TLorentzVector     wajp,
				const TLorentzVector     wbjp,
				const TLorentzVector     topajp,
				const TLorentzVector     topbjp,
				Float_t            & fit_chi2,
				Int_t              & fit_NDF,
				Int_t              & fit_status);
				*/
		virtual void     calculateAngles( TLorentzVector& thep4M11, TLorentzVector& thep4M12, TLorentzVector& thep4M21, TLorentzVector& thep4M22, double& costheta1, double& costheta2, double& phi, double& costhetastar, double& phistar1, double& phistar2);
		virtual void     InitCounters(const char* input_file_name, TH1F* h_events, TH1F* h_events_weighted);

};

#endif //kanamuon_ZJets_h

#ifdef kanamuon_ZJets_cxx
kanamuon::kanamuon(TTree *tree)
{
	// if parameter tree is not specified (or zero), connect the file
	// used to generate this class and read the Tree.
	if (tree == 0) {
		//TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/uscms/home/zixu/work/ZPlusJet/CMSSW_5_3_2_patch4/src/ElectroWeakAnalysis/VPlusJets/test/ZPlusJets_jetsubstructure/ZmumuJetsAnalysisntuple.root");
		TFile *f = new TFile("/uscms/home/zixu/work/ZPlusJet/CMSSW_5_3_2_patch4/src/ElectroWeakAnalysis/VPlusJets/test/ZPlusJets_jetsubstructure/ZmumuJetsAnalysisntuple.root");
		if (!f) { cout<<"can't find the root file"<<endl; f = new TFile("tmp.root"); }
		tree = (TTree*)gDirectory->Get("ZJet");

	}
	Init(tree);
}

kanamuon::~kanamuon()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

Int_t kanamuon::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
Long64_t kanamuon::LoadTree(Long64_t entry)
{
	// Set the environment to read one entry
	if (!fChain) return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if (centry < 0) return centry;
	if (!fChain->InheritsFrom(TChain::Class()))  return centry;
	TChain *chain = (TChain*)fChain;
	if (chain->GetTreeNumber() != fCurrent) {
		fCurrent = chain->GetTreeNumber();
		Notify();
	}
	return centry;
}

void kanamuon::Init(TTree *tree)
{
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).

	// Set branch addresses and branch pointers
	if (!tree) return;
	fChain = tree;
	fCurrent = -1;
	fChain->SetMakeClass(1);

	fChain->SetBranchAddress("numPFCorJets", &numPFCorJets, &b_numPFCorJets);
	fChain->SetBranchAddress("numPFCorJetBTags", &numPFCorJetBTags, &b_numPFCorJetBTags);
	fChain->SetBranchAddress("JetPFCor_Et", JetPFCor_Et, &b_JetPFCor_Et);
	fChain->SetBranchAddress("JetPFCor_Pt", JetPFCor_Pt, &b_JetPFCor_Pt);
	fChain->SetBranchAddress("JetPFCor_Eta", JetPFCor_Eta, &b_JetPFCor_Eta);
	fChain->SetBranchAddress("JetPFCor_Phi", JetPFCor_Phi, &b_JetPFCor_Phi);
	fChain->SetBranchAddress("JetPFCor_Theta", JetPFCor_Theta, &b_JetPFCor_Theta);
	fChain->SetBranchAddress("JetPFCor_Px", JetPFCor_Px, &b_JetPFCor_Px);
	fChain->SetBranchAddress("JetPFCor_Py", JetPFCor_Py, &b_JetPFCor_Py);
	fChain->SetBranchAddress("JetPFCor_Pz", JetPFCor_Pz, &b_JetPFCor_Pz);
	fChain->SetBranchAddress("JetPFCor_E", JetPFCor_E, &b_JetPFCor_E);
	fChain->SetBranchAddress("JetPFCor_Y", JetPFCor_Y, &b_JetPFCor_Y);
	fChain->SetBranchAddress("JetPFCor_Mass", JetPFCor_Mass, &b_JetPFCor_Mass);
	fChain->SetBranchAddress("JetPFCor_etaetaMoment", JetPFCor_etaetaMoment, &b_JetPFCor_etaetaMoment);
	fChain->SetBranchAddress("JetPFCor_phiphiMoment", JetPFCor_phiphiMoment, &b_JetPFCor_phiphiMoment);
	fChain->SetBranchAddress("JetPFCor_etaphiMoment", JetPFCor_etaphiMoment, &b_JetPFCor_etaphiMoment);
	fChain->SetBranchAddress("JetPFCor_maxDistance", JetPFCor_maxDistance, &b_JetPFCor_maxDistance);
	fChain->SetBranchAddress("JetPFCor_nConstituents", JetPFCor_nConstituents, &b_JetPFCor_nConstituents);
	fChain->SetBranchAddress("JetPFCor_Area", JetPFCor_Area, &b_JetPFCor_Area);
	fChain->SetBranchAddress("VplusPFCorJet_Mass", VplusPFCorJet_Mass, &b_VplusPFCorJet_Mass);
	fChain->SetBranchAddress("JetPFCor_dphiBoson", JetPFCor_dphiBoson, &b_JetPFCor_dphiBoson);
	fChain->SetBranchAddress("JetPFCor_detaBoson", JetPFCor_detaBoson, &b_JetPFCor_detaBoson);
	fChain->SetBranchAddress("JetPFCor_dRBoson", JetPFCor_dRBoson, &b_JetPFCor_dRBoson);
	fChain->SetBranchAddress("JetPFCor_dphiMET", JetPFCor_dphiMET, &b_JetPFCor_dphiMET);
	fChain->SetBranchAddress("JetPFCor_bDiscriminator", JetPFCor_bDiscriminator, &b_JetPFCor_bDiscriminator);
	fChain->SetBranchAddress("JetPFCor_bDiscriminatorSSVHE", JetPFCor_bDiscriminatorSSVHE, &b_JetPFCor_bDiscriminatorSSVHE);
	fChain->SetBranchAddress("JetPFCor_bDiscriminatorTCHE", JetPFCor_bDiscriminatorTCHE, &b_JetPFCor_bDiscriminatorTCHE);
	fChain->SetBranchAddress("JetPFCor_bDiscriminatorCSV", JetPFCor_bDiscriminatorCSV, &b_JetPFCor_bDiscriminatorCSV);
	fChain->SetBranchAddress("JetPFCor_bDiscriminatorJP", JetPFCor_bDiscriminatorJP, &b_JetPFCor_bDiscriminatorJP);
	fChain->SetBranchAddress("JetPFCor_bDiscriminatorSSVHP", JetPFCor_bDiscriminatorSSVHP, &b_JetPFCor_bDiscriminatorSSVHP);
	fChain->SetBranchAddress("JetPFCor_bDiscriminatorTCHP", JetPFCor_bDiscriminatorTCHP, &b_JetPFCor_bDiscriminatorTCHP);
	fChain->SetBranchAddress("JetPFCor_secVertexMass", JetPFCor_secVertexMass, &b_JetPFCor_secVertexMass);
	fChain->SetBranchAddress("JetPFCor_ChargedHadronEnergy", JetPFCor_ChargedHadronEnergy, &b_JetPFCor_ChargedHadronEnergy);
	fChain->SetBranchAddress("JetPFCor_ChargedHadronEnergyFrac", JetPFCor_ChargedHadronEnergyFrac, &b_JetPFCor_ChargedHadronEnergyFrac);
	fChain->SetBranchAddress("JetPFCor_NeutralHadronEnergy", JetPFCor_NeutralHadronEnergy, &b_JetPFCor_NeutralHadronEnergy);
	fChain->SetBranchAddress("JetPFCor_NeutralHadronEnergyFrac", JetPFCor_NeutralHadronEnergyFrac, &b_JetPFCor_NeutralHadronEnergyFrac);
	fChain->SetBranchAddress("JetPFCor_ChargedEmEnergy", JetPFCor_ChargedEmEnergy, &b_JetPFCor_ChargedEmEnergy);
	fChain->SetBranchAddress("JetPFCor_ChargedEmEnergyFrac", JetPFCor_ChargedEmEnergyFrac, &b_JetPFCor_ChargedEmEnergyFrac);
	fChain->SetBranchAddress("JetPFCor_ChargedMuEnergy", JetPFCor_ChargedMuEnergy, &b_JetPFCor_ChargedMuEnergy);
	fChain->SetBranchAddress("JetPFCor_ChargedMuEnergyFrac", JetPFCor_ChargedMuEnergyFrac, &b_JetPFCor_ChargedMuEnergyFrac);
	fChain->SetBranchAddress("JetPFCor_NeutralEmEnergy", JetPFCor_NeutralEmEnergy, &b_JetPFCor_NeutralEmEnergy);
	fChain->SetBranchAddress("JetPFCor_NeutralEmEnergyFrac", JetPFCor_NeutralEmEnergyFrac, &b_JetPFCor_NeutralEmEnergyFrac);
	fChain->SetBranchAddress("JetPFCor_ChargedMultiplicity", JetPFCor_ChargedMultiplicity, &b_JetPFCor_ChargedMultiplicity);
	fChain->SetBranchAddress("JetPFCor_NeutralMultiplicity", JetPFCor_NeutralMultiplicity, &b_JetPFCor_NeutralMultiplicity);
	fChain->SetBranchAddress("JetPFCor_MuonMultiplicity", JetPFCor_MuonMultiplicity, &b_JetPFCor_MuonMultiplicity);
	fChain->SetBranchAddress("JetPFCor_PhotonEnergy", JetPFCor_PhotonEnergy, &b_JetPFCor_PhotonEnergy);
	fChain->SetBranchAddress("JetPFCor_PhotonEnergyFraction", JetPFCor_PhotonEnergyFraction, &b_JetPFCor_PhotonEnergyFraction);
	fChain->SetBranchAddress("JetPFCor_ElectronEnergy", JetPFCor_ElectronEnergy, &b_JetPFCor_ElectronEnergy);
	fChain->SetBranchAddress("JetPFCor_ElectronEnergyFraction", JetPFCor_ElectronEnergyFraction, &b_JetPFCor_ElectronEnergyFraction);
	fChain->SetBranchAddress("JetPFCor_MuonEnergy", JetPFCor_MuonEnergy, &b_JetPFCor_MuonEnergy);
	fChain->SetBranchAddress("JetPFCor_MuonEnergyFraction", JetPFCor_MuonEnergyFraction, &b_JetPFCor_MuonEnergyFraction);
	fChain->SetBranchAddress("JetPFCor_HFHadronEnergy", JetPFCor_HFHadronEnergy, &b_JetPFCor_HFHadronEnergy);
	fChain->SetBranchAddress("JetPFCor_HFHadronEnergyFraction", JetPFCor_HFHadronEnergyFraction, &b_JetPFCor_HFHadronEnergyFraction);
	fChain->SetBranchAddress("JetPFCor_HFEMEnergy", JetPFCor_HFEMEnergy, &b_JetPFCor_HFEMEnergy);
	fChain->SetBranchAddress("JetPFCor_HFEMEnergyFraction", JetPFCor_HFEMEnergyFraction, &b_JetPFCor_HFEMEnergyFraction);
	fChain->SetBranchAddress("JetPFCor_ChargedHadronMultiplicity", JetPFCor_ChargedHadronMultiplicity, &b_JetPFCor_ChargedHadronMultiplicity);
	fChain->SetBranchAddress("JetPFCor_NeutralHadronMultiplicity", JetPFCor_NeutralHadronMultiplicity, &b_JetPFCor_NeutralHadronMultiplicity);
	fChain->SetBranchAddress("JetPFCor_PhotonMultiplicity", JetPFCor_PhotonMultiplicity, &b_JetPFCor_PhotonMultiplicity);
	fChain->SetBranchAddress("JetPFCor_ElectronMultiplicity", JetPFCor_ElectronMultiplicity, &b_JetPFCor_ElectronMultiplicity);
	fChain->SetBranchAddress("JetPFCor_HFHadronMultiplicity", JetPFCor_HFHadronMultiplicity, &b_JetPFCor_HFHadronMultiplicity);
	fChain->SetBranchAddress("JetPFCor_HFEMMultiplicity", JetPFCor_HFEMMultiplicity, &b_JetPFCor_HFEMMultiplicity);
	fChain->SetBranchAddress("JetPFCor_SumPtCands", JetPFCor_SumPtCands, &b_JetPFCor_SumPtCands);
	fChain->SetBranchAddress("JetPFCor_SumPt2Cands", JetPFCor_SumPt2Cands, &b_JetPFCor_SumPt2Cands);
	fChain->SetBranchAddress("JetPFCor_rmsCands", JetPFCor_rmsCands, &b_JetPFCor_rmsCands);
	fChain->SetBranchAddress("JetPFCor_PtD", JetPFCor_PtD, &b_JetPFCor_PtD);
	fChain->SetBranchAddress("JetPFCor_QGLikelihood", JetPFCor_QGLikelihood, &b_JetPFCor_QGLikelihood);
	fChain->SetBranchAddress("MassV2j_PFCor_MVAMET", &MassV2j_PFCor_MVAMET, &b_MassV2j_PFCor_MVAMET);
	fChain->SetBranchAddress("MassV2j_PFCor", &MassV2j_PFCor, &b_MassV2j_PFCor);
	fChain->SetBranchAddress("MassV3j_PFCor", &MassV3j_PFCor, &b_MassV3j_PFCor);
	fChain->SetBranchAddress("MassV4j_PFCor", &MassV4j_PFCor, &b_MassV4j_PFCor);
	fChain->SetBranchAddress("MassV5j_PFCor", &MassV5j_PFCor, &b_MassV5j_PFCor);
	fChain->SetBranchAddress("MassV6j_PFCor", &MassV6j_PFCor, &b_MassV6j_PFCor);
	fChain->SetBranchAddress("Mass2j_PFCor", &Mass2j_PFCor, &b_Mass2j_PFCor);
	fChain->SetBranchAddress("Mass3j_PFCor", &Mass3j_PFCor, &b_Mass3j_PFCor);
	fChain->SetBranchAddress("Mass4j_PFCor", &Mass4j_PFCor, &b_Mass4j_PFCor);
	fChain->SetBranchAddress("Mass5j_PFCor", &Mass5j_PFCor, &b_Mass5j_PFCor);
	fChain->SetBranchAddress("Mass6j_PFCor", &Mass6j_PFCor, &b_Mass6j_PFCor);
	fChain->SetBranchAddress("cosJacksonAngleV2j_PFCor", &cosJacksonAngleV2j_PFCor, &b_cosJacksonAngleV2j_PFCor);
	fChain->SetBranchAddress("cosJacksonAngle2j_PFCor", &cosJacksonAngle2j_PFCor, &b_cosJacksonAngle2j_PFCor);
	fChain->SetBranchAddress("cosJacksonAngleV3j_PFCor", &cosJacksonAngleV3j_PFCor, &b_cosJacksonAngleV3j_PFCor);
	fChain->SetBranchAddress("cosJacksonAngle3j12_PFCor", &cosJacksonAngle3j12_PFCor, &b_cosJacksonAngle3j12_PFCor);
	fChain->SetBranchAddress("cosJacksonAngle3j23_PFCor", &cosJacksonAngle3j23_PFCor, &b_cosJacksonAngle3j23_PFCor);
	fChain->SetBranchAddress("cosJacksonAngle3j31_PFCor", &cosJacksonAngle3j31_PFCor, &b_cosJacksonAngle3j31_PFCor);
	fChain->SetBranchAddress("cosphiDecayPlane_PFCor", &cosphiDecayPlane_PFCor, &b_cosphiDecayPlane_PFCor);
	fChain->SetBranchAddress("cosThetaLnu_PFCor", &cosThetaLnu_PFCor, &b_cosThetaLnu_PFCor);
	fChain->SetBranchAddress("cosThetaJJ_PFCor", &cosThetaJJ_PFCor, &b_cosThetaJJ_PFCor);
	fChain->SetBranchAddress("colorCorrPull01PFCor", &colorCorrPull01PFCor, &b_colorCorrPull01PFCor);
	fChain->SetBranchAddress("colorCorrPull02PFCor", &colorCorrPull02PFCor, &b_colorCorrPull02PFCor);
	fChain->SetBranchAddress("colorCorrPull12PFCor", &colorCorrPull12PFCor, &b_colorCorrPull12PFCor);
	fChain->SetBranchAddress("colorCorrPull03PFCor", &colorCorrPull03PFCor, &b_colorCorrPull03PFCor);
	fChain->SetBranchAddress("colorCorrPull13PFCor", &colorCorrPull13PFCor, &b_colorCorrPull13PFCor);
	fChain->SetBranchAddress("colorCorrPull23PFCor", &colorCorrPull23PFCor, &b_colorCorrPull23PFCor);
	fChain->SetBranchAddress("colorCorrPull04PFCor", &colorCorrPull04PFCor, &b_colorCorrPull04PFCor);
	fChain->SetBranchAddress("colorCorrPull14PFCor", &colorCorrPull14PFCor, &b_colorCorrPull14PFCor);
	fChain->SetBranchAddress("colorCorrPull24PFCor", &colorCorrPull24PFCor, &b_colorCorrPull24PFCor);
	fChain->SetBranchAddress("colorCorrPull34PFCor", &colorCorrPull34PFCor, &b_colorCorrPull34PFCor);
	fChain->SetBranchAddress("colorCorrPull05PFCor", &colorCorrPull05PFCor, &b_colorCorrPull05PFCor);
	fChain->SetBranchAddress("colorCorrPull15PFCor", &colorCorrPull15PFCor, &b_colorCorrPull15PFCor);
	fChain->SetBranchAddress("colorCorrPull25PFCor", &colorCorrPull25PFCor, &b_colorCorrPull25PFCor);
	fChain->SetBranchAddress("colorCorrPull35PFCor", &colorCorrPull35PFCor, &b_colorCorrPull35PFCor);
	fChain->SetBranchAddress("colorCorrPull45PFCor", &colorCorrPull45PFCor, &b_colorCorrPull45PFCor);
	fChain->SetBranchAddress("cosThetaJ1HiggsCM_PFCor", &cosThetaJ1HiggsCM_PFCor, &b_cosThetaJ1HiggsCM_PFCor);
	fChain->SetBranchAddress("cosThetaJ2HiggsCM_PFCor", &cosThetaJ2HiggsCM_PFCor, &b_cosThetaJ2HiggsCM_PFCor);
	fChain->SetBranchAddress("cosThetaL1HiggsCM_PFCor", &cosThetaL1HiggsCM_PFCor, &b_cosThetaL1HiggsCM_PFCor);
	fChain->SetBranchAddress("cosThetaL2HiggsCM_PFCor", &cosThetaL2HiggsCM_PFCor, &b_cosThetaL2HiggsCM_PFCor);
	fChain->SetBranchAddress("cosThetaV1HiggsCM_PFCor", &cosThetaV1HiggsCM_PFCor, &b_cosThetaV1HiggsCM_PFCor);
	fChain->SetBranchAddress("cosThetaV2HiggsCM_PFCor", &cosThetaV2HiggsCM_PFCor, &b_cosThetaV2HiggsCM_PFCor);
	fChain->SetBranchAddress("JetPFCor_isPileUpJetLoose", JetPFCor_isPileUpJetLoose, &b_JetPFCor_isPileUpJetLoose);
	fChain->SetBranchAddress("JetPFCor_isPileUpJetMedium", JetPFCor_isPileUpJetMedium, &b_JetPFCor_isPileUpJetMedium);
	fChain->SetBranchAddress("JetPFCor_isPileUpJetTight", JetPFCor_isPileUpJetTight, &b_JetPFCor_isPileUpJetTight);
	fChain->SetBranchAddress("GroomedJet_CA8_pt_uncorr", GroomedJet_CA8_pt_uncorr, &b_GroomedJet_CA8_pt_uncorr);
	fChain->SetBranchAddress("GroomedJet_CA8_mass_uncorr", GroomedJet_CA8_mass_uncorr, &b_GroomedJet_CA8_mass_uncorr);
	fChain->SetBranchAddress("GroomedJet_CA8_mass_tr_uncorr", GroomedJet_CA8_mass_tr_uncorr, &b_GroomedJet_CA8_mass_tr_uncorr);
	fChain->SetBranchAddress("GroomedJet_CA8_mass_ft_uncorr", GroomedJet_CA8_mass_ft_uncorr, &b_GroomedJet_CA8_mass_ft_uncorr);
	fChain->SetBranchAddress("GroomedJet_CA8_mass_pr_uncorr", GroomedJet_CA8_mass_pr_uncorr, &b_GroomedJet_CA8_mass_pr_uncorr);
	fChain->SetBranchAddress("GroomedJet_CA8_tau2tau1", GroomedJet_CA8_tau2tau1, &b_GroomedJet_CA8_tau2tau1);
	fChain->SetBranchAddress("GroomedJet_CA8_tau1", GroomedJet_CA8_tau1, &b_GroomedJet_CA8_tau1);
	fChain->SetBranchAddress("GroomedJet_CA8_tau2", GroomedJet_CA8_tau2, &b_GroomedJet_CA8_tau2);
	fChain->SetBranchAddress("GroomedJet_CA8_tau3", GroomedJet_CA8_tau3, &b_GroomedJet_CA8_tau3);
	fChain->SetBranchAddress("GroomedJet_CA8_tau4", GroomedJet_CA8_tau4, &b_GroomedJet_CA8_tau4);
	fChain->SetBranchAddress("GroomedJet_CA8_massdrop_pr_uncorr", GroomedJet_CA8_massdrop_pr_uncorr, &b_GroomedJet_CA8_massdrop_pr_uncorr);
	fChain->SetBranchAddress("GroomedJet_CA8_pt", GroomedJet_CA8_pt, &b_GroomedJet_CA8_pt);
	fChain->SetBranchAddress("GroomedJet_CA8_eta", GroomedJet_CA8_eta, &b_GroomedJet_CA8_eta);
	fChain->SetBranchAddress("GroomedJet_CA8_phi", GroomedJet_CA8_phi, &b_GroomedJet_CA8_phi);
	fChain->SetBranchAddress("GroomedJet_CA8_e", GroomedJet_CA8_e, &b_GroomedJet_CA8_e);
	fChain->SetBranchAddress("GroomedJet_CA8_pt_tr_uncorr", GroomedJet_CA8_pt_tr_uncorr, &b_GroomedJet_CA8_pt_tr_uncorr);
	fChain->SetBranchAddress("GroomedJet_CA8_pt_tr", GroomedJet_CA8_pt_tr, &b_GroomedJet_CA8_pt_tr);
	fChain->SetBranchAddress("GroomedJet_CA8_eta_tr", GroomedJet_CA8_eta_tr, &b_GroomedJet_CA8_eta_tr);
	fChain->SetBranchAddress("GroomedJet_CA8_phi_tr", GroomedJet_CA8_phi_tr, &b_GroomedJet_CA8_phi_tr);
	fChain->SetBranchAddress("GroomedJet_CA8_e_tr", GroomedJet_CA8_e_tr, &b_GroomedJet_CA8_e_tr);
	fChain->SetBranchAddress("GroomedJet_CA8_pt_ft_uncorr", GroomedJet_CA8_pt_ft_uncorr, &b_GroomedJet_CA8_pt_ft_uncorr);
	fChain->SetBranchAddress("GroomedJet_CA8_pt_ft", GroomedJet_CA8_pt_ft, &b_GroomedJet_CA8_pt_ft);
	fChain->SetBranchAddress("GroomedJet_CA8_eta_ft", GroomedJet_CA8_eta_ft, &b_GroomedJet_CA8_eta_ft);
	fChain->SetBranchAddress("GroomedJet_CA8_phi_ft", GroomedJet_CA8_phi_ft, &b_GroomedJet_CA8_phi_ft);
	fChain->SetBranchAddress("GroomedJet_CA8_e_ft", GroomedJet_CA8_e_ft, &b_GroomedJet_CA8_e_ft);
	fChain->SetBranchAddress("GroomedJet_CA8_pt_pr_uncorr", GroomedJet_CA8_pt_pr_uncorr, &b_GroomedJet_CA8_pt_pr_uncorr);
	fChain->SetBranchAddress("GroomedJet_CA8_pt_pr", GroomedJet_CA8_pt_pr, &b_GroomedJet_CA8_pt_pr);
	fChain->SetBranchAddress("GroomedJet_CA8_eta_pr", GroomedJet_CA8_eta_pr, &b_GroomedJet_CA8_eta_pr);
	fChain->SetBranchAddress("GroomedJet_CA8_phi_pr", GroomedJet_CA8_phi_pr, &b_GroomedJet_CA8_phi_pr);
	fChain->SetBranchAddress("GroomedJet_CA8_e_pr", GroomedJet_CA8_e_pr, &b_GroomedJet_CA8_e_pr);
	fChain->SetBranchAddress("GroomedJet_CA8_prsubjet1_px", GroomedJet_CA8_prsubjet1_px, &b_GroomedJet_CA8_prsubjet1_px);
	fChain->SetBranchAddress("GroomedJet_CA8_prsubjet1_py", GroomedJet_CA8_prsubjet1_py, &b_GroomedJet_CA8_prsubjet1_py);
	fChain->SetBranchAddress("GroomedJet_CA8_prsubjet1_pz", GroomedJet_CA8_prsubjet1_pz, &b_GroomedJet_CA8_prsubjet1_pz);
	fChain->SetBranchAddress("GroomedJet_CA8_prsubjet1_e", GroomedJet_CA8_prsubjet1_e, &b_GroomedJet_CA8_prsubjet1_e);
	fChain->SetBranchAddress("GroomedJet_CA8_prsubjet2_px", GroomedJet_CA8_prsubjet2_px, &b_GroomedJet_CA8_prsubjet2_px);
	fChain->SetBranchAddress("GroomedJet_CA8_prsubjet2_py", GroomedJet_CA8_prsubjet2_py, &b_GroomedJet_CA8_prsubjet2_py);
	fChain->SetBranchAddress("GroomedJet_CA8_prsubjet2_pz", GroomedJet_CA8_prsubjet2_pz, &b_GroomedJet_CA8_prsubjet2_pz);
	fChain->SetBranchAddress("GroomedJet_CA8_prsubjet2_e", GroomedJet_CA8_prsubjet2_e, &b_GroomedJet_CA8_prsubjet2_e);
	fChain->SetBranchAddress("GroomedJet_CA8_mass", GroomedJet_CA8_mass, &b_GroomedJet_CA8_mass);
	fChain->SetBranchAddress("GroomedJet_CA8_mass_tr", GroomedJet_CA8_mass_tr, &b_GroomedJet_CA8_mass_tr);
	fChain->SetBranchAddress("GroomedJet_CA8_mass_ft", GroomedJet_CA8_mass_ft, &b_GroomedJet_CA8_mass_ft);
	fChain->SetBranchAddress("GroomedJet_CA8_mass_pr", GroomedJet_CA8_mass_pr, &b_GroomedJet_CA8_mass_pr);
	fChain->SetBranchAddress("GroomedJet_CA8_massdrop_pr", GroomedJet_CA8_massdrop_pr, &b_GroomedJet_CA8_massdrop_pr);
	fChain->SetBranchAddress("GroomedJet_CA8_area", GroomedJet_CA8_area, &b_GroomedJet_CA8_area);
	fChain->SetBranchAddress("GroomedJet_CA8_area_tr", GroomedJet_CA8_area_tr, &b_GroomedJet_CA8_area_tr);
	fChain->SetBranchAddress("GroomedJet_CA8_area_ft", GroomedJet_CA8_area_ft, &b_GroomedJet_CA8_area_ft);
	fChain->SetBranchAddress("GroomedJet_CA8_area_pr", GroomedJet_CA8_area_pr, &b_GroomedJet_CA8_area_pr);
	fChain->SetBranchAddress("GroomedJet_CA8_jetconstituents", GroomedJet_CA8_jetconstituents, &b_GroomedJet_CA8_jetconstituents);
	fChain->SetBranchAddress("GroomedJet_CA8_jetcharge", GroomedJet_CA8_jetcharge, &b_GroomedJet_CA8_jetcharge);
	fChain->SetBranchAddress("GroomedJet_CA8_rcores", GroomedJet_CA8_rcores, &b_GroomedJet_CA8_rcores);
	fChain->SetBranchAddress("GroomedJet_CA8_ptcores", GroomedJet_CA8_ptcores, &b_GroomedJet_CA8_ptcores);
	fChain->SetBranchAddress("GroomedJet_CA8_planarflow", GroomedJet_CA8_planarflow, &b_GroomedJet_CA8_planarflow);
	fChain->SetBranchAddress("GroomedJet_CA8_qjetmass", GroomedJet_CA8_qjetmass, &b_GroomedJet_CA8_qjetmass);
	fChain->SetBranchAddress("GroomedJet_CA8_qjetmassdrop", GroomedJet_CA8_qjetmassdrop, &b_GroomedJet_CA8_qjetmassdrop);
	fChain->SetBranchAddress("GroomedJet_CA8_constituents0_eta", GroomedJet_CA8_constituents0_eta, &b_GroomedJet_CA8_constituents0_eta);
	fChain->SetBranchAddress("GroomedJet_CA8_constituents0_phi", GroomedJet_CA8_constituents0_phi, &b_GroomedJet_CA8_constituents0_phi);
	fChain->SetBranchAddress("GroomedJet_CA8_constituents0_e", GroomedJet_CA8_constituents0_e, &b_GroomedJet_CA8_constituents0_e);
	fChain->SetBranchAddress("GroomedJet_CA8_nconstituents0", &GroomedJet_CA8_nconstituents0, &b_GroomedJet_CA8_nconstituents0);
	fChain->SetBranchAddress("GroomedJet_CA8_constituents0pr_eta", GroomedJet_CA8_constituents0pr_eta, &b_GroomedJet_CA8_constituents0pr_eta);
	fChain->SetBranchAddress("GroomedJet_CA8_constituents0pr_phi", GroomedJet_CA8_constituents0pr_phi, &b_GroomedJet_CA8_constituents0pr_phi);
	fChain->SetBranchAddress("GroomedJet_CA8_constituents0pr_e", GroomedJet_CA8_constituents0pr_e, &b_GroomedJet_CA8_constituents0pr_e);
	fChain->SetBranchAddress("GroomedJet_CA8_nconstituents0pr", &GroomedJet_CA8_nconstituents0pr, &b_GroomedJet_CA8_nconstituents0pr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_pt_uncorr", GenGroomedJet_CA8_pt_uncorr, &b_GenGroomedJet_CA8_pt_uncorr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_mass_uncorr", GenGroomedJet_CA8_mass_uncorr, &b_GenGroomedJet_CA8_mass_uncorr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_mass_tr_uncorr", GenGroomedJet_CA8_mass_tr_uncorr, &b_GenGroomedJet_CA8_mass_tr_uncorr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_mass_ft_uncorr", GenGroomedJet_CA8_mass_ft_uncorr, &b_GenGroomedJet_CA8_mass_ft_uncorr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_mass_pr_uncorr", GenGroomedJet_CA8_mass_pr_uncorr, &b_GenGroomedJet_CA8_mass_pr_uncorr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_tau2tau1", GenGroomedJet_CA8_tau2tau1, &b_GenGroomedJet_CA8_tau2tau1);
	fChain->SetBranchAddress("GenGroomedJet_CA8_tau1", GenGroomedJet_CA8_tau1, &b_GenGroomedJet_CA8_tau1);
	fChain->SetBranchAddress("GenGroomedJet_CA8_tau2", GenGroomedJet_CA8_tau2, &b_GenGroomedJet_CA8_tau2);
	fChain->SetBranchAddress("GenGroomedJet_CA8_tau3", GenGroomedJet_CA8_tau3, &b_GenGroomedJet_CA8_tau3);
	fChain->SetBranchAddress("GenGroomedJet_CA8_tau4", GenGroomedJet_CA8_tau4, &b_GenGroomedJet_CA8_tau4);
	fChain->SetBranchAddress("GenGroomedJet_CA8_massdrop_pr_uncorr", GenGroomedJet_CA8_massdrop_pr_uncorr, &b_GenGroomedJet_CA8_massdrop_pr_uncorr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_pt", GenGroomedJet_CA8_pt, &b_GenGroomedJet_CA8_pt);
	fChain->SetBranchAddress("GenGroomedJet_CA8_eta", GenGroomedJet_CA8_eta, &b_GenGroomedJet_CA8_eta);
	fChain->SetBranchAddress("GenGroomedJet_CA8_phi", GenGroomedJet_CA8_phi, &b_GenGroomedJet_CA8_phi);
	fChain->SetBranchAddress("GenGroomedJet_CA8_e", GenGroomedJet_CA8_e, &b_GenGroomedJet_CA8_e);
	fChain->SetBranchAddress("GenGroomedJet_CA8_pt_tr_uncorr", GenGroomedJet_CA8_pt_tr_uncorr, &b_GenGroomedJet_CA8_pt_tr_uncorr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_pt_tr", GenGroomedJet_CA8_pt_tr, &b_GenGroomedJet_CA8_pt_tr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_eta_tr", GenGroomedJet_CA8_eta_tr, &b_GenGroomedJet_CA8_eta_tr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_phi_tr", GenGroomedJet_CA8_phi_tr, &b_GenGroomedJet_CA8_phi_tr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_e_tr", GenGroomedJet_CA8_e_tr, &b_GenGroomedJet_CA8_e_tr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_pt_ft_uncorr", GenGroomedJet_CA8_pt_ft_uncorr, &b_GenGroomedJet_CA8_pt_ft_uncorr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_pt_ft", GenGroomedJet_CA8_pt_ft, &b_GenGroomedJet_CA8_pt_ft);
	fChain->SetBranchAddress("GenGroomedJet_CA8_eta_ft", GenGroomedJet_CA8_eta_ft, &b_GenGroomedJet_CA8_eta_ft);
	fChain->SetBranchAddress("GenGroomedJet_CA8_phi_ft", GenGroomedJet_CA8_phi_ft, &b_GenGroomedJet_CA8_phi_ft);
	fChain->SetBranchAddress("GenGroomedJet_CA8_e_ft", GenGroomedJet_CA8_e_ft, &b_GenGroomedJet_CA8_e_ft);
	fChain->SetBranchAddress("GenGroomedJet_CA8_pt_pr_uncorr", GenGroomedJet_CA8_pt_pr_uncorr, &b_GenGroomedJet_CA8_pt_pr_uncorr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_pt_pr", GenGroomedJet_CA8_pt_pr, &b_GenGroomedJet_CA8_pt_pr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_eta_pr", GenGroomedJet_CA8_eta_pr, &b_GenGroomedJet_CA8_eta_pr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_phi_pr", GenGroomedJet_CA8_phi_pr, &b_GenGroomedJet_CA8_phi_pr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_e_pr", GenGroomedJet_CA8_e_pr, &b_GenGroomedJet_CA8_e_pr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_prsubjet1_px", GenGroomedJet_CA8_prsubjet1_px, &b_GenGroomedJet_CA8_prsubjet1_px);
	fChain->SetBranchAddress("GenGroomedJet_CA8_prsubjet1_py", GenGroomedJet_CA8_prsubjet1_py, &b_GenGroomedJet_CA8_prsubjet1_py);
	fChain->SetBranchAddress("GenGroomedJet_CA8_prsubjet1_pz", GenGroomedJet_CA8_prsubjet1_pz, &b_GenGroomedJet_CA8_prsubjet1_pz);
	fChain->SetBranchAddress("GenGroomedJet_CA8_prsubjet1_e", GenGroomedJet_CA8_prsubjet1_e, &b_GenGroomedJet_CA8_prsubjet1_e);
	fChain->SetBranchAddress("GenGroomedJet_CA8_prsubjet2_px", GenGroomedJet_CA8_prsubjet2_px, &b_GenGroomedJet_CA8_prsubjet2_px);
	fChain->SetBranchAddress("GenGroomedJet_CA8_prsubjet2_py", GenGroomedJet_CA8_prsubjet2_py, &b_GenGroomedJet_CA8_prsubjet2_py);
	fChain->SetBranchAddress("GenGroomedJet_CA8_prsubjet2_pz", GenGroomedJet_CA8_prsubjet2_pz, &b_GenGroomedJet_CA8_prsubjet2_pz);
	fChain->SetBranchAddress("GenGroomedJet_CA8_prsubjet2_e", GenGroomedJet_CA8_prsubjet2_e, &b_GenGroomedJet_CA8_prsubjet2_e);
	fChain->SetBranchAddress("GenGroomedJet_CA8_mass", GenGroomedJet_CA8_mass, &b_GenGroomedJet_CA8_mass);
	fChain->SetBranchAddress("GenGroomedJet_CA8_mass_tr", GenGroomedJet_CA8_mass_tr, &b_GenGroomedJet_CA8_mass_tr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_mass_ft", GenGroomedJet_CA8_mass_ft, &b_GenGroomedJet_CA8_mass_ft);
	fChain->SetBranchAddress("GenGroomedJet_CA8_mass_pr", GenGroomedJet_CA8_mass_pr, &b_GenGroomedJet_CA8_mass_pr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_massdrop_pr", GenGroomedJet_CA8_massdrop_pr, &b_GenGroomedJet_CA8_massdrop_pr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_area", GenGroomedJet_CA8_area, &b_GenGroomedJet_CA8_area);
	fChain->SetBranchAddress("GenGroomedJet_CA8_area_tr", GenGroomedJet_CA8_area_tr, &b_GenGroomedJet_CA8_area_tr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_area_ft", GenGroomedJet_CA8_area_ft, &b_GenGroomedJet_CA8_area_ft);
	fChain->SetBranchAddress("GenGroomedJet_CA8_area_pr", GenGroomedJet_CA8_area_pr, &b_GenGroomedJet_CA8_area_pr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_jetconstituents", GenGroomedJet_CA8_jetconstituents, &b_GenGroomedJet_CA8_jetconstituents);
	fChain->SetBranchAddress("GenGroomedJet_CA8_jetcharge", GenGroomedJet_CA8_jetcharge, &b_GenGroomedJet_CA8_jetcharge);
	fChain->SetBranchAddress("GenGroomedJet_CA8_rcores", GenGroomedJet_CA8_rcores, &b_GenGroomedJet_CA8_rcores);
	fChain->SetBranchAddress("GenGroomedJet_CA8_ptcores", GenGroomedJet_CA8_ptcores, &b_GenGroomedJet_CA8_ptcores);
	fChain->SetBranchAddress("GenGroomedJet_CA8_planarflow", GenGroomedJet_CA8_planarflow, &b_GenGroomedJet_CA8_planarflow);
	fChain->SetBranchAddress("GenGroomedJet_CA8_qjetmass", GenGroomedJet_CA8_qjetmass, &b_GenGroomedJet_CA8_qjetmass);
	fChain->SetBranchAddress("GenGroomedJet_CA8_qjetmassdrop", GenGroomedJet_CA8_qjetmassdrop, &b_GenGroomedJet_CA8_qjetmassdrop);
	fChain->SetBranchAddress("GenGroomedJet_CA8_constituents0_eta", GenGroomedJet_CA8_constituents0_eta, &b_GenGroomedJet_CA8_constituents0_eta);
	fChain->SetBranchAddress("GenGroomedJet_CA8_constituents0_phi", GenGroomedJet_CA8_constituents0_phi, &b_GenGroomedJet_CA8_constituents0_phi);
	fChain->SetBranchAddress("GenGroomedJet_CA8_constituents0_e", GenGroomedJet_CA8_constituents0_e, &b_GenGroomedJet_CA8_constituents0_e);
	fChain->SetBranchAddress("GenGroomedJet_CA8_nconstituents0", &GenGroomedJet_CA8_nconstituents0, &b_GenGroomedJet_CA8_nconstituents0);
	fChain->SetBranchAddress("GenGroomedJet_CA8_constituents0pr_eta", GenGroomedJet_CA8_constituents0pr_eta, &b_GenGroomedJet_CA8_constituents0pr_eta);
	fChain->SetBranchAddress("GenGroomedJet_CA8_constituents0pr_phi", GenGroomedJet_CA8_constituents0pr_phi, &b_GenGroomedJet_CA8_constituents0pr_phi);
	fChain->SetBranchAddress("GenGroomedJet_CA8_constituents0pr_e", GenGroomedJet_CA8_constituents0pr_e, &b_GenGroomedJet_CA8_constituents0pr_e);
	fChain->SetBranchAddress("GenGroomedJet_CA8_nconstituents0pr", &GenGroomedJet_CA8_nconstituents0pr, &b_GenGroomedJet_CA8_nconstituents0pr);
	fChain->SetBranchAddress("numGenJets", &numGenJets, &b_numGenJets);
	fChain->SetBranchAddress("numGenJetBTags", &numGenJetBTags, &b_numGenJetBTags);
	fChain->SetBranchAddress("JetGen_Et", JetGen_Et, &b_JetGen_Et);
	fChain->SetBranchAddress("JetGen_Pt", JetGen_Pt, &b_JetGen_Pt);
	fChain->SetBranchAddress("JetGen_Eta", JetGen_Eta, &b_JetGen_Eta);
	fChain->SetBranchAddress("JetGen_Phi", JetGen_Phi, &b_JetGen_Phi);
	fChain->SetBranchAddress("JetGen_Theta", JetGen_Theta, &b_JetGen_Theta);
	fChain->SetBranchAddress("JetGen_Px", JetGen_Px, &b_JetGen_Px);
	fChain->SetBranchAddress("JetGen_Py", JetGen_Py, &b_JetGen_Py);
	fChain->SetBranchAddress("JetGen_Pz", JetGen_Pz, &b_JetGen_Pz);
	fChain->SetBranchAddress("JetGen_E", JetGen_E, &b_JetGen_E);
	fChain->SetBranchAddress("JetGen_Y", JetGen_Y, &b_JetGen_Y);
	fChain->SetBranchAddress("JetGen_Mass", JetGen_Mass, &b_JetGen_Mass);
	fChain->SetBranchAddress("JetGen_etaetaMoment", JetGen_etaetaMoment, &b_JetGen_etaetaMoment);
	fChain->SetBranchAddress("JetGen_phiphiMoment", JetGen_phiphiMoment, &b_JetGen_phiphiMoment);
	fChain->SetBranchAddress("JetGen_etaphiMoment", JetGen_etaphiMoment, &b_JetGen_etaphiMoment);
	fChain->SetBranchAddress("JetGen_maxDistance", JetGen_maxDistance, &b_JetGen_maxDistance);
	fChain->SetBranchAddress("JetGen_nConstituents", JetGen_nConstituents, &b_JetGen_nConstituents);
	fChain->SetBranchAddress("JetGen_Area", JetGen_Area, &b_JetGen_Area);
	fChain->SetBranchAddress("VplusGenJet_Mass", VplusGenJet_Mass, &b_VplusGenJet_Mass);
	fChain->SetBranchAddress("JetGen_dphiBoson", JetGen_dphiBoson, &b_JetGen_dphiBoson);
	fChain->SetBranchAddress("JetGen_detaBoson", JetGen_detaBoson, &b_JetGen_detaBoson);
	fChain->SetBranchAddress("JetGen_dRBoson", JetGen_dRBoson, &b_JetGen_dRBoson);
	fChain->SetBranchAddress("JetGen_dphiMET", JetGen_dphiMET, &b_JetGen_dphiMET);
	fChain->SetBranchAddress("JetGen_bDiscriminator", JetGen_bDiscriminator, &b_JetGen_bDiscriminator);
	fChain->SetBranchAddress("JetGen_bDiscriminatorSSVHE", JetGen_bDiscriminatorSSVHE, &b_JetGen_bDiscriminatorSSVHE);
	fChain->SetBranchAddress("JetGen_bDiscriminatorTCHE", JetGen_bDiscriminatorTCHE, &b_JetGen_bDiscriminatorTCHE);
	fChain->SetBranchAddress("JetGen_bDiscriminatorCSV", JetGen_bDiscriminatorCSV, &b_JetGen_bDiscriminatorCSV);
	fChain->SetBranchAddress("JetGen_bDiscriminatorJP", JetGen_bDiscriminatorJP, &b_JetGen_bDiscriminatorJP);
	fChain->SetBranchAddress("JetGen_bDiscriminatorSSVHP", JetGen_bDiscriminatorSSVHP, &b_JetGen_bDiscriminatorSSVHP);
	fChain->SetBranchAddress("JetGen_bDiscriminatorTCHP", JetGen_bDiscriminatorTCHP, &b_JetGen_bDiscriminatorTCHP);
	fChain->SetBranchAddress("JetGen_secVertexMass", JetGen_secVertexMass, &b_JetGen_secVertexMass);
	fChain->SetBranchAddress("MassV2j_Gen_MVAMET", &MassV2j_Gen_MVAMET, &b_MassV2j_Gen_MVAMET);
	fChain->SetBranchAddress("MassV2j_Gen", &MassV2j_Gen, &b_MassV2j_Gen);
	fChain->SetBranchAddress("MassV3j_Gen", &MassV3j_Gen, &b_MassV3j_Gen);
	fChain->SetBranchAddress("MassV4j_Gen", &MassV4j_Gen, &b_MassV4j_Gen);
	fChain->SetBranchAddress("MassV5j_Gen", &MassV5j_Gen, &b_MassV5j_Gen);
	fChain->SetBranchAddress("MassV6j_Gen", &MassV6j_Gen, &b_MassV6j_Gen);
	fChain->SetBranchAddress("Mass2j_Gen", &Mass2j_Gen, &b_Mass2j_Gen);
	fChain->SetBranchAddress("Mass3j_Gen", &Mass3j_Gen, &b_Mass3j_Gen);
	fChain->SetBranchAddress("Mass4j_Gen", &Mass4j_Gen, &b_Mass4j_Gen);
	fChain->SetBranchAddress("Mass5j_Gen", &Mass5j_Gen, &b_Mass5j_Gen);
	fChain->SetBranchAddress("Mass6j_Gen", &Mass6j_Gen, &b_Mass6j_Gen);
	fChain->SetBranchAddress("cosJacksonAngleV2j_Gen", &cosJacksonAngleV2j_Gen, &b_cosJacksonAngleV2j_Gen);
	fChain->SetBranchAddress("cosJacksonAngle2j_Gen", &cosJacksonAngle2j_Gen, &b_cosJacksonAngle2j_Gen);
	fChain->SetBranchAddress("cosJacksonAngleV3j_Gen", &cosJacksonAngleV3j_Gen, &b_cosJacksonAngleV3j_Gen);
	fChain->SetBranchAddress("cosJacksonAngle3j12_Gen", &cosJacksonAngle3j12_Gen, &b_cosJacksonAngle3j12_Gen);
	fChain->SetBranchAddress("cosJacksonAngle3j23_Gen", &cosJacksonAngle3j23_Gen, &b_cosJacksonAngle3j23_Gen);
	fChain->SetBranchAddress("cosJacksonAngle3j31_Gen", &cosJacksonAngle3j31_Gen, &b_cosJacksonAngle3j31_Gen);
	fChain->SetBranchAddress("cosphiDecayPlane_Gen", &cosphiDecayPlane_Gen, &b_cosphiDecayPlane_Gen);
	fChain->SetBranchAddress("cosThetaLnu_Gen", &cosThetaLnu_Gen, &b_cosThetaLnu_Gen);
	fChain->SetBranchAddress("cosThetaJJ_Gen", &cosThetaJJ_Gen, &b_cosThetaJJ_Gen);
	fChain->SetBranchAddress("NumPhotons", &NumPhotons, &b_NumPhotons);
	fChain->SetBranchAddress("Photon_Et", Photon_Et, &b_Photon_Et);
	fChain->SetBranchAddress("Photon_E", Photon_E, &b_Photon_E);
	fChain->SetBranchAddress("Photon_Eta", Photon_Eta, &b_Photon_Eta);
	fChain->SetBranchAddress("Photon_Phi", Photon_Phi, &b_Photon_Phi);
	fChain->SetBranchAddress("Photon_Theta", Photon_Theta, &b_Photon_Theta);
	fChain->SetBranchAddress("Photon_Px", Photon_Px, &b_Photon_Px);
	fChain->SetBranchAddress("Photon_Py", Photon_Py, &b_Photon_Py);
	fChain->SetBranchAddress("Photon_Pz", Photon_Pz, &b_Photon_Pz);
	fChain->SetBranchAddress("Photon_Vx", Photon_Vx, &b_Photon_Vx);
	fChain->SetBranchAddress("Photon_Vy", Photon_Vy, &b_Photon_Vy);
	fChain->SetBranchAddress("Photon_Vz", Photon_Vz, &b_Photon_Vz);
	fChain->SetBranchAddress("Photon_SC_Et", Photon_SC_Et, &b_Photon_SC_Et);
	fChain->SetBranchAddress("Photon_SC_E", Photon_SC_E, &b_Photon_SC_E);
	fChain->SetBranchAddress("Photon_SC_Eta", Photon_SC_Eta, &b_Photon_SC_Eta);
	fChain->SetBranchAddress("Photon_SC_Phi", Photon_SC_Phi, &b_Photon_SC_Phi);
	fChain->SetBranchAddress("Photon_SC_Theta", Photon_SC_Theta, &b_Photon_SC_Theta);
	fChain->SetBranchAddress("Photon_SC_x", Photon_SC_x, &b_Photon_SC_x);
	fChain->SetBranchAddress("Photon_SC_y", Photon_SC_y, &b_Photon_SC_y);
	fChain->SetBranchAddress("Photon_SC_z", Photon_SC_z, &b_Photon_SC_z);
	fChain->SetBranchAddress("PFisocharged03", PFisocharged03, &b_PFisocharged03);
	fChain->SetBranchAddress("PFisophoton03", PFisophoton03, &b_PFisophoton03);
	fChain->SetBranchAddress("PFisoneutral03", PFisoneutral03, &b_PFisoneutral03);
	fChain->SetBranchAddress("trkSumPtHollowConeDR04_Photon11", trkSumPtHollowConeDR04_Photon11, &b_trkSumPtHollowConeDR04_Photon11);
	fChain->SetBranchAddress("ecalRecHitSumEtConeDR04_Photon11", ecalRecHitSumEtConeDR04_Photon11, &b_ecalRecHitSumEtConeDR04_Photon11);
	fChain->SetBranchAddress("hcalTowerSumEtConeDR04_Photon11", hcalTowerSumEtConeDR04_Photon11, &b_hcalTowerSumEtConeDR04_Photon11);
	fChain->SetBranchAddress("Photon_HoverE", Photon_HoverE, &b_Photon_HoverE);
	fChain->SetBranchAddress("Photon_HoverE2011", Photon_HoverE2011, &b_Photon_HoverE2011);
	fChain->SetBranchAddress("Photon_SigmaIetaIeta", Photon_SigmaIetaIeta, &b_Photon_SigmaIetaIeta);
	fChain->SetBranchAddress("Photon_hasPixelSeed", Photon_hasPixelSeed, &b_Photon_hasPixelSeed);
	fChain->SetBranchAddress("Photon_passElecVeto", Photon_passElecVeto, &b_Photon_passElecVeto);
	fChain->SetBranchAddress("Photon_Id2011", Photon_Id2011, &b_Photon_Id2011);
	fChain->SetBranchAddress("Photon_Id2012", Photon_Id2012, &b_Photon_Id2012);
	fChain->SetBranchAddress("Z_mass", &Z_mass, &b_Z_mass);
	fChain->SetBranchAddress("Z_mt", &Z_mt, &b_Z_mt);
	fChain->SetBranchAddress("Z_mtMVA", &Z_mtMVA, &b_Z_mtMVA);
	fChain->SetBranchAddress("Z_px", &Z_px, &b_Z_px);
	fChain->SetBranchAddress("Z_py", &Z_py, &b_Z_py);
	fChain->SetBranchAddress("Z_pz", &Z_pz, &b_Z_pz);
	fChain->SetBranchAddress("Z_e", &Z_e, &b_Z_e);
	fChain->SetBranchAddress("Z_pt", &Z_pt, &b_Z_pt);
	fChain->SetBranchAddress("Z_et", &Z_et, &b_Z_et);
	fChain->SetBranchAddress("Z_eta", &Z_eta, &b_Z_eta);
	fChain->SetBranchAddress("Z_phi", &Z_phi, &b_Z_phi);
	fChain->SetBranchAddress("Z_vx", &Z_vx, &b_Z_vx);
	fChain->SetBranchAddress("Z_vy", &Z_vy, &b_Z_vy);
	fChain->SetBranchAddress("Z_vz", &Z_vz, &b_Z_vz);
	fChain->SetBranchAddress("Z_y", &Z_y, &b_Z_y);
	fChain->SetBranchAddress("Z_muplus_px", &Z_muplus_px, &b_Z_muplus_px);
	fChain->SetBranchAddress("Z_muplus_py", &Z_muplus_py, &b_Z_muplus_py);
	fChain->SetBranchAddress("Z_muplus_pz", &Z_muplus_pz, &b_Z_muplus_pz);
	fChain->SetBranchAddress("Z_muplus_e", &Z_muplus_e, &b_Z_muplus_e);
	fChain->SetBranchAddress("Z_muplus_pt", &Z_muplus_pt, &b_Z_muplus_pt);
	fChain->SetBranchAddress("Z_muplus_et", &Z_muplus_et, &b_Z_muplus_et);
	fChain->SetBranchAddress("Z_muplus_eta", &Z_muplus_eta, &b_Z_muplus_eta);
	fChain->SetBranchAddress("Z_muplus_theta", &Z_muplus_theta, &b_Z_muplus_theta);
	fChain->SetBranchAddress("Z_muplus_phi", &Z_muplus_phi, &b_Z_muplus_phi);
	fChain->SetBranchAddress("Z_muplus_charge", &Z_muplus_charge, &b_Z_muplus_charge);
	fChain->SetBranchAddress("Z_muplus_vx", &Z_muplus_vx, &b_Z_muplus_vx);
	fChain->SetBranchAddress("Z_muplus_vy", &Z_muplus_vy, &b_Z_muplus_vy);
	fChain->SetBranchAddress("Z_muplus_vz", &Z_muplus_vz, &b_Z_muplus_vz);
	fChain->SetBranchAddress("Z_muplus_y", &Z_muplus_y, &b_Z_muplus_y);
	fChain->SetBranchAddress("Z_muplus_trackiso", &Z_muplus_trackiso, &b_Z_muplus_trackiso);
	fChain->SetBranchAddress("Z_muplus_hcaliso", &Z_muplus_hcaliso, &b_Z_muplus_hcaliso);
	fChain->SetBranchAddress("Z_muplus_ecaliso", &Z_muplus_ecaliso, &b_Z_muplus_ecaliso);
	fChain->SetBranchAddress("Z_muplus_type", &Z_muplus_type, &b_Z_muplus_type);
	fChain->SetBranchAddress("Z_muplus_numberOfChambers", &Z_muplus_numberOfChambers, &b_Z_muplus_numberOfChambers);
	fChain->SetBranchAddress("Z_muplus_numberOfMatches", &Z_muplus_numberOfMatches, &b_Z_muplus_numberOfMatches);
	fChain->SetBranchAddress("Z_muplus_d0bsp", &Z_muplus_d0bsp, &b_Z_muplus_d0bsp);
	fChain->SetBranchAddress("Z_muplus_dz000", &Z_muplus_dz000, &b_Z_muplus_dz000);
	fChain->SetBranchAddress("Z_muplus_dzPV", &Z_muplus_dzPV, &b_Z_muplus_dzPV);
	fChain->SetBranchAddress("Z_muplus_pfiso_sumChargedHadronPt", &Z_muplus_pfiso_sumChargedHadronPt, &b_Z_muplus_pfiso_sumChargedHadronPt);
	fChain->SetBranchAddress("Z_muplus_pfiso_sumChargedParticlePt", &Z_muplus_pfiso_sumChargedParticlePt, &b_Z_muplus_pfiso_sumChargedParticlePt);
	fChain->SetBranchAddress("Z_muplus_pfiso_sumNeutralHadronEt", &Z_muplus_pfiso_sumNeutralHadronEt, &b_Z_muplus_pfiso_sumNeutralHadronEt);
	fChain->SetBranchAddress("Z_muplus_pfiso_sumPhotonEt", &Z_muplus_pfiso_sumPhotonEt, &b_Z_muplus_pfiso_sumPhotonEt);
	fChain->SetBranchAddress("Z_muplus_pfiso_sumPUPt", &Z_muplus_pfiso_sumPUPt, &b_Z_muplus_pfiso_sumPUPt);

	fChain->SetBranchAddress("Z_muminus_px",               &Z_muminus_px,               &b_Z_muminus_px);
	fChain->SetBranchAddress("Z_muminus_py",               &Z_muminus_py,               &b_Z_muminus_py);
	fChain->SetBranchAddress("Z_muminus_pz",               &Z_muminus_pz,               &b_Z_muminus_pz);
	fChain->SetBranchAddress("Z_muminus_e" ,               &Z_muminus_e,                &b_Z_muminus_e);
	fChain->SetBranchAddress("Z_muminus_pt",               &Z_muminus_pt,               &b_Z_muminus_pt);
	fChain->SetBranchAddress("Z_muminus_et",               &Z_muminus_et,               &b_Z_muminus_et);
	fChain->SetBranchAddress("Z_muminus_eta",              &Z_muminus_eta,              &b_Z_muminus_eta);
	fChain->SetBranchAddress("Z_muminus_theta",            &Z_muminus_theta,            &b_Z_muminus_theta);
	fChain->SetBranchAddress("Z_muminus_phi",              &Z_muminus_phi,              &b_Z_muminus_phi);
	fChain->SetBranchAddress("Z_muminus_charge",           &Z_muminus_charge,           &b_Z_muminus_charge);
	fChain->SetBranchAddress("Z_muminus_vx",               &Z_muminus_vx,               &b_Z_muminus_vx);
	fChain->SetBranchAddress("Z_muminus_vy",               &Z_muminus_vy,               &b_Z_muminus_vy);
	fChain->SetBranchAddress("Z_muminus_vz",               &Z_muminus_vz,               &b_Z_muminus_vz);
	fChain->SetBranchAddress("Z_muminus_y",                &Z_muminus_y,                &b_Z_muminus_y);
	fChain->SetBranchAddress("Z_muminus_trackiso",         &Z_muminus_trackiso,         &b_Z_muminus_trackiso);
	fChain->SetBranchAddress("Z_muminus_hcaliso",          &Z_muminus_hcaliso,          &b_Z_muminus_hcaliso);
	fChain->SetBranchAddress("Z_muminus_ecaliso",          &Z_muminus_ecaliso,          &b_Z_muminus_ecaliso);
	fChain->SetBranchAddress("Z_muminus_type",             &Z_muminus_type,             &b_Z_muminus_type);
	fChain->SetBranchAddress("Z_muminus_numberOfChambers", &Z_muminus_numberOfChambers, &b_Z_muminus_numberOfChambers);
	fChain->SetBranchAddress("Z_muminus_numberOfMatches",  &Z_muminus_numberOfMatches,  &b_Z_muminus_numberOfMatches);
	//fChain->SetBranchAddress("Z_muminus_d0bsp",            &Z_muminus_d0bsp,            &b_Z_muminus_d0bsp);
	//fChain->SetBranchAddress("Z_muminus_dz000",            &Z_muminus_dz000,            &b_Z_muminus_dz000);
	//fChain->SetBranchAddress("Z_muminus_dzPV",             &Z_muminus_dzPV,             &b_Z_muminus_dzPV);
	//fChain->SetBranchAddress("Z_muminus_pfiso_sumChargedHadronPt",   &Z_muminus_pfiso_sumChargedHadronPt,   &b_Z_muminus_pfiso_sumChargedHadronPt);
	//fChain->SetBranchAddress("Z_muminus_pfiso_sumChargedParticlePt", &Z_muminus_pfiso_sumChargedParticlePt, &b_Z_muminus_pfiso_sumChargedParticlePt);
	//fChain->SetBranchAddress("Z_muminus_pfiso_sumNeutralHadronEt",   &Z_muminus_pfiso_sumNeutralHadronEt,   &b_Z_muminus_pfiso_sumNeutralHadronEt);
	//fChain->SetBranchAddress("Z_muminus_pfiso_sumPhotonEt",          &Z_muminus_pfiso_sumPhotonEt,          &b_Z_muminus_pfiso_sumPhotonEt);
	//fChain->SetBranchAddress("Z_muminus_pfiso_sumPUPt",              &Z_muminus_pfiso_sumPUPt,              &b_Z_muminus_pfiso_sumPUPt);
	fChain->SetBranchAddress("Z_Photon_pt_gen", &Z_Photon_pt_gen, &b_Z_Photon_pt_gen);
	fChain->SetBranchAddress("Z_muplus_px_gen", &Z_muplus_px_gen, &b_Z_muplus_px_gen);
	fChain->SetBranchAddress("Z_muplus_py_gen", &Z_muplus_py_gen, &b_Z_muplus_py_gen);
	fChain->SetBranchAddress("Z_muplus_pz_gen", &Z_muplus_pz_gen, &b_Z_muplus_pz_gen);
	fChain->SetBranchAddress("Z_muplus_e_gen", &Z_muplus_e_gen, &b_Z_muplus_e_gen);
	fChain->SetBranchAddress("Z_muplus_pt_gen", &Z_muplus_pt_gen, &b_Z_muplus_pt_gen);
	fChain->SetBranchAddress("Z_muplus_et_gen", &Z_muplus_et_gen, &b_Z_muplus_et_gen);
	fChain->SetBranchAddress("Z_muplus_eta_gen", &Z_muplus_eta_gen, &b_Z_muplus_eta_gen);
	fChain->SetBranchAddress("Z_muplus_theta_gen", &Z_muplus_theta_gen, &b_Z_muplus_theta_gen);
	fChain->SetBranchAddress("Z_muplus_phi_gen", &Z_muplus_phi_gen, &b_Z_muplus_phi_gen);
	fChain->SetBranchAddress("Z_muplus_charge_gen", &Z_muplus_charge_gen, &b_Z_muplus_charge_gen);
	fChain->SetBranchAddress("Z_muplus_vx_gen", &Z_muplus_vx_gen, &b_Z_muplus_vx_gen);
	fChain->SetBranchAddress("Z_muplus_vy_gen", &Z_muplus_vy_gen, &b_Z_muplus_vy_gen);
	fChain->SetBranchAddress("Z_muplus_vz_gen", &Z_muplus_vz_gen, &b_Z_muplus_vz_gen);
	fChain->SetBranchAddress("Z_muplus_y_gen", &Z_muplus_y_gen, &b_Z_muplus_y_gen);
	fChain->SetBranchAddress("Z_muminus_px_gen", &Z_muminus_px_gen, &b_Z_muminus_px_gen);
	fChain->SetBranchAddress("Z_muminus_py_gen", &Z_muminus_py_gen, &b_Z_muminus_py_gen);
	fChain->SetBranchAddress("Z_muminus_pz_gen", &Z_muminus_pz_gen, &b_Z_muminus_pz_gen);
	fChain->SetBranchAddress("Z_muminus_e_gen", &Z_muminus_e_gen, &b_Z_muminus_e_gen);
	fChain->SetBranchAddress("Z_muminus_pt_gen", &Z_muminus_pt_gen, &b_Z_muminus_pt_gen);
	fChain->SetBranchAddress("Z_muminus_et_gen", &Z_muminus_et_gen, &b_Z_muminus_et_gen);
	fChain->SetBranchAddress("Z_muminus_eta_gen", &Z_muminus_eta_gen, &b_Z_muminus_eta_gen);
	fChain->SetBranchAddress("Z_muminus_theta_gen", &Z_muminus_theta_gen, &b_Z_muminus_theta_gen);
	fChain->SetBranchAddress("Z_muminus_phi_gen", &Z_muminus_phi_gen, &b_Z_muminus_phi_gen);
	fChain->SetBranchAddress("Z_muminus_charge_gen", &Z_muminus_charge_gen, &b_Z_muminus_charge_gen);
	fChain->SetBranchAddress("Z_muminus_vx_gen", &Z_muminus_vx_gen, &b_Z_muminus_vx_gen);
	fChain->SetBranchAddress("Z_muminus_vy_gen", &Z_muminus_vy_gen, &b_Z_muminus_vy_gen);
	fChain->SetBranchAddress("Z_muminus_vz_gen", &Z_muminus_vz_gen, &b_Z_muminus_vz_gen);
	fChain->SetBranchAddress("Z_muminus_y_gen", &Z_muminus_y_gen, &b_Z_muminus_y_gen);
	fChain->SetBranchAddress("event_runNo", &event_runNo, &b_event_runNo);
	fChain->SetBranchAddress("event_evtNo", &event_evtNo, &b_event_evtNo);
	fChain->SetBranchAddress("event_lumi", &event_lumi, &b_event_lumi);
	fChain->SetBranchAddress("event_bunch", &event_bunch, &b_event_bunch);
	fChain->SetBranchAddress("event_nPV", &event_nPV, &b_event_nPV);
	fChain->SetBranchAddress("event_met_pfmet", &event_met_pfmet, &b_event_met_pfmet);
	fChain->SetBranchAddress("event_met_pfsumet", &event_met_pfsumet, &b_event_met_pfsumet);
	fChain->SetBranchAddress("event_met_pfmetsignificance", &event_met_pfmetsignificance, &b_event_met_pfmetsignificance);
	fChain->SetBranchAddress("event_met_pfmetPhi", &event_met_pfmetPhi, &b_event_met_pfmetPhi);
	fChain->SetBranchAddress("event_metMVA_met", &event_metMVA_met, &b_event_metMVA_met);
	fChain->SetBranchAddress("event_metMVA_sumet", &event_metMVA_sumet, &b_event_metMVA_sumet);
	fChain->SetBranchAddress("event_metMVA_metsignificance", &event_metMVA_metsignificance, &b_event_metMVA_metsignificance);
	fChain->SetBranchAddress("event_metMVA_metPhi", &event_metMVA_metPhi, &b_event_metMVA_metPhi);
	fChain->SetBranchAddress("event_fastJetRho", &event_fastJetRho, &b_event_fastJetRho);
	fChain->SetBranchAddress("event_met_genmet", &event_met_genmet, &b_event_met_genmet);
	fChain->SetBranchAddress("event_met_gensumet", &event_met_gensumet, &b_event_met_gensumet);
	fChain->SetBranchAddress("event_met_genmetsignificance", &event_met_genmetsignificance, &b_event_met_genmetsignificance);
	fChain->SetBranchAddress("event_met_genmetPhi", &event_met_genmetPhi, &b_event_met_genmetPhi);
	fChain->SetBranchAddress("event_mcPU_totnvtx", &event_mcPU_totnvtx, &b_event_mcPU_totnvtx);
	fChain->SetBranchAddress("event_mcPU_trueInteractions", &event_mcPU_trueInteractions, &b_event_mcPU_trueInteractions);
	fChain->SetBranchAddress("event_mcPU_bx", event_mcPU_bx, &b_event_mcPU_bx);
	fChain->SetBranchAddress("event_mcPU_nvtx", event_mcPU_nvtx, &b_event_mcPU_nvtx);
	Notify();
}

Bool_t kanamuon::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}

void kanamuon::Show(Long64_t entry)
{
	// Print contents of entry.
	// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}
Int_t kanamuon::Cut(Long64_t entry)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	Long64_t tmp; tmp=entry;
	return 1;
}
#endif // #ifdef kanamuon_ZJets_cxx
