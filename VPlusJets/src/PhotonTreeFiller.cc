/*****************************************************************************
* Storing the Photon information 
 *****************************************************************************/
#include "ElectroWeakAnalysis/VPlusJets/interface/PhotonTreeFiller.h"

// CMS includes
#include "TMath.h" 
#include <TLorentzVector.h>
#include "DataFormats/Common/interface/View.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidatePhotonExtraFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidatePhotonExtra.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidatePhotonExtraFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidatePhotonExtra.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

// Monte Carlo stuff

// Header file

using namespace std;
using namespace edm;
using namespace pat;
using namespace reco;
//using namespace isodeposit;

ewk::PhotonTreeFiller::PhotonTreeFiller(const char *name, TTree* tree, 
				  const edm::ParameterSet iConfig):
    inputTagIsoValPhotonsPFId_(iConfig.getParameter< std::vector<edm::InputTag> >("IsoValPhoton"))
{

  //*********************  Run Over AOD or PAT  ***********//
  if( iConfig.existsAs<bool>("runningOverAOD"))
    runoverAOD = iConfig.getParameter<bool>("runningOverAOD");
    		
  tree_     = tree;

  if( !(tree==0) ) SetBranches();

}



void ewk::PhotonTreeFiller::SetBranches()
{


  // Declare photon branches

  SetBranchSingle( &NumPhotons, "NumPhotons");
  SetBranch( Et, "Photon_Et");
  SetBranch( E, "Photon_E");
  SetBranch( Eta, "Photon_Eta");
  SetBranch( Phi, "Photon_Phi");
  SetBranch( Theta, "Photon_Theta");

  SetBranch( Px, "Photon_Px");
  SetBranch( Py, "Photon_Py");
  SetBranch( Pz, "Photon_Pz");
  SetBranch( Vx, "Photon_Vx");
  SetBranch( Vy, "Photon_Vy");
  SetBranch( Vz, "Photon_Vz");
  
  SetBranch( SC_Et, "Photon_SC_Et");
  SetBranch( SC_E, "Photon_SC_E");
  SetBranch( SC_Eta, "Photon_SC_Eta");
  SetBranch( SC_Phi, "Photon_SC_Phi");
  SetBranch( SC_Theta, "Photon_SC_Theta");
  SetBranch( SC_x, "Photon_SC_x");
  SetBranch( SC_y, "Photon_SC_y");
  SetBranch( SC_z, "Photon_SC_z");

  SetBranch( PFisocharged03, "PFisocharged03");
  SetBranch( PFisophoton03, "PFisophoton03");
  SetBranch( PFisoneutral03, "PFisoneutral03");
  SetBranch( trkSumPtHollConeDR04_Pho11, "trkSumPtHollowConeDR04_Photon11");
  SetBranch( ecalRecHitSumEtConeDR04_Pho11, "ecalRecHitSumEtConeDR04_Photon11");
  SetBranch( hcalTowSumEtConeDR04_Pho11, "hcalTowerSumEtConeDR04_Photon11");

  SetBranch( HoverE, "Photon_HoverE");
  SetBranch( HoverE2011, "Photon_HoverE2011");
  SetBranch( SigmaIetaIeta, "Photon_SigmaIetaIeta");

  SetBranch( hasPixelSeed, "Photon_hasPixelSeed");
  SetBranch( passElecVeto, "Photon_passElecVeto");
  SetBranch( Id2011, "Photon_Id2011");
  SetBranch( Id2012, "Photon_Id2012");

}


//////////////////////////////////////////////////////////////////
/////// Helper for Effective Areas ///////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
float ewk::PhotonTreeFiller::EAch( float x){
 float EA = 0.012;
 if(x>1.0)   EA = 0.010;
 if(x>1.479) EA = 0.014;
 if(x>2.0)   EA = 0.012;
 if(x>2.2)   EA = 0.016;
 if(x>2.3)   EA = 0.020;
 if(x>2.4)   EA = 0.012;
 return EA;
}

float ewk::PhotonTreeFiller::EAnh( float x){
 float EA = 0.030;
 if(x>1.0)   EA = 0.057;
 if(x>1.479) EA = 0.039;
 if(x>2.0)   EA = 0.015;
 if(x>2.2)   EA = 0.024;
 if(x>2.3)   EA = 0.039;
 if(x>2.4)   EA = 0.072;
 return EA;
}

float ewk::PhotonTreeFiller::EApho( float x){
 float EA = 0.148;
 if(x>1.0)   EA = 0.130;
 if(x>1.479) EA = 0.112;
 if(x>2.0)   EA = 0.216;
 if(x>2.2)   EA = 0.262;
 if(x>2.3)   EA = 0.260;
 if(x>2.4)   EA = 0.266;
 return EA;
}

//////////////////////////////////////////////////////////////////
/////// Helper for above function ////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

void ewk::PhotonTreeFiller::SetBranchSingle( float* x, std::string name)
{
  tree_->Branch( name.c_str(), x, ( name+"/F").c_str() );
  bnames.push_back( name );
}

void ewk::PhotonTreeFiller::SetBranchSingle( int* x, std::string name)
{
  tree_->Branch( name.c_str(), x, ( name+"/I").c_str() );
  bnames.push_back( name );
}



void ewk::PhotonTreeFiller::SetBranch( float* x, std::string name)
{
  tree_->Branch( name.c_str(), x, ( name+"[NumPhotons]/F").c_str() );
  bnames.push_back( name );
}


void ewk::PhotonTreeFiller::SetBranch( int* x, std::string name)
{
  tree_->Branch( name.c_str(), x, ( name+"[NumPhotons]/I").c_str() );
  bnames.push_back( name );
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////



void ewk::PhotonTreeFiller::FillBranches() const 
{

  for(std::vector<std::string>::iterator it = bnames.begin();
      it != bnames.end(); ++it) {
    if(TBranch *br = tree_->GetBranch( (*it).c_str() ) ) br->Fill();
  }

}




void ewk::PhotonTreeFiller::init()   
{

  // initialize private data members
  NumPhotons = 0; 

  for (int j =0; j< NUM_PHO_MAX; ++j) {
    Et[j] = -1.0;
    E[j] = -1.0;
    Eta[j] = -10.0;
    Phi[j] = -10.0;
    Theta[j] = -10.0;

    Px[j] = -999999.9;
    Py[j] = -999999.9;
    Pz[j] = -999999.9;
    Vx[j] = -999999.9;
    Vy[j] = -999999.9;
    Vz[j] = -999999.9;

    SC_Et[j] = -1.0;
    SC_E[j] = -1.0;
    SC_Eta[j] = -10.0;
    SC_Phi[j] = -10.0;
    SC_Theta[j] = -10.0;

    SC_x[j] = -999999.9;
    SC_y[j] = -999999.9;
    SC_z[j] = -999999.9;
       
    PFisocharged03[j] = -1;
    PFisophoton03[j]  = -1;
    PFisoneutral03[j] = -1;
    trkSumPtHollConeDR04_Pho11[j] = -1;
    ecalRecHitSumEtConeDR04_Pho11[j] = -1;
    hcalTowSumEtConeDR04_Pho11[j] = -1;

    HoverE[j] = -1.; 
    SigmaIetaIeta[j] = -1.; 

    hasPixelSeed[j] = 0.;
    passElecVeto[j]=0.;

    Id2011[j]=0;
    Id2012[j]=0;

  }
  // initialization done

   
}



void ewk::PhotonTreeFiller::fill(const edm::Event& iEvent)
{
  // first initialize to the default values
  init();


   edm::Handle<reco::BeamSpot> bsHandle;
   iEvent.getByLabel("offlineBeamSpot", bsHandle);
   const reco::BeamSpot &beamspot = *bsHandle.product();

   edm::Handle<reco::ConversionCollection> hConversions;
   iEvent.getByLabel("allConversions", hConversions);

   edm::Handle<reco::GsfElectronCollection> hElectrons;
   iEvent.getByLabel("gsfElectrons", hElectrons);

   // Photons - from reco
    edm::Handle<reco::PhotonCollection> photonH;
    iEvent.getByLabel("photons",photonH);
    unsigned nrecopho=photonH->size();

/*
     edm::Handle<edm::View<pat::Photon> > phoHandle;
     iEvent.getByLabel("selectedPatPhotons",phoHandle);
     edm::View<pat::Photon>::const_iterator photon;
*/

     unsigned nTypes=3;
     IsoDepositVals photonIsoValPFId(nTypes);
     for (size_t j = 0; j<inputTagIsoValPhotonsPFId_.size(); ++j) {
       iEvent.getByLabel(inputTagIsoValPhotonsPFId_[j], photonIsoValPFId[j]);
     }

     const IsoDepositVals * photonIsoVals = &photonIsoValPFId;

     edm::Handle<double> rho;
     const edm::InputTag eventrho("kt6PFJetsPFlow", "rho");
     iEvent.getByLabel(eventrho,rho);
     double fastJetRho = *rho;


      for(unsigned ipho=0; ipho<nrecopho;++ipho) {
	 reco::PhotonRef myPhotonRef(photonH,ipho);

	 if (!myPhotonRef->isEB()&&!myPhotonRef->isEE()) continue;
         if ( myPhotonRef->et()<20.) continue;

	 passElecVeto[NumPhotons] = !ConversionTools::hasMatchedPromptElectron(myPhotonRef->superCluster(), hElectrons, hConversions, beamspot.position());
         hasPixelSeed[NumPhotons] = myPhotonRef->hasPixelSeed();

         Et[NumPhotons] = myPhotonRef->et();
         E[NumPhotons] = myPhotonRef->energy();
         Eta[NumPhotons] = myPhotonRef->eta();
         Phi[NumPhotons] = myPhotonRef->phi();
         Theta[NumPhotons] = myPhotonRef->theta();
 
         Px[NumPhotons] = myPhotonRef->px();
         Py[NumPhotons] = myPhotonRef->py();
         Pz[NumPhotons] = myPhotonRef->pz();
         Vx[NumPhotons] = myPhotonRef->vx();
         Vy[NumPhotons] = myPhotonRef->vy();
         Vz[NumPhotons] = myPhotonRef->vz();

         PFisocharged03[NumPhotons] =  max((*(*photonIsoVals)[0])[myPhotonRef] - fastJetRho*EAch(fabs(myPhotonRef->eta())),0.);
         PFisophoton03[NumPhotons]  = max((*(*photonIsoVals)[1])[myPhotonRef] - fastJetRho*EAnh(fabs(myPhotonRef->eta())),0.);
         PFisoneutral03[NumPhotons] = max((*(*photonIsoVals)[2])[myPhotonRef] - fastJetRho*EApho(fabs(myPhotonRef->eta())),0.);

         trkSumPtHollConeDR04_Pho11[NumPhotons] = myPhotonRef->trkSumPtHollowConeDR04();
         ecalRecHitSumEtConeDR04_Pho11[NumPhotons] = myPhotonRef->ecalRecHitSumEtConeDR04();
         hcalTowSumEtConeDR04_Pho11[NumPhotons] = myPhotonRef->hcalTowerSumEtConeDR04();
         
         SC_Et[NumPhotons] = (myPhotonRef->superCluster()->energy())/(cosh(myPhotonRef->superCluster()->position().eta()));
         SC_E[NumPhotons] = myPhotonRef->superCluster()->energy();
         SC_Eta[NumPhotons] = myPhotonRef->superCluster()->position().eta();
         SC_Phi[NumPhotons] = myPhotonRef->superCluster()->position().phi();
         SC_Theta[NumPhotons] = myPhotonRef->superCluster()->position().theta();

         SC_x[NumPhotons] = (myPhotonRef->superCluster())->x();
         SC_y[NumPhotons] = (myPhotonRef->superCluster())->y();
         SC_z[NumPhotons] = (myPhotonRef->superCluster())->z();

         HoverE[NumPhotons] = myPhotonRef->hadTowOverEm();
         HoverE2011[NumPhotons] = myPhotonRef->hadronicOverEm();

         SigmaIetaIeta[NumPhotons]= myPhotonRef->sigmaIetaIeta();
	 
         // Photon ID 2012 - Medium (EB 80%)
         if ( myPhotonRef->hadTowOverEm()<0.05 
              && ((myPhotonRef->isEB()&& myPhotonRef->sigmaIetaIeta()<0.011)||(myPhotonRef->isEE()&& myPhotonRef->sigmaIetaIeta()<0.033))
              && passElecVeto[NumPhotons]
//              && ((myPhotonRef->isEB()&& PFisocharged03[NumPhotons] <1.5 ) || (myPhotonRef->isEE()&& PFisocharged03[NumPhotons] <1.2 ))
//              && ( (myPhotonRef->isEB()&& PFisophoton03[NumPhotons] <(0.7+0.005*(myPhotonRef->et())) )|| (myPhotonRef->isEE()&& PFisophoton03[NumPhotons] <(1.0+0.005*(myPhotonRef->et()))) )
//              && ( (myPhotonRef->isEB()&& PFisoneutral03[NumPhotons] <(1.0+0.04*(myPhotonRef->et())) )||(myPhotonRef->isEE()&& PFisoneutral03[NumPhotons] <(1.5+0.04*(myPhotonRef->et()))) )
//       Tight (EB 70%)
              && ((myPhotonRef->isEB()&& PFisocharged03[NumPhotons] <0.7 ) || (myPhotonRef->isEE()&& PFisocharged03[NumPhotons] <0.5 ))
              && ( (myPhotonRef->isEB()&& PFisophoton03[NumPhotons] <(0.5+0.005*(myPhotonRef->et())) )|| (myPhotonRef->isEE()&& PFisophoton03[NumPhotons] <(1.0+0.005*(myPhotonRef->et()))) )
              && ( (myPhotonRef->isEB()&& PFisoneutral03[NumPhotons] <(0.4+0.04*(myPhotonRef->et())) )||(myPhotonRef->isEE()&& PFisoneutral03[NumPhotons] <(1.5+0.04*(myPhotonRef->et()))) )
         ) Id2012[NumPhotons] =1;

         // Photon ID 2011 - Vgamma
	 if (  myPhotonRef->hadronicOverEm() <0.05
               && !hasPixelSeed[NumPhotons]
               && ((myPhotonRef->isEB()&& myPhotonRef->sigmaIetaIeta()<0.011)||(myPhotonRef->isEE()&& myPhotonRef->sigmaIetaIeta()<0.03))
               && ((myPhotonRef->isEB()&& myPhotonRef->sigmaIetaIeta()>0.001)||(myPhotonRef->isEE()&& myPhotonRef->sigmaIetaIeta()>-0.0))
               && ( (myPhotonRef->isEB()&&(myPhotonRef->trkSumPtHollowConeDR04()<(2.0 + 0.001*(myPhotonRef->et()) + 0.0167*fastJetRho))) 
			|| (myPhotonRef->isEE()&&(myPhotonRef->trkSumPtHollowConeDR04()<(2.0 + 0.001 *(myPhotonRef->et()) + 0.032*fastJetRho) )) )
               && ( (myPhotonRef->isEB()&&(myPhotonRef->ecalRecHitSumEtConeDR04()<(4.2 + 0.006 *(myPhotonRef->et()) + 0.183 *fastJetRho) )) 
			|| (myPhotonRef->isEE()&&(myPhotonRef->ecalRecHitSumEtConeDR04()<(4.2 + 0.006 *(myPhotonRef->et()) + 0.090 *fastJetRho) )) )
               && ( (myPhotonRef->isEB()&&(myPhotonRef->hcalTowerSumEtConeDR04()<(2.2 + 0.0025 *(myPhotonRef->et()) + 0.062 *fastJetRho) )) 
			|| (myPhotonRef->isEE()&&(myPhotonRef->hcalTowerSumEtConeDR04()<(2.2 + 0.0025 *(myPhotonRef->et()) + 0.180 *fastJetRho) )) )
         ) Id2011[NumPhotons] =1;

         NumPhotons++;
       }


}


