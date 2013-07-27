// system include files
#include <memory>
#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <TMath.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
//
// Code to calculate the pull(color correlation) between W jets
//   and the helicity angle (Higgs system).
//

//TVector2 getTvect( reco::PFJet* patJet ){

template <typename T> TVector2 getTvect(const T&  patJet, bool isReco){

  TVector2 t_Vect(0,0);
  TVector2 null(0,0);
  TVector2 ci(0,0);
  TLorentzVector pi(0,0,0,0);
  TLorentzVector J(0,0,0,0);
  TVector2 r(0,0);
  double patJetpfcPt = 1e10;
  double r_mag = 1e10;
  unsigned int nOfconst = 0;

//Re-reconstruct the jet direction with the charged tracks
  //patJetpfc = patJet->getPFConstituents();
  std::vector<reco::PFCandidatePtr> patJetpfc = getPFConstituents(patJet, isReco);

  for(size_t idx = 0; idx < patJetpfc.size(); idx++){
    if( patJetpfc.at(idx)->charge() != 0 ){
      pi.SetPtEtaPhiE( patJetpfc.at(idx)->pt(), 
                       patJetpfc.at(idx)->eta(),
                       patJetpfc.at(idx)->phi(),
                       patJetpfc.at(idx)->energy() );
      J += pi;
      nOfconst++;
    }
  }

// if there are less than two charged tracks do not calculate the pull
//  (there is not enough info). It returns a null vector
  if( nOfconst < 2 ) return null;

  TVector2 v_J( J.Rapidity(), J.Phi() );

//calculate TVector using only charged tracks
  for(size_t idx = 0; idx < patJetpfc.size(); idx++){
    if( patJetpfc.at(idx)->charge() != 0  ){
      patJetpfcPt = patJetpfc.at(idx)->pt();
      pi.SetPtEtaPhiE( patJetpfc.at(idx)->pt(),
                       patJetpfc.at(idx)->eta(), 
                       patJetpfc.at(idx)->phi(),
                       patJetpfc.at(idx)->energy() );
      r.Set( pi.Rapidity() - J.Rapidity(), 
             Geom::deltaPhi(patJetpfc.at(idx)->phi(), J.Phi() ) );
      r_mag = r.Mod();
      t_Vect += ( patJetpfcPt / J.Pt() ) * r_mag * r;
    }
  }

 return t_Vect;

}



//double getDeltaTheta(reco::PFJet* j1, reco::PFJet* j2 ){
template <typename T1,typename T2> 
double getDeltaTheta(const T1&  j1,const T2& j2, bool isReco=false){

  //the pull vector is going to be calculated using j1 & j2
  double deltaTheta = 1e10;
  TLorentzVector pi(0,0,0,0);
  TLorentzVector v_j1(0,0,0,0);
  TLorentzVector v_j2(0,0,0,0);

  //re-reconstruct the jet direction with the charged tracks
  std::vector<reco::PFCandidatePtr> j1pfc = getPFConstituents(j1, isReco);
  for(size_t idx = 0; idx < j1pfc.size(); idx++){
    if( j1pfc.at(idx)->charge() != 0 ){
      pi.SetPtEtaPhiE( j1pfc.at(idx)->pt(), 
		       j1pfc.at(idx)->eta(),
		       j1pfc.at(idx)->phi(), 
		       j1pfc.at(idx)->energy() );
      v_j1 += pi;
    }
  }

  //re-reconstruct the jet direction with the charged tracks
  std::vector<reco::PFCandidatePtr> j2pfc = getPFConstituents(j2, isReco);
  for(size_t idx = 0; idx < j2pfc.size(); idx++){
    if( j2pfc.at(idx)->charge() != 0 ){
      pi.SetPtEtaPhiE( j2pfc.at(idx)->pt(), 
                       j2pfc.at(idx)->eta(),
                       j2pfc.at(idx)->phi(), 
                       j2pfc.at(idx)->energy() );
      v_j2 += pi;
    }
  }

  //  if( v_j2.Mag() < == 0 or v_j1.Mag() < == 0 ) return deltaTheta = 1e10;
  if( v_j2.Mag() <= 0 || v_j1.Mag() <= 0 ) return deltaTheta = 1e10;


  TVector2 v2_j1( v_j1.Rapidity(), v_j1.Phi());
  TVector2 v2_j2( v_j2.Rapidity(), v_j2.Phi());

  //use j1 to calculate the pull vector
  TVector2 t = getTvect(j1, isReco);

  if( t.Mod() == 0 ) return deltaTheta = 1e10;

  Double_t deltaphi = Geom::deltaPhi( v_j2.Phi(), v_j1.Phi() );
  Double_t deltaeta = v_j2.Rapidity() - v_j1.Rapidity();
  TVector2 BBdir( deltaeta, deltaphi );

  deltaTheta = t.DeltaPhi(BBdir);

  return deltaTheta;
} // end color correlation calculation 






//
// The code for cos(theta*) in the Higgs rest frame
//

double getHelicity( TLorentzVector p4 , TVector3 boost ){
  double hel = 1e10;
  p4.Boost( -boost );
  hel = TMath::Cos( p4.Vect().Angle( boost ) );
  return hel;
}





// Get constituents of a jet
template <typename T> 
std::vector<reco::PFCandidatePtr> getPFConstituents(const T&  jet, bool isReco) {

  std::vector<reco::PFCandidatePtr> jetpfc ;

  if(isReco)  {
    const reco::PFJet* pfjet  = static_cast<const reco::PFJet *>(jet);
    jetpfc = pfjet->getPFConstituents(); 
  }
  else {
    const pat::Jet* pfjet  = static_cast<const pat::Jet *>(jet);
    if(pfjet->isPFJet()) jetpfc = pfjet->getPFConstituents();
  } 
  return jetpfc;
}


