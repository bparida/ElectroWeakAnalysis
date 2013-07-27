// -*- C++ -*-
//
// Package:    JetsORboostedV
// Class:      JetsORboostedV
// 
/**\class JetsORboostedV JetsORboostedV.cc ElectroWeakAnalysis/VPlusJets/src/JetsORboostedV.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jordan Damgov
//         Created:  Tue Sep 11 19:16:16 CDT 2012
// $Id: JetsORboostedV.cc,v 1.2 2013/02/06 21:05:11 jdamgov Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/PatCandidates/interface/Photon.h"

//
// class declaration
//
namespace ewk {

class JetsORboostedV : public edm::EDFilter {
   public:
      explicit JetsORboostedV(const edm::ParameterSet&);
      ~JetsORboostedV();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
      edm::InputTag mInputBoson;
      edm::InputTag mInputJets;
      unsigned int minNumberJets;
      unsigned int maxNumberJets;
      double minVPt;
      unsigned int minNumberPhotons;
};
}
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

ewk::JetsORboostedV::JetsORboostedV(const edm::ParameterSet& iConfig):
minNumberJets (iConfig.getUntrackedParameter<int>("minNumber")),
maxNumberJets (iConfig.getUntrackedParameter<int>("maxNumber")),
minVPt (iConfig.getUntrackedParameter<double>("minVpt")),
minNumberPhotons (iConfig.getUntrackedParameter<int>("minNumberPhotons"))
{
   //now do what ever initialization is needed
  if(  iConfig.existsAs<edm::InputTag>("srcVectorBoson"))
    mInputBoson = iConfig.getParameter<edm::InputTag>("srcVectorBoson") ; 
  if(  iConfig.existsAs<edm::InputTag>("srcJets"))
    mInputJets = iConfig.getParameter<edm::InputTag>("srcJets");

}


ewk::JetsORboostedV::~JetsORboostedV()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called on each new Event  ------------
bool ewk::JetsORboostedV::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
//  using namespace edm;

  bool result = false;

////////////////////////////////////
// V boson PT cut
  edm::Handle <reco::CandidateView> boson;
  iEvent.getByLabel( mInputBoson, boson);
  if( boson->size()!=1 ) return false; // Nothing to analyze ...

  const reco::Candidate *Vboson = &((*boson)[0]); 
  if( Vboson == 0) return false;

  if(Vboson->pt()>minVPt) result = true;

////////////////////////////////////
// Jets counting

  edm::Handle<edm::View<reco::Jet> > jets;
  iEvent.getByLabel( mInputJets, jets );

  if(jets->size() >= minNumberJets && jets->size() <=maxNumberJets) result = true;
  
////////////////////////////////////
// Photon counting
  edm::Handle<reco::PhotonCollection> photonH;
  iEvent.getByLabel("photons",photonH);

  if(photonH->size() < minNumberPhotons ) result = false;
  return result;
}

// ------------ method called once each job just before starting event loop  ------------
void 
ewk::JetsORboostedV::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ewk::JetsORboostedV::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
ewk::JetsORboostedV::beginRun(edm::Run&, edm::EventSetup const&)
{ 
return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
ewk::JetsORboostedV::endRun(edm::Run&, edm::EventSetup const&)
{

return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
ewk::JetsORboostedV::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
return true;

}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
ewk::JetsORboostedV::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
return true;

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ewk::JetsORboostedV::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
using ewk::JetsORboostedV;
DEFINE_FWK_MODULE(JetsORboostedV);
