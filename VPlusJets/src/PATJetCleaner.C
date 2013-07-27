////////////////////////////////////////////////////////////////////////////////
//
// JetCleaner
// --------------
//
////////////////////////////////////////////////////////////////////////////////

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
 
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
 
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"  
#include "DataFormats/JetReco/interface/CaloJetCollection.h"  
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/JPTJet.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <memory>
#include <vector>
#include <sstream>


////////////////////////////////////////////////////////////////////////////////
// class definition
////////////////////////////////////////////////////////////////////////////////
template<typename T>
class PATJetCleaner : public edm::EDProducer
{
public:
  typedef std::vector<T> JetCollection;
  // construction/destruction
  explicit PATJetCleaner(const edm::ParameterSet& iConfig);
  virtual ~PATJetCleaner();
  
  // member functions
  void produce(edm::Event& iEvent,const edm::EventSetup& iSetup);
  void endJob();


private:  
  // member data
  edm::InputTag              srcJets_;
  std::vector<edm::InputTag> srcObjects_;
  double                     deltaRMin_;
  std::string  moduleLabel_;
  unsigned int nJetsTot_;
  unsigned int nJetsClean_;
};


using namespace std;


////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////
// idLevel: 0==NoId, 1==Loose, 2==Medium, 3==Tight.
//______________________________________________________________________________
template<typename T>
PATJetCleaner<T>::PATJetCleaner(const edm::ParameterSet& iConfig)
  : srcJets_    (iConfig.getParameter<edm::InputTag>         ("srcJets"))
  , srcObjects_ (iConfig.getParameter<vector<edm::InputTag> >("srcObjects"))
  , deltaRMin_  (iConfig.getParameter<double>                ("deltaRMin"))
  , moduleLabel_(iConfig.getParameter<string>                ("@module_label"))
  , nJetsTot_(0)
  , nJetsClean_(0)
{
  produces<JetCollection>();
}


//______________________________________________________________________________
template<typename T>
PATJetCleaner<T>::~PATJetCleaner(){}



////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
template<typename T>
void PATJetCleaner<T>::produce(edm::Event& iEvent,const edm::EventSetup& iSetup)
{
  auto_ptr<JetCollection> cleanJets(new JetCollection);
 edm::Handle<edm::View<reco::Jet> > jets;

//  edm::Handle<reco::JetView> jets;
  iEvent.getByLabel(srcJets_,jets);


  bool* isClean = new bool[jets->size()];
  for (unsigned int iJet=0;iJet<jets->size();iJet++) isClean[iJet] = true;
  
  for (unsigned int iSrc=0;iSrc<srcObjects_.size();iSrc++) {
    edm::Handle<reco::CandidateView> objects;
    iEvent.getByLabel(srcObjects_[iSrc],objects);
    
    for (unsigned int iJet=0;iJet<jets->size();iJet++) {
      const reco::Jet& jet = jets->at(iJet);
      for (unsigned int iObj=0;iObj<objects->size();iObj++) {
	const reco::Candidate& obj = objects->at(iObj);
	double deltaR = reco::deltaR(jet,obj);
	if (deltaR<deltaRMin_)  isClean[iJet] = false;
	// if (deltaR<deltaRMin_)  cout<<"deltaR="<<deltaR<<"iJet="<<iJet<<endl;
      }
    }
  }
  
  for (unsigned int iJet=0;iJet<jets->size();iJet++)
    if (isClean[iJet]) {      
      const T& goodJet = static_cast<const T&>((*jets)[iJet]);
      cleanJets->push_back( goodJet );
    }

  nJetsTot_  +=jets->size();
  nJetsClean_+=cleanJets->size();

  delete [] isClean;  
  iEvent.put(cleanJets);
}




//______________________________________________________________________________
template<typename T>
void PATJetCleaner<T>::endJob()
{
  stringstream ss;
  ss<<"nJetsTot="<<nJetsTot_<<" nJetsClean="<<nJetsClean_
    <<" fJetsClean="<<100.*(nJetsClean_/(double)nJetsTot_)<<"%\n";
  cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++"
      <<"\n"<<moduleLabel_<<"(PATJetCleaner) SUMMARY:\n"<<ss.str()
      <<"++++++++++++++++++++++++++++++++++++++++++++++++++"
      <<endl;
}


////////////////////////////////////////////////////////////////////////////////
// plugin definition
////////////////////////////////////////////////////////////////////////////////

typedef PATJetCleaner<reco::CaloJet> CaloPATJetCleaner;
//typedef PATJetCleaner<reco::PFJet>   PFPATJetCleaner;
typedef PATJetCleaner<pat::Jet>   PFPATJetCleaner;
typedef PATJetCleaner<reco::JPTJet>  JPTPATJetCleaner;
typedef PATJetCleaner<reco::GenJet>  GenPATJetCleaner;

DEFINE_FWK_MODULE(CaloPATJetCleaner);
DEFINE_FWK_MODULE(PFPATJetCleaner);
DEFINE_FWK_MODULE(JPTPATJetCleaner);
DEFINE_FWK_MODULE(GenPATJetCleaner);
