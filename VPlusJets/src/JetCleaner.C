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
class JetCleaner : public edm::EDProducer
{
public:
  typedef std::vector<T> JetCollection;
  // construction/destruction
  explicit JetCleaner(const edm::ParameterSet& iConfig);
  virtual ~JetCleaner();
  
  // member functions
  void produce(edm::Event& iEvent,const edm::EventSetup& iSetup);
  void endJob();


private:  
  // member data
  edm::InputTag              srcJets_;
  std::vector<edm::InputTag> srcObjects_;
  double                     deltaRMin_;


  std::string  moduleLabel_;
  int idLevel_;
  double etaMax_;
  double etaMin_;
  double ptMin_;

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
JetCleaner<T>::JetCleaner(const edm::ParameterSet& iConfig)
  : srcJets_    (iConfig.getParameter<edm::InputTag>         ("srcJets"))
  , srcObjects_ (iConfig.getParameter<vector<edm::InputTag> >("srcObjects"))
  , deltaRMin_  (iConfig.getParameter<double>                ("deltaRMin"))
  , moduleLabel_(iConfig.getParameter<string>                ("@module_label"))
  , idLevel_    (iConfig.getParameter<int>                   ("idLevel"))
  , etaMax_     (iConfig.getParameter<double>                ("etaMax"))
  , etaMin_     (iConfig.getParameter<double>                ("etaMin"))
  , ptMin_      (iConfig.getParameter<double>                ("ptMin"))
  , nJetsTot_(0)
  , nJetsClean_(0)
{
  produces<JetCollection>();
}


//______________________________________________________________________________
template<typename T>
JetCleaner<T>::~JetCleaner(){}



////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
template<typename T>
void JetCleaner<T>::produce(edm::Event& iEvent,const edm::EventSetup& iSetup)
{
  auto_ptr<JetCollection> cleanJets(new JetCollection);
  edm::Handle<reco::JetView> jets;
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
      }
    }
  }
  
  for (unsigned int iJet=0;iJet<jets->size();iJet++)
    if (isClean[iJet]) {
      
      //calculate the Calo jetID
      bool passedId=false;
      bool ThisIsClean=false;

      const std::type_info & type = typeid((*jets)[iJet]); 
      if( type == typeid(reco::CaloJet) ) {
	passedId = true;
      }
      //calculate the PF jetID
      if ( type == typeid(reco::PFJet) ) {
	const reco::PFJet pfjet = static_cast<const reco::PFJet &>((*jets)[iJet]);
	ThisIsClean=true;
	//apply following only if |eta|<2.4: CHF>0, CEMF<0.99, chargedMultiplicity>0   
	if(( pfjet.chargedHadronEnergy()/ pfjet.energy())<= 0.0  
	   && fabs(pfjet.eta())<2.4) ThisIsClean=false; 
	if( (pfjet.chargedEmEnergy()/pfjet.energy())>= 0.99 
	    && fabs(pfjet.eta())<2.4 ) ThisIsClean=false;
	if( pfjet.chargedMultiplicity()<=0 && fabs(pfjet.eta())<2.4 ) 
	  ThisIsClean=false;
	
	// always require #Constituents > 1
	if( pfjet.nConstituents() <=1 ) ThisIsClean=false;

	if(ThisIsClean && 
	   (pfjet.neutralHadronEnergy()/pfjet.energy())< 0.99 
	   && (pfjet.neutralEmEnergy()/pfjet.energy())<0.99) 
	  passedId=true;	
      }
      // in case of GenJet apply no jet ID
      if ( type == typeid(reco::GenJet) ) passedId=true;

      
      bool isPassing = false;
      if(idLevel_==0) isPassing = true;
      if(idLevel_==1) isPassing = passedId;

      const T& goodJet = static_cast<const T&>((*jets)[iJet]);
      double pt = goodJet.pt();
      double eta = goodJet.eta();
      if(isPassing && fabs(eta)<etaMax_ && fabs(eta)>=etaMin_ && pt>ptMin_) cleanJets->push_back( goodJet );
    }

  nJetsTot_  +=jets->size();
  nJetsClean_+=cleanJets->size();

  delete [] isClean;  
  iEvent.put(cleanJets);
}




//______________________________________________________________________________
template<typename T>
void JetCleaner<T>::endJob()
{
  stringstream ss;
  ss<<"nJetsTot="<<nJetsTot_<<" nJetsClean="<<nJetsClean_
    <<" fJetsClean="<<100.*(nJetsClean_/(double)nJetsTot_)<<"%\n";
  cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++"
      <<"\n"<<moduleLabel_<<"(JetCleaner) SUMMARY:\n"<<ss.str()
      <<"++++++++++++++++++++++++++++++++++++++++++++++++++"
      <<endl;
}


////////////////////////////////////////////////////////////////////////////////
// plugin definition
////////////////////////////////////////////////////////////////////////////////

typedef JetCleaner<reco::CaloJet> CaloJetCleaner;
typedef JetCleaner<reco::PFJet>   PFJetCleaner;
typedef JetCleaner<reco::JPTJet>  JPTJetCleaner;
typedef JetCleaner<reco::GenJet>  GenJetCleaner;

DEFINE_FWK_MODULE(CaloJetCleaner);
DEFINE_FWK_MODULE(PFJetCleaner);
DEFINE_FWK_MODULE(JPTJetCleaner);
DEFINE_FWK_MODULE(GenJetCleaner);
