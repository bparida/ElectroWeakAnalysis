#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "ElectroWeakAnalysis/VPlusJets/interface/ZMuMuFilter.h"

using namespace std;
using namespace reco;
using namespace edm;

ZMuMuFilter::ZMuMuFilter(ParameterSet const& cfg) 
{
  muons_    = cfg.getParameter<edm::InputTag> ("muons");
  minPt_    = cfg.getParameter<double>        ("minPt");
  maxEta_   = cfg.getParameter<double>        ("maxEta");
  minMass_  = cfg.getParameter<double>        ("minMass");
  debug_    = cfg.getUntrackedParameter<int>  ("debug");
}
////////////////////////////////////////////////////////////////////
void ZMuMuFilter::beginJob() 
{
  NTotal_ = 0;
  NPass_  = 0;
  hMass_  = fs_->make<TH1F>("mass","mass",200,0,200);
  hPt_    = fs_->make<TH1F>("muonPt","muonPt",200,0,200);
  hEta_   = fs_->make<TH1F>("muonEta","muonEta",120,-3,3);
}
////////////////////////////////////////////////////////////////////
bool ZMuMuFilter::filter(edm::Event &event, const edm::EventSetup &iSetup)
{
  edm::Handle<View<Muon> > muons;
  event.getByLabel(muons_,muons); 
  bool check(false),cut_pt(true),cut_eta(true),cut_mass(true);
  NTotal_++;
  if (muons->size() > 1) {
    for(int i=0;i<2;i++) {
      cut_pt  *= ((*muons)[i].pt() > minPt_);
      cut_eta *= (fabs((*muons)[i].eta()) < maxEta_);
      if (debug_ == 2) {
        cout<<i<<" "<<(*muons)[i].pt()<<" "<<(*muons)[i].eta()<<" ";
      }
    }
    double m = ((*muons)[0].p4()+(*muons)[1].p4()).mass();
    if (debug_ == 2) { 
      cout<<m<<" ";
    }
    cut_mass = (m > minMass_);
    if (cut_pt && cut_eta && cut_mass) {
      check = true;
      hMass_->Fill(m);
      for(int i=0;i<2;i++) {
        hPt_->Fill((*muons)[i].pt());
        hEta_->Fill((*muons)[i].eta());
      }
      NPass_++;
    }
    if (debug_ == 2) {
      cout<<check<<endl;
    }
  }
  return check;
}
//////////////////////////////////////////////////////////////////////////////////////////
void ZMuMuFilter::endJob() 
{
  if (debug_ > 0) {
    cout<<"Number of events read:   "<<NTotal_<<endl;
    cout<<"Number of events passed: "<<NPass_<<endl;
  }
}
DEFINE_FWK_MODULE(ZMuMuFilter);

