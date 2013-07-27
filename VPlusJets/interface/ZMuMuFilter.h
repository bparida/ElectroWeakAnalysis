#ifndef ZMUMUFILTER_H
#define ZMUMUFILTER_H
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"

class ZMuMuFilter : public edm::EDFilter 
{
  public:
    explicit ZMuMuFilter(const edm::ParameterSet & cfg);
    virtual void beginJob();
    bool filter(edm::Event &event, const edm::EventSetup &iSetup);
    virtual void endJob();  
  private:
    edm::Service<TFileService> fs_;
    TH1F *hMass_,*hPt_,*hEta_;
    edm::InputTag muons_;
    double minPt_;
    double maxEta_;
    double minMass_;
    int    debug_; 
    int    NTotal_;
    int    NPass_;
};
#endif

