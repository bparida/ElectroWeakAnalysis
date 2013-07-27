//#include "ElectroWeakAnalysis/VPlusJets/plugins/AllPassFilter.h"
#include "AllPassFilter.h"






//! ctor
AllPassFilter::AllPassFilter(const edm::ParameterSet& iConfig)
{
  edm::Service<TFileService> fs;
  
  m_totalEvents = fs -> make<TH1F>("totalEvents", "totalEvents", 1,  0., 1.);
}

// ----------------------------------------------------------------






//! dtor
AllPassFilter::~AllPassFilter()
{}

// ----------------------------------------------------------------






bool AllPassFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
  m_totalEvents -> Fill(0.5);
  
  return true;
}

// ----------------------------------------------------------------

DEFINE_FWK_MODULE(AllPassFilter);
