#include "EffTableLoader.h"
#include "EffTableReader.h"

#include <iostream>

namespace {
  const unsigned nParameters = 2;
}

EffTableLoader::EffTableLoader () 
  : mParameters (0) 
{}

EffTableLoader::EffTableLoader (const std::string& fDataFile) 
  : mParameters (new EffTableReader (fDataFile)) 
{
  std::cout << fDataFile << '\n';
}

EffTableLoader::~EffTableLoader () {
  delete mParameters;
}

int EffTableLoader::GetBandIndex(float fEt, float fEta)const {
  int index=mParameters->bandIndex(fEt, fEta);
  return index;
}
std::vector<float> EffTableLoader::GetEfficiencyAndError (float fEt,float fEta) const {
  int index=mParameters->bandIndex(fEt, fEta);
  EffTableReader::Record rec=mParameters->record(index);
  std::vector<float> param=rec.parameters();
  return param;
}
std::vector<float> EffTableLoader::GetEfficiencyAndError (int index) const {
  EffTableReader::Record rec=mParameters->record(index);
  std::vector<float> param=rec.parameters();
  return param;
}


float EffTableLoader::GetEfficiency (float fEt,float fEta) const {
   std::vector<float> param = GetEfficiencyAndError (fEt, fEta);
   return param[0];
}

float EffTableLoader::GetError (float fEt,float fEta) const {
   std::vector<float> param = GetEfficiencyAndError (fEt, fEta);
   return param[1];
}

float EffTableLoader::GetEfficiency (int index) const {
   std::vector<float> param = GetEfficiencyAndError (index);
   return param[0];
}

float EffTableLoader::GetError (int index) const {
   std::vector<float> param = GetEfficiencyAndError (index);
   return param[1];
}



std::vector<std::pair<float, float> > EffTableLoader::GetCellInfo(int index)const {
  EffTableReader::Record rec=mParameters->record(index);
  std::pair<float, float> PtBin;
  PtBin.first = rec.EtMin();
  PtBin.second= rec.EtMax();
   
  std::pair<float, float> EtaBin;
  EtaBin.first = rec.etaMin();
  EtaBin.second= rec.etaMax();
   
  std::vector<std::pair<float, float> > BinInfo;
  BinInfo.push_back(PtBin);
  BinInfo.push_back(EtaBin);
  return BinInfo ;
}



std::pair<float, float> EffTableLoader::GetCellCenter(int index )const {
  EffTableReader::Record rec=mParameters->record(index);
  std::pair<float, float> BinCenter;
  BinCenter.first = rec.EtMiddle();
  BinCenter.second= rec.etaMiddle();
  return BinCenter ;
}





std::vector<std::pair<float, float> > EffTableLoader::GetCellInfo(float fEt, float fEta)const {
  int index=mParameters->bandIndex(fEt, fEta);
  return (this->GetCellInfo(index)) ;
}





std::pair<float, float> EffTableLoader::GetCellCenter(float fEt, float fEta )const {
  int index=mParameters->bandIndex(fEt, fEta);
  return (this->GetCellCenter(index)); 
}

int EffTableLoader::size(void) {
  return mParameters->size();
}


