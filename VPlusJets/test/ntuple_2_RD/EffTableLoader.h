// -*- mode: C++ -*-
#ifndef EffTableLoader_h
#define EffTableLoader_h

#include <string>
#include <vector>
#include <utility>

class EffTableReader;
class EffTableLoader {
public:
  EffTableLoader ();
  EffTableLoader (const std::string& fDataFile);
  virtual ~EffTableLoader ();
  std::vector<float> GetEfficiencyAndError (float fEt, float fEta) const; 
  std::vector<float> GetEfficiencyAndError (int index) const;
  float GetEfficiency (float fEt,float fEta) const;
  float GetError (float fEt,float fEta) const;
  float GetEfficiency (int index) const;
  float GetError (int index) const;

  int GetBandIndex(float fEt, float fEta) const;
  std::vector<std::pair<float, float> > GetCellInfo(int index)const;
  std::vector<std::pair<float, float> > GetCellInfo(float fEt, float fEta)const;
  std::pair<float, float> GetCellCenter(int index )const;
  std::pair<float, float> GetCellCenter(float fEt, float fEta )const;
  int size(void);

private:
  EffTableReader* mParameters;
};

#endif
