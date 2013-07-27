//---------------------------------------------------------------------------
// Adapted From Following Description:
// Original Author:  Fedor Ratnikov Nov 9, 2007
// $Id: EffTableReader.cc,v 1.2 2012/10/15 19:44:06 andersj Exp $
// Generic parameters for Jet corrections
//----------------------------------------------------------------------------
#include "EffTableReader.h"
#include <iostream>
#include <stdio.h>
#include <ctype.h>
#include <fstream>
#include <stdlib.h>

namespace {
  float getFloat (const std::string& token) {
    char* endptr;
    float result = strtod (token.c_str(), &endptr);
    if (endptr == token.c_str()) {
       std:: cout << "EffTableReadter error: " 
                  << "can not convert token " 
                  << token << " to float value" << std::endl;
    }
    return result;
  }
  unsigned getUnsigned (const std::string& token) {
    char* endptr;
    unsigned result = strtoul (token.c_str(), &endptr, 0);
    if (endptr == token.c_str()) {
       std:: cout << "EffTableReadter error: " 
                  << "can not convert token " 
                  << token << " to float value" << std::endl;
    }
    return result;
  }
}


EffTableReader::Record::Record (const std::string& fLine) 
  : mEtaMin (0), mEtaMax (0), mEtMax(0), mEtMin(0) 
{
  // quckly parse the line
  std::vector<std::string> tokens;
  std::string currentToken;
  for (unsigned ipos = 0; ipos < fLine.length (); ++ipos) {
    char c = fLine[ipos];

    if (c == '#') break; // ignore comments
    else if (isspace(c)) { // flush current token if any
      if (!currentToken.empty()) {
	tokens.push_back (currentToken);
	currentToken.clear();
      }
    }
    else {
      currentToken += c;
    }
  }
  if (!currentToken.empty()) tokens.push_back (currentToken); // flush end
  if (!tokens.empty ()) { // pure comment line
     unsigned nParam = 0;

     if (tokens.size() < 4) {
        std::cout << "EffTableReader error: " 
                  << "at least 4 tokens are expected: minPt, maxPt, eff, error. " 
                  << "Line ->>" << fLine << "<<-" << std::endl;  
     }
     else if (tokens.size() == 4) {
        mEtMin = getFloat (tokens[0]);
        mEtMax = getFloat (tokens[1]);
        mEtaMin = -5.0;
        mEtaMax = 5.0;
        nParam = 2;
     }
     else if (tokens.size() != 6) {
        std::cout << "EffTableReader error: " 
                  << "exactly 6 tokens are expected: minPt, maxPt, minEta, maxEta, eff, error. " 
                  << "Line ->>" << fLine << "<<-" << std::endl;  
     }
     else {
        // get parameters
        mEtMin = getFloat (tokens[0]);
        mEtMax = getFloat (tokens[1]);
        mEtaMin = getFloat (tokens[2]);
        mEtaMax = getFloat (tokens[3]);
        nParam = 4;
     }
  
    for (unsigned i = nParam; i < tokens.size(); ++i) {
      mParameters.push_back (getFloat (tokens[i]));
    }
  }
}


EffTableReader::EffTableReader (const std::string& fFile) {
  std::ifstream input (fFile.c_str());
  std::string line;
  while (std::getline (input, line)) {
    Record newrecord (line);
    if (!(newrecord.nParameters() == 0)) {
      mRecords.push_back (newrecord);
    }
  }
  if (mRecords.empty()) mRecords.push_back (Record ());

}


int EffTableReader::bandIndex (float fEt, float fEta) const{
  int bandInd =0;
  // std::cout << '\n';
  // std::cout << "finding band for: " << fEt << ", " << fEta << '\n';
  // std::cout << &mRecords << ":" << mRecords.size();
  for (unsigned i = 0; i < mRecords.size(); ++i) {
    // std::cout << "record " << i ;
    // std::cout << " Et: " << mRecords[i].EtMin() << '-' << mRecords[i].EtMax()
    // 	      << " eta: " << mRecords[i].etaMin() << '-' << mRecords[i].etaMax()
    // 	      << '\n';
    if(fEt>=mRecords[i].EtMin() && fEt<mRecords[i].EtMax()){
      if(fEta>=mRecords[i].etaMin() && fEta<mRecords[i].etaMax()){
	bandInd=i;
	break;                        
      }
    }
    
  }
  return bandInd;
}



