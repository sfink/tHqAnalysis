#ifndef SELECTION_HPP
#define SELECTION_HPP

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "tHqAnalysis/tHqAnalyzer/interface/InputCollections.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/Cutflow.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/tHqUtils.hpp"

class Selection{

public:
  virtual void Init(const edm::ParameterSet& iConfig, Cutflow& cutflow)=0;
  virtual bool IsSelected(const InputCollections& input,Cutflow& cutflow)=0;

protected:
  bool initialized;


};


#endif
