#ifndef THQANALYSIS_THQANALYZER_THQGENVARPROCESSOR_HPP
#define THQANALYSIS_THQANALYZER_THQGENVARPROCESSOR_HPP

#include "tHqAnalysis/tHqAnalyzer/interface/TreeProcessor.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/tHqUtils.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/tHqEvent.hpp"

class tHqGenVarProcessor: public TreeProcessor{

public:

  tHqGenVarProcessor();
  ~tHqGenVarProcessor();

  void Init(const InputCollections& input,VariableContainer& var);
  void Process(const InputCollections& input,VariableContainer& var);

};

#endif
