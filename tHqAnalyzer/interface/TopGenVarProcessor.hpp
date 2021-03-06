#ifndef THQANALYSIS_THQANALYZER_MCMATCHVARPROCESSOR_HPP
#define THQANALYSIS_THQANALYZER_MCMATCHVARPROCESSOR_HPP

#include <vector>

#include "tHqAnalysis/tHqAnalyzer/interface/TreeProcessor.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/tHqUtils.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/GenTopEvent.hpp"

class TopGenVarProcessor: public TreeProcessor{
  
public:
  
  TopGenVarProcessor();
  ~TopGenVarProcessor();
    
  void Init(const InputCollections& input,VariableContainer& var);
  void Process(const InputCollections& input,VariableContainer& var);
  
};

#endif
