#ifndef THQANALYSIS_THQANALYZER_BASEVARPROCESSOR_HPP
#define THQANALYSIS_THQANALYZER_BASEVARPROCESSOR_HPP

#include "tHqAnalysis/tHqAnalyzer/interface/TreeProcessor.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/tHqUtils.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/tHqEvent.hpp"

class BaseVarProcessor: public TreeProcessor{
  
  public:
    
    BaseVarProcessor();
    ~BaseVarProcessor();
    
    void Init(const InputCollections& input,VariableContainer& var);
    void Process(const InputCollections& input,VariableContainer& var);
  
};

#endif
