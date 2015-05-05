#ifndef THQANALYSIS_THQANALYZER_RECOVARPROCESSOR_HPP
#define THQANALYSIS_THQANALYZER_RECOVARPROCESSOR_HPP

#include "tHqAnalysis/tHqAnalyzer/interface/TreeProcessor.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/tHqUtils.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/tHqEvent.hpp"

class RecoVarProcessor: public TreeProcessor{
  
  public:
    
    RecoVarProcessor();
    ~RecoVarProcessor();
    
    void Init(const InputCollections& input,VariableContainer& var);
    void Process(const InputCollections& input,VariableContainer& var);
  
};

#endif
