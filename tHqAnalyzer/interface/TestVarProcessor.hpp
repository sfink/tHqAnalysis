#ifndef THQANALYSIS_THQANALYZER_MVAVARPROCESSOR_HPP
#define THQANALYSIS_THQANALYZER_MVAVARPROCESSOR_HPP

#include "tHqAnalysis/tHqAnalyzer/interface/TreeProcessor.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/tHqUtils.hpp"

class TestVarProcessor: public TreeProcessor{
  
  public:
    
    TestVarProcessor();
    ~TestVarProcessor();
    
    void Init(const InputCollections& input,VariableContainer& var);
    void Process(const InputCollections& input,VariableContainer& var);
  
};

#endif
