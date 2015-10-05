#ifndef THQANALYSIS_THQANALYZER_WEIGHTPROCESSOR_HPP
#define THQANALYSIS_THQANALYZER_WEIGHTPROCESSOR_HPP

#include <vector>

#include "tHqAnalysis/tHqAnalyzer/interface/TreeProcessor.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/tHqUtils.hpp"


class WeightProcessor: public TreeProcessor{
  
public:
  
  WeightProcessor();
  ~WeightProcessor();
    
  void Init(const InputCollections& input,VariableContainer& var);
  void Process(const InputCollections& input,VariableContainer& var);
   
};

#endif
