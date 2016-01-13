#ifndef THQANALYSIS_THQANALYZER_TRIGGERVARPROCESSOR_HPP
#define THQANALYSIS_THQANALYZER_TRIGGERVARPROCESSOR_HPP

#include <vector>

#include "tHqAnalysis/tHqAnalyzer/interface/TreeProcessor.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/tHqUtils.hpp"


class TriggerVarProcessor: public TreeProcessor{
  
public:
  
  TriggerVarProcessor(std::vector<std::string> relevantTriggers);
  ~TriggerVarProcessor();
    
  void Init(const InputCollections& input,VariableContainer& var);
  void Process(const InputCollections& input,VariableContainer& var);

private:
  const std::vector<std::string> relevantTriggers;
  std::string replaceAsterix(std::string triggername);
};

#endif
