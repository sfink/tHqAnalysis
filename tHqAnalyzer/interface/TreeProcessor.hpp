#ifndef TREEPROCESSOR_HPP
#define TREEPROCESSOR_HPP

#include "tHqAnalysis/tHqAnalyzer/interface/InputCollections.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/VariableContainer.hpp"

class TreeProcessor{
  public:
  
    virtual void Init(const InputCollections& input,VariableContainer& var)=0;
    virtual void Process(const InputCollections& input,VariableContainer& var) =0;


  protected:
  
    bool initialized;
};

#endif
