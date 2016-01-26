#include "tHqAnalysis/tHqAnalyzer/interface/TriggerVarProcessor.hpp"

using namespace std;

TriggerVarProcessor::TriggerVarProcessor(const std::vector<std::string> relevantTriggers_):relevantTriggers(relevantTriggers_){}
TriggerVarProcessor::~TriggerVarProcessor(){}


void TriggerVarProcessor::Init(const InputCollections& input,VariableContainer& vars){
  for (auto it=relevantTriggers.begin(); it!=relevantTriggers.end(); ++it){
    //    if(input.triggerInfo.Exists(*it)){
      cout << "Initializing " << replaceAsterix(*it) << endl;
      vars.InitVar(replaceAsterix(*it),"I");  
      // }
      //else cout << "TRIGGER DOES NOT EXIST !!!!!!!!!!!!!!!!!!!" << endl << endl;
  }
  
  initialized=true;
}

void TriggerVarProcessor::Process(const InputCollections& input,VariableContainer& vars){
  if(!initialized) cerr << "tree processor not initialized" << endl;
  for (auto it=relevantTriggers.begin(); it!=relevantTriggers.end(); ++it){
    vars.FillVar(replaceAsterix(*it),int(input.triggerInfo.IsTriggered(*it)));  
    cout << replaceAsterix(*it) << " : " << int(input.triggerInfo.IsTriggered(*it)) << endl;
  }  
}

std::string TriggerVarProcessor::replaceAsterix(std::string triggername){
  int asterix = triggername.find("*");
  if(triggername.find("*")==std::string::npos){    
    return triggername;
  }
  else{
    return triggername.replace(asterix,1,"X");
  }
}
