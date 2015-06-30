#include  "tHqAnalysis/tHqAnalyzer/interface/TriggerInfo.hpp"

TriggerInfo::TriggerInfo(std::map<std::string,bool> triggers_):triggers(triggers_){
      
}
bool TriggerInfo::IsTriggered(std::string triggername) const {
  if(triggers.count(triggername)==0){ 
    //    std::cerr << "trigger "  << triggername << " not existing (you might need to add it in the python cfg)" << std::endl;
    return false;
  }
  std::cout << "trigger "  << triggername << " with response " <<triggers.at(triggername) << std::endl;
  return triggers.at(triggername);
}

bool TriggerInfo::IsAnyTriggered(std::vector< std::string > triggernames) const {
  for(auto name=triggernames.begin(); name!=triggernames.end();name++){
    if(IsTriggered(*name)) return true;
  }
  return false;
}

