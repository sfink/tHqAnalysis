#include "tHqAnalysis/tHqAnalyzer/interface/RecoVarProcessor.hpp"

using namespace std;

RecoVarProcessor::RecoVarProcessor(){}
RecoVarProcessor::~RecoVarProcessor(){}


void RecoVarProcessor::Init(const InputCollections& input,VariableContainer& vars){

  vars.InitVar("best_recbdtout");
  vars.InitVar("hyp_posbdt", "I");
  vars.InitVar("hyp_posdR", "I");
  vars.InitVar("hbbm", "F");
  vars.InitVar("hbbpt", "F");
  vars.InitVar("hbbphi", "F");
  vars.InitVar("hbbeta", "F");
  vars.InitVar("hbbdr", "F");
  vars.InitVars("hbbjtidx", 3);
  vars.InitVar("topm", "F");
  vars.InitVar("toppt", "F");
  vars.InitVar("topphi", "F");
  vars.InitVar("topeta", "F");
  vars.InitVar("topjtidx", "I");
  vars.InitVar("ljtidx", "I");
  vars.InitVar("coststh_rec", "F");

  initialized=true;
}

void RecoVarProcessor::Process(const InputCollections& input,VariableContainer& vars){

}
