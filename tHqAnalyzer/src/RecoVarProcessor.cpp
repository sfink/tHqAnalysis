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
  vars.InitIntVars("hbbjtidx", -1, 3);
  vars.InitVar("topm", "F");
  vars.InitVar("toppt", "F");
  vars.InitVar("topphi", "F");
  vars.InitVar("topeta", "F");
  vars.InitVar("topjtidx", "I");
  vars.InitVar("ljtidx", "I");
  vars.InitVar("coststh_rec", "F");


  vars.InitVar("top_best_recbdtout", "F");
  vars.InitVar("top_hyp_posbdt", "I");
  vars.InitVar("top_hyp_posdR", "I");
  vars.InitVar("tophadwm", "F");
  vars.InitVar("tophadwpt", "F");
  vars.InitVar("tophadwphi", "F");
  vars.InitVar("tophadweta", "F");
  vars.InitVar("tophadm", "F");
  vars.InitVar("tophadpt", "F");
  vars.InitVar("tophadphi", "F");
  vars.InitVar("tophadeta", "F");
  vars.InitVar("tophaddr", "F");
  vars.InitIntVars("tophadjtidx", -1, 3);
  vars.InitVar("toplepm", "F");
  vars.InitVar("topleppt", "F");
  vars.InitVar("toplepphi", "F");
  vars.InitVar("toplepeta", "F");
  vars.InitVar("toplepjtidx", "I");

  vars.InitVar("mlpout", "F");

  initialized=true;
}

void RecoVarProcessor::Process(const InputCollections& input,VariableContainer& vars){

}
