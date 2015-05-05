#include "tHqAnalysis/tHqAnalyzer/interface/RecoVarProcessor.hpp"

using namespace std;

RecoVarProcessor::RecoVarProcessor(){}
RecoVarProcessor::~RecoVarProcessor(){}


void RecoVarProcessor::Init(const InputCollections& input,VariableContainer& vars){

  vars.InitVar("best_recbdtout");
  vars.InitVar("hyp_posbdt");
  vars.InitVar("hyp_posdR");
  vars.InitVar("hbbm");
  vars.InitVar("hbbpt");
  vars.InitVar("hbbphi");
  vars.InitVar("hbbeta");
  vars.InitVar("hbbdr");
  vars.InitVar("hbbjtidx");
  vars.InitVar("topm");
  vars.InitVar("toppt");
  vars.InitVar("topphi");
  vars.InitVar("topeta");
  vars.InitVar("topjtidx");
  vars.InitVar("ljtidx");
  vars.InitVar("coststh_rec");

  initialized=true;
}

void RecoVarProcessor::Process(const InputCollections& input,VariableContainer& vars){

}
