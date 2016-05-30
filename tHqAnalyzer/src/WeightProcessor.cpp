#include "tHqAnalysis/tHqAnalyzer/interface/WeightProcessor.hpp"

using namespace std;

WeightProcessor::WeightProcessor(){}
WeightProcessor::~WeightProcessor(){}


void WeightProcessor::Init(const InputCollections& input,VariableContainer& vars){
  for (auto it=input.weights.begin(); it!=input.weights.end(); ++it){
    vars.InitVar(it->first);  
  }
  
  cout << "Init new Weights at the WeightProcessor" << endl;
  vars.InitVar( "Weight_orig",1.,"F" );
  vars.InitVar( "nweights",-9, "I");
  vars.InitVars( "weights_syst","nweights" );
  vars.InitStrings( "weights_syst_id","999","nweights" );
  
  initialized=true;
}

void WeightProcessor::Process(const InputCollections& input,VariableContainer& vars){
  if(!initialized) cerr << "tree processor not initialized" << endl;
    
  for (auto it=input.weights.begin(); it!=input.weights.end(); ++it){
    vars.FillVar( it->first,it->second); 
  }
  
  vars.FillVar("Weight_orig",input.Weight_orig);
  vars.FillVar("nweights",input.syst_weights.size());
  for(std::vector<float>::const_iterator itWeight = input.syst_weights.begin() ; itWeight != input.syst_weights.end(); ++itWeight){
    int iWeight = itWeight - input.syst_weights.begin();
    vars.FillVars( "weights_syst",iWeight,input.syst_weights[iWeight] );
  }


  for(std::vector<string>::const_iterator itWeight = input.syst_weights_id.begin() ; itWeight != input.syst_weights_id.end(); ++itWeight){
    int iWeight = itWeight - input.syst_weights_id.begin();
    vars.FillStrings("weights_syst_id",iWeight,input.syst_weights_id[iWeight]);
  }

  //  vars.DumpBasic();
}

