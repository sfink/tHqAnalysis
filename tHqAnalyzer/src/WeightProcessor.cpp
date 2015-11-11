#include "tHqAnalysis/tHqAnalyzer/interface/WeightProcessor.hpp"

using namespace std;

WeightProcessor::WeightProcessor(){}
WeightProcessor::~WeightProcessor(){}


void WeightProcessor::Init(const InputCollections& input,VariableContainer& vars){

  for (auto it=input.weights.begin(); it!=input.weights.end(); ++it){
    vars.InitVar(it->first);  
  }


  vars.InitVar( "Weight_orig",1.,"F" );

  initialized=true;
}

void WeightProcessor::Process(const InputCollections& input,VariableContainer& vars){
  if(!initialized) cerr << "tree processor not initialized" << endl;
  
  for (auto it=input.weights.begin(); it!=input.weights.end(); ++it){
    vars.FillVar( it->first,it->second);  
  }
  
  vars.FillVar("Weight_orig",input.Weight_orig);

  cout << "%%%%%%%%%%%%%%%%%%%%%%%%% /n Filling tha weights" << endl;
}

