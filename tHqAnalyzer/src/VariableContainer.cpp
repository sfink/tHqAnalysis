#include "tHqAnalysis/tHqAnalyzer/interface/VariableContainer.hpp"

using namespace std;

VariableContainer::VariableContainer(){

}


VariableContainer::~VariableContainer(){

}


void VariableContainer::InitVar( TString name,float defaultValue, std::string type ) {
  if(intMap.count(name)>0||floatMap.count(name)>0||arrayMap.count(name)>0||arrayIntMap.count(name)>0){
    cerr << name << " already initialized!" << endl;
  }
  if(type=="F"){
    floatMap[name] = 0;
    floatMapDefaults[name] = defaultValue;
    floatMapFilled[name] = false;
  }
  else if(type=="I"){
    intMap[name] = 0;
    intMapDefaults[name] = defaultValue;
    intMapFilled[name] = false;
  }
  else if(type=="L"){
    longMap[name] = 0;
    longMapDefaults[name] = defaultValue;
    longMapFilled[name] = false;
  }

  else
    cout << "unknown type " << type << endl;
}


void VariableContainer::InitVar( TString name, std::string type ) {
  if(intMap.count(name)>0||floatMap.count(name)>0||arrayMap.count(name)>0|| arrayIntMap.count(name)>0){
    cerr << name << " already initialized!" << endl;
  }

  InitVar(name,-1.,type);
}

void VariableContainer::InitString( TString name,TString defaultValue ) {
  if(stringMap.count(name)>0||arrayStringMap.count(name)>0){
    cerr << name << " already initialized!" << endl;
  }
  stringMap[name] = "0";
  stringMapDefaults[name] = defaultValue;
  stringMapFilled[name] = false;
}


void VariableContainer::InitVars( TString name, float defaultValue, TString nEntryVariable, int maxentries ){
  if(intMap.count(name)>0 || floatMap.count(name)>0 || arrayMap.count(name)>0 || arrayIntMap.count(name)>0){
    cerr << name << " already initialized!" << endl;
  }

  arrayMap[name] = new float[maxentries];
  arrayMapDefaults[name] = defaultValue;
  arrayMapFilled[name] = false;
  maxEntriesArrays[name] = maxentries;
  entryVariableOf[name] = nEntryVariable;
}

void VariableContainer::InitVars( TString name, TString nEntryVariable, int maxentries ){
  InitVars(name,-1.,nEntryVariable,maxentries);
}

void VariableContainer::InitStrings( TString name, TString defaultValue, TString nEntryVariable, int maxentries ){
  if(stringMap.count(name)>0 || arrayStringMap.count(name)>0){
    cerr << name << " already initialized!" << endl;
  }

  arrayStringMap[name] = new TString[maxentries];
  arrayStringMapDefaults[name] = defaultValue;
  arrayStringMapFilled[name] = false;
  maxEntriesArraysString[name] = maxentries;
  entryVariableOf[name] = nEntryVariable;
  

}


void VariableContainer::InitIntVars( TString name, int defaultValue, int nEntries ){
  if(intMap.count(name)>0 || floatMap.count(name)>0 || arrayMap.count(name)>0 || arrayIntMap.count(name)>0){
    cerr << name << " already initialized!" << endl;
  }

  arrayIntMap[name] = new int[nEntries];
  arrayIntMapDefaults[name] = defaultValue;
  arrayIntMapFilled[name] = false;
  maxEntriesArraysInt[name] = nEntries;
  entryVariableOf[name] = name;
  
}


void VariableContainer::FillVar( TString name, float value ) {
  if(intMap.count(name)==0&&floatMap.count(name)==0&&longMap.count(name)==0){
    cerr << name << " does not exist!" << endl;
  }
  if(floatMap.count(name)!=0){
    if(floatMapFilled[name]){
      cerr << name << " already filled!" << endl;
    }
    floatMap[name] = value;
    floatMapFilled[name] = true;
  }
  if(intMap.count(name)!=0){
    if(intMapFilled[name]){
      cerr << name << " already filled!" << endl;
    }
    intMap[name] = value;
    intMapFilled[name] = true;
  }

  if(longMap.count(name)!=0){
    if(longMapFilled[name]){
      cerr << name << " already filled!" << endl;
    }
    longMap[name] = value;
    longMapFilled[name] = true;
  }
//     std::cout << "filling variable " << name << " with " << value << std::endl;
}



void VariableContainer::FillVars( TString name, int index, float value ) {
  if(arrayMap.count(name)==0&&arrayIntMap.count(name)==0){
    cerr << name << " does not exist!" << endl;
  }
  else if(arrayMap.count(name)==1){
    if(maxEntriesArrays[name]<index){
      cerr << "array " << name << " is shorter than " << index << endl;
    }
    else arrayMap[name][index]=value;
  }
  else if(arrayIntMap.count(name)==1){
    if(maxEntriesArraysInt[name]<index){
      cerr << "array " << name << " is shorter than " << index << endl;
    }
    else arrayIntMap[name][index]=value;
  }
}

void VariableContainer::FillStrings( TString name, int index, TString value ) {
  if(arrayStringMap.count(name)==0){
    cerr << name << " does not exist!" << endl;
  }
  else if(arrayStringMap.count(name)==1){
    if(maxEntriesArraysString[name]<index){
      cerr << "array " << name << " is shorter than " << index << endl;
    }
    else {
      //      cout << "Assigning " << value << " to arrayStringMap["<<name<< "][" << index << "]" << endl;
      arrayStringMap[name][index]=value;
    }
  }
}


void VariableContainer::SetDefaultValues(){
  
  auto itF= floatMap.begin();
  auto itFdefault = floatMapDefaults.begin();
  auto itFfilled = floatMapFilled.begin();
  while (itF != floatMap.end()) {
    itF->second = itFdefault->second;
    itFfilled->second=false;
    ++itF;
    ++itFdefault;
    ++itFfilled;
  }
  
  auto itI= intMap.begin();
  auto itIdefault = intMapDefaults.begin();
  auto itIfilled = intMapFilled.begin();
  while (itI != intMap.end()) {
    itI->second = itIdefault->second;
    itIfilled->second=false;
    ++itI;
    ++itIdefault;
    ++itIfilled;
  }

  auto itL= longMap.begin();
  auto itLdefault = longMapDefaults.begin();
  auto itLfilled = longMapFilled.begin();
  while (itL != longMap.end()) {
    itL->second = itLdefault->second;
    itLfilled->second=false;
    ++itL;
    ++itLdefault;
    ++itLfilled;
  }

  auto itS= stringMap.begin();
  auto itSdefault = stringMapDefaults.begin();
  auto itSfilled = stringMapFilled.begin();
  while (itS != stringMap.end()) {
    itS->second = itSdefault->second;
    itSfilled->second=false;
    ++itS;
    ++itSdefault;
    ++itSfilled;
  }

  auto itA = arrayMap.begin();
  auto itAdefault = arrayMapDefaults.begin();
  auto itAmaxEntriesArrays = maxEntriesArrays.begin();
  auto itAfilled = arrayMapFilled.begin();
  while (itA != arrayMap.end()) {    
    for(int i=0;i<itAmaxEntriesArrays->second;++i)
      itA->second[i] = itAdefault->second;
    
    itAfilled->second=false;
    ++itA;
    ++itAdefault;
    ++itAmaxEntriesArrays;
    ++itAfilled;
  }

  auto itIA = arrayIntMap.begin();
  auto itIAdefault = arrayIntMapDefaults.begin();
  auto itIAmaxEntriesArrays = maxEntriesArraysInt.begin();
  auto itIAfilled = arrayIntMapFilled.begin();
  while (itIA != arrayIntMap.end()) {    
    for(int i=0;i<itIAmaxEntriesArrays->second;++i)
      itIA->second[i] = itIAdefault->second;
    
    itIAfilled->second=false;
    ++itIA;
    ++itIAdefault;
    ++itIAmaxEntriesArrays;
    ++itIAfilled;
  }

  auto itSA = arrayStringMap.begin();
  auto itSAdefault = arrayStringMapDefaults.begin();
  auto itSAmaxEntriesArrays = maxEntriesArraysString.begin();
  auto itSAfilled = arrayStringMapFilled.begin();
  while (itSA != arrayStringMap.end()) {    
    for(int i=0;i<itSAmaxEntriesArrays->second;++i)
      itSA->second[i] = itSAdefault->second;
    
    itSAfilled->second=false;
    ++itSA;
    ++itSAdefault;
    ++itSAmaxEntriesArrays;
    ++itSAfilled;
  }
  
}


void VariableContainer::ConnectTree(TTree* tree){
  auto itF= floatMap.begin();
  while (itF != floatMap.end()) {
    tree->Branch(itF->first, &(itF->second), itF->first+"/F" );
    itF++;
  }
  auto itI= intMap.begin();
  while (itI != intMap.end()) {
    tree->Branch(itI->first, &(itI->second), itI->first+"/I" );
    itI++;
  }
  auto itL= longMap.begin();
  while (itL != longMap.end()) {
    tree->Branch(itL->first, &(itL->second), itL->first+"/L" );
    itL++;
  }
  auto itS= stringMap.begin();
  while (itS != stringMap.end()) {
    tree->Branch(itS->first, &(itS->second), itS->first+"/C" );
    itS++;
  }
  auto itA= arrayMap.begin();
  while (itA != arrayMap.end()) {
    tree->Branch(itA->first,itA->second , itA->first+"["+entryVariableOf[itA->first]+"]/F" );
    itA++;
  }
  auto itIA= arrayIntMap.begin();
  while (itIA != arrayIntMap.end()) {
    tree->Branch(itIA->first,itIA->second , itIA->first+"["+to_string(maxEntriesArraysInt[itIA->first])+"]/I" );
    itIA++;
  }
  auto itSA= arrayStringMap.begin();
  while (itSA != arrayStringMap.end()) {
    tree->Branch(itSA->first,itSA->second , itSA->first+"["+to_string(maxEntriesArraysString[itSA->first])+"]/C" );
    itSA++;
  }
}

void VariableContainer::PrintArrayValue(TString name){
  
  float* printArray = arrayMap[name];
  int nEntries = intMap[entryVariableOf[name]];
  
  
  for(int i=0;i<nEntries;++i){
    cout << name+"[" << i << "] : " << printArray[i] << std::endl;
  }
}


void VariableContainer::Dump(){
  auto itF= floatMap.begin();
  auto itFdefault = floatMapDefaults.begin();
  cout << "floats: " << endl;
  while (itF != floatMap.end()) {
    cout << itF->first << " : " << itF->second << " : " << itFdefault->second << endl;
    ++itF;
    ++itFdefault;
  }
  auto itI= intMap.begin();
  auto itIdefault = intMapDefaults.begin();
  cout << "ints: " << endl;
  while (itI != intMap.end()) {
    cout << itI->first << " : " << itI->second << " : " << itIdefault->second << endl;
    ++itI;
    ++itIdefault;
  }
  auto itS= stringMap.begin();
  auto itSdefault = stringMapDefaults.begin();
  cout << "strings: " << endl;
  while (itS != stringMap.end()) {
    cout << itS->first << " : " << itS->second << " : " << itSdefault->second << endl;
    ++itS;
    ++itSdefault;
  }
}


float* VariableContainer::GetFloatVarPointer(TString name){
  if(floatMap.count(name)==0){
    cerr << name << " does not exist!" << endl;
    return 0;
  }
  return &(floatMap[name]);
}

float* VariableContainer::GetArrayVarPointer(TString name, int entry){
  if(arrayMap.count(name)==0){
    cerr << name << " does not exist!" << endl;
    return 0;
  }
  return &(arrayMap[name][entry]);
}


int* VariableContainer::GetIntVarPointer(TString name){
  if(intMap.count(name)==0&&floatMap.count(name)==0&&arrayMap.count(name)==0){
    cerr << name << " does not exist!" << endl;
    return 0;
  }
  return &(intMap[name]);
}

TString* VariableContainer::GetStringVarPointer(TString name){
  if(stringMap.count(name)==0){
    cerr << name << " does not exist!" << endl;
    return 0;
  }
  return &(stringMap[name]);
}

TString* VariableContainer::GetArrayStringVarPointer(TString name, int entry){
  if(arrayStringMap.count(name)==0){
    cerr << name << " does not exist!" << endl;
    return 0;
  }
  return &(arrayStringMap[name][entry]);
}


float VariableContainer::GetFloatVar(TString name){
  float* x=GetFloatVarPointer(name);
  if(x!=0)
    return *x;
  else 
    return -999;
}


float VariableContainer::GetArrayVar(TString name, int entry){
  float* x=GetArrayVarPointer(name,entry);
  if(x!=0)
    return *x;
  else 
    return -999;
}


int VariableContainer::GetIntVar(TString name){
  int* x=GetIntVarPointer(name);
  if(x!=0)
    return *x;
  else 
    return -999;
}

/*
TString VariableContainer::GetStringVar(TString name){
  TString* x=GetStringVarPointer(name);
  if(x!="0")
    return *x;
  else 
    return "-999";
}
*/
