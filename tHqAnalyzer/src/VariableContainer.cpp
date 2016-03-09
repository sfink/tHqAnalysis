#include "tHqAnalysis/tHqAnalyzer/interface/VariableContainer.hpp"
#include <iostream>
#include <iomanip>

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
    else {
      arrayMap[name][index]=value;
      arrayMapFilled[name]=true;
    }
  }
  else if(arrayIntMap.count(name)==1){
    if(maxEntriesArraysInt[name]<index){
      cerr << "array " << name << " is shorter than " << index << endl;
    }
    else{
      arrayIntMap[name][index]=value;
      arrayIntMapFilled[name]=true;
    }
    
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


void VariableContainer::DumpBasic(){
  cout << "#######################################" << endl;
  cout << "#########     DUMP BASIC  #############" << endl;
  cout << "#######################################" << endl << endl;
  
  cout << "Eventinfo : " << endl;
  if(intMapFilled["run"]) cout << setw(25) << "Run" << setw(10) << intMap["run"] << endl;
  if(intMapFilled["lbn"])  cout << setw(25) << "Lumi" << setw(10) << intMap["lbn"] << endl;
  if(longMapFilled["evt"])  cout << setw(25) << "Evt" << setw(10) << longMap["evt"] << endl<< endl;
  
  if(intMapFilled["nmu"] && intMapFilled["nel"] ) cout << setw(25) << "# of Leptons"  << setw(10) << intMap["nmu"] << endl;
  if(intMapFilled["njt"])  cout << setw(25) << "# of Jets"  << setw(10) << intMap["njt"] << endl;

  if(intMapFilled["nbtagl"]) cout << setw(25) << "# of btagged(CSVL) Jets"  << setw(10) << intMap["nbtagl"] << endl;
  if(intMapFilled["nbtagm"]) cout << setw(25) << "# of btagged(CSVM) Jets"  << setw(10) << intMap["nbtagm"] << endl;
  if(intMapFilled["nbtagt"]) cout << setw(25) << "# of btagged(CSVT) Jets"  << setw(10) << intMap["nbtagt"] << endl;
  
  if(intMapFilled["nbtagl_mva"]) cout << setw(25) << "# of btagged(MVAL) Jets"  << setw(10) << intMap["nbtagl_mva"] << endl;
  if(intMapFilled["nbtagm_mva"]) cout << setw(25) << "# of btagged(MVAM) Jets"  << setw(10) << intMap["nbtagm_mva"] << endl;  
  if(intMapFilled["nbtagt_mva"]) cout << setw(25) << "# of btagged(MVAT) Jets"  << setw(10) << intMap["nbtagt_mva"] << endl;
  
  if(intMapFilled["nfwdjt"])   cout << setw(25) << "# of forward Jets"  << setw(10) << intMap["nfwdjt"] << endl << endl;

  cout << "Leptons : " << endl;
  if(intMapFilled["nel"]) cout << setw(25) << "# of Electrons" << setw(10) << intMap["nel"] << endl;
  if(intMapFilled["nmu"]) cout << setw(25) << "# of Muons" << setw(10) << intMap["nmu"] << endl << endl;
  
  if(arrayMapFilled["leppt"]) cout << setw(25) << "Lep1 Pt" << setw(10) << arrayMap["leppt"][0] << endl;
  if(arrayMapFilled["lepeta"]) cout << setw(25) << "Lep1 Eta" << setw(10) << arrayMap["lepeta"][0] << endl;
  if(arrayMapFilled["lepphi"]) cout << setw(25) << "Lep1 Phi" << setw(10) << arrayMap["lepphi"][0] << endl<< endl;
  if(arrayMapFilled["leppdg"]) cout << setw(25) << "Lep1 PDG" << setw(10) << arrayMap["leppdg"][0] << endl<< endl;

  
  if(intMap["nmu"]+intMap["nel"]>1){
    if(arrayMapFilled["leppt"]) cout << setw(25) << "Lep2 Pt" << setw(10) << arrayMap["leppt"][1] << endl;
    if(arrayMapFilled["lepeta"]) cout << setw(25) << "Lep2 Eta" << setw(10) << arrayMap["lepeta"][1] << endl;
    if(arrayMapFilled["lepphi"]) cout << setw(25) << "Lep2 Phi" << setw(10) << arrayMap["lepphi"][1] << endl<< endl;
    if(arrayMapFilled["leppdg"]) cout << setw(25) << "Lep2 PDG" << setw(10) << arrayMap["leppdg"][1] << endl<< endl;
  }

  cout << "Jets : " << endl;
  if(intMap["njt"]>0){
    if(arrayMapFilled["jtpt"])  cout << setw(25) << "Jet1 Pt" << setw(10) << arrayMap["jtpt"][0] << endl;
    if(arrayMapFilled["jteta"])  cout << setw(25) << "Jet1 Eta"<< setw(10) << arrayMap["jteta"][0] << endl;
    if(arrayMapFilled["jtphi"])  cout << setw(25) << "Jet1 Phi" << setw(10) << arrayMap["jtphi"][0]<< endl;
    if(arrayMapFilled["jtcsvt"])  cout << setw(25) << "Jet1 CSV" << setw(10) << arrayMap["jtcsvt"][0]<<  endl;
    if(arrayMapFilled["jtCvsL"])  cout << setw(25) << "Jet1 CvsL" << setw(10) << arrayMap["jtCvsL"][0]<< endl;
    if(arrayMapFilled["jtCvsB"])  cout << setw(25) << "Jet1 CvsB" << setw(10) << arrayMap["jtCvsB"][0]<< endl;
    if(arrayMapFilled["jtcostheta_l"])  cout << setw(25) << "Jet1 CosTheta_L" << setw(10) << arrayMap["jtcostheta_l"][0]<< endl;
    if(arrayMapFilled["jtcostheta_j1"])  cout << setw(25) << "Jet1 CosTheta_j1" << setw(10) << arrayMap["jtcostheta_j1"][0]<< endl;
    if(arrayMapFilled["jtcostheta_cm"])  cout << setw(25) << "Jet1 CosTheta_CM" << setw(10) << arrayMap["jtcostheta_cm"][0]<< endl << endl;

  }

  if(intMap["njt"]>1){
    if(arrayMapFilled["jtpt"])  cout << setw(25) << "Jet2 Pt" << setw(10) << arrayMap["jtpt"][1] << endl;
    if(arrayMapFilled["jteta"])  cout << setw(25) << "Jet2 Eta" << setw(10) << arrayMap["jteta"][1] << endl;
    if(arrayMapFilled["jtphi"])  cout << setw(25) << "Jet2 Phi" << setw(10) << arrayMap["jtphi"][1] << endl;
    if(arrayMapFilled["jtcsvt"])  cout << setw(25) << "Jet2 CSV" << setw(10) << arrayMap["jtcsvt"][1] << endl;
    if(arrayMapFilled["jtCvsL"])  cout << setw(25) << "Jet2 CvsL" << setw(10) << arrayMap["jtCvsL"][1]<< endl;
    if(arrayMapFilled["jtCvsB"])  cout << setw(25) << "Jet2 CvsB" << setw(10) << arrayMap["jtCvsB"][1]<< endl;
    if(arrayMapFilled["jtcostheta_l"])  cout << setw(25) << "Jet2 CosTheta_L" << setw(10) << arrayMap["jtcostheta_l"][1]<< endl;
    if(arrayMapFilled["jtcostheta_j1"])  cout << setw(25) << "Jet2 CosTheta_j1" << setw(10) << arrayMap["jtcostheta_j1"][1]<< endl;
    if(arrayMapFilled["jtcostheta_cm"])  cout << setw(25) << "Jet2 CosTheta_CM" << setw(10) << arrayMap["jtcostheta_cm"][1]<< endl << endl;

  }

  if(intMap["njt"]>2){
    if(arrayMapFilled["jtpt"])  cout << setw(25) << "Jet3 Pt" << setw(10) << arrayMap["jtpt"][2] << endl;
    if(arrayMapFilled["jteta"])  cout << setw(25) << "Jet3 Eta" << setw(10) << arrayMap["jteta"][2] << endl;
    if(arrayMapFilled["jtphi"])  cout << setw(25) << "Jet3 Phi" << setw(10) << arrayMap["jtphi"][2] << endl;
    if(arrayMapFilled["jtcsvt"])  cout << setw(25) << "Jet3 CSV" << setw(10) << arrayMap["jtcsvt"][2] << endl;
    if(arrayMapFilled["jtCvsL"])  cout << setw(25) << "Jet3 CvsL" << setw(10) << arrayMap["jtCvsL"][2]<< endl;
    if(arrayMapFilled["jtCvsB"])  cout << setw(25) << "Jet3 CvsB" << setw(10) << arrayMap["jtCvsB"][2]<< endl;    
    if(arrayMapFilled["jtcostheta_l"])  cout << setw(25) << "Jet3 CosTheta_L" << setw(10) << arrayMap["jtcostheta_l"][2]<< endl;
    if(arrayMapFilled["jtcostheta_j1"])  cout << setw(25) << "Jet3 CosTheta_j1" << setw(10) << arrayMap["jtcostheta_j1"][2]<< endl;
    if(arrayMapFilled["jtcostheta_cm"])  cout << setw(25) << "Jet3 CosTheta_CM" << setw(10) << arrayMap["jtcostheta_cm"][2]<< endl << endl;

  }
  
    cout << "MET : " << endl;
  if(floatMapFilled["met"])  cout << setw(25) << "MET" << setw(10) << floatMap["met"] << endl;
  if(floatMapFilled["metphi"])  cout << setw(25) << "METphi" <<  setw(10) <<  floatMap["metphi"] <<endl << endl;

  
  cout << "Weights : " << endl;
  if(floatMapFilled["Weight_CSV"])  cout << setw(25) << "CSV Weight" << setw(10) << floatMap["Weight_CSV"] << endl;
  if(floatMapFilled["Weight_TopPt"])  cout << setw(25) << "TopPt Weight" << setw(10) << floatMap["Weight_TopPt"] << endl;
  if(floatMapFilled["Weight_PU"])  cout << setw(25) << "Pileup Weight" << setw(10) << floatMap["Weight_PU"] << endl;
  if(floatMapFilled["Weight_CSVLFup"])  cout << setw(25) << "CSV LFup" << setw(10) << floatMap["Weight_CSVLFup"] << endl;
  if(floatMapFilled["Weight_CSVLFdown"])  cout << setw(25) << "CSV LFdown" << setw(10) << floatMap["Weight_CSVLFdown"] << endl;
  if(floatMapFilled["Weight_CSVHFup"])  cout << setw(25) << "CSV HFup" << setw(10) << floatMap["Weight_CSVHFup"] << endl;
  if(floatMapFilled["Weight_CSVHFdown"])  cout << setw(25) << "CSV HFdown" << setw(10) << floatMap["Weight_CSVHFdown"] << endl;
  if(arrayMapFilled["weight_syst"])  cout << setw(25) << "syst weight[0]" << setw(10) << arrayMap["weight_syst"][0] << endl;
  if(floatMapFilled["Weight_NJets"])  cout << setw(25) << "NJet Weight" << setw(10) << floatMap["Weight_NJets"] << endl;
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
