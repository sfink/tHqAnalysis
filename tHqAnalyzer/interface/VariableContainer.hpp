#ifndef VARIABLECONTAINER_HPP
#define VARIABLECONTAINER_HPP

#include <map>
#include <string>
#include <iostream>

#include "TString.h"
#include "TTree.h"

class VariableContainer{
  public:
  
    VariableContainer();
    ~VariableContainer();

    void InitVar( TString name, std::string type="F" );
    void InitVar( TString name, float defaultValue, std::string type="F" );
    void FillVar( TString name, float value );
    void InitVars( TString name, float defaultValue, TString nEntryVariable, int maxentries =999 );
    void InitVars( TString name, TString nEntryVariable, int maxentries=999 );
    void InitIntVars( TString name, int defaulValue, int nEntries=999 );
    void FillVars( TString name, int index, float value ); 
    float* GetFloatVarPointer( TString name); 
    float* GetArrayVarPointer( TString name,int entry); 
    int* GetIntVarPointer( TString name); 
    int GetIntVar( TString name); 
    float GetFloatVar( TString name); 
    float GetArrayVar( TString name,int index ); 
    void ConnectTree(TTree* tree);
    
    void SetDefaultValues();
    void PrintArrayValue( TString name );
    void Dump();


  private:
  
    std::map<TString,Float_t> floatMap;
    std::map<TString,bool> floatMapFilled;
    std::map<TString,Float_t> floatMapDefaults;
    std::map<TString,int> intMap;
    std::map<TString,int> intMapDefaults;
    std::map<TString,bool> intMapFilled;
    std::map<TString,Float_t*> arrayMap;
    std::map<TString,Float_t> arrayMapDefaults;
    std::map<TString,bool> arrayMapFilled;
    std::map<TString,int> maxEntriesArrays;
    std::map<TString,TString > entryVariableOf;
    std::map<TString,Int_t*> arrayIntMap;
    std::map<TString,Int_t> arrayIntMapDefaults;
    std::map<TString,bool> arrayIntMapFilled;
    std::map<TString,int> maxEntriesArraysInt;
};

#endif
