#include "tHqAnalysis/tHqAnalyzer/interface/TreeWriter.hpp"

using namespace std;

TreeWriter::TreeWriter(){
  tree=0;
  initialized=false;
  outFile=0;
}

TreeWriter::~TreeWriter(){

  if(outFile!=0){
    outFile->cd();
  }
  if(tree!=0){
    tree->Write(); 
  }  
  if(outFile!=0){
    outFile->Write();
    outFile->Close();
  }
  if(tree!=0){
    cout << "Tree Written to " << outFile->GetPath() << endl;
  }

}

void TreeWriter::Init( std::string fileName){
  cout << treename << " Tree initialized." << endl;
    
  outFile = new TFile( (fileName+"_Tree.root").c_str(), "RECREATE" );
  outFile->cd();
  cout << "Created File named " << fileName<<"_Tree.root" << endl;
  
  dir = (TDirectory*) outFile->mkdir("utm");
  dir->cd();
  tree = new TTree("t","t");
  vars = VariableContainer();
}

std::vector<TreeProcessor*> TreeWriter::GetTreeProcessors() const{
  return processors;
}

std::vector<std::string> TreeWriter::GetTreeProcessorNames() const{
  return processorNames;
}

void TreeWriter::SetTreeName(std::string name){
  treename = name;
}

std::string TreeWriter::GetTreeName(){
  return treename;
}

void TreeWriter::AddTreeProcessor(TreeProcessor* processor){
  processors.push_back(processor);
  stopwatches.push_back(TStopwatch());
}

void TreeWriter::AddTreeProcessor(TreeProcessor* processor, string name){
  processorNames.push_back(name);
  processors.push_back(processor);
  stopwatches.push_back(TStopwatch());
}

void TreeWriter::FillProcessorName(string name){
  processorNames.push_back(name);
}




void TreeWriter::FillProcessorMap(){
  assert(processors.size() == processorNames.size());
  for (size_t i = 0; i < processorNames.size(); ++i){
    ProcessorMap[processorNames[i]] = processors[i];
  }
}

void TreeWriter::RemoveTreeProcessor(string name){
  auto it = ProcessorMap.find(name);
  if (it != ProcessorMap.end()){
    cout << "Found Processor and removing it: " << name << endl;
    ProcessorMap.erase (it);
    std::vector<TreeProcessor*> processors_temp;
    std::vector<std::string> processorNames_temp;

    auto it_map = ProcessorMap.begin();
    while (it_map != ProcessorMap.end()) {
      processors_temp.push_back(it_map->second);
      processorNames_temp.push_back(it_map->first);
      ++it_map;
    }
    processors=processors_temp;
    processorNames=processorNames_temp;
  }
}



bool TreeWriter::Process(const InputCollections& input) {  
  
  if(!initialized){
    for(uint i=0; i<processors.size(); i++){
      processors[i]->Init(input,vars);
    }
    
    vars.ConnectTree(tree);
    
    initialized=true;
  }
  
  vars.SetDefaultValues();
  
  for(uint i=0; i<processors.size(); i++){
    stopwatches[i].Start(false);
    processors[i]->Process(input,vars);
    stopwatches[i].Stop();
  }

  FillTree();

  return true;
}


void TreeWriter::FillTree(){
  tree->Fill();
}


void TreeWriter::AddSampleInformation(){
  double ntotal = tree->GetEntries();
  double efficiency = 1.0;
  double kfac = 1.0;
  double xsec = 1.0;

  stringstream ss; ss<<"efficiency:"<<efficiency;
  stringstream ss1; ss1<<"ntotal:"<<ntotal;
  stringstream ss2; ss2<<"xsection:"<<xsec;
  stringstream ss3; ss3<<"kfactor:"<<kfac;

  
  std::cout << "SampleInformation is called. Tree contains " << ntotal << " Entries."  << std::endl;

  tree->GetUserInfo()->Add(new TObjString(ss.str().c_str()));
  tree->GetUserInfo()->Add(new TObjString(ss1.str().c_str()));
  tree->GetUserInfo()->Add(new TObjString(ss2.str().c_str()));
  tree->GetUserInfo()->Add(new TObjString(ss3.str().c_str()));

  
}
