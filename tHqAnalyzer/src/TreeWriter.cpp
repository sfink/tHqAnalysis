#include "tHqAnalysis/tHqAnalyzer/interface/TreeWriter.hpp"

using namespace std;

TreeWriter::TreeWriter(){
  initialized=false;
  cout << "Tree initialized." << endl;
}

TreeWriter::~TreeWriter(){
  for(uint i=0; i<stopwatches.size(); i++){
    cout << "time spent in " << processorNames[i] << " -- real time: " << stopwatches[i].RealTime() << ", cpu time: " << stopwatches[i].CpuTime() << endl;
  }
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
/*TreeWriter::~TreeWriter(){
  outFile->cd();

  tree->Write(); 
  outFile->Write();
  outFile->Close();
  cout << "Tree Written to " << outFile->GetPath() << endl;
}
*/

void TreeWriter::Init( std::string fileName){
  
  outFile = new TFile( (fileName+"_Tree.root").c_str(), "RECREATE" );
  //  outFile->cd();
      
  dir = (TDirectory*) outFile->mkdir("utm");
  dir->cd();
  tree = new TTree("t","t");
  vars = VariableContainer();
}


void TreeWriter::AddTreeProcessor(TreeProcessor* processor){
  processors.push_back(processor);
}

void TreeWriter::FillProcessorName(string name){
  processorNames.push_back(name);
}

void TreeWriter::FillProcessorMap(){
  //  cout << "FILLING THE PROCESSOR MAP"<< processors.size() << " &" << processorNames.size() << endl;
  assert(processors.size() == processorNames.size());
  for (size_t i = 0; i < processorNames.size(); ++i)
    ProcessorMap[processorNames[i]] = processors[i];
}

void TreeWriter::RemoveTreeProcessor(string name){
  //  cout << "TRYING TO REMOVE " << name << endl;
  auto it = ProcessorMap.find(name);
  if (it != ProcessorMap.end()){
    //cout << "IM GOING IN!" << endl;
    ProcessorMap.erase (it);
    std::vector<TreeProcessor*> processors_temp;
    std::vector<std::string> processorNames_temp;

    auto it_map = ProcessorMap.begin();
    while (it_map != ProcessorMap.end()) {
      processors_temp.push_back(it_map->second);
      processorNames_temp.push_back(it_map->first);
      cout << it_map->first << endl;
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
    processors[i]->Process(input,vars);
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
