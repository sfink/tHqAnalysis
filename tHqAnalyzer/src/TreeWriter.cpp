#include "tHqAnalysis/tHqAnalyzer/interface/TreeWriter.hpp"

using namespace std;

TreeWriter::TreeWriter(){
  initialized=false;
  cout << "Tree initialized." << endl;
}


TreeWriter::~TreeWriter(){
  outFile->cd();

  tree->Write(); 
  outFile->Write();
  outFile->Close();
  cout << "Tree Written to " << outFile->GetPath() << endl;
}


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
