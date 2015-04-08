#include "tHqAnalysis/tHqAnalyzer/interface/tHqEvent.hpp"


tHqEvent::tHqEvent(const InputCollections& input_):input(input_),verbose(false),btagger("combinedSecondaryVertexBJetTags"){
  ResetEvent();
}


tHqEvent::~tHqEvent(){

}


void tHqEvent::ResetEvent(){

  // Charged Lepton
  lepVecCand = math::XYZTLorentzVector(0.,0.,0.,0.);

  // Neutrino
  nuVecCand = math::XYZTLorentzVector(0.,0.,0.,0.);

  // Anti-kt 5 Jets
  selectedJets.clear();
  BTagL.clear();
  BTagM.clear();
  BTagT.clear();
  nJets = 0;
  nBTagL = 0;
  nBTagM = 0;
  nBTagT = 0;
  
  cleanedak5Jets.clear();
  nCleanedak5Jets = 0;
  nCleanedBTagL = 0;
  nCleanedBTagM = 0;
  nCleanedBTagT = 0;
}


void tHqEvent::LeptonRec(){
  lepVecCand = BoostedUtils::GetPrimLepVec(input.selectedElectronsLoose,input.selectedMuonsLoose);
}


void tHqEvent::NeutrinoRec(){
  TVector2 metvec(input.pfMets[0].px(),input.pfMets[0].py());
  
  if(lepVecCand.Pt()<=0.001) return;
  
  BoostedUtils::GetNuVecs(lepVecCand,metvec,nu1VecCand,nu2VecCand);
}


// Anti-kt 5 Jets Methods    

void tHqEvent::ak5JetsRec(){

  selectedJets.clear();
  BTagL.clear();
  BTagM.clear();
  BTagT.clear();
  nJets = 0;
  nBTagL = 0;
  nBTagM = 0;
  nBTagT = 0;
  
  cleanedak5Jets.clear();
  nCleanedak5Jets = 0;
  nCleanedBTagL = 0;
  nCleanedBTagM = 0;
  nCleanedBTagT = 0;
  
  for(std::vector<pat::Jet>::const_iterator itJet=input.selectedJets.begin();itJet!=input.selectedJets.end();++itJet){
    if(lepVecCand.Pt()>0.001 && BoostedUtils::DeltaR(lepVecCand,itJet->p4())<.5) continue;
    
    nJets++;
    selectedJets.push_back(*itJet);
    cleanedak5Jets.push_back(*itJet);
  }
  for(std::vector<pat::Jet>::const_iterator itJet=selectedJets.begin();itJet!=selectedJets.end();++itJet){
    if(BoostedUtils::PassesCSV(*itJet,'L')) nBTagL++;
    BTagL.push_back(BoostedUtils::PassesCSV(*itJet,'L'));
    if(BoostedUtils::PassesCSV(*itJet,'M')) nBTagM++;
    BTagM.push_back(BoostedUtils::PassesCSV(*itJet,'M'));
    if(BoostedUtils::PassesCSV(*itJet,'T')) nBTagT++;
    BTagT.push_back(BoostedUtils::PassesCSV(*itJet,'T'));
  }
  
  nCleanedak5Jets = nJets;
  nCleanedBTagL = nBTagL;
  nCleanedBTagM = nBTagM;
  nCleanedBTagT = nBTagT;
}


void tHqEvent::ak5JetsClean(bool cleanHiggsCand, bool cleanTopHadCand, bool cleanTopLepCand){
  
  cleanedak5Jets.clear();
  nCleanedak5Jets = 0;
  nCleanedBTagL = 0;
  nCleanedBTagM = 0;
  nCleanedBTagT = 0;

  for(std::vector<pat::Jet>::const_iterator itJet=selectedJets.begin();itJet!=selectedJets.end();++itJet){
    int iJet = itJet-selectedJets.begin();

    nCleanedak5Jets++;
    if(BTagL[iJet]) nCleanedBTagL++;
    if(BTagM[iJet]) nCleanedBTagM++;
    if(BTagT[iJet]) nCleanedBTagT++;
    
    cleanedak5Jets.push_back(*itJet);
  }
}




const InputCollections& tHqEvent::GetInput(){
  return input;
} 


math::XYZTLorentzVector tHqEvent::GetLeptonVec(){
  return lepVecCand;
} 

    
math::XYZTLorentzVector tHqEvent::GetNeutrinoVec(){
  return nuVecCand;
}


std::vector<pat::Jet> tHqEvent::Getak5JetsAll(){
  return selectedJets;
}


int tHqEvent::GetNJets(){
  return nJets;
}


int tHqEvent::GetNBTagL(){
  return nBTagL;
}


int tHqEvent::GetNBTagM(){
  return nBTagM;
}

int tHqEvent::GetNBTagT(){
  return nBTagT;
}

float tHqEvent::GetAverageCSV(){
  
  float avgCSV = 0.;
  for(std::vector<pat::Jet>::const_iterator itJet=selectedJets.begin();itJet!=selectedJets.end();++itJet)
    avgCSV += fmax(itJet->bDiscriminator(btagger),0.);
    
  return avgCSV/nJets;
}


std::vector<pat::Jet> tHqEvent::Getak5JetsCleaned(){
  return cleanedak5Jets;
}


int tHqEvent::GetNCleanedak5Jets(){
  return nCleanedak5Jets;
}


int tHqEvent::GetNCleanedBTagL(){
  return nCleanedBTagL;
}


int tHqEvent::GetNCleanedBTagM(){
  return nCleanedBTagM;
}


int tHqEvent::GetNCleanedBTagT(){
  return nCleanedBTagT;
}


float tHqEvent::GetAverageCSVClean(){
  
  float avgCSV = 0.;
  
  for(size_t iJet=0;iJet<cleanedak5Jets.size();iJet++)
    avgCSV += fmax(cleanedak5Jets[iJet].bDiscriminator(btagger),0.);
    
  return avgCSV/nCleanedak5Jets;
}
