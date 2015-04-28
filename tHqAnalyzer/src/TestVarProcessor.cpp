#include "tHqAnalysis/tHqAnalyzer/interface/TestVarProcessor.hpp"

using namespace std;

TestVarProcessor::TestVarProcessor(){}
TestVarProcessor::~TestVarProcessor(){}


void TestVarProcessor::Init(const InputCollections& input,VariableContainer& vars){





  vars.InitVar( "njt", "I");
  vars.InitVar( "njt15", "I"); 
  vars.InitVar( "nlmu", "I");
  vars.InitVar( "nlel", "I");
  vars.InitVar( "nmu", "I");
  vars.InitVar( "nel", "I");
  vars.InitVar( "nvetomu", "I");
  vars.InitVar( "nvetoel", "I");

  vars.InitVar( "nbtagl", "I");  //New btag multiplicity variables
  vars.InitVar( "nbtagm", "I");  //
  vars.InitVar( "nbtagt", "I");  //

  vars.InitVar( "npv", "I");

  vars.InitVars( "jte","njt" );
  vars.InitVars( "jtpt","njt" );
  vars.InitVars( "jtphi","njt" );
  vars.InitVars( "jteta","njt" );
  vars.InitVars( "jtcsvt","njt" );

  vars.InitVars( "jtntracks","njt" );
  vars.InitVars( "jtarea","njt" );
  vars.InitVars( "jtpull","njt" );
  vars.InitVars( "jtcharge","njt" );
  vars.InitVars( "jtid","njt" );
  vars.InitVars( "jtchhadmult","njt" ); 

  vars.InitVars( "jtndaughters","njt" );
  vars.InitVars( "jtchmult","njt" );
  vars.InitVars( "jtnhadronfrac","njt" );
  vars.InitVars( "jthfhadronfrac","njt" );
  vars.InitVars( "jtjer","njt" );

  vars.InitVars( "jtgenflv","njt" );
  vars.InitVars( "jtgenpt","njt" );
  vars.InitVars( "jtgenphi","njt" );
  vars.InitVars( "jtgeneta","njt" );
  vars.InitVars( "jtgene","njt" );





  /*
  vars.InitVar( "Evt_E_PrimaryLepton" );
  vars.InitVar( "Evt_M_PrimaryLepton" );
  vars.InitVar( "Evt_Pt_PrimaryLepton" );
  vars.InitVar( "Evt_Eta_PrimaryLepton" );
  vars.InitVar( "Evt_Phi_PrimaryLepton" );
  
  vars.InitVars( "LooseLepton_E","N_LooseLeptons" );
  vars.InitVars( "LooseLepton_M","N_LooseLeptons" );
  vars.InitVars( "LooseLepton_Pt","N_LooseLeptons" );
  vars.InitVars( "LooseLepton_Eta","N_LooseLeptons" );
  vars.InitVars( "LooseLepton_Phi","N_LooseLeptons" );
  */



  vars.InitVars( "lmue","nlmu" );
  vars.InitVars( "lmupt","nlmu" );
  vars.InitVars( "lmueta","nlmu" );
  vars.InitVars( "lmuphi","nlmu" );
  vars.InitVars( "lmuiso","nlmu" );
  vars.InitVars( "lmucharge","nlmu" );

  vars.InitVars( "mue","nmu" );
  vars.InitVars( "mupt","nmu" );
  vars.InitVars( "mueta","nmu" );
  vars.InitVars( "muphi","nmu" );
  vars.InitVars( "muiso","nmu" );
  vars.InitVars( "mucharge","nmu" );


  vars.InitVars( "ele","nlel" );
  vars.InitVars( "elpt","nlel" );
  vars.InitVars( "eleta","nlel" );
  vars.InitVars( "elphi","nlel" );
  vars.InitVars( "eliso","nlel" );
  vars.InitVars( "elcharge","nlel" );
  
  vars.InitVar( "met" );
  vars.InitVar( "metphi" );
  vars.InitVar( "meteta" );

  vars.InitVar("m3"); //new var
  vars.InitVar("mtw"); //new var
  vars.InitVar("sumHtTotal");
  vars.InitVar("sumHt");
 
  vars.InitVars( "jtcsvt","njt" );
 
  vars.InitVar( "aplanarity" );
  vars.InitVar( "sphericity" );
   
  vars.InitVar( "wolframh0" ); 
  vars.InitVar( "wolframh1" );
  vars.InitVar( "wolframh2" );
  vars.InitVar( "wolframh3" );
  vars.InitVar( "wolframh4" );


  vars.InitVar( "coststh" );
  vars.InitVar( "costst" );

  initialized=true;
}

void TestVarProcessor::Process(const InputCollections& input,VariableContainer& vars){
  if(!initialized) cerr << "tree processor not initialized" << endl;

  const char* btagger="combinedInclusiveSecondaryVertexV2BJetTags";
  std::vector<pat::Jet> selectedTaggedJets;
  std::vector<pat::Jet> selectedTaggedJetsT;
  std::vector<pat::Jet> selectedTaggedJetsL;
  std::vector<pat::Jet> selectedUntaggedJets;
  for(std::vector<pat::Jet>::const_iterator itJet = input.selectedJets.begin(); itJet != input.selectedJets.end(); ++itJet){
    if(tHqUtils::PassesCSV(*itJet, 'M')){
      selectedTaggedJets.push_back(*itJet);
    }
    else{
      selectedUntaggedJets.push_back(*itJet);
    }
    if(tHqUtils::PassesCSV(*itJet, 'L')){
      selectedTaggedJetsL.push_back(*itJet);
    }
    if(tHqUtils::PassesCSV(*itJet, 'T')){
      selectedTaggedJetsT.push_back(*itJet);
    }
  }
  
  // Fill Multiplicity Variables
  vars.FillVar( "npv",input.selectedPVs.size());  
  vars.FillVar( "njt",input.selectedJets.size());
  vars.FillVar( "njt15",input.selectedJetsLoose.size());
  vars.FillVar( "nel",input.selectedElectrons.size());  
  vars.FillVar( "nlel",input.selectedElectronsLoose.size());  
  vars.FillVar( "nmu",input.selectedMuons.size());  
  vars.FillVar( "nlmu",input.selectedMuonsLoose.size());  
  
  // Fill Jet Variables
  // All Jets
  for(std::vector<pat::Jet>::const_iterator itJet = input.selectedJets.begin() ; itJet != input.selectedJets.end(); ++itJet){
    int iJet = itJet - input.selectedJets.begin();
    vars.FillVars( "jte",iJet,itJet->energy() );
    vars.FillVars( "jtpt",iJet,itJet->pt() );
    vars.FillVars( "jteta",iJet,itJet->eta() );
    vars.FillVars( "jtphi",iJet,itJet->phi() );
    vars.FillVars( "jtcsvt",iJet,fmax(itJet->bDiscriminator(btagger),-.1) );        

    vars.FillVars( "jtntracks",iJet,itJet->associatedTracks().size() );
    vars.FillVars( "jtarea",iJet,itJet->jetArea() );
    vars.FillVars( "jtpull",iJet,itJet->associatedTracks().size() );  // to implement 
    vars.FillVars( "jtcharge",iJet,itJet->jetCharge() );
    vars.FillVars( "jtid",iJet,itJet->associatedTracks().size() ); //to implement
    //    vars.FillVars( "jtndaughters",iJet,itJet.numberOfDaughters());
    //    vars.FillVars( "jtjer",iJet,itJet.pfSpecific().mChargedHadronEnergy ); //to implement
    vars.FillVars( "jtgenflv",iJet,itJet->partonFlavour() );        
    /*
    if(itJet.isPFJet()){
      vars.FillVars( "jtchhadmult",iJet,itJet.pfSpecific().mChargedHadronEnergy );
      vars.FillVars( "jtchmult",iJet,itJet.pfSpecific().mChargedMultiplicity );
      vars.FillVars( "jtnhadronfrac",iJet,itJet.pfSpecific().mChargedHadronEnergy ); //to implement /missing raw e
      vars.FillVars( "jthfhadronfrac",iJet,itJet.pfSpecific().mChargedHadronEnergy ); //to implement
      } */
  }

  
  for(std::vector<reco::GenJet>::const_iterator itGenJet = input.selectedGenJets.begin() ; itGenJet != input.selectedGenJets.end(); ++itGenJet){
    int iGenJet = itGenJet - input.selectedGenJets.begin();
    
    vars.FillVars( "jtgene",iGenJet,itGenJet->energy() );
    vars.FillVars( "jtgenpt",iGenJet,itGenJet->pt() );
    vars.FillVars( "jtgeneta",iGenJet,itGenJet->eta() );
    vars.FillVars( "jtgenphi",iGenJet,itGenJet->phi() );
  }


  math::XYZTLorentzVector primLepVec = math::XYZTLorentzVector();
  if(input.selectedElectrons.size()>0 || input.selectedMuons.size()>0){
    primLepVec = tHqUtils::GetPrimLepVec(input.selectedElectrons,input.selectedMuons);
  }

  
  // Fill Lepton Variables
  for(std::vector<pat::Electron>::const_iterator itEle = input.selectedElectronsLoose.begin(); itEle != input.selectedElectronsLoose.end(); ++itEle){
    int iEle = itEle - input.selectedElectronsLoose.begin();
    vars.FillVars( "ele",iEle,itEle->energy() );
    vars.FillVars( "elpt",iEle,itEle->pt() );
    vars.FillVars( "eleta",iEle,itEle->eta() );
    vars.FillVars( "elphi",iEle,itEle->phi() ); 
    vars.FillVars( "eliso",iEle,itEle->phi() ); //FIXME
    vars.FillVars( "elcharge",iEle,itEle->charge() );
  }

  for(std::vector<pat::Muon>::const_iterator itMu = input.selectedMuonsLoose.begin(); itMu != input.selectedMuonsLoose.end(); ++itMu){
    int iMu = itMu - input.selectedMuonsLoose.begin();
    vars.FillVars( "lmue",iMu,itMu->energy() );
    vars.FillVars( "lmupt",iMu,itMu->pt() );
    vars.FillVars( "lmueta",iMu,itMu->eta() );
    vars.FillVars( "lmuphi",iMu,itMu->phi() );
    vars.FillVars( "lmuiso",iMu,tHqUtils::GetMuondBetaIso(*itMu) );  //FIXME
    vars.FillVars( "lmucharge",iMu,itMu->charge() ); 
  }

  for(std::vector<pat::Muon>::const_iterator itMu = input.selectedMuons.begin(); itMu != input.selectedMuons.end(); ++itMu){
    int iMu = itMu - input.selectedMuons.begin();
    vars.FillVars( "mue",iMu,itMu->energy() );
    vars.FillVars( "mupt",iMu,itMu->pt() );
    vars.FillVars( "mueta",iMu,itMu->eta() );
    vars.FillVars( "muphi",iMu,itMu->phi() );
    vars.FillVars( "muiso",iMu,tHqUtils::GetMuondBetaIso(*itMu) );  //FIXME
    vars.FillVars( "mucharge",iMu,itMu->charge() ); 
  }
  
  vars.FillVar( "met",input.pfMets[0].pt() );
  vars.FillVar( "metphi",input.pfMets[0].phi() );
  vars.FillVar( "meteta",input.pfMets[0].eta() );
  
  std::vector<math::XYZTLorentzVector> jetvecs = tHqUtils::GetJetVecs(input.selectedJets);
  math::XYZTLorentzVector metvec = input.pfMets[0].p4();
  
  // Fill M3 Variables
  float m3_helper = -1.;
  float maxpt=-1;
  for(std::vector<math::XYZTLorentzVector>::iterator itJetVec1 = jetvecs.begin() ; itJetVec1 != jetvecs.end(); ++itJetVec1){
    for(std::vector<math::XYZTLorentzVector>::iterator itJetVec2 = itJetVec1+1 ; itJetVec2 != jetvecs.end(); ++itJetVec2){
      for(std::vector<math::XYZTLorentzVector>::iterator itJetVec3 = itJetVec2+1 ; itJetVec3 != jetvecs.end(); ++itJetVec3){
    
        math::XYZTLorentzVector m3vec = *itJetVec1 + *itJetVec2 + *itJetVec3;
        
	      if(m3vec.Pt() > maxpt){
	        maxpt = m3vec.Pt();
	        m3_helper = m3vec.M();
	      }
      } 
    }
  }
  vars.FillVar("m3",m3_helper);
  
  // Fill MTW
  float mtw_helper = -1.;
  if(input.selectedElectrons.size()>0 || input.selectedMuons.size()>0){
    mtw_helper = sqrt(2*(primLepVec.Pt()*input.pfMets[0].pt() - primLepVec.Px()*input.pfMets[0].px() - primLepVec.Py()*input.pfMets[0].py()));
  }
  vars.FillVar("mtw",mtw_helper);
  
  // Fill Ht Variables
  float ht = 0.;
  float htjets = 0.;
  for(std::vector<pat::Jet>::const_iterator itJet = input.selectedJets.begin() ; itJet != input.selectedJets.end(); ++itJet){
    ht += itJet->pt();
    htjets += itJet->pt();
  }
  for(std::vector<pat::Electron>::const_iterator itEle = input.selectedElectronsLoose.begin(); itEle != input.selectedElectronsLoose.end(); ++itEle){
    ht += itEle->pt();
  }
  for(std::vector<pat::Muon>::const_iterator itMu = input.selectedMuonsLoose.begin(); itMu != input.selectedMuonsLoose.end(); ++itMu){
    ht += itMu->pt();
  }
  ht += input.pfMets[0].pt();
  
  vars.FillVar("sumHtTotal",ht);
  vars.FillVar("sumHt",htjets);
   
  // Fill Number of b Tags
  
  vars.FillVar( "nbtagm",selectedTaggedJets.size() );  
  vars.FillVar( "nbtagl",selectedTaggedJetsL.size() );  
  vars.FillVar( "nbtagt",selectedTaggedJetsT.size() );
  
  // Fill CSV Variables
  // All Jets
  std::vector<float> csvJets;
  for(std::vector<pat::Jet>::const_iterator itJet = input.selectedJets.begin() ; itJet != input.selectedJets.end(); ++itJet){
    csvJets.push_back(fmax(itJet->bDiscriminator(btagger),-.1));
  }
  
  std::sort(csvJets.begin(),csvJets.end(),tHqUtils::FirstIsLarger);
  //  vars.FillVar("Evt_CSV_Min",csvJets.size()>0 ? csvJets.back() : -.1);
    
  // Event Shape Variables
  // Fox Wolfram Moments
  float h0,h1,h2,h3,h4;
  h0=-9;
  h1=-9;
  h2=-9;
  h3=-9;
  h4=-9;
  tHqUtils::GetFoxWolframMoments(jetvecs, h0,h1,h2,h3,h4);
  vars.FillVar( "wolframh0", h0 );
  vars.FillVar( "wolframh1", h1 );
  vars.FillVar( "wolframh2", h2 );
  vars.FillVar( "wolframh3", h3 );
  vars.FillVar( "wolframh4", h4 );
  
  // Aplanarity and Sphericity;
  float aplanarity,sphericity;
  tHqUtils::GetAplanaritySphericity(primLepVec, metvec, jetvecs, aplanarity, sphericity) ;
  vars.FillVar( "aplanarity", aplanarity );
  vars.FillVar( "sphericity", sphericity );
  /*
  // Event Angle Variables 
  float drmax_lj=-1;
  float detamax_lj=-1;
  float drmax_j1j=-1;
  float drmax_j2j=-1;
  float drmax_j3j=-1;
  float drmax_j4j=-1;
  float detamax_j1j=-1;
  float detamax_j2j=-1;
  float detamax_j3j=-1;
  float detamax_j4j=-1;
  float costhetamax_jcm=-1;
  for(std::vector<math::XYZTLorentzVector>::iterator itJetVec = jetvecs.begin();itJetVec != jetvecs.end();++itJetVec){
    int iJetVec = itJetVec - jetvecs.begin();
    
    float c_lj=-1.5;
    float c_j1j=-1.5;
    float c_j2j=-1.5;
    float c_j3j=-1.5;
    float c_j4j=-1.5;
    float deta_lj=-1.;
    float deta_j1j=-1.;
    float deta_j2j=-1.;
    float deta_j3j=-1.;
    float deta_j4j=-1.;
    float dr_lj=-1.;
    float dr_j1j=-1.;
    float dr_j2j=-1.;
    float dr_j3j=-1.;
    float dr_j4j=-1.;
    float dkt_lj=-50.;
    float dkt_j1j=-50.;
    float dkt_j2j=-50.;
    float dkt_j3j=-50.;
    float dkt_j4j=-50.;
    
    if(primLepVec.Pt()>1){
      deta_lj = tHqUtils::DeltaEta(*itJetVec,primLepVec);
      dr_lj = tHqUtils::DeltaR(*itJetVec,primLepVec);
      dkt_lj = tHqUtils::DeltaKt(*itJetVec,primLepVec);
      c_lj = tHqUtils::CosThetaStar(*itJetVec,primLepVec);
    }
    if(jetvecs.size()>0){
      deta_j1j = tHqUtils::DeltaEta(*itJetVec,jetvecs[0]);
      dr_j1j = tHqUtils::DeltaR(*itJetVec,jetvecs[0]);
      dkt_j1j = tHqUtils::DeltaKt(*itJetVec,jetvecs[0]);
      c_j1j = tHqUtils::CosThetaStar(*itJetVec,jetvecs[0]);
    }
    if(jetvecs.size()>1){
      deta_j2j = tHqUtils::DeltaEta(*itJetVec,jetvecs[1]);
      dr_j2j = tHqUtils::DeltaR(*itJetVec,jetvecs[1]);
      dkt_j2j = tHqUtils::DeltaKt(*itJetVec,jetvecs[1]);
      c_j2j = tHqUtils::CosThetaStar(*itJetVec,jetvecs[1]);
    }
    if(jetvecs.size()>2){
      deta_j3j = tHqUtils::DeltaEta(*itJetVec,jetvecs[2]);
      dr_j3j = tHqUtils::DeltaR(*itJetVec,jetvecs[2]);
      dkt_j3j = tHqUtils::DeltaKt(*itJetVec,jetvecs[2]);
      c_j3j = tHqUtils::CosThetaStar(*itJetVec,jetvecs[2]);
    }
    if(jetvecs.size()>3){
      deta_j4j = tHqUtils::DeltaEta(*itJetVec,jetvecs[3]);
      dr_j4j = tHqUtils::DeltaR(*itJetVec,jetvecs[3]);
      dkt_j4j = tHqUtils::DeltaKt(*itJetVec,jetvecs[3]);
      c_j4j = tHqUtils::CosThetaStar(*itJetVec,jetvecs[3]);
    }
    
    vars.FillVars("Jet_Deta_Lepton",iJetVec,deta_lj);
    vars.FillVars("Jet_Deta_Jet1",iJetVec,deta_j1j);
    vars.FillVars("Jet_Deta_Jet2",iJetVec,deta_j2j);
    vars.FillVars("Jet_Deta_Jet3",iJetVec,deta_j3j);
    vars.FillVars("Jet_Deta_Jet4",iJetVec,deta_j4j);
    
    vars.FillVars("Jet_Dr_Lepton",iJetVec,dr_lj);
    vars.FillVars("Jet_Dr_Jet1",iJetVec,dr_j1j);
    vars.FillVars("Jet_Dr_Jet2",iJetVec,dr_j2j);
    vars.FillVars("Jet_Dr_Jet3",iJetVec,dr_j3j);
    vars.FillVars("Jet_Dr_Jet4",iJetVec,dr_j4j);
    
    vars.FillVars("Jet_Dkt_Lepton",iJetVec,dkt_lj);
    vars.FillVars("Jet_Dkt_Jet1",iJetVec,dkt_j1j);
    vars.FillVars("Jet_Dkt_Jet2",iJetVec,dkt_j2j);
    vars.FillVars("Jet_Dkt_Jet3",iJetVec,dkt_j3j);
    vars.FillVars("Jet_Dkt_Jet4",iJetVec,dkt_j4j);
    
    vars.FillVars("Jet_CosThetaStar_Lepton",iJetVec,c_lj);
    vars.FillVars("Jet_CosThetaStar_Jet1",iJetVec,c_j1j);
    vars.FillVars("Jet_CosThetaStar_Jet2",iJetVec,c_j2j);
    vars.FillVars("Jet_CosThetaStar_Jet3",iJetVec,c_j3j);
    vars.FillVars("Jet_CosThetaStar_Jet4",iJetVec,c_j4j);
    
    if(drmax_lj < dr_lj){
      drmax_lj = dr_lj;
    }
    if(detamax_lj < deta_lj){
      detamax_lj = deta_lj;
    }
    if(drmax_j1j < dr_j1j){
      drmax_j1j = dr_j1j;
    }
    if(drmax_j2j < dr_j2j){
      drmax_j2j = dr_j2j;
    }
    if(drmax_j3j < dr_j3j){
      drmax_j3j = dr_j3j;
    }
    if(drmax_j4j < dr_j4j){
      drmax_j4j = dr_j4j;
    }
    if(detamax_j1j < deta_j1j){
      detamax_j1j = deta_j1j;
    }
    if(detamax_j2j < deta_j2j){
      detamax_j2j = deta_j2j;
    }
    if(detamax_j3j < deta_j3j){
      detamax_j3j = deta_j3j;
    }
    if(detamax_j4j < deta_j4j){
      detamax_j4j = deta_j4j;
    }
    
    float costheta_jcm= tHqUtils::CosThetaCM(*itJetVec,p4all);
    vars.FillVars("Jet_CosTheta_cm",iJetVec,costheta_jcm  );
    if(costhetamax_jcm<fabs(costheta_jcm)){
      costhetamax_jcm=fabs(costheta_jcm);
    }
  }
  vars.FillVar("Evt_Jet_Drmax_Lepton",drmax_lj);
  vars.FillVar("Evt_Jet_Detamax_Lepton",detamax_lj);
  vars.FillVar("Evt_Jet_Drmax_Jet1",drmax_j1j);
  vars.FillVar("Evt_Jet_Detamax_Jet1",detamax_j1j);
  vars.FillVar("Evt_Jet_CosThetamax_cm",costhetamax_jcm );
  
  // Ohio Variables
  std::vector<pat::Jet> selectedJetsLooseExclusive;
  for(std::vector<pat::Jet>::const_iterator itJet = input.selectedJetsLoose.begin() ; itJet != input.selectedJetsLoose.end(); ++itJet){
    if(itJet->pt() >= 30) continue;
    selectedJetsLooseExclusive.push_back(*itJet);
  }
  */

}
