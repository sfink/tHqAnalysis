#include "tHqAnalysis/tHqAnalyzer/interface/tHqGenVarProcessor.hpp"

using namespace std;

tHqGenVarProcessor::tHqGenVarProcessor (){}
tHqGenVarProcessor::~tHqGenVarProcessor (){}


void tHqGenVarProcessor::Init(const InputCollections& input,VariableContainer& vars){
 
  vars.InitVar( "Hpt",-9. );
  vars.InitVar( "Heta",-9. );
  vars.InitVar( "Hphi",-9. );
  vars.InitVar( "Hm",-9. );

  vars.InitVar( "Wpt",-9. );
  vars.InitVar( "Weta",-9. );
  vars.InitVar( "Wphi",-9. );
  vars.InitVar( "Wm",-9. );

  vars.InitVar( "tpt",-9. );
  vars.InitVar( "teta",-9. );
  vars.InitVar( "tphi",-9. );
  vars.InitVar( "tm",-9. );

  vars.InitVar( "btoppt",-9. );
  vars.InitVar( "btopeta",-9. );
  vars.InitVar( "btopphi",-9. );
  vars.InitVar( "btopm",-9. );

  vars.InitVar( "sbpt",-9. );
  vars.InitVar( "sbeta",-9. );
  vars.InitVar( "sbphi",-9. );
  vars.InitVar( "sbm",-9. );

  vars.InitVar( "lqpt",-9. );
  vars.InitVar( "lqeta",-9. );
  vars.InitVar( "lqphi",-9. );
  vars.InitVar( "lqm",-9. );

  vars.InitVar( "nHdau", -9); 

  vars.InitVars( "Hdaupt",-9. , "nHdau" );
  vars.InitVars( "Hdaueta",-9.,  "nHdau");
  vars.InitVars( "Hdauphi",-9., "nHdau" );
  vars.InitVars( "Hdaum",-9.,  "nHdau");

  vars.InitVar( "nWdau", -9); 

  vars.InitVars( "Wdaupt",-9. , "nWdau" );
  vars.InitVars( "Wdaueta",-9.,  "nWdau");
  vars.InitVars( "Wdauphi",-9., "nWdau" );
  vars.InitVars( "Wdaum",-9.,  "nWdau");
 
  vars.InitVar( "Hb1pt",-9. );
  vars.InitVar( "Hb1eta",-9. );
  vars.InitVar( "Hb1phi",-9. );
  vars.InitVar( "Hb1idx",-9. );
 
  vars.InitVar( "Hb2pt",-9. );
  vars.InitVar( "Hb2eta",-9. );
  vars.InitVar( "Hb2phi",-9. );
  vars.InitVar( "Hb2idx",-9. );
 

  initialized = true;
}

void tHqGenVarProcessor::Process(const InputCollections& input,VariableContainer& vars){
  
  if(!initialized) cerr << "tree processor not initialized" << endl;

  
  std::vector<reco::GenParticle> H;
  std::vector<reco::GenParticle> W;
  std::vector<reco::GenParticle> t;
  std::vector<reco::GenParticle> btop;
  std::vector<reco::GenParticle> sb;
  std::vector<reco::GenParticle> lq;
  std::vector<reco::GenParticle> Hdau;
  std::vector<reco::GenParticle> Wdau;
  if(input.gentHqEvt.IsFilled()){
    H=input.gentHqEvt.GetAllHiggs();
    W=input.gentHqEvt.GetAllWs();
    t=input.gentHqEvt.GetAllTops();
    btop=input.gentHqEvt.GetAllBtops();
    sb=input.gentHqEvt.GetAllSecondBs();
    lq=input.gentHqEvt.GetAllLightQuarks();
    Hdau=input.gentHqEvt.GetAllHdaus();
    Wdau=input.gentHqEvt.GetAllWdaus();
  }

  std::cout << "#H: " << H.size() << endl;
  std::cout << "#W: " << W.size() << endl;
  std::cout << "#t: " << t.size() << endl;
  std::cout << "#btops: " << btop.size() << endl;
  std::cout << "#2ndBs: " << sb.size() << endl;
  std::cout << "#LightQuarks: " << lq.size() << endl;
  std::cout << "#Hdaus: " << Hdau.size() << endl;
  std::cout << "#Wdaus: " << Wdau.size() << endl;

  reco::GenParticle b1;
  reco::GenParticle b2;
  reco::GenParticle decProd1;
  reco::GenParticle decProd2;

  if(Hdau.size()>2)std::cout<<"MORE THAN TWO HIGGS PRODUCTS"<<std::endl;
  bool dfirst=true;
  for(auto p =Hdau.begin(); p!=Hdau.end(); p++){
    if(p->pdgId()==5) b1=*p;
    if(p->pdgId()==-5) b2=*p;
    if(dfirst){
      decProd1=*p;
      dfirst=false;
    }
    else{
      decProd2=*p;
    }
  }
  
  if(decProd1.pt()>0.){
    vars.FillVar("Hb1pt",decProd1.pt());
    vars.FillVar("Hb2pt",decProd2.pt());
    vars.FillVar("Hb1eta",decProd1.eta());
    vars.FillVar("Hb2eta",decProd2.eta());
    vars.FillVar("Hb1phi",decProd1.phi());
    vars.FillVar("Hb2phi",decProd2.phi());
    vars.FillVar("Hb1id",decProd1.pdgId());
    vars.FillVar("Hb2id",decProd2.pdgId());
    
  }
  
  vector<math::XYZTLorentzVector> jetvecs = tHqUtils::GetJetVecs(input.selectedJets);

  for(size_t i=0;i<t.size();i++){
    vars.FillVars( "tpt",i,t[i].pt());
    vars.FillVars( "teta",i,t[i].eta());
    vars.FillVars( "tphi",i,t[i].phi());
    vars.FillVars( "tm",i,t[i].M());
    vars.FillVars( "Wpt",i,W[i].pt());
    vars.FillVars( "Weta",i,W[i].eta());
    vars.FillVars( "Wphi",i,W[i].phi());
    vars.FillVars( "Wm",i,W[i].M());
    vars.FillVars( "btoppt",i,btop[i].pt());
    vars.FillVars( "btopeta",i,btop[i].eta());
    vars.FillVars( "btopphi",i,btop[i].phi());
    vars.FillVars( "btopm",i,btop[i].M());
    vars.FillVars( "Wdaupt",i,Wdau[i].pt());
    vars.FillVars( "Wdaueta",i,Wdau[i].eta());
    vars.FillVars( "Wdauphi",i,Wdau[i].phi());
    vars.FillVars( "Wdaum",i,Wdau[i].M());
    
    int idxblep = -1;
    double minDrTopLep = 999;
    for(std::vector<math::XYZTLorentzVector>::iterator itJetVec = jetvecs.begin() ; itJetVec != jetvecs.end(); ++itJetVec){
      if(tHqUtils::DeltaR(*itJetVec,btop[i].p4())<minDrTopLep){
        idxblep = itJetVec-jetvecs.begin();
        minDrTopLep = tHqUtils::DeltaR(*itJetVec,btop[i].p4());
      }
    }
    
    if(minDrTopLep<.25){
      vars.FillVars( "btop_idx",i,idxblep);
    }
    std::cout << "Found Jet to b from top: #" << idxblep << "  with dR= " << minDrTopLep << endl; 
  }
  
  if(H.pt()>0.){
    vars.FillVar( "Hpt",H.pt());
    vars.FillVar( "Heta",H.eta());
    vars.FillVar( "Hphi",H.phi());
    vars.FillVar( "Hm",H.M());
  }
  if(b1.pt()>0.){
    vars.FillVar("Hb1pt",b1.pt());
    vars.FillVar("Hb2pt",b2.pt());
    vars.FillVar("Hb1eta",b1.eta());
    vars.FillVar("Hb2eta",b2.eta());
    vars.FillVar("Hb1phi",b1.phi());
    vars.FillVar("Hb2phi",b2.phi());
    
    int idxb1=-1;
    int idxb2=-1;
    
    double minDrB1 = 999;
    double minDrB2 = 999;

    for(std::vector<math::XYZTLorentzVector>::iterator itJetVec = jetvecs.begin() ; itJetVec != jetvecs.end(); ++itJetVec){
      assert(itJetVec->pt()>0);
      assert(b1.pt()>0);
      assert(b2.pt()>0);
      if(tHqUtils::DeltaR(*itJetVec,b1.p4())<minDrB1){
        idxb1 = itJetVec-jetvecs.begin();
        minDrB1 = tHqUtils::DeltaR(*itJetVec,b1.p4());
      }
      if(tHqUtils::DeltaR(*itJetVec,b2.p4())<minDrB2){
        idxb2 = itJetVec-jetvecs.begin();
        minDrB2 = tHqUtils::DeltaR(*itJetVec,b2.p4());
      }
    }
    
    if(minDrB1<.25){
      vars.FillVar( "Hb1_idx",idxb1);
    }
    if(minDrB2<.25){
      vars.FillVar( "Hb2_idx",idxb2);
    }
  }

  if(sb.pt()>0.){
    vars.FillVar( "sbpt",sb.pt());
    vars.FillVar( "sbeta",sb.eta());
    vars.FillVar( "sbphi",sb.phi());
    vars.FillVar( "sbm",sb.M());
  }
  if(lq.pt()>0.){
    vars.FillVar("lqpt",lq.pt());
    vars.FillVar("lqeta",lq.eta());
    vars.FillVar("lqphi",lq.phi());
    
    int idxsb=-1;
    int idxlq=-1;
    
    double minDrSB = 999;
    double minDrLQ = 999;

    for(std::vector<math::XYZTLorentzVector>::iterator itJetVec = jetvecs.begin() ; itJetVec != jetvecs.end(); ++itJetVec){
      assert(itJetVec->pt()>0);
      assert(sb.pt()>0);
      assert(lq.pt()>0);
      if(tHqUtils::DeltaR(*itJetVec,sb.p4())<minDrSB){
        idxsb = itJetVec-jetvecs.begin();
        minDrSB = tHqUtils::DeltaR(*itJetVec,sb.p4());
      }
      if(tHqUtils::DeltaR(*itJetVec,lq.p4())<minDrLQ){
        idxLQ = itJetVec-jetvecs.begin();
        minDrLQ = tHqUtils::DeltaR(*itJetVec,lq.p4());
      }
    }
    
    if(minDrSB<.25){
      vars.FillVar( "sb_idx",idxsb);
    }
    if(minDrLQ<.25){
      vars.FillVar( "lq_idx",idxlq);
    }
  }

  
  std::cout << "IsFilled: "<< input.genTopEvt.IsFilled() << " | TTxIsFilled: " << input.genTopEvt.TTxIsFilled() << " | IsSemiLepton: " << input.genTopEvt.IsSemiLepton() << endl;

  if(input.genTopEvt.IsFilled()&&input.genTopEvt.TTxIsFilled()&&input.genTopEvt.IsSemiLepton()){
    std::vector<reco::GenJet> bhad_genjet=input.genTopEvt.GetAllTopHadBGenJets();
    std::vector<reco::GenJet> blep_genjet=input.genTopEvt.GetAllTopLepBGenJets();
    reco::GenJet b1_genjet=input.genTopEvt.GetHiggsBBarGenJet();
    reco::GenJet b2_genjet=input.genTopEvt.GetHiggsBGenJet();

    std::vector<reco::GenParticle> bhad_hadron=input.genTopEvt.GetAllTopHadBHadrons();
    std::vector<reco::GenParticle> blep_hadron=input.genTopEvt.GetAllTopLepBHadrons();
    reco::GenParticle b1_hadron=input.genTopEvt.GetHiggsBBarHadron();
    reco::GenParticle b2_hadron=input.genTopEvt.GetHiggsBHadron();

    vars.FillVar( "GenHiggs_B1_GenJet_Pt",b1_genjet.pt() );
    vars.FillVar( "GenHiggs_B2_GenJet_Pt",b2_genjet.pt() );
    vars.FillVar( "GenHiggs_B1_GenJet_Eta",b1_genjet.eta() );
    vars.FillVar( "GenHiggs_B2_GenJet_Eta",b2_genjet.eta());
    vars.FillVar( "GenHiggs_B1_GenJet_Phi",b1_genjet.phi());
    vars.FillVar( "GenHiggs_B2_GenJet_Phi",b2_genjet.phi() );

    vars.FillVar( "GenHiggs_B1_Hadron_Pt",b1_hadron.pt() );
    vars.FillVar( "GenHiggs_B2_Hadron_Pt",b2_hadron.pt() );
    vars.FillVar( "GenHiggs_B1_Hadron_Eta",b1_hadron.eta() );
    vars.FillVar( "GenHiggs_B2_Hadron_Eta",b2_hadron.eta());
    vars.FillVar( "GenHiggs_B1_Hadron_Phi",b1_hadron.phi());
    vars.FillVar( "GenHiggs_B2_Hadron_Phi",b2_hadron.phi() );
 
    for(uint i=0;i<bhad_genjet.size();i++){
      if(bhad_genjet[i].pt()>1){
	vars.FillVars( "GenTopHad_B_GenJet_Pt",i,bhad_genjet[i].pt() );
	vars.FillVars( "GenTopHad_B_GenJet_Eta",i,bhad_genjet[i].eta() );
	vars.FillVars( "GenTopHad_B_GenJet_Phi",i,bhad_genjet[i].phi());
      }
      if(bhad_hadron[i].pt()>1){
	vars.FillVars( "GenTopHad_B_Hadron_Pt",i,bhad_hadron[i].pt() );
	vars.FillVars( "GenTopHad_B_Hadron_Eta",i,bhad_hadron[i].eta() );
	vars.FillVars( "GenTopHad_B_Hadron_Phi",i,bhad_hadron[i].phi());
      }
    }
    
    for(uint i=0;i<blep_genjet.size();i++){
      if(blep_genjet[i].pt()>1){
	vars.FillVars( "GenTopLep_B_GenJet_Phi",i,blep_genjet[i].phi());
	vars.FillVars( "GenTopLep_B_GenJet_Pt",i,blep_genjet[i].pt() );
	vars.FillVars( "GenTopLep_B_GenJet_Eta",i,blep_genjet[i].eta());
      }
      if(blep_hadron[i].pt()>1){
	vars.FillVars( "GenTopLep_B_Hadron_Pt",i,blep_hadron[i].pt() );
	vars.FillVars( "GenTopLep_B_Hadron_Eta",i,blep_hadron[i].eta());
	vars.FillVars( "GenTopLep_B_Hadron_Phi",i,blep_hadron[i].phi());
      }
    }
  }
}
