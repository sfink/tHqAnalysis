#include "tHqAnalysis/tHqAnalyzer/interface/tHqGenVarProcessor.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/GentHqEvent.hpp"

using namespace std;

tHqGenVarProcessor::tHqGenVarProcessor (){}
tHqGenVarProcessor::~tHqGenVarProcessor (){}


void tHqGenVarProcessor::Init(const InputCollections& input,VariableContainer& vars){
 
  vars.InitVar( "Hpt",-9.,"F" );
  vars.InitVar( "Heta",-9.,"F" );
  vars.InitVar( "Hphi",-9., "F" );
  vars.InitVar( "Hm",-9., "F" );

  vars.InitVar( "Wpt",-9.,"F" );
  vars.InitVar( "Weta",-9., "F" );
  vars.InitVar( "Wphi",-9., "F" );
  vars.InitVar( "Wm",-9., "F" );

  vars.InitVar( "tpt",-9., "F" );
  vars.InitVar( "teta",-9., "F" );
  vars.InitVar( "tphi",-9., "F" );
  vars.InitVar( "tm",-9., "F" );

  vars.InitVar( "btoppt",-9., "F" );
  vars.InitVar( "btopeta",-9., "F" );
  vars.InitVar( "btopphi",-9., "F" );
  vars.InitVar( "btopm",-9., "F" );

  vars.InitVar( "sbpt",-9., "F" );
  vars.InitVar( "sbeta",-9., "F" );
  vars.InitVar( "sbphi",-9., "F" );
  vars.InitVar( "sbm",-9., "F" );

  vars.InitVar( "lqpt",-9., "F" );
  vars.InitVar( "lqeta",-9., "F" );
  vars.InitVar( "lqphi",-9.,"F" );
  vars.InitVar( "lqm",-9., "F" );

  vars.InitVar( "nHdau", -9, "I"); 

  vars.InitVars( "Hdaupt" , -9., "nHdau", 2 );
  vars.InitVars( "Hdaueta", -9., "nHdau", 2 );
  vars.InitVars( "Hdauphi", -9., "nHdau", 2 );
  vars.InitVars( "Hdaum",   -9., "nHdau", 2 );
  vars.InitVars( "Hdauid",  -9., "nHdau", 2 );

  vars.InitVar( "nWdau", -9, "I"); 

  vars.InitVars( "Wdaupt",-9. , "nWdau",2 );
  vars.InitVars( "Wdaueta",-9.,  "nWdau",2);
  vars.InitVars( "Wdauphi",-9., "nWdau",2 );
  vars.InitVars( "Wdaum",-9.,  "nWdau",2);
  vars.InitVars( "Wdauid",-9.,  "nWdau",2);
 
  vars.InitVars( "Hjtidx", -9, "nHdau",2); 
  vars.InitVars( "Wjtidx", -9, "nWdau",2); 
  vars.InitVar( "tjtidx", -9, "I");
  vars.InitVar( "lqjtidx", -9, "I");
  vars.InitVar( "sbjtidx", -9, "I");

  initialized = true;
}

void tHqGenVarProcessor::Process(const InputCollections& input,VariableContainer& vars){
  
  if(!initialized) cerr << "tree processor not initialized" << endl;
    
  reco::GenParticle H;
  reco::GenParticle W;
  reco::GenParticle t;
  reco::GenParticle btop;
  reco::GenParticle sb;
  reco::GenParticle lq;
  std::vector<reco::GenParticle> Hdau;
  std::vector<reco::GenParticle> Wdau;
  if(input.gentHqEvt.IsFilled()){
    H=input.gentHqEvt.GetHiggs();
    W=input.gentHqEvt.GetW();
    t=input.gentHqEvt.GetTop();
    btop=input.gentHqEvt.GetTopDecayQuark();
    sb=input.gentHqEvt.GetSecondb();
    lq=input.gentHqEvt.GetLightQuark();
    Hdau=input.gentHqEvt.GetHiggsDecayProducts();
    Wdau=input.gentHqEvt.GetWDecayProducts();
  }

  assert(Hdau.size()!=0);

  
  vector<math::XYZTLorentzVector> jetvecs = tHqUtils::GetJetVecs(input.selectedJets);

  vars.FillVar( "tpt",t.pt());
  vars.FillVar( "teta",t.eta());
  vars.FillVar( "tphi",t.phi());
  vars.FillVar( "tm",t.mass());
  vars.FillVar( "Wpt",W.pt());
  vars.FillVar( "Weta",W.eta());
  vars.FillVar( "Wphi",W.phi());
  vars.FillVar( "Wm",W.mass());
  vars.FillVar( "Hpt",H.pt());
  vars.FillVar( "Heta",H.eta());
  vars.FillVar( "Hphi",H.phi());
  vars.FillVar( "Hm",H.mass());
  vars.FillVar( "btoppt",btop.pt());
  vars.FillVar( "btopeta",btop.eta());
  vars.FillVar( "btopphi",btop.phi());
  vars.FillVar( "btopm",btop.mass());
  vars.FillVar( "sbpt",sb.pt());
  vars.FillVar( "sbeta",sb.eta());
  vars.FillVar( "sbphi",sb.phi());
  vars.FillVar( "sbm",sb.mass());
  vars.FillVar( "lqpt",lq.pt());
  vars.FillVar( "lqeta",lq.eta());
  vars.FillVar( "lqphi",lq.phi());
  vars.FillVar( "lqm",lq.mass());
  vars.FillVar( "nHdau",Hdau.size());
  vars.FillVar( "nWdau",Wdau.size());
  for(size_t i=0;i<2;i++){
    vars.FillVars( "Wdaupt",i,Wdau[i].pt());
    vars.FillVars( "Wdaueta",i,Wdau[i].eta());
    vars.FillVars( "Wdauphi",i,Wdau[i].phi());
    vars.FillVars( "Wdaum",i,Wdau[i].mass());
    vars.FillVars( "Wdauid",i,Wdau[i].pdgId());
    vars.FillVars( "Hdaupt",i,Hdau[i].pt());
    vars.FillVars( "Hdaueta",i,Hdau[i].eta());
    vars.FillVars( "Hdauphi",i,Hdau[i].phi());
    vars.FillVars( "Hdaum",i,Hdau[i].mass());
    vars.FillVars( "Hdauid",i,Hdau[i].pdgId());
  }

  int idxblep = -1;
  double minDrTopLep = 999;
  for(std::vector<math::XYZTLorentzVector>::iterator itJetVec = jetvecs.begin() ; itJetVec != jetvecs.end(); ++itJetVec){
    if(tHqUtils::DeltaR(*itJetVec,btop.p4())<minDrTopLep){
      idxblep = itJetVec-jetvecs.begin();
      minDrTopLep = tHqUtils::DeltaR(*itJetVec,btop.p4());
    }
  }
  
  if(minDrTopLep>4)
    idxblep=-99;
  
  
  vars.FillVar("tjtidx",idxblep);
    
    //    std::cout << "Found Jet to b from top: #" << idxblep << "  with dR= " << minDrTopLep << endl; 
  
  int idxhbb[2];
  
  double minDrB1 = 999;
  double minDrB2 = 999;
  
  for(std::vector<math::XYZTLorentzVector>::iterator itJetVec = jetvecs.begin() ; itJetVec != jetvecs.end(); ++itJetVec){
    if (itJetVec->pt()<0) continue;
    if (Hdau[0].pt()<0) continue;
    if (Hdau[1].pt()<0) continue;
    
    if(tHqUtils::DeltaR(*itJetVec,Hdau[0].p4())<minDrB1){
      idxhbb[0] = itJetVec-jetvecs.begin();
      minDrB1 = tHqUtils::DeltaR(*itJetVec,Hdau[0].p4());
    }
    if(tHqUtils::DeltaR(*itJetVec,Hdau[1].p4())<minDrB2){
      idxhbb[1] = itJetVec-jetvecs.begin();
      minDrB2 = tHqUtils::DeltaR(*itJetVec,Hdau[1].p4());
    }
  }
    
  if(minDrB1>.4){
    idxhbb[0]=-99.;
  }
  if(minDrB2>.4){
    idxhbb[1]=-99.;
  }
  
  for (int i =0; i<2; i++)
    vars.FillVars("Hjtidx",i,idxhbb[i]);

  int idxsb=-1;
  int idxlq=-1;
    
  double minDrSB = 999;
  double minDrLQ = 999;

  for(std::vector<math::XYZTLorentzVector>::iterator itJetVec = jetvecs.begin() ; itJetVec != jetvecs.end(); ++itJetVec){
    if (itJetVec->pt()>0) continue;
    if (sb.pt()>0) continue;
    if (lq.pt()>0) continue;
    if(tHqUtils::DeltaR(*itJetVec,sb.p4())<minDrSB){
      idxsb = itJetVec-jetvecs.begin();
      minDrSB = tHqUtils::DeltaR(*itJetVec,sb.p4());
    }
    if(tHqUtils::DeltaR(*itJetVec,lq.p4())<minDrLQ){
      idxlq = itJetVec-jetvecs.begin();
      minDrLQ = tHqUtils::DeltaR(*itJetVec,lq.p4());
    }
  }

  if(minDrSB>.4)
    idxsb=-99;
  
  vars.FillVar("sbjtidx",idxsb);

  if(minDrLQ>.4)
    idxlq=-99;
  
  vars.FillVar("lqjtidx",idxlq);

}
