#include "tHqAnalysis/tHqAnalyzer/interface/tHqGenVarProcessor.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/GentHqEvent.hpp"

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
  vars.InitVars( "Hdauid",-9.,  "nHdau");

  vars.InitVar( "nWdau", -9); 

  vars.InitVars( "Wdaupt",-9. , "nWdau" );
  vars.InitVars( "Wdaueta",-9.,  "nWdau");
  vars.InitVars( "Wdauphi",-9., "nWdau" );
  vars.InitVars( "Wdaum",-9.,  "nWdau");
  vars.InitVars( "Wdauid",-9.,  "nWdau");
 
  vars.InitVars( "hbbjtidx", -9, "nHdau"); 
  vars.InitVars( "Wjtidx", -9, "nWdau"); 
  vars.InitVar( "tjtidx", -9);
  vars.InitVar( "lqjtidx", -9);
  vars.InitVar( "sbjtidx", -9);

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
  if(input.GentHqEvent.IsFilled()){
    H=input.GentHqEvent.GetHiggs();
    W=input.GentHqEvent.GetW();
    t=input.GentHqEvent.GetTop();
    btop=input.GentHqEvent.GetTopDecayQuark();
    sb=input.GentHqEvent.GetSecondb();
    lq=input.GentHqEvent.GetLightQuark();
    Hdau=input.GentHqEvent.GetHiggsDecayProducts();
    Wdau=input.GentHqEvent.GetWDecayProducts();
  }

  std::cout << "#Hdaus: " << Hdau.size() << endl;
  std::cout << "#Wdaus: " << Wdau.size() << endl;

  
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
  
  
  vars.FillVar("btop_idx",idxblep);
    
    //    std::cout << "Found Jet to b from top: #" << idxblep << "  with dR= " << minDrTopLep << endl; 
  
  int idxhbb[2];
  
  double minDrB1 = 999;
  double minDrB2 = 999;
  
  for(std::vector<math::XYZTLorentzVector>::iterator itJetVec = jetvecs.begin() ; itJetVec != jetvecs.end(); ++itJetVec){
    assert(itJetVec->pt()>0);
    assert(Hdau[0].pt()>0);
    assert(Hdau[1].pt()>0);
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
    vars.FillVars("hbbjtidx",i,idxhbb[i]);

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
