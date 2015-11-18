#include "tHqAnalysis/tHqAnalyzer/interface/TopGenVarProcessor.hpp"

using namespace std;

TopGenVarProcessor::TopGenVarProcessor (){}
TopGenVarProcessor::~TopGenVarProcessor (){}


void TopGenVarProcessor::Init(const InputCollections& input,VariableContainer& vars){

 
  vars.InitVar( "genevt_ttcc",-1,"I" );
  vars.InitVar( "genevt_ttbb",-1,"I" );
  vars.InitVar( "genevt_ttxid",-1,"I" );
  
  vars.InitVar( "n_tophad", 0, "I");
  vars.InitVar( "n_toplep", 0, "I");

  vars.InitVar( "n_top", 0, "I" );
  vars.InitVar( "n_topwdau", 0, "I" );
  vars.InitVar( "top_tpt",-9.,"F" );
  vars.InitVar( "top_twpt",-9.,"F" );
  vars.InitVar( "top_tbpt",-9.,"F" );
  vars.InitVars( "top_twdaupt",-9.,"n_topwdau" );
 
  vars.InitVar( "top_teta",-9.,"F" );
  vars.InitVar( "top_tweta",-9.,"F" );
  vars.InitVar( "top_tbeta",-9.,"F" );
  vars.InitVars( "top_twdaueta",-9.,"n_topwdau" );
  
  vars.InitVar( "top_tphi",-9.,"F" );
  vars.InitVar( "top_twphi",-9.,"F" );
  vars.InitVar( "top_tbphi",-9.,"F" );
  vars.InitVars( "top_twdauphi",-9.,"n_topwdau" );
  
  vars.InitVar( "top_tm",-9.,"F" );
  vars.InitVar( "top_twm",-9.,"F" );
  vars.InitVar( "top_tbm",-9.,"F" );
  vars.InitVars( "top_twdaum",-9.,"n_topwdau" );

 
  vars.InitVar( "n_topbar",0,"I" );
  vars.InitVar( "n_topbarwdau",0,"I" );
  vars.InitVar( "top_tbarpt",-9.,"F" );
  vars.InitVar( "top_tbarwpt",-9.,"F" );
  vars.InitVar( "top_tbarbpt",-9.,"F" );
  vars.InitVars( "top_tbarwdaupt",-9.,"n_topbarwdau" );
  
  vars.InitVar( "top_tbareta",-9.,"F" );
  vars.InitVar( "top_tbarweta",-9.,"F" );
  vars.InitVar( "top_tbarbeta",-9.,"F" );
  vars.InitVars( "top_tbarwdaueta",-9.,"n_topbarwdau" );
  
  vars.InitVar( "top_tbarphi",-9.,"F" );
  vars.InitVar( "top_tbarwphi",-9.,"F" );
  vars.InitVar( "top_tbarbphi",-9.,"F" );
  vars.InitVars( "top_tbarwdauphi",-9.,"n_topbarwdau" );
  
  vars.InitVar( "top_tbarm",-9.,"F" );
  vars.InitVar( "top_tbarwm",-9.,"F" );
  vars.InitVar( "top_tbarbm",-9.,"F" );
  vars.InitVars( "top_tbarwdaum",-9.,"n_topbarwdau" );

  vars.InitVar( "top_lepcharge",-9,"F" );
  
  vars.InitVar( "top_tbidx",-1,"I" );
  vars.InitVar( "tophadq1idx",-1,"I" );
  vars.InitVar( "tophadq2idx",-1,"I" );
  vars.InitVar( "top_tbarbidx",-1,"I" );


  vars.InitVars( "GenTopLep_B_GenJet_Pt",-9., "n_toplep" );
  vars.InitVars( "GenTopHad_B_GenJet_Pt",-9., "n_tophad");
  vars.InitVars( "GenTopLep_B_GenJet_Eta",-9., "n_toplep" );
  vars.InitVars( "GenTopHad_B_GenJet_Eta",-9., "n_tophad");
  vars.InitVars( "GenTopLep_B_GenJet_Phi",-9., "n_toplep" );
  vars.InitVars( "GenTopHad_B_GenJet_Phi",-9., "n_tophad");
  vars.InitVars( "GenTopLep_B_GenJet_M",-9., "n_toplep" );
  vars.InitVars( "GenTopHad_B_GenJet_M",-9., "n_tophad");


  vars.InitVars( "GenTopLep_B_Hadron_Pt",-9., "n_toplep" );
  vars.InitVars( "GenTopHad_B_Hadron_Pt",-9., "n_tophad");
  vars.InitVars( "GenTopLep_B_Hadron_Eta",-9., "n_toplep" );
  vars.InitVars( "GenTopHad_B_Hadron_Eta",-9., "n_tophad");
  vars.InitVars( "GenTopLep_B_Hadron_Phi",-9., "n_toplep" );
  vars.InitVars( "GenTopHad_B_Hadron_Phi",-9., "n_tophad");
  vars.InitVars( "GenTopLep_B_Hadron_M",-9., "n_toplep" );
  vars.InitVars( "GenTopHad_B_Hadron_M",-9., "n_tophad");


  initialized = true;
}

void TopGenVarProcessor::Process(const InputCollections& input,VariableContainer& vars){
  
  if(!initialized) cerr << "tree processor not initialized" << endl;
  cout << "Starting the TopGenVarProcessor for ttbar samples..." << endl;  
  
  int iBB = 0;
  int iCC = 0;
  
  std::cout << "Before filling ttx" << std::endl;

  if(input.sampleType == SampleType::ttbb) iBB = 3;
  if(input.sampleType == SampleType::ttb) iBB = 1;
  if(input.sampleType == SampleType::tt2b) iBB = 2;
  if(input.sampleType == SampleType::ttcc) iCC = 1;
  
  vars.FillVar( "genevt_ttcc",iCC );
  vars.FillVar( "genevt_ttbb",iBB );
  vars.FillVar( "genevt_ttxid",input.genTopEvt.GetTTxIdFromHelper());
  
  
  reco::GenParticle top;
  reco::GenParticle topw;
  reco::GenParticle topb;
  std::vector<reco::GenParticle> topwdau;
  reco::GenParticle topbar;
  reco::GenParticle topbarw;
  reco::GenParticle topbarb;
  std::vector<reco::GenParticle> topbarwdau;

  std::vector<reco::GenParticle> tophad;
  std::vector<reco::GenParticle> toplep;
  reco::GenParticle toplep_lep;
  
  std::vector<reco::GenParticle> q1;
  std::vector<reco::GenParticle> q2;

  if(input.genTopEvt.IsFilled()){
    top=input.genTopEvt.GetTop();
    topw=input.genTopEvt.GetWplus();
    topb=input.genTopEvt.GetTopDecayQuark();
    topwdau=input.genTopEvt.GetWplusDecayProducts();
    topbar=input.genTopEvt.GetTopBar();
    topbarw=input.genTopEvt.GetWminus();
    topbarb=input.genTopEvt.GetTopBarDecayQuark();
    topbarwdau=input.genTopEvt.GetWminusDecayProducts();


    tophad=input.genTopEvt.GetAllTopHads();
    toplep=input.genTopEvt.GetAllTopLeps();
    toplep_lep=input.genTopEvt.GetLepton();
    q1=input.genTopEvt.GetAllWQuarks();
    q2=input.genTopEvt.GetAllWAntiQuarks();


   vars.FillVar("n_top", 1);
   vars.FillVar("n_topbar", 1);
   vars.FillVar("n_topwdau", topwdau.size());
   vars.FillVar("n_topbarwdau", topbarwdau.size());
   
  }

  std::cout << "#W from top daughters: " << topwdau.size() << endl;
  std::cout << "#W from topbar daughters: " << topbarwdau.size() << endl;

  vars.FillVar("n_tophad", tophad.size());
  vars.FillVar("n_toplep", toplep.size());
  
  vars.FillVar("top_lepcharge", toplep_lep.charge());
  cout << "LepCharge: " << toplep_lep.charge() << endl;

  vector<math::XYZTLorentzVector> jetvecs = tHqUtils::GetJetVecs(input.selectedJets);
  
  vars.FillVar( "top_tpt",top.pt());
  vars.FillVar( "top_twpt",topw.pt());
  vars.FillVar( "top_tbpt",topb.pt());
  vars.FillVars( "top_twdaupt",0,topwdau[0].pt());
  vars.FillVars( "top_twdaupt",1,topwdau[1].pt()); //Fix 
  
  cout << "Pt of TopDaughters: " << topwdau[0].pt() << "     "   << topwdau[1].pt() << endl;

  vars.FillVar( "top_teta",top.eta());  
  vars.FillVar( "top_tweta",topw.eta());
  vars.FillVar( "top_tbeta",topb.eta());
  vars.FillVars( "top_twdaueta",0,topwdau[0].eta());
  vars.FillVars( "top_twdaueta",1,topwdau[1].eta());
  
  vars.FillVar( "top_tphi",top.phi());
  vars.FillVar( "top_twphi",topw.phi());
  vars.FillVar( "top_tbphi",topb.phi());
  vars.FillVars( "top_twdauphi",0,topwdau[0].phi());
  vars.FillVars( "top_twdauphi",1,topwdau[1].phi());
  
  vars.FillVar( "top_tm",top.mass());
  vars.FillVar( "top_twm",topw.mass());
  vars.FillVar( "top_tbm",topb.mass());
  vars.FillVars( "top_twdaum",0,topwdau[0].mass());
  vars.FillVars( "top_twdaum",1,topwdau[1].mass());
  
  
  int idxtopb = -1;
  double minDrTop = 999;
  for(std::vector<math::XYZTLorentzVector>::iterator itJetVec = jetvecs.begin() ; itJetVec != jetvecs.end(); ++itJetVec){
    if(tHqUtils::DeltaR(*itJetVec,topb.p4())<minDrTop){
      idxtopb = itJetVec-jetvecs.begin();
      minDrTop = tHqUtils::DeltaR(*itJetVec,topb.p4());
    }
  }
  
  if(minDrTop<.25){
    vars.FillVar( "top_tbidx",idxtopb);
  }
  

  vars.FillVar( "top_tbarpt",topbar.pt());
  vars.FillVar( "top_tbarwpt",topbarw.pt());
  vars.FillVar( "top_tbarbpt",topbarb.pt());
  vars.FillVars( "top_tbarwdaupt",0,topbarwdau[0].pt());
  vars.FillVars( "top_tbarwdaupt",1,topbarwdau[1].pt());

  vars.FillVar( "top_tbareta",topbar.eta()); 
  vars.FillVar( "top_tbarweta",topbarw.eta());
  vars.FillVar( "top_tbarbeta",topbarb.eta());
  vars.FillVars( "top_tbarwdaueta",0,topbarwdau[0].eta());
  vars.FillVars( "top_tbarwdaueta",1,topbarwdau[1].eta());

  vars.FillVar( "top_tbarphi",topbar.phi());
  vars.FillVar( "top_tbarwphi",topbarw.phi());
  vars.FillVar( "top_tbarbphi",topbarb.phi());
  vars.FillVars( "top_tbarwdauphi",0,topbarwdau[0].phi());
  vars.FillVars( "top_tbarwdauphi",1,topbarwdau[1].phi());

  vars.FillVar( "top_tbarm",topbar.mass());
  vars.FillVar( "top_tbarwm",topbarw.mass());
  vars.FillVar( "top_tbarbm",topbarb.mass());
  vars.FillVars( "top_tbarwdaum",0,topbarwdau[0].mass());
  vars.FillVars( "top_tbarwdaum",1,topbarwdau[1].mass());

  
  int idxtopbarb=-1;
  
  double minDrTopBarB = 999;
  for(size_t j=0; j<jetvecs.size(); j++){
    if(tHqUtils::DeltaR(jetvecs[j],topbarb.p4())<minDrTopBarB){
      idxtopbarb = j;
      minDrTopBarB = tHqUtils::DeltaR(jetvecs[j],topbarb.p4());
    }
  }
  
  if(minDrTopBarB<.25){
    vars.FillVar( "top_tbarbidx",idxtopbarb);
  }
  
  
  int idxq1=-1;
  int idxq2=-1;
  double minDrTopHadQ1 = 999;
  double minDrTopHadQ2 = 999;
  
  for(size_t i=0;i<tophad.size();i++){
    for(size_t j=0; j<jetvecs.size(); j++){
      if(tHqUtils::DeltaR(jetvecs[j],q1[i].p4())<minDrTopHadQ1){
	idxq1 = j;
	minDrTopHadQ1 = tHqUtils::DeltaR(jetvecs[j],q1[i].p4());
      }
      if(tHqUtils::DeltaR(jetvecs[j],q2[i].p4())<minDrTopHadQ2){
	idxq2 = j;
	minDrTopHadQ2 = tHqUtils::DeltaR(jetvecs[j],q2[i].p4());
      }
    }
        
    if(minDrTopHadQ1<.25){
      vars.FillVar( "tophadq1idx",idxq1);
    }
    if(minDrTopHadQ2<.25){
      vars.FillVar( "tophadq2dx",idxq2);
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
    vars.FillVar( "GenHiggs_B1_GenJet_M",b1_genjet.mass());
    vars.FillVar( "GenHiggs_B2_GenJet_M",b2_genjet.mass() );

    vars.FillVar( "GenHiggs_B1_Hadron_Pt",b1_hadron.pt() );
    vars.FillVar( "GenHiggs_B2_Hadron_Pt",b2_hadron.pt() );
    vars.FillVar( "GenHiggs_B1_Hadron_Eta",b1_hadron.eta() );
    vars.FillVar( "GenHiggs_B2_Hadron_Eta",b2_hadron.eta());
    vars.FillVar( "GenHiggs_B1_Hadron_Phi",b1_hadron.phi());
    vars.FillVar( "GenHiggs_B2_Hadron_Phi",b2_hadron.phi() );
    vars.FillVar( "GenHiggs_B1_Hadron_M",b1_hadron.mass());
    vars.FillVar( "GenHiggs_B2_Hadron_M",b2_hadron.mass() );
 

    for(uint i=0;i<bhad_genjet.size();i++){
      if(bhad_genjet[i].pt()>1){
	vars.FillVars( "GenTopHad_B_GenJet_Pt",i,bhad_genjet[i].pt() );
	vars.FillVars( "GenTopHad_B_GenJet_Eta",i,bhad_genjet[i].eta() );
	vars.FillVars( "GenTopHad_B_GenJet_Phi",i,bhad_genjet[i].phi());
	vars.FillVars( "GenTopHad_B_GenJet_M",i,bhad_genjet[i].mass());
      }
      if(bhad_hadron[i].pt()>1){
	vars.FillVars( "GenTopHad_B_Hadron_Pt",i,bhad_hadron[i].pt() );
	vars.FillVars( "GenTopHad_B_Hadron_Eta",i,bhad_hadron[i].eta() );
	vars.FillVars( "GenTopHad_B_Hadron_Phi",i,bhad_hadron[i].phi());
	vars.FillVars( "GenTopHad_B_Hadron_M",i,bhad_hadron[i].mass());
      }
    }
    
    for(uint i=0;i<blep_genjet.size();i++){
      if(blep_genjet[i].pt()>1){
	vars.FillVars( "GenTopLep_B_GenJet_Pt",i,blep_genjet[i].pt() );
	vars.FillVars( "GenTopLep_B_GenJet_Eta",i,blep_genjet[i].eta());
	vars.FillVars( "GenTopLep_B_GenJet_Phi",i,blep_genjet[i].phi());
	vars.FillVars( "GenTopLep_B_GenJet_M",i,blep_genjet[i].mass());

      }
      if(blep_hadron[i].pt()>1){
	vars.FillVars( "GenTopLep_B_Hadron_Pt",i,blep_hadron[i].pt() );
	vars.FillVars( "GenTopLep_B_Hadron_Eta",i,blep_hadron[i].eta());
	vars.FillVars( "GenTopLep_B_Hadron_Phi",i,blep_hadron[i].phi());
	vars.FillVars( "GenTopLep_B_Hadron_M",i,blep_hadron[i].mass());
      }
    }
  }
}
