#include "tHqAnalysis/tHqAnalyzer/interface/tHqUtils.hpp"

using namespace std;


std::string tHqUtils::GetAnalyzerPath(){
  char* cpath;
  cpath = getenv("CMSSW_BASE");
  std::string path = cpath;
  
  if (path.size() == 0){
    cerr << "could not find path of exe" << endl;
    return path;
  }
  else
    return path+"/src/tHqAnalysis/tHqAnalyzer";
}

TLorentzVector tHqUtils::GetTLorentzVector(const math::XYZTLorentzVector& vec){
  
  TLorentzVector result(vec.Px(),vec.Py(),vec.Pz(),vec.E());
  
  return result;
  
}


math::XYZTLorentzVector tHqUtils::GetXYZTLorentzVector(const TLorentzVector& vec){
  
  math::XYZTLorentzVector result(vec.Px(),vec.Py(),vec.Pz(),vec.E());
  
  return result;
  
}


bool tHqUtils::FirstIsLarger(float val1,float val2){
  return val1 > val2;
}


bool tHqUtils::FirstIsHarder(math::XYZTLorentzVector vec1,math::XYZTLorentzVector vec2){
  return vec1.Pt()>vec2.Pt();
}


bool tHqUtils::FirstJetIsHarder(pat::Jet jet1, pat::Jet jet2){
  return jet1.pt()>jet2.pt();
}


bool tHqUtils::FirstHasHigherCSV(pat::Jet jet1,pat::Jet jet2){
  return jet1.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags") > jet2.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
}

bool tHqUtils::FirstHasHigherCSVold(pat::Jet jet1,pat::Jet jet2){
  return jet1.bDiscriminator("combinedSecondaryVertexBJetTags") > jet2.bDiscriminator("combinedSecondaryVertexBJetTags");
}


float tHqUtils::DeltaEta(const math::XYZTLorentzVector& vec1,const math::XYZTLorentzVector& vec2){
  if(vec1.Pt()<0.001||vec2.Pt()<0.001) return -2;
  
  float deta = fabs(vec1.Eta()-vec2.Eta());
  
  return deta;
}


float tHqUtils::DeltaEta(const pat::Jet& jet1,const pat::Jet& jet2){
  if(jet1.pt()<0.001||jet2.pt()<0.001) return -2;
  
  math::XYZTLorentzVector vec1 = jet1.p4();
  math::XYZTLorentzVector vec2 = jet2.p4();
  
  float deta = fabs(vec1.Eta()-vec2.Eta());
  
  return deta;
}


float tHqUtils::DeltaPhi(const math::XYZTLorentzVector& vec1,const math::XYZTLorentzVector& vec2){
  if(vec1.Pt()<0.001||vec2.Pt()<0.001) return -2;
  
  float dphi = ROOT::Math::VectorUtil::DeltaPhi(vec1,vec2);
  
  return fabs(dphi);
}


float tHqUtils::DeltaPhi(const pat::Jet& jet1,const pat::Jet& jet2){
  if(jet1.pt()<0.001||jet2.pt()<0.001) return -2;
  
  math::XYZTLorentzVector vec1 = jet1.p4();
  math::XYZTLorentzVector vec2 = jet2.p4();
  
  float dphi = ROOT::Math::VectorUtil::DeltaPhi(vec1,vec2);
  
  return fabs(dphi);
}


float tHqUtils::DeltaR(const math::XYZTLorentzVector& vec1,const math::XYZTLorentzVector& vec2){
  if(vec1.Pt()<0.001||vec2.Pt()<0.001) return -2;
  
  float dr = ROOT::Math::VectorUtil::DeltaR(vec1,vec2);
  
  return dr;
}


float tHqUtils::DeltaR(const pat::Jet& jet1,const pat::Jet& jet2){
  if(jet1.pt()<0.001||jet2.pt()<0.001) return -2;
  
  math::XYZTLorentzVector vec1 = jet1.p4();
  math::XYZTLorentzVector vec2 = jet2.p4();
  
  float dr = ROOT::Math::VectorUtil::DeltaR(vec1,vec2);
  
  return dr;
}


float tHqUtils::DeltaKt(const math::XYZTLorentzVector& vec1,const math::XYZTLorentzVector& vec2){
  if(vec1.Pt()<0.001 || vec2.Pt()<0.001) return -2;
  
  float dr = ROOT::Math::VectorUtil::DeltaR(vec1,vec2);
  float ptmin=min(vec1.Pt(),vec2.Pt());
  
  return sqrt(dr*dr*ptmin*ptmin);
}


float tHqUtils::DeltaKt(const pat::Jet& jet1,const pat::Jet& jet2){
  if(jet1.pt()<0.001||jet2.pt()<0.001) return -2;
  
  math::XYZTLorentzVector vec1 = jet1.p4();
  math::XYZTLorentzVector vec2 = jet2.p4();
  
  float dr = ROOT::Math::VectorUtil::DeltaR(vec1,vec2);
  float ptmin=min(vec1.Pt(),vec2.Pt());
  
  return sqrt(dr*dr*ptmin*ptmin);
}


float tHqUtils::CosThetaStar(const math::XYZTLorentzVector& vec1, const math::XYZTLorentzVector& vec2){
  if(vec1.Pt()<0.001||vec2.Pt()<0.001) return -2;
  
  TLorentzVector sumVec = tHqUtils::GetTLorentzVector(vec1+vec2);
  TVector3 cmboost = -sumVec.BoostVector();
  
  TLorentzVector boostedvec1 = GetTLorentzVector(vec1);
  boostedvec1.Boost(cmboost);
  
  return cos( sumVec.Angle(boostedvec1.Vect()) );
}


float tHqUtils::CosThetaCM(const math::XYZTLorentzVector& vec,const math::XYZTLorentzVector& boostVec){
  if(vec.Pt()<0.001||boostVec.Pt()<0.001) return -2;
  
  TLorentzVector vec_ = GetTLorentzVector(vec);
  TLorentzVector boostVec_ = GetTLorentzVector(boostVec);
  
  TVector3 zBoost = TVector3(0,0,-boostVec_.BoostVector().Pz());
  TLorentzVector boostedVec = vec_;
  boostedVec.Boost(zBoost);
  float theta = boostedVec.Theta();
  return TMath::Cos(theta);
}


float tHqUtils::GetMuondBetaIso(const pat::Muon& muon){
  float iso_photon = muon.photonIso();
  float iso_neutrals = muon.neutralHadronIso();
  float iso_sumPUPt = muon.pfIsolationR04().sumPUPt;
  float iso_charged = muon.chargedHadronIso();
  float pt = muon.pt();

  float iso = std::max(0.,(iso_photon+iso_neutrals-0.5*iso_sumPUPt)+iso_charged)/pt;

  return iso;
}


std::vector<math::XYZTLorentzVector> tHqUtils::GetGenParticleVecs(const std::vector<reco::GenParticle>& genParticles){
  std::vector<math::XYZTLorentzVector> genParticleVecs;
  
  for(std::vector<reco::GenParticle>::const_iterator itPart=genParticles.begin();itPart!=genParticles.end();++itPart){
    genParticleVecs.push_back(itPart->p4());
  }
  
  return genParticleVecs;
}


bool tHqUtils::MCContainsTTbar(const std::vector<reco::GenParticle>& genParticles){
  bool foundT=false;
  bool foundTbar=false;
  for(size_t i=0; i<genParticles.size();i++){
    if(genParticles[i].pdgId()==6) foundT=true;
    if(genParticles[i].pdgId()==-6) foundTbar=true;
  }
  return foundT&&foundTbar;
}

bool tHqUtils::MCContainsST(const std::vector<reco::GenParticle>& genParticles){
  bool foundT=false;
  bool foundTbar=false;
  for(size_t i=0; i<genParticles.size();i++){
    if(genParticles[i].pdgId()==6) foundT=true;
    if(genParticles[i].pdgId()==-6) foundTbar=true;
  }
  return (foundT&&!foundTbar)||(!foundT&&foundTbar);
}



bool tHqUtils::MCContainsHiggs(const std::vector<reco::GenParticle>& genParticles){
  for(size_t i=0; i<genParticles.size();i++){
    if(genParticles[i].pdgId()==25) return true;
  }
  return false;
}


void tHqUtils::GetttHMCParticles(const std::vector<reco::GenParticle>& genParticles, std::vector<reco::GenParticle>& tophad, std::vector<reco::GenParticle>& whad, std::vector<reco::GenParticle>& bhad, std::vector<reco::GenParticle>& q1, std::vector<reco::GenParticle>& q2, std::vector<reco::GenParticle>& toplep, std::vector<reco::GenParticle>& wlep, std::vector<reco::GenParticle>& blep, std::vector<reco::GenParticle>& lep, std::vector<reco::GenParticle>& nu, reco::GenParticle& higgs, reco::GenParticle& b1, reco::GenParticle& b2){
  std::vector<reco::GenParticle> leptonsFromWplus;
  std::vector<reco::GenParticle> leptonsFromWminus;
  std::vector<reco::GenParticle> quarksFromWplus;
  std::vector<reco::GenParticle> quarksFromWminus;
  std::vector<reco::GenParticle> quarksFromT;
  std::vector<reco::GenParticle> quarksFromTbar;
  std::vector<reco::GenParticle> bquarksFromHiggs;
  reco::GenParticle Wplus;
  reco::GenParticle Wminus;
  reco::GenParticle top;
  reco::GenParticle topbar;
  
  for(std::vector<reco::GenParticle>::const_iterator itPart = genParticles.begin();itPart != genParticles.end();itPart++) {
    if(itPart->numberOfMothers()==0) continue;
    if(itPart->mother()->pdgId()==6    && abs(itPart->pdgId())<6) quarksFromT.push_back(*itPart);
    if(itPart->mother()->pdgId()==-6   && abs(itPart->pdgId())<6) quarksFromTbar.push_back(*itPart);
    if(itPart->mother()->pdgId()==24   && abs(itPart->pdgId())<6) quarksFromWplus.push_back(*itPart);
    if(itPart->mother()->pdgId()==-24  && abs(itPart->pdgId())<6) quarksFromWminus.push_back(*itPart);
    if(itPart->mother()->pdgId()==24   && abs(itPart->pdgId())>=11&&abs(itPart->pdgId())<=16) leptonsFromWplus.push_back(*itPart);
    if(itPart->mother()->pdgId()==-24  && abs(itPart->pdgId())>=11&&abs(itPart->pdgId())<=16) leptonsFromWminus.push_back(*itPart);
    if(itPart->mother()->pdgId()==25   && abs(itPart->pdgId())==5) bquarksFromHiggs.push_back(*itPart);
    
    if(itPart->mother()->pdgId()==6    && itPart->pdgId()==24)    Wplus = *itPart;
    if(itPart->mother()->pdgId()==-6   && itPart->pdgId()==-24)   Wminus = *itPart;
    
    if(itPart->pdgId()==25 || fabs(itPart->pdgId())==6){
      bool lastHiggs = true;
      bool lastTop = true;
      bool lastTopBar = true;

      for(size_t iDaughter = 0;iDaughter<itPart->numberOfDaughters();++iDaughter){
        if(itPart->daughter(iDaughter)->pdgId()==25) lastHiggs = false;
        if(itPart->daughter(iDaughter)->pdgId()==6) lastTop = false;
        if(itPart->daughter(iDaughter)->pdgId()==-6) lastTopBar = false;
      }

      if(lastHiggs  && itPart->pdgId()==25) higgs = *itPart;
      if(lastTop    && itPart->pdgId()==6)  top = *itPart;
      if(lastTopBar && itPart->pdgId()==-6) topbar = *itPart; 
    }
  }
  
  if(leptonsFromWplus.size()==2 && quarksFromT.size()>=1){
    toplep.push_back(top);
    wlep.push_back(Wplus);
    blep.push_back(quarksFromT[0]);
    if(leptonsFromWplus[0].pdgId()%2!=0){
      lep.push_back(leptonsFromWplus[0]);
      nu.push_back(leptonsFromWplus[1]);
    }
    else{
      lep.push_back(leptonsFromWplus[1]);
      nu.push_back(leptonsFromWplus[0]);
    }
  }
  if(leptonsFromWminus.size()==2 && quarksFromTbar.size()>=1){
    toplep.push_back(topbar);
    wlep.push_back(Wminus);
    blep.push_back(quarksFromTbar[0]);
    if(leptonsFromWminus[0].pdgId()%2!=0){
      lep.push_back(leptonsFromWminus[0]);
      nu.push_back(leptonsFromWminus[1]);
    }
    else{
      lep.push_back(leptonsFromWminus[1]);
      nu.push_back(leptonsFromWminus[0]);
    }
  }
  if(quarksFromWplus.size()==2 && quarksFromT.size()>=1){
    tophad.push_back(top);
    whad.push_back(Wplus);
    bhad.push_back(quarksFromT[0]);
    q1.push_back(quarksFromWplus[0]);
    q2.push_back(quarksFromWplus[1]);
  }
  if(quarksFromWminus.size()==2 && quarksFromTbar.size()>=1){
    tophad.push_back(topbar);
    whad.push_back(Wminus);
    bhad.push_back(quarksFromTbar[0]);
    q1.push_back(quarksFromWminus[0]);
    q2.push_back(quarksFromWminus[1]);
  }
  if(bquarksFromHiggs.size()==2){
    b1 = bquarksFromHiggs[0];
    b2 = bquarksFromHiggs[1];
  }  
}


void tHqUtils::GetttHMCVecs(const std::vector<reco::GenParticle>& genParticles, std::vector<math::XYZTLorentzVector>& tophadvecs, std::vector<math::XYZTLorentzVector>& whadvecs, std::vector<math::XYZTLorentzVector>& bhadvecs, std::vector<math::XYZTLorentzVector>& q1vecs, std::vector<math::XYZTLorentzVector>& q2vecs, std::vector<math::XYZTLorentzVector>& toplepvecs, std::vector<math::XYZTLorentzVector>& wlepvecs, std::vector<math::XYZTLorentzVector>& blepvecs, std::vector<math::XYZTLorentzVector>& lepvecs, std::vector<math::XYZTLorentzVector>& nuvecs, math::XYZTLorentzVector& higgsvec, math::XYZTLorentzVector& b1vec, math::XYZTLorentzVector& b2vec){
  
  std::vector<reco::GenParticle> tophad;
  std::vector<reco::GenParticle> whad;
  std::vector<reco::GenParticle> bhad;
  std::vector<reco::GenParticle> q1;
  std::vector<reco::GenParticle> q2;
  std::vector<reco::GenParticle> toplep;
  std::vector<reco::GenParticle> wlep;
  std::vector<reco::GenParticle> blep;
  std::vector<reco::GenParticle> lep;
  std::vector<reco::GenParticle> nu;
  reco::GenParticle higgs;
  reco::GenParticle b1;
  reco::GenParticle b2;
  
  GetttHMCParticles(genParticles,tophad,whad,bhad,q1,q2,toplep,wlep,blep,lep,nu,higgs,b1,b2);
  
  tophadvecs = GetGenParticleVecs(tophad);
  whadvecs = GetGenParticleVecs(whad);
  bhadvecs = GetGenParticleVecs(bhad);
  q1vecs = GetGenParticleVecs(q1);
  q2vecs = GetGenParticleVecs(q2);
  toplepvecs = GetGenParticleVecs(toplep);
  wlepvecs = GetGenParticleVecs(wlep);
  blepvecs = GetGenParticleVecs(blep);
  lepvecs = GetGenParticleVecs(lep);
  nuvecs = GetGenParticleVecs(nu);
  higgsvec = higgs.p4();
  b1vec = b1.p4();
  b2vec = b2.p4();
}


bool tHqUtils::IsAnyTriggerBitFired(const std::vector<string> & targetTriggers, const edm::TriggerResults& triggerResults){
  
  // check to see if you passed the trigger by looping over the bits
  // looking for your bit
  // and see if it is 1
  
  for(vector<string>::const_iterator iTarget = targetTriggers.begin();iTarget != targetTriggers.end();iTarget++) {
    
    if(*iTarget == "None") return true;
    
    // if this is the right name and the bit is set to one
    int TriggerID = triggerResults.find(*iTarget);
    
    if(triggerResults.accept(TriggerID)==1) return true;
    
  }// end for each target
  
  return false;
  
}


std::vector<math::XYZTLorentzVector> tHqUtils::GetLepVecs(const std::vector<pat::Electron>& selectedElectrons, const std::vector<pat::Muon> selectedMuons){
  
  std::vector<math::XYZTLorentzVector> leptonVecs;
  
  for(std::vector<pat::Electron>::const_iterator itEle = selectedElectrons.begin(); itEle != selectedElectrons.end(); ++itEle){
    leptonVecs.push_back(itEle->p4());
  }
  for(std::vector<pat::Muon>::const_iterator itMu = selectedMuons.begin(); itMu != selectedMuons.end(); ++itMu){
    leptonVecs.push_back(itMu->p4());
  }
  
  std::sort(leptonVecs.begin(), leptonVecs.end(),tHqUtils::FirstIsHarder);
  
  return leptonVecs;
  
} 


math::XYZTLorentzVector tHqUtils::GetPrimLepVec(const std::vector<pat::Electron>& selectedElectrons,const std::vector<pat::Muon> selectedMuons){
  
  std::vector<math::XYZTLorentzVector> leptonVecs = GetLepVecs(selectedElectrons,selectedMuons);
  
  if(leptonVecs.size()==0){
    std::cerr<< "No PrimLep Found!" << std::endl;
    return math::XYZTLorentzVector();
  }
  
  return leptonVecs[0];
}


void tHqUtils::GetNuVec(const math::XYZTLorentzVector& lepvec, const TVector2& metvec, math::XYZTLorentzVector& nuvec){
  
  float nu_e  = std::sqrt(metvec.Mod2());
  float nu_px = metvec.Px();
  float nu_py = metvec.Py();

  float lep_e  = lepvec.E();
  float lep_pt = lepvec.Pt();
  float lep_px = lepvec.Px();
  float lep_py = lepvec.Py();
  float lep_pz = lepvec.Pz();

  float mw  = 80.43;

  //  float nupt_temp  = std::sqrt((lep_px+nu_px)*(lep_px+nu_px)+
  //(lep_py+nu_py)*(lep_py+nu_py));

  float mtsq= ((lep_pt+nu_e)*(lep_pt+nu_e) -
	       (lep_px+nu_px)*(lep_px+nu_px) -
	       (lep_py+nu_py)*(lep_py+nu_py));

  float mt  = (mtsq>0.0) ? std::sqrt(mtsq) : 0.0;

  float A, B, Csq, C, S1, S2;
  float scf(1.0);

  if (mt < mw) {
    A = mw*mw/2.0;
  }
  else {
    A = mt*mt/2.0;
    float k = nu_e*lep_pt - nu_px*lep_px - nu_py*lep_py;
    if (k <0.0001) k = 0.0001;
    scf = 0.5*(mw*mw)/k;
    nu_px *= scf;
    nu_py *= scf;
    nu_e  = sqrt(nu_px*nu_px + nu_py*nu_py);
  }

  B  = nu_px*lep_px + nu_py*lep_py;
  Csq= 1.0 + nu_e*nu_e*(lep_pz*lep_pz-lep_e*lep_e)/(A+B)/(A+B);
  C  = (Csq>0.0) ? std::sqrt(Csq) : 0.0;
  S1 = (-(A+B)*lep_pz + (A+B)*lep_e*C)/(lep_pz*lep_pz-lep_e*lep_e);
  S2 = (-(A+B)*lep_pz - (A+B)*lep_e*C)/(lep_pz*lep_pz-lep_e*lep_e);

  float nu_pz = (std::abs(S1) < std::abs(S2)) ? S1 : S2;
  nu_e = sqrt(metvec.Mod2()+nu_pz*nu_pz);
//  nuvec.SetPtEtaPhiM(nupt,nueta,nuphi,0);
  nuvec.SetPxPyPzE(nu_px,nu_py,nu_pz,nu_e);

  return;
}


vector<math::XYZTLorentzVector> tHqUtils::GetJetVecs(const std::vector<pat::Jet>& jets){
  std::vector<math::XYZTLorentzVector> jetVecs;
  
  for(std::vector<pat::Jet>::const_iterator itJet=jets.begin();itJet!=jets.end();++itJet){
    jetVecs.push_back(itJet->p4());
  }
  
  return jetVecs;
}


bool tHqUtils::PassesCSV(const pat::Jet& jet, const char workingPoint){
  
  float CSVLv2wp = 0.423;
  float CSVMv2wp = 0.814;
  float CSVTv2wp = 0.941;
  
  float csvValue = jet.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
  
  switch(workingPoint){
    case 'L': if(csvValue > CSVLv2wp){ return true; } break;
    case 'M': if(csvValue > CSVMv2wp){ return true; } break;
    case 'T': if(csvValue > CSVTv2wp){ return true; } break;
    case '-': return true; break;
  }
  
  return false;
}


float tHqUtils::GetClosestJetIDs(int& idJet1, int& idJet2, const std::vector<pat::Jet>& jets){
  
  float minDr = 9.;
  for(std::vector<pat::Jet>::const_iterator itJet1 = jets.begin();itJet1 != jets.end();++itJet1){
    for(std::vector<pat::Jet>::const_iterator itJet2 = itJet1+1;itJet2 != jets.end();++itJet2){
      math::XYZTLorentzVector jetVec1 = itJet1->p4();
      math::XYZTLorentzVector jetVec2 = itJet2->p4();
      
      if(tHqUtils::DeltaR(jetVec1,jetVec2)<minDr){
	      idJet1 = itJet1 - jets.begin();
	      idJet2 = itJet2 - jets.begin();
	      minDr = tHqUtils::DeltaR(jetVec1,jetVec2);
      }
    }
  }
  
  return minDr;
}


float tHqUtils::GetClosestLepJetID(int& idJet, const math::XYZTLorentzVector& lepVec, const std::vector<pat::Jet>& jets){
  
  float minDr = 9.;
  for(std::vector<pat::Jet>::const_iterator itJet = jets.begin(); itJet != jets.end(); ++itJet){
    math::XYZTLorentzVector jetVec = itJet->p4();
    if(DeltaR(lepVec,jetVec)<minDr){
	    idJet = itJet - jets.begin();
	    minDr = DeltaR(lepVec,jetVec);
    }
  }

  return minDr;
}


float tHqUtils::GetJetAverageJetEtaMax(const std::vector<pat::Jet>& jets1, const std::vector<pat::Jet>& jets2){
  int count=0;
  float avgval=0.;
  for(std::vector<pat::Jet>::const_iterator itJet1 = jets1.begin(); itJet1 != jets1.end(); ++itJet1){
    for(std::vector<pat::Jet>::const_iterator itJet2 = itJet1+1; itJet2 != jets1.end(); ++itJet2){
      avgval += itJet1->eta()-itJet2->eta();
      count++;
    }
  }
  avgval /= count;
  
  float imax = 0.;
  float etamax=0.;
  for(std::vector<pat::Jet>::const_iterator itJet = jets2.begin(); itJet != jets2.end(); ++itJet){
    imax = abs(itJet->eta()-avgval);
    if(imax>etamax) {
      etamax = imax;
    }
  }
  return etamax;
}


void tHqUtils::GetFoxWolframMoments(std::vector<math::XYZTLorentzVector> jetVecs, float &h0, float &h1, float &h2, float &h3, float &h4) {
  //
  // Aplanarity and sphericity
  //
  
  float eVis = 0.0;
  for(std::vector<math::XYZTLorentzVector>::iterator itJetVec = jetVecs.begin();itJetVec != jetVecs.end();++itJetVec){
    eVis += itJetVec->energy();
  }
  h0 = 0.0;
  h1 = 0.0;
  h2 = 0.0;
  h3 = 0.0;
  h4 = 0.0;
  for(std::vector<math::XYZTLorentzVector>::iterator itJetVec1 = jetVecs.begin();itJetVec1 != jetVecs.end();++itJetVec1){
    for(std::vector<math::XYZTLorentzVector>::iterator itJetVec2 = itJetVec1+1;itJetVec2 != jetVecs.end();++itJetVec2){
      float costh = ROOT::Math::VectorUtil::CosTheta(*itJetVec1,*itJetVec2);
      float p0 = 1.0;
      float p1 = costh;
      float p2 = 0.5*(3.0*costh*costh - 1.0);
      float p3 = 0.5*(5.0*costh*costh - 3.0*costh);
      float p4 = 0.125*(35.0*costh*costh*costh*costh - 30.0*costh*costh + 3.0);
      float pipj = itJetVec1->P()*itJetVec2->P();
      h0 += (pipj/(eVis*eVis))*p0;
      h1 += (pipj/(eVis*eVis))*p1;
      h2 += (pipj/(eVis*eVis))*p2;
      h3 += (pipj/(eVis*eVis))*p3;
      h4 += (pipj/(eVis*eVis))*p4;
    }
  }
  return;
}


void tHqUtils::GetAplanaritySphericity(math::XYZTLorentzVector leptonVec, math::XYZTLorentzVector metVec, std::vector<math::XYZTLorentzVector> jetVecs, float &aplanarity, float &sphericity){
  //
  // Aplanarity and sphericity
  //
  float mxx = leptonVec.Px()*leptonVec.Px() + metVec.Px()*metVec.Px();
  float myy = leptonVec.Py()*leptonVec.Py() + metVec.Py()*metVec.Py();
  float mzz = leptonVec.Pz()*leptonVec.Pz() + metVec.Pz()*metVec.Pz();
  float mxy = leptonVec.Px()*leptonVec.Py() + metVec.Px()*metVec.Py();
  float mxz = leptonVec.Px()*leptonVec.Pz() + metVec.Px()*metVec.Pz();
  float myz = leptonVec.Py()*leptonVec.Pz() + metVec.Px()*metVec.Pz();
  
  for(std::vector<math::XYZTLorentzVector>::iterator itJetVec=jetVecs.begin();itJetVec!=jetVecs.end();++itJetVec){
    mxx += itJetVec->Px()*itJetVec->Px();
    myy += itJetVec->Py()*itJetVec->Py();
    mzz += itJetVec->Pz()*itJetVec->Pz();
    mxy += itJetVec->Px()*itJetVec->Py();
    mxz += itJetVec->Px()*itJetVec->Pz();
    myz += itJetVec->Py()*itJetVec->Pz();
  }
  float sum = mxx + myy + mzz;
  mxx /= sum;
  myy /= sum;
  mzz /= sum;
  mxy /= sum;
  mxz /= sum;
  myz /= sum;
  
  TMatrix tensor(3,3);
  tensor(0,0) = mxx;
  tensor(1,1) = myy;
  tensor(2,2) = mzz;
  tensor(0,1) = mxy;
  tensor(1,0) = mxy;
  tensor(0,2) = mxz;
  tensor(2,0) = mxz;
  tensor(1,2) = myz;
  tensor(2,1) = myz;
  TVector eigenval(3);
  tensor.EigenVectors(eigenval);
  
  sphericity = 3.0*(eigenval(1)+eigenval(2))/2.0;
  aplanarity = 3.0*eigenval(2)/2.0;
  
  return;
}


