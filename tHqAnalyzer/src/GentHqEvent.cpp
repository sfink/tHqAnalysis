#include "tHqAnalysis/tHqAnalyzer/interface/GentHqEvent.hpp"

GentHqEvent::GentHqEvent (){
  isFilled=false;
}
GentHqEvent::~GentHqEvent(){}

bool GentHqEvent::IsFilled() const{
  return isFilled;
}

void GentHqEvent::Fill(const std::vector<reco::GenParticle>& prunedGenParticles){

  //  for(auto p=prunedGenParticles.begin(); p!=prunedGenParticles.end(); p++){
  // GentHqEvent::PrintParticle(*p); 
  //}
  
  //std::cout << " " << std::endl;
  //std::cout << " " << std::endl;

  for(auto p=prunedGenParticles.begin(); p!=prunedGenParticles.end(); p++){

    
    //    GentHqEvent::PrintParticle(*p);
    // top quark
    if (abs(p->pdgId())==6 && p->isLastCopy()){
      top=*p;
    }
    // Higgs
    if (abs(p->pdgId())==25 && p->isLastCopy()){
      higgs=*p;
    }
    // Wboson - tricky: can also come from the Higgs decay in an inclusive sample
    if (abs(p->pdgId())==24 && p->isLastCopy()){

      const reco::Candidate* temp=&(*p);
      const reco::Candidate* mother=&(*p);

      bool hassamemother=true;
      while (hassamemother) {
	hassamemother=false;
	for (uint m=0; m<mother->numberOfMothers();m++){
	  if (mother->pdgId()==mother->mother(m)->pdgId()){
	    mother = mother->mother(m);
	    hassamemother=true;
	    break;
	  }
	}
      } 
      bool hasnohiggsmother = true;
      for (uint m=0; m<mother->numberOfMothers();m++){
	if (abs(mother->mother(m)->pdgId())==25)
	  hasnohiggsmother = false;
      }
      if (hasnohiggsmother)
	wboson=*(reco::GenParticle*)temp;
    }
  }

  // (b) quark from top
  for(auto p=prunedGenParticles.begin(); p!=prunedGenParticles.end(); p++){
    
    if (abs(p->pdgId())<6 && p->statusFlags().isLastCopy() && p->statusFlags().fromHardProcess()){

      const reco::Candidate* temp=&(*p);
      const reco::Candidate* mother=&(*p);
     
      bool hassamemother=true;
      while (hassamemother) {
	hassamemother=false;
	for (uint m=0; m<mother->numberOfMothers();m++){
	  if (mother->pdgId()==mother->mother(m)->pdgId()){
	    mother = mother->mother(m);
	    hassamemother=true;
	    break;
	  }
	}
      } 
      for (uint m=0; m<mother->numberOfMothers();m++){
	if (abs(mother->mother(m)->pdgId())==6)
	  top_decay_quark=*(reco::GenParticle*)temp;
      }
    }
  }

  // Higgs decay objects
  for(auto p=prunedGenParticles.begin(); p!=prunedGenParticles.end(); p++){
    
    if (p->statusFlags().isLastCopy() && p->statusFlags().fromHardProcess()){

      const reco::Candidate* temp=&(*p);
      const reco::Candidate* mother=&(*p);

      bool hassamemother=true;
      while (hassamemother) {
	hassamemother=false;
	for (uint m=0; m<mother->numberOfMothers();m++){
          if (mother->pdgId()==mother->mother(m)->pdgId()){
	    mother = mother->mother(m);
            hassamemother=true;
            break;
	  }
	}
      }
      for (uint m=0; m<mother->numberOfMothers();m++){
        if (abs(mother->mother(m)->pdgId())==25) {
	    if ( higgs_decay_products.size()==2 && higgs_decay_products[0].pdgId()==21)
	      break;
	    higgs_decay_products.push_back(*(reco::GenParticle*)temp);
	}
      }
    }
  }

  // W decay objects
  for(auto p=prunedGenParticles.begin(); p!=prunedGenParticles.end(); p++){
    
    if (p->statusFlags().isLastCopy() && p->statusFlags().fromHardProcess()){

      const reco::Candidate* temp=&(*p);
      const reco::Candidate* mother=&(*p);

      bool hassamemother=true;
      while (hassamemother) {
	hassamemother=false;
	for (uint m=0; m<mother->numberOfMothers();m++){
          if (mother->pdgId()==mother->mother(m)->pdgId()){
	    mother = mother->mother(m);
            hassamemother=true;
            break;
	  }
	}
      }
      for (uint m=0; m<mother->numberOfMothers();m++){

        if (abs(mother->mother(m)->pdgId())==24 ){
	  
	  
	  const reco::Candidate* mother_W = &(*(reco::GenParticle*)mother->mother(m));
	  
	  bool hassamemother=true;
	  while (hassamemother) {
	    hassamemother=false;
	    for (uint m=0; m<mother_W->numberOfMothers();m++){
	      if (mother_W->pdgId()==mother_W->mother(m)->pdgId()){
		mother_W = mother_W->mother(m);
		hassamemother=true;
		break;
	      }
	    }
	  } 
      
	  bool isnotfromhiggstoww = true;
	  for (uint m=0; m<mother_W->numberOfMothers();m++){
	    if (abs(mother_W->mother(m)->pdgId())==25)
	      isnotfromhiggstoww = false;
	  }
	  
	  if (isnotfromhiggstoww)
	    w_decay_products.push_back(*(reco::GenParticle*)temp);
	}
      }
    }
  }

  // light quark
  for(auto p=prunedGenParticles.begin(); p!=prunedGenParticles.end(); p++){
    
    if (std::abs(p->pdgId())<5 && p->statusFlags().isLastCopy() && p->statusFlags().fromHardProcess()){

      if (p->pt()!=top_decay_quark.pt())
	lightquark=*p;
      
    }
  }

  // second b quark
  for(auto p=prunedGenParticles.begin(); p!=prunedGenParticles.end(); p++){

    
    //    if (p->pdgId()*(-1)==top_decay_quark.pdgId() && p->statusFlags().isLastCopy() && p->statusFlags().fromHardProcess()){
    if (p->pdgId()*(-1)==top_decay_quark.pdgId() ){//&& p->statusFlags().isLastCopy() && p->statusFlags().fromHardProcess()){

      //      std::cout << "Status of second b: " << p->status() << std::endl;
      if (p->pt()!=higgs_decay_products[0].pt() && p->pt()!=higgs_decay_products[1].pt())
	secondb=*p;
      
    }
  }

  /*  std::cout << "Top decay quark id: " << top_decay_quark.pdgId() << std::endl;
  std::cout << "Second b pt: " << secondb.pt() << std::endl;
  std::cout << "Higgs energy: " << higgs.energy() << std::endl;
  std::cout << "W energy: " << wboson.energy() << std::endl;
  std::cout << "W id: " << wboson.pdgId() << std::endl; */

  if(w_decay_products.size()!=2) std::cerr << "GentHqEvent: error 2"<<std::endl;
  if(top.energy()<1||wboson.energy()<1) std::cerr << "GentHqEvent: error 3"<<std::endl;

  isFilled=true;

}

void GentHqEvent::Print() const{
  std::cout << "top" << std::endl;
  PrintParticle(GetTop());
  std::cout << "============================" << std::endl;
}

void GentHqEvent::PrintParticle(reco::GenParticle p) const{
  std::cout << "pdgId: " << p.pdgId() << ", pt: " << p.pt() << ", eta: " << p.eta() << ", phi: " << p.phi() << ", status:" << p.status() << std::endl;
}
void GentHqEvent::PrintParticles(std::vector<reco::GenParticle> ps) const{
  for(auto p=ps.begin();p!=ps.end();p++){
    PrintParticle(*p);
  }
}
reco::GenParticle GentHqEvent::GetHiggs() const{
  if(!isFilled) std::cerr << "Trying to access GentHqEvent but it is not filled" << std::endl;
  return higgs;
}
reco::GenParticle GentHqEvent::GetTop() const{
  if(!isFilled) std::cerr << "Trying to access GentHqEvent but it is not filled" << std::endl;
  return top;
}
reco::GenParticle GentHqEvent::GetW() const{
  if(!isFilled) std::cerr << "Trying to access GentHqEvent but it is not filled" << std::endl;
  return wboson;
}
reco::GenParticle GentHqEvent::GetLightQuark() const{
  if(!isFilled) std::cerr << "Trying to access GentHqEvent but it is not filled" << std::endl;
  return lightquark;
}
reco::GenParticle GentHqEvent::GetSecondb() const{
  if(!isFilled) std::cerr << "Trying to access GentHqEvent but it is not filled" << std::endl;
  return secondb;
}
std::vector<reco::GenParticle> GentHqEvent::GetHiggsDecayProducts() const{
  if(!isFilled) std::cerr << "Trying to access GentHqEvent but it is not filled" << std::endl;
  return higgs_decay_products;
}
std::vector<reco::GenParticle> GentHqEvent::GetWDecayProducts() const{
  if(!isFilled) std::cerr << "Trying to access GentHqEvent but it is not filled" << std::endl;
  return w_decay_products;
}
reco::GenParticle GentHqEvent::GetTopDecayQuark() const{
  if(!isFilled) std::cerr << "Trying to access GentHqEvent but it is not filled" << std::endl;
  return top_decay_quark;
}
