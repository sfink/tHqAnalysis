#include "tHqAnalysis/tHqAnalyzer/interface/GentHqEvent.hpp"

GentHqEvent::GentHqEvent (){
  isFilled=false;
}
GentHqEvent::~GentHqEvent(){}

bool GentHqEvent::IsFilled() const{
  return isFilled;
}

void GentHqEvent::Fill(const std::vector<reco::GenParticle>& prunedGenParticles){

  for(auto p=prunedGenParticles.begin(); p!=prunedGenParticles.end(); p++){
    

    // Search for top quark

    if (abs(p->pdgId())==6 && p->isLastCopy()){
      top=*p;
      for(uint i=0;i<p->numberOfDaughters();i++){
	if(abs(p->daughter(i)->pdgId())<6){
	  top_decay_quark=*(reco::GenParticle*)p->daughter(i);

	  bool goon = true;

          while (goon == true){
            goon = false;
            for (unsigned i=0; i<top_decay_quark.numberOfDaughters();i++){
              if (top_decay_quark.pdgId()==top_decay_quark.daughter(i)->pdgId()){
                top_decay_quark = *(reco::GenParticle*)top_decay_quark.daughter(i);
                goon = true;
                break;
              }
            }
          }
          //top_decay_quark=*(reco::GenParticle*)top_decay_quark;
	}
      }
    }

    // Search for W boson

    if (abs(p->pdgId())==24 && p->isLastCopy()){
      wboson=*p;
      for(uint i=0;i<p->numberOfDaughters();i++){
	if(p->pdgId()==24 && abs(p->daughter(i)->pdgId())<=16){

	  W_decay_product = *(reco::GenParticle*)p->daughter(i);
	  
	  bool goon = true;

	  while (goon == true){
	    goon = false;
	    for (unsigned i=0; i<W_decay_product.numberOfDaughters();i++){
	      if (W_decay_product.pdgId()==W_decay_product.daughter(i)->pdgId()){
		W_decay_product = *(reco::GenParticle*)W_decay_product.daughter(i);
		goon = true;
		break;
	      }
	    }
	  }
	  w_decay_products.push_back((reco::GenParticle)W_decay_product);   
	}
      }
    }

    // Search for Higgs boson

    if (abs(p->pdgId())==25 && p->isLastCopy()){
      higgs=*p;
      for(uint i=0;i<p->numberOfDaughters();i++){
	if(p->pdgId()==25 && abs(p->daughter(i)->pdgId())!=25){
	  
	  H_decay_product = *(reco::GenParticle*)p->daughter(i);

          bool goon = true;

          while (goon == true){
            goon = false;
            for (unsigned i=0; i<H_decay_product.numberOfDaughters();i++){
              if (H_decay_product.pdgId()==H_decay_product.daughter(i)->pdgId()){
                H_decay_product = *(reco::GenParticle*)H_decay_product.daughter(i);
                goon = true;
                break;
              }
            }
          }
          higgs_decay_products.push_back((reco::GenParticle)H_decay_product);
	}
      }
    }

    // Search for light quark

    if (abs(p->pdgId())<5 && p->status()==23){
      
      bool hastopmother = false;
      for(uint i=0;i<p->numberOfMothers();i++){
        if(p->mother(i)->pdgId()==6){
          hastopmother = true;
        }
      }
      
      if (!hastopmother)
	lightquark = *p;
	
      bool goon = true;

      while (goon == true){
	goon = false;
	for (unsigned i=0; i<lightquark.numberOfDaughters();i++){
	  if (lightquark.pdgId()==lightquark.daughter(i)->pdgId()){
	    lightquark = *(reco::GenParticle*)lightquark.daughter(i);
	    goon = true;
	    break;
	  }
	}
      }
    }

    // Search for additional b quark

    if ( (p->status()== 23 ||  p->status()== 43) && abs(p->pdgId())==5 && p->pdgId()*top.pdgId()<0){
      secondb=*p;
	
      bool goon = true;

      while (goon == true){
	goon = false;
	for (unsigned i=0; i<secondb.numberOfDaughters();i++){
	  if (secondb.pdgId()==secondb.daughter(i)->pdgId()){
	    secondb = *(reco::GenParticle*)secondb.daughter(i);
	    goon = true;
	    break;
	  }
	}
      }
    }
    //    if (abs(p->pdgId())<5 && p->fromHardProcess() && isLastCopy())
  }  if(w_decay_products.size()!=2) std::cerr << "GentHqEvent: error 2"<<std::endl;
  if(top.energy()<1||wboson.energy()<1) std::cerr << "GentHqEvent: error 3"<<std::endl;

  /*
  int nquarks_from_w=0;
  for(auto p=w_decay_products.begin(); p!=w_decay_products.end();p++){
    if(abs(p->pdgId())<6) n_from_wplus++;
  }
  topIsHadronic=nquarks_from_w==2;
  isFilled=true;
  */
}

void GentHqEvent::Print() const{
  std::cout << "top" << std::endl;
  PrintParticle(GetTop());
  std::cout << "============================" << std::endl;
}

void GentHqEvent::PrintParticle(reco::GenParticle p) const{
  std::cout << "pdgId: " << p.pdgId() << ", pt: " << p.pt() << ", eta: " << p.eta() << ", phi: " << p.phi() << std::endl;
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
