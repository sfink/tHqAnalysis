#ifndef THQANALYSIS_THQANALYZER_GENTHQEVENT_HPP
#define THQANALYSIS_THQANALYZER_GENTHQEVENT_HPP
#include <vector>
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"

class GentHqEvent{
  
public:
  GentHqEvent();
  ~GentHqEvent();
  void Fill(const std::vector<reco::GenParticle>& prunedGenParticles);
  

  std::vector<reco::GenParticle> GetHiggsDecayProducts() const;
  bool IsFilled() const;
  reco::GenParticle GetTop() const;
  reco::GenParticle GetHiggs() const;
  reco::GenParticle GetTopDecayQuark() const;
  reco::GenParticle GetW() const;
  
  std::vector<reco::GenParticle> GetWDecayProducts() const;

  reco::GenParticle GetSecondb() const;
  reco::GenParticle GetLightQuark() const;

  void Print() const;
private:
  math::XYZTLorentzVector GetLV(const reco::GenParticle& p) const;
  std::vector<math::XYZTLorentzVector> GetLVs(const std::vector<reco::GenParticle>& ps) const;
  void PrintParticle(reco::GenParticle) const;
  void PrintParticles(std::vector<reco::GenParticle>) const;

  reco::GenParticle higgs;
  reco::GenParticle top;
  reco::GenParticle wboson;
  reco::GenParticle lightquark;
  reco::GenParticle secondb;
  reco::GenParticle top_decay_quark;
  reco::GenParticle W_decay_product;
  reco::GenParticle H_decay_product;
  std::vector<reco::GenParticle> w_decay_products;
  std::vector<reco::GenParticle> higgs_decay_products;

  bool isFilled;

};

#endif
