#ifndef THQANALYSIS_THQANALYZER_BOOSTEDUTILS_HPP
#define THQANALYSIS_THQANALYZER_BOOSTEDUTILS_HPP

#include <vector>

#include "TMath.h"
#include "TVector.h"
#include "TLorentzVector.h"

#include "Math/VectorUtil.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

class tHqUtils{
  
  public:
    
    static std::string GetAnalyzerPath();
    
    static TLorentzVector GetTLorentzVector(const math::XYZTLorentzVector& vec);
    static math::XYZTLorentzVector GetXYZTLorentzVector(const TLorentzVector& vec);
    
    static bool FirstIsLarger(float val1,float val2);
    static bool FirstIsHarder(math::XYZTLorentzVector vec1,math::XYZTLorentzVector vec2);
    static bool FirstHasHigherCSV(pat::Jet jet1,pat::Jet jet2);
    static bool FirstHasHigherCSVold(pat::Jet jet1,pat::Jet jet2);
    static bool FirstJetIsHarder(pat::Jet jet1, pat::Jet jet2);
    template<typename boostedJetType>
    static bool FirstFatJetIsHarder(boostedJetType jet1, boostedJetType jet2){
      return tHqUtils::FirstJetIsHarder(jet1.fatjet,jet2.fatjet);
    }
    
    static float DeltaEta(const math::XYZTLorentzVector& vec1,const math::XYZTLorentzVector& vec2);
    static float DeltaEta(const pat::Jet& jet1,const pat::Jet& jet2);
    static float DeltaPhi(const math::XYZTLorentzVector& vec1,const math::XYZTLorentzVector& vec2);
    static float DeltaPhi(const pat::Jet& jet1,const pat::Jet& jet2);
    static float DeltaR(const math::XYZTLorentzVector& vec1,const math::XYZTLorentzVector& vec2);
    static float DeltaR(const pat::Jet& jet1,const pat::Jet& jet2);
    static float DeltaKt(const math::XYZTLorentzVector& vec1,const math::XYZTLorentzVector& vec2);
    static float DeltaKt(const pat::Jet& jet1,const pat::Jet& jet2);
    
    static float CosThetaStar(const math::XYZTLorentzVector& vec1, const math::XYZTLorentzVector& vec2);
    static float CosThetaCM(const math::XYZTLorentzVector& vec,const math::XYZTLorentzVector& boostVec);
    static float GetMuondBetaIso(const pat::Muon& muon);
    static float GetElectronIso(const pat::Electron& elec);

    
    static std::vector<math::XYZTLorentzVector> GetGenParticleVecs(const std::vector<reco::GenParticle>& genParticles);
    static bool MCContainsTTbar(const std::vector<reco::GenParticle>& genParticles);
    static bool MCContainsST(const std::vector<reco::GenParticle>& genParticles);
    static bool MCContainsHiggs(const std::vector<reco::GenParticle>& genParticles);
    static void GetttHMCParticles(const std::vector<reco::GenParticle>& genParticles, std::vector<reco::GenParticle>& tophad, std::vector<reco::GenParticle>& whad, std::vector<reco::GenParticle>& bhad, std::vector<reco::GenParticle>& q1, std::vector<reco::GenParticle>& q2, std::vector<reco::GenParticle>& toplep, std::vector<reco::GenParticle>& wlep, std::vector<reco::GenParticle>& blep, std::vector<reco::GenParticle>& lep, std::vector<reco::GenParticle>& nu, reco::GenParticle& higgs, reco::GenParticle& b1, reco::GenParticle& b2);
    static void GetttHMCVecs(const std::vector<reco::GenParticle>& genParticles, std::vector<math::XYZTLorentzVector>& tophadvecs, std::vector<math::XYZTLorentzVector>& whadvecs, std::vector<math::XYZTLorentzVector>& bhadvecs, std::vector<math::XYZTLorentzVector>& q1vecs, std::vector<math::XYZTLorentzVector>& q2vecs, std::vector<math::XYZTLorentzVector>& toplepvecs, std::vector<math::XYZTLorentzVector>& wlepvecs, std::vector<math::XYZTLorentzVector>& blepvecs, std::vector<math::XYZTLorentzVector>& lepvecs, std::vector<math::XYZTLorentzVector>& nuvecs, math::XYZTLorentzVector& higgsvec, math::XYZTLorentzVector& b1vec, math::XYZTLorentzVector& b2vec);
    
    static bool IsAnyTriggerBitFired(const std::vector<std::string>& targetTriggers, const edm::TriggerResults& triggerResults);
    
    static std::vector<math::XYZTLorentzVector> GetLepVecs(const std::vector<pat::Electron>& selectedElectrons, const std::vector<pat::Muon> selectedMuons);
    static math::XYZTLorentzVector GetPrimLepVec(const std::vector<pat::Electron>& selectedElectrons, const std::vector<pat::Muon> selectedMuons);
    
    static void GetNuVec(const math::XYZTLorentzVector& lepvec, const TVector2& metvec, math::XYZTLorentzVector& nuvec, math::XYZTLorentzVector& lepwvec);
    
    static std::vector<math::XYZTLorentzVector> GetJetVecs(const std::vector<pat::Jet>& jets);
    
    static bool PassesCSV(const pat::Jet& jet, const char workingPoint);
    
    static float GetClosestJetIDs(int& idJet1, int& idJet2, const std::vector<pat::Jet>& jets);
    static float GetClosestLepJetID(int& idJet, const math::XYZTLorentzVector& lepVec, const std::vector<pat::Jet>& jets);
    
    static float GetJetAverageJetEtaMax(const std::vector<pat::Jet>& jets1, const std::vector<pat::Jet>& jets2);
    
    static void GetFoxWolframMoments(std::vector<math::XYZTLorentzVector> jetVecs, float &h0, float &h1, float &h2, float &h3, float &h4);
    static void GetAplanaritySphericity(math::XYZTLorentzVector leptonVec, math::XYZTLorentzVector metVec, std::vector<math::XYZTLorentzVector> jetVecs, float &aplanarity, float &sphericity);
    
  template <typename T, typename S> static std::vector<T> RemoveOverlaps( const std::vector<S>&, const std::vector<T>& );
  template <typename T, typename S> static T RemoveOverlap( const std::vector<S>&, const T& );
  
  private:
  
};


template <typename PATObj1, typename PATObj2> 
PATObj1 tHqUtils::RemoveOverlap( const std::vector<PATObj2>& other, const PATObj1& unclean ){

  unsigned int nSources1 = unclean.numberOfSourceCandidatePtrs();
  bool hasOverlaps = false;

  std::vector<reco::CandidatePtr> overlaps;

  for( typename std::vector<PATObj2>::const_iterator iobj2 = other.begin(); iobj2!=other.end(); ++iobj2 ){

    unsigned int nSources2 = iobj2->numberOfSourceCandidatePtrs();

    for( unsigned int i1=0; i1<nSources1; i1++ ){
      reco::CandidatePtr source1 = unclean.sourceCandidatePtr(i1);

      if( !(source1.isNonnull() && source1.isAvailable()) ) continue;

      for( unsigned int i2=0; i2<nSources2; i2++ ){
	reco::CandidatePtr source2 = iobj2->sourceCandidatePtr(i2);

	if( !(source2.isNonnull() && source2.isAvailable()) ) continue;

	if( source1==source2 ){
	  hasOverlaps = true;
	  overlaps.push_back(source2);
	}
      }

    }
  }// end loop over iobj22


  PATObj1 cleaned = unclean;
  if( hasOverlaps ){
    math::XYZTLorentzVector original = cleaned.p4();

    for( int iOverlap=0; iOverlap<int(overlaps.size()); iOverlap++ ){

      const reco::Candidate & cOverlap = *(overlaps[iOverlap]);
      math::XYZTLorentzVector overlaper = cOverlap.p4();

      original -= overlaper;
    }

    cleaned.setP4( original );
  }

  return cleaned;
}


template <typename PATObj1, typename PATObj2> 
std::vector<PATObj1> tHqUtils::RemoveOverlaps( const std::vector<PATObj2>& other, const std::vector<PATObj1>& unclean ){

  std::vector<PATObj1> cleaned;
  
  for( typename std::vector<PATObj1>::const_iterator iobj1 = unclean.begin(); iobj1!=unclean.end(); ++iobj1 ){

    PATObj1 myobj = (*iobj1);
    PATObj1 clean = RemoveOverlap(other, myobj);

    cleaned.push_back(clean);
  }

  return cleaned;
}
    


#endif
