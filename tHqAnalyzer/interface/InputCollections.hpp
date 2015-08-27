#ifndef THQ_ANALYSIS_THQANALYZER_INPUTCOLLECTIONS_HPP
#define THQ_ANALYSIS_THQANALYZER_INPUTCOLLECTIONS_HPP

#include <vector>
#include <map>


#include "tHqAnalysis/tHqObjects/interface/Event.h"
#include "tHqAnalysis/tHqAnalyzer/interface/TriggerInfo.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/GenTopEvent.hpp"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"


enum SampleType{data,tth,tt,nonttbkg};
namespace HiggsDecay{enum HiggsDecay{NA,bb,nonbb};};

struct InputCollections{
  InputCollections( const boosted::Event&                         event_,
		    const TriggerInfo&                            triggerInfo_,
		    //                    const pat::TriggerObjectStandAloneCollection& selectedTrigger_,
		    //                    const edm::TriggerResults&                    triggerResults_,
		    //                    const HLTConfigProvider&                      hlt_config_,
                    const reco::VertexCollection&                 selectedPVs_,
                    const std::vector<pat::Muon>&                 selectedMuons_,
                    const std::vector<pat::Muon>&                 selectedMuonsLoose_,
                    const std::vector<pat::Electron>&             selectedElectrons_,
                    const std::vector<pat::Electron>&             selectedElectronsLoose_,
                    const std::vector<pat::Jet>&                  selectedJets_,
                    const std::vector<pat::Jet>&                  selectedJetsLoose_,
                    const std::vector<pat::Jet>&                  selectedPuppiJets_,
                    const std::vector<pat::MET>&                  pfMets_,
                    const std::vector<reco::GenParticle>&         genParticles_,
                    const std::vector<reco::GenJet>&              selectedGenJets_,
		    const GenTopEvent&                            genTopEvt_,
                    const SampleType                              sampleType_,
                    const std::map<std::string,float>&            weights_,
		    const edm::EventSetup&	                  setup_,
        	    const edm::Event&	                          edmevent_

                  ):
                    event(event_),
		    triggerInfo(triggerInfo_),
		    //		    selectedTrigger(selectedTrigger_),
		    //                    triggerResults(triggerResults_),
		    //                    hlt_config(hlt_config_),
                    selectedPVs(selectedPVs_),
                    selectedMuons(selectedMuons_),
                    selectedMuonsLoose(selectedMuonsLoose_),
                    selectedElectrons(selectedElectrons_),
                    selectedElectronsLoose(selectedElectronsLoose_),
                    selectedJets(selectedJets_),
                    selectedJetsLoose(selectedJetsLoose_),
                    selectedPuppiJets(selectedPuppiJets_),
                    pfMets(pfMets_),
                    genParticles(genParticles_),
                    selectedGenJets(selectedGenJets_),
		    genTopEvt(genTopEvt_),
                    sampleType(sampleType_),
                    weights(weights_),
		    setup(setup_),
                    edmevent(edmevent_){}
  
  const boosted::Event&                         event;
  const TriggerInfo&                            triggerInfo;
  //  const pat::TriggerObjectStandAloneCollection& selectedTrigger;
  //  const edm::TriggerResults&                    triggerResults;
  //  const HLTConfigProvider&                       hlt_config;
  const reco::VertexCollection&                 selectedPVs;
  const std::vector<pat::Muon>&                 selectedMuons;
  const std::vector<pat::Muon>&                 selectedMuonsLoose;
  const std::vector<pat::Electron>&             selectedElectrons;
  const std::vector<pat::Electron>&             selectedElectronsLoose;
  const std::vector<pat::Jet>&                  selectedJets;
  const std::vector<pat::Jet>&                  selectedJetsLoose;
  const std::vector<pat::Jet>&                  selectedPuppiJets;
  const std::vector<pat::MET>&                  pfMets;
  const std::vector<reco::GenParticle>&         genParticles;
  const std::vector<reco::GenJet>&              selectedGenJets;
  const GenTopEvent&                            genTopEvt;
  const SampleType                              sampleType;
  const std::map<std::string,float>&            weights;
  const edm::EventSetup& 			setup;
  const edm::Event& 			        edmevent;

  void Dump();
};

#endif
