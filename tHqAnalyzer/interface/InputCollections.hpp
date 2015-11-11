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
#include "tHqAnalysis/tHqAnalyzer/interface/GenTopEvent.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/GentHqEvent.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/EventInfo.hpp"
//#include "tHqAnalysis/tHqAnalyzer/interface/BoostedUtils.hpp"
#include "MiniAOD/MiniAODHelper/interface/MiniAODHelper.h"


enum SampleType{data,tth,ttl,ttbb,ttb,tt2b,ttcc,ttc,nonttbkg,thq};
namespace HiggsDecay{enum HiggsDecay{NA,bb,nonbb};};

struct InputCollections{
  InputCollections( const EventInfo&                              eventInfo_,
		    const TriggerInfo&                            triggerInfo_,
                    const std::vector<reco::Vertex>&              selectedPVs_,
		    const std::vector<pat::Muon>&                 rawMuons_,
                    const std::vector<pat::Muon>&                 selectedMuons_,
                    const std::vector<pat::Muon>&                 selectedMuonsLoose_,
                    const std::vector<pat::Electron>&             rawElectrons_,
                    const std::vector<pat::Electron>&             selectedElectrons_,
                    const std::vector<pat::Electron>&             selectedElectronsLoose_,
                    const std::vector<pat::Jet>&                  rawJets_,
                    const std::vector<pat::Jet>&                  selectedJets_,
                    const std::vector<pat::Jet>&                  rawPuppiJets_,
                    const std::vector<pat::Jet>&                  selectedJetsLoose_,
                    const std::vector<pat::Jet>&                  selectedPuppiJets_,
                    const pat::MET&                               pfMET_,
		    const GenTopEvent&                            genTopEvt_,
                    const GentHqEvent&                            gentHqEvt_,
                    const std::vector<reco::GenJet>&              selectedGenJets_,
                    const SampleType                              sampleType_,
                    const std::map<std::string,float>&            weights_,
		    const float                                   Weight_orig_,
		    const std::vector<float>&                     syst_weights_,
		    const std::vector<string>&                    syst_weights_id_
                  ):
                    eventInfo(eventInfo_),
		    triggerInfo(triggerInfo_),
                    selectedPVs(selectedPVs_),
		    rawMuons(rawMuons_),
                    selectedMuons(selectedMuons_),
                    selectedMuonsLoose(selectedMuonsLoose_),
		    rawElectrons(rawElectrons_),
                    selectedElectrons(selectedElectrons_),
                    selectedElectronsLoose(selectedElectronsLoose_),
		    rawJets(rawJets_),
                    selectedJets(selectedJets_),
                    rawPuppiJets(rawPuppiJets_),
                    selectedJetsLoose(selectedJetsLoose_),
                    selectedPuppiJets(selectedPuppiJets_),
                    pfMET(pfMET_),
		    genTopEvt(genTopEvt_),
		    gentHqEvt(gentHqEvt_),
                    selectedGenJets(selectedGenJets_),
                    sampleType(sampleType_),
                    weights(weights_),
		    Weight_orig(Weight_orig_),
		    syst_weights(syst_weights_),
		    syst_weights_id(syst_weights_id_)
                    {}
  
  const EventInfo&                              eventInfo;
  const TriggerInfo&                            triggerInfo;
  const std::vector<reco::Vertex>&              selectedPVs;
  const std::vector<pat::Muon>&                 rawMuons;
  const std::vector<pat::Muon>&                 selectedMuons;
  const std::vector<pat::Muon>&                 selectedMuonsLoose;
  const std::vector<pat::Electron>&             rawElectrons;
  const std::vector<pat::Electron>&             selectedElectrons;
  const std::vector<pat::Electron>&             selectedElectronsLoose;
  const std::vector<pat::Jet>&                  rawJets;
  const std::vector<pat::Jet>&                  selectedJets;
  const std::vector<pat::Jet>&                  rawPuppiJets;
  const std::vector<pat::Jet>&                  selectedJetsLoose;
  const std::vector<pat::Jet>&                  selectedPuppiJets;
  const pat::MET&                               pfMET;
  const GenTopEvent&                            genTopEvt;
  const GentHqEvent&                            gentHqEvt;
  const std::vector<reco::GenJet>&              selectedGenJets;
  const SampleType                              sampleType;
  const std::map<std::string,float>&            weights;
  const float                                   Weight_orig;
  const vector<float>&                          syst_weights;
  const vector<string>&                         syst_weights_id;
};

#endif
