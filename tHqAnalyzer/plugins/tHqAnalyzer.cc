// -*- C++ -*-
//
// Package:    tHqAnalysis/tHqAnalyzer
// Class:      BoostedAnalyzer
// 
/**\class BoostedAnalyzer BoostedAnalyzer.cc tHqAnalysis/tHqAnalyzer/plugins/BoostedAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Shawn Williamson, Hannes Mildner
//         
//
//


// system include files
#include <memory>
#include <vector>
#include <map>

#include "TStopwatch.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "MiniAOD/MiniAODHelper/interface/MiniAODHelper.h"

#include "tHqAnalysis/tHqAnalyzer/interface/tHqUtils.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/InputCollections.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/Cutflow.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/TreeWriter.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/TriggerInfo.hpp"

#include "tHqAnalysis/tHqAnalyzer/interface/Selection.hpp"
//#include "tHqAnalysis/tHqAnalyzer/interface/LeptonSelection.hpp"
//#include "tHqAnalysis/tHqAnalyzer/interface/JetTagSelection.hpp"
//#include "tHqAnalysis/tHqAnalyzer/interface/SynchSelection.hpp"

//#include "tHqAnalysis/tHqAnalyzer/interface/WeightProcessor.hpp"
//#include "tHqAnalysis/tHqAnalyzer/interface/MCMatchVarProcessor.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/MVAVarProcessor.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/BaseVarProcessor.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/RecoVarProcessor.hpp"

//#include "tHqAnalysis/tHqAnalyzer/interface/BDTVarProcessor.hpp"
//#include "tHqAnalysis/tHqAnalyzer/interface/BoostedJetVarProcessor.hpp"
//#include "tHqAnalysis/tHqAnalyzer/interface/ttHVarProcessor.hpp"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

//
// class declaration
//

class tHqAnalyzer : public edm::EDAnalyzer {
   public:
      explicit tHqAnalyzer(const edm::ParameterSet&);
      ~tHqAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override;
  
      boosted::Event FillEvent(const edm::Event& iEvent, const edm::Handle<GenEventInfoProduct>& genEvtInfo, const edm::Handle<reco::BeamSpot>& beamSpot, const edm::Handle<HcalNoiseSummary>& hcalNoiseSummary, const edm::Handle< std::vector<PileupSummaryInfo> >& puSummaryInfo);
      map<string,float> GetWeights(const boosted::Event& event, const reco::VertexCollection& selectedPVs, const std::vector<pat::Jet>& selectedJets, const std::vector<pat::Electron>& selectedElectrons, const std::vector<pat::Muon>& selectedMuons, const std::vector<reco::GenParticle>& genParticles);
      std::vector<pat::Electron> ElectronSelection( std::vector<pat::Electron> selectedElectrons, const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> pvposition);
      std::vector<pat::Muon> MuonSelection( std::vector<pat::Muon> selectedMuons, const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> pvposition);      
  std::vector<pat::Jet> JetSelection( std::vector<pat::Jet> selectedJets, std::vector<pat::Electron> selectedElectrons, std::vector<pat::Muon> selectedMuons, InputCollections input);
      
      
      
      // ----------member data ---------------------------
      
      /** the beanhelper is used for selections and reweighting */
      MiniAODHelper helper;
      
      /** writes flat trees for MVA analysis */
      TreeWriter treewriter;
      
      /** stores cutflow*/
      Cutflow cutflow;
  
      /** selections that are applied */
      vector<Selection*> selections;
      
      /** sample ID */
      int sampleID;
      
      /** stops the time needed */
      TStopwatch watch;
      
      /** events to be analyzed */
      int maxEvents;
      
      /** data luminosity */
      float luminosity;
      
      /** sample cross section*/
      float xs;
      
      /** total number of events in input file(s) */
      int totalMCevents;
     
      /** Event counter */
      int eventcount;
      
       /** is analyzed sample data? */
      bool isData;

       /** triggers that are checked */
       vector<std::string> relevantTriggers;

      /** disable some object selections for synch exe? */
      bool disableObjectSelections;

      /** calculate and store systematic weights? */
      bool doSystematics;
      
      /** jet systematic that is applied (the outher systematics are done at a different place with reweighting)*/
      sysType::sysType jsystype;
      
      
      /** pu summary data access token **/
      edm::EDGetTokenT< std::vector<PileupSummaryInfo> > EDMPUInfoToken;
      
      /** hcal noise data access token **/
      edm::EDGetTokenT< HcalNoiseSummary > EDMHcalNoiseToken;
      
      /** selected trigger data access token **/
      edm::EDGetTokenT< pat::TriggerObjectStandAloneCollection > EDMSelectedTriggerToken;
      
      /** trigger results data access token **/
      edm::EDGetTokenT< edm::TriggerResults > EDMTriggerResultToken;
      HLTConfigProvider hlt_config;

      /** beam spot data access token **/
      edm::EDGetTokenT< reco::BeamSpot > EDMBeamSpotToken;
      
      /** vertex data access token **/
      edm::EDGetTokenT< reco::VertexCollection > EDMVertexToken;
      
      /** muons data access token **/
      edm::EDGetTokenT< std::vector<pat::Muon> > EDMMuonsToken;
      
      /** electrons data access token **/
      edm::EDGetTokenT< std::vector<pat::Electron> > EDMElectronsToken;
      
      /** jets data access token **/
      edm::EDGetTokenT< std::vector<pat::Jet> > EDMJetsToken;

      /** puppi jets data access token **/
      edm::EDGetTokenT< std::vector<pat::Jet> > EDMPuppiJetsToken;


      /** mets data access token **/
      edm::EDGetTokenT< std::vector<pat::MET> > EDMMETsToken;
      
      /** hep top jets data access token **/
  //      edm::EDGetTokenT< boosted::HEPTopJetCollection > EDMHEPTopJetsToken;
      
      /** subjet filterjets data access token **/
  //      edm::EDGetTokenT< boosted::SubFilterJetCollection > EDMSubFilterJetsToken;
      
      /** gen info data access token **/
      edm::EDGetTokenT< GenEventInfoProduct > EDMGenInfoToken;
      
      /** gen particles data access token **/
      edm::EDGetTokenT< std::vector<reco::GenParticle> > EDMGenParticlesToken;
      
      /** gen jets data access token **/
      edm::EDGetTokenT< std::vector<reco::GenJet> > EDMGenJetsToken;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
tHqAnalyzer::tHqAnalyzer(const edm::ParameterSet& iConfig)
{

  std::string era = iConfig.getParameter<std::string>("era");
  string analysisType = iConfig.getParameter<std::string>("analysisType");
  analysisType::analysisType iAnalysisType = analysisType::LJ;
  if(analysisType == "LJ") iAnalysisType = analysisType::LJ;
  else if(analysisType == "DIL") iAnalysisType = analysisType::DIL;
  else if(analysisType == "TauLJ") iAnalysisType = analysisType::TauLJ;
  else if(analysisType == "TauDIL") iAnalysisType = analysisType::TauDIL;
  else cerr << "No matching analysis type found for: " << analysisType << endl;
  
  luminosity = iConfig.getParameter<double>("luminostiy");
  sampleID = iConfig.getParameter<int>("sampleID");
  xs = iConfig.getParameter<double>("xs");
  totalMCevents = iConfig.getParameter<int>("nMCEvents");
  isData = iConfig.getParameter<bool>("isData");
  
  //  useFatJets = iConfig.getParameter<bool>("useFatJets");
  disableObjectSelections = iConfig.getParameter<bool>("disableObjectSelections");

  string outfileName = iConfig.getParameter<std::string>("outfileName");
  
  std::cout << "Outfile Name: " << outfileName << std::endl;
  
  // REGISTER DATA ACCESS
  EDMPUInfoToken          = consumes< std::vector<PileupSummaryInfo> >(edm::InputTag("addPileupInfo","","HLT"));
  EDMHcalNoiseToken       = consumes< HcalNoiseSummary >(edm::InputTag("hcalnoise","","RECO"));
  EDMSelectedTriggerToken = consumes< pat::TriggerObjectStandAloneCollection > (edm::InputTag("selectedPatTrigger","","PAT"));
  EDMTriggerResultToken   = consumes< edm::TriggerResults > (edm::InputTag("TriggerResults","","HLT"));
  EDMBeamSpotToken        = consumes< reco::BeamSpot > (edm::InputTag("offlineBeamSpot","","RECO"));
  EDMVertexToken          = consumes< reco::VertexCollection > (edm::InputTag("offlineSlimmedPrimaryVertices","","PAT"));
  EDMMuonsToken           = consumes< std::vector<pat::Muon> >(edm::InputTag("slimmedMuons","","PAT"));
  EDMElectronsToken       = consumes< std::vector<pat::Electron> >(edm::InputTag("slimmedElectrons","","PAT"));
  EDMJetsToken            = consumes< std::vector<pat::Jet> >(edm::InputTag("slimmedJets","","PAT"));
  EDMPuppiJetsToken       = consumes< std::vector<pat::Jet> >(edm::InputTag("slimmedJetsPuppi","","PAT"));
  EDMMETsToken            = consumes< std::vector<pat::MET> >(edm::InputTag("slimmedMETs","","PAT"));
  //EDMHEPTopJetsToken      = consumes< boosted::HEPTopJetCollection >(edm::InputTag("HEPTopJetsPFMatcher","heptopjets","p"));
  // EDMSubFilterJetsToken   = consumes< boosted::SubFilterJetCollection >(edm::InputTag("CA12JetsCA3FilterjetsPFMatcher","subfilterjets","p"));
  EDMGenInfoToken         = consumes< GenEventInfoProduct >(edm::InputTag("generator","","SIM"));
  EDMGenParticlesToken    = consumes< std::vector<reco::GenParticle> >(edm::InputTag("prunedGenParticles","","PAT"));
  EDMGenJetsToken         = consumes< std::vector<reco::GenJet> >(edm::InputTag("slimmedGenJets","","PAT"));
  
  // INITIALIZE MINIAOD HELPER
  helper.SetUp(era, sampleID, iAnalysisType, isData);
  
  // INITIALIZE SELECTION & CUTFLOW
  cutflow.Init((outfileName+"_Cutflow.txt").c_str());
  cutflow.AddStep("all");
  
  std::vector<std::string> selectionNames = iConfig.getParameter< std::vector<std::string> >("selectionNames");
  /*  for(vector<string>::const_iterator itSel = selectionNames.begin();itSel != selectionNames.end();itSel++) {
    
    if(*itSel == "LeptonSelection") selections.push_back(new LeptonSelection());
    else if(*itSel == "JetTagSelection") selections.push_back(new JetTagSelection());
    else if(*itSel == "SynchSelection") selections.push_back(new SynchSelection());
    else cout << "No matching selection found for: " << *itSel << endl;
    
    selections.back()->Init(iConfig,cutflow);
    } */
  
  
  relevantTriggers = iConfig.getParameter< std::vector<std::string> >("relevantTriggers");

  // INITIALIZE TREEWRITER
  treewriter.Init(outfileName);
  std::vector<std::string> processorNames = iConfig.getParameter< std::vector<std::string> >("processorNames");
  for(vector<string>::const_iterator itPro = processorNames.begin();itPro != processorNames.end();++itPro) {
    
    /* if(*itPro == "WeightProcessor") treewriter.AddTreeProcessor(new WeightProcessor());
    else if(*itPro == "MCMatchVarProcessor") treewriter.AddTreeProcessor(new MCMatchVarProcessor());
    else if(*itPro == "MVAVarProcessor") treewriter.AddTreeProcessor(new MVAVarProcessor());
    else if(*itPro == "tHqJetVarProcessor") treewriter.AddTreeProcessor(new tHqJetVarProcessor());
    else if(*itPro == "tHqTopHiggsVarProcessor") treewriter.AddTreeProcessor(new ttHVarProcessor(tHqRecoType::tHqTopHiggs,"TopLikelihood","HiggsCSV","tHqTopHiggs_"));
    else if(*itPro == "tHqTopVarProcessor") treewriter.AddTreeProcessor(new ttHVarProcessor(tHqRecoType::tHqTop,"TopLikelihood","HiggsCSV","tHqTop_"));
    else if(*itPro == "tHqHiggsVarProcessor") treewriter.AddTreeProcessor(new ttHVarProcessor(tHqRecoType::tHqHiggs,"TopLikelihood","HiggsCSV","tHqHiggs_"));
    // the BDT processor rely on the variables filled py the other producers and should be added at the end
    else if(*itPro == "BDTVarProcessor") treewriter.AddTreeProcessor(new BDTVarProcessor());
    */
    if(*itPro == "BaseVarProcessor") treewriter.AddTreeProcessor(new BaseVarProcessor());
    if(*itPro == "RecoVarProcessor") treewriter.AddTreeProcessor(new RecoVarProcessor());
    if(*itPro == "MVAVarProcessor") treewriter.AddTreeProcessor(new MVAVarProcessor());
    else cout << "No matching processor found for: " << *itPro << endl;    
    } 
}


tHqAnalyzer::~tHqAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
tHqAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if(eventcount<10||eventcount%1000==0){
    cout << "Analyzing event " << eventcount << endl;
    watch.Print();
    watch.Continue();
  }
  
  eventcount++;
  
  /**** GET PILEUPSUMMARYINFO ****/
  edm::Handle< std::vector<PileupSummaryInfo> >  h_puinfosummary;
  iEvent.getByToken( EDMPUInfoToken, h_puinfosummary);
  
  /**** GET HCALNOISESUMMARY ****/
  edm::Handle<HcalNoiseSummary> h_hcalnoisesummary;
  iEvent.getByToken( EDMHcalNoiseToken,h_hcalnoisesummary );
  
  /**** GET TRIGGER ****/
  //edm::Handle<pat::TriggerObjectStandAloneCollection> h_selectedtrigger;
  //iEvent.getByToken( EDMSelectedTriggerToken,h_selectedtrigger );
  //  pat::TriggerObjectStandAloneCollection const &selectedTrigger = *h_selectedtrigger;
  
  edm::Handle<edm::TriggerResults> h_triggerresults;
  iEvent.getByToken( EDMTriggerResultToken,h_triggerresults );
  edm::TriggerResults const &triggerResults = *h_triggerresults;
  
  /**** GET BEAMSPOT ****/
  edm::Handle<reco::BeamSpot> h_beamspot;
  iEvent.getByToken( EDMBeamSpotToken,h_beamspot );
  
  /**** GET PRIMARY VERTICES ****/
  edm::Handle< reco::VertexCollection > h_vtxs;
  iEvent.getByToken( EDMVertexToken,h_vtxs );
  reco::VertexCollection const &vtxs = *h_vtxs;
  
  // Primary Vertex Selection
  reco::VertexCollection selectedPVs;
  if( h_vtxs.isValid() ){
    for( reco::VertexCollection::const_iterator itvtx = vtxs.begin(); itvtx!=vtxs.end(); ++itvtx ){
      
      bool isGood = ( !(itvtx->isFake()) &&
		                  (itvtx->ndof() >= 4.0) &&
		                  (abs(itvtx->z()) <= 24.0) &&
		                  (abs(itvtx->position().Rho()) <= 2.0) 
		                );

      if( !isGood ) continue;
      
      selectedPVs.push_back(*itvtx);
    }
  }

  if( selectedPVs.size()>0 ) helper.SetVertex( selectedPVs.at(0) );
  if( selectedPVs.size() == 0 ) return;
  
  /**** GET LEPTONS ****/
  // MUONS
  edm::Handle< std::vector<pat::Muon> > h_muons;
  iEvent.getByToken( EDMMuonsToken,h_muons );
  std::vector<pat::Muon> const &muons = *h_muons; 
  std::vector<pat::Muon> selectedMuons = helper.GetSelectedMuons( muons, 20., muonID::muonTight );
  std::vector<pat::Muon> selectedMuonsLoose = helper.GetSelectedMuons( muons, 10., muonID::muonLoose );

  // ELECTRONS
  edm::Handle< std::vector<pat::Electron> > h_electrons;
  iEvent.getByToken( EDMElectronsToken,h_electrons );
  std::vector<pat::Electron> const &electrons = *h_electrons;
  std::vector<pat::Electron> selectedElectrons = helper.GetSelectedElectrons( electrons, 20., electronID::electronTight );
  std::vector<pat::Electron> selectedElectronsLoose = helper.GetSelectedElectrons( electrons, 10., electronID::electronLoose );

  // Leptons
  /*
  std::vector<pat::Pat>BNleptonCollection selectedLeptons_loose;
  for(size_t i=0; i<selectedElectronsLoose.size();i++){
    selectedLeptons_loose.push_back(&(selectedElectronsLoose[i]));
  }
  for(size_t i=0; i<selectedMuonsLoose.size();i++){
    selectedLeptons_loose.push_back(&(selectedMuonsLoose[i]));
  }    
  */

  /**** GET JETS ****/
  edm::Handle< std::vector<pat::Jet> > h_pfjets;
  iEvent.getByToken( EDMJetsToken,h_pfjets );
  std::vector<pat::Jet> const &pfjets = *h_pfjets;

  /**** GET PUPPI JETS (BECAUSE PUPPI PERFORMS BEST) ****/
  edm::Handle< std::vector<pat::Jet> > h_pfpuppijets;
  iEvent.getByToken( EDMPuppiJetsToken,h_pfpuppijets );
  std::vector<pat::Jet> const &pfpuppijets = *h_pfpuppijets;

  
  //  const JetCorrector* corrector = JetCorrector::getJetCorrector( "ak4PFchsL1L2L3", iSetup );   
  //helper.SetJetCorrector(corrector);
  
  // Get raw jets
  std::vector<pat::Jet> rawJets = helper.GetUncorrectedJets(pfjets);
  // Clean muons from jets
  std::vector<pat::Jet> jetsNoMu = helper.RemoveOverlaps(selectedMuonsLoose, rawJets);
  // Clean electrons from jets
  std::vector<pat::Jet> jetsNoEle = helper.RemoveOverlaps(selectedElectronsLoose, jetsNoMu);
  // Apply jet corrections
  //  std::vector<pat::Jet> correctedJets = helper.GetCorrectedJets(jetsNoEle, iEvent, iSetup);
  // Get jet Collection which pass selection
  std::vector<pat::Jet> selectedJets = helper.GetSelectedJets(jetsNoEle, 30., 2.4, jetID::jetLoose, '-' );
  // Get jet Collection which pass loose selection
  std::vector<pat::Jet> selectedJetsLoose = helper.GetSelectedJets(jetsNoEle, 20., 2.5, jetID::jetLoose, '-' );




  // Get raw puppi jets
  std::vector<pat::Jet> rawPuppiJets = helper.GetUncorrectedJets(pfpuppijets);
  // Clean muons from jets
  std::vector<pat::Jet> puppiJetsNoMu = helper.RemoveOverlaps(selectedMuonsLoose, rawPuppiJets);
  // Clean electrons from jets
  std::vector<pat::Jet> puppiJetsNoEle = helper.RemoveOverlaps(selectedElectronsLoose, puppiJetsNoMu);

  // Get puppi jet Collection which pass selection
  std::vector<pat::Jet> selectedPuppiJets = helper.GetSelectedJets(jetsNoEle, 30., 2.4, jetID::jetLoose, '-' );

  /**** GET MET ****/
  edm::Handle< std::vector<pat::MET> > h_pfmet;
  iEvent.getByToken( EDMMETsToken,h_pfmet );
  std::vector<pat::MET> const &pfMETs = *h_pfmet;

  /**** GET TOPJETS ****/
  /* edm::Handle<boosted::HEPTopJetCollection> h_heptopjet;
  boosted::HEPTopJetCollection heptopjets;
  if(useFatJets){
    iEvent.getByToken( EDMHEPTopJetsToken,h_heptopjet);
    boosted::HEPTopJetCollection const &heptopjets_unsorted = *h_heptopjet;
    heptopjets = tHqUtils::GetSortedByPt(heptopjets_unsorted);
    }*/
  
  /**** GET SUBFILTERJETS ****/
  /*  edm::Handle<boosted::SubFilterJetCollection> h_subfilterjet;                   
  boosted::SubFilterJetCollection subfilterjets;
  if(useFatJets){
    iEvent.getByToken( EDMSubFilterJetsToken,h_subfilterjet );
    boosted::SubFilterJetCollection const &subfilterjets_unsorted = *h_subfilterjet;
    subfilterjets = tHqUtils::GetSortedByPt(subfilterjets_unsorted);
    }*/
  
  /**** GET GENEVENTINFO ****/
  edm::Handle<GenEventInfoProduct> h_geneventinfo;
  iEvent.getByToken( EDMGenInfoToken, h_geneventinfo );
  
  /**** GET GENPARTICLES ****/
  edm::Handle< std::vector<reco::GenParticle> > h_genParticles;
  iEvent.getByToken( EDMGenParticlesToken,h_genParticles );
  std::vector<reco::GenParticle> const &genParticles = *h_genParticles;
  
  /**** GET GENJETS ****/
  edm::Handle< std::vector<reco::GenJet> > h_genjets;
  iEvent.getByToken( EDMGenJetsToken,h_genjets );
  std::vector<reco::GenJet> const &genjets = *h_genjets;
  std::vector<reco::GenJet> selectedGenJets;
  for(size_t i=0; i<genjets.size();i++){
    if(genjets[i].pt()>30&&fabs(genjets[i].eta())<2.5)
      selectedGenJets.push_back(genjets[i]);
  }
  
  // Fill tHq Event Object
  boosted::Event event = FillEvent(iEvent,h_geneventinfo,h_beamspot,h_hcalnoisesummary,h_puinfosummary);
  
  // Fill Trigger Info
  map<string,bool> triggerMap;
  for(auto name=relevantTriggers.begin(); name!=relevantTriggers.end();name++){
    unsigned int TriggerID =  hlt_config.triggerIndex(*name);
    if( TriggerID >= triggerResults.size() ) { 
      //      std::cerr <<"triggerID > trigger results.size: "<<TriggerID<<" > "<<triggerResults.size()<<std::endl; 
      triggerMap[*name]=false;
    }
    else{
      triggerMap[*name]=triggerResults.accept(TriggerID);
      std::cout << "Jo: Trigger ID: " << TriggerID << "  Name: " << *name << endl;
    }
  }

  for(auto name=triggerResults.begin(); name!=triggerResults.end();name++){
    //    unsigned int TriggerID =  hlt_config.triggerIndex(*name);
    std::cout << "Jo: Name: " << *name << endl;
  }

  TriggerInfo triggerInfo(triggerMap);



  // FIGURE OUT SAMPLE
  SampleType sampleType;
  if(isData)
    sampleType = SampleType::data;
  else if(tHqUtils::MCContainsTTbar(genParticles) && tHqUtils::MCContainsHiggs(genParticles)){
    sampleType = SampleType::tth;
  }
  else if(tHqUtils::MCContainsTTbar(genParticles)){
    sampleType = SampleType::tt;
  }
  else{
    sampleType = SampleType::nonttbkg;
  }

  // DO REWEIGHTING
  map<string,float> weights = GetWeights(event,selectedPVs,selectedJets,selectedElectrons,selectedMuons,genParticles);

  // DEFINE INPUT
  InputCollections input( event,
			  //selectedTrigger,
			  //triggerResults,
                          triggerInfo,
			  //			  hlt_config_,
			  selectedPVs,
			  selectedMuons,
			  selectedMuonsLoose,
			  selectedElectrons,
                          selectedElectronsLoose,
                          selectedJets,
			  selectedPuppiJets,
                          selectedJetsLoose,
                          pfMETs,
			  // heptopjets,
                          //subfilterjets,
                          genParticles,
                          selectedGenJets,
                          sampleType,
                          weights,
			  iSetup,
                          iEvent
			  );
  InputCollections unselected_input( event,
				     //			  selectedTrigger,
				     //			  triggerResults,
			  triggerInfo,
				     //			  hlt_config_,
			  vtxs,
			  muons,
			  muons,
			  electrons,
                          electrons,
                          pfjets,
                          pfjets,
			  pfpuppijets,
                          pfMETs,
				     //heptopjets,
				     //subfilterjets,
                          genParticles,
                          selectedGenJets,
                          sampleType,
                          weights,
			  iSetup,
                          iEvent
			  );
        
  
  // DO SELECTION
  cutflow.EventSurvivedStep("all");
  bool selected=true;
  for(size_t i=0; i<selections.size() && selected; i++){
    if(disableObjectSelections){
      if(!selections.at(i)->IsSelected(unselected_input,cutflow))
	selected=false;
    }
    else {
      if(!selections.at(i)->IsSelected(input,cutflow))
	selected=false;
    }
  }

  if(!selected) return;    

  selectedElectrons = ElectronSelection(selectedElectrons,input.selectedPVs[0].position());
  selectedMuons = MuonSelection(selectedMuons,input.selectedPVs[0].position());
  selectedJets = JetSelection(selectedJets,selectedElectrons, selectedMuons, input);

  // WRITE TREE
  if(disableObjectSelections)
    treewriter.Process(unselected_input);  
  else
    treewriter.Process(input);  
}

std::vector<pat::Muon> tHqAnalyzer::MuonSelection( std::vector<pat::Muon> selectedMuons, const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> pvposition){


  std::vector<pat::Muon> newselectedMuons;
  for( std::vector<pat::Muon>::const_iterator it = selectedMuons.begin(), ed = selectedMuons.end(); it != ed; ++it ){
    pat::Muon iMuon = *it;
    bool passesKinematics=false;
    bool isPFandGlobal=false;
    bool passesISO=false;
    bool passesGlobalTrackID=false;
    bool passesBestTrackID=false;
    bool passesinnerTrackID=false;
    bool passesTrackID=false;
    bool passesID=false;
    
    if(iMuon.globalTrack().isAvailable()){
    double chi2ndof = iMuon.globalTrack()->normalizedChi2();
    int nValidHits = iMuon.globalTrack()->hitPattern().numberOfValidMuonHits();
    passesGlobalTrackID=(chi2ndof<10.0 && nValidHits>0);
    } 

    if(iMuon.muonBestTrack().isAvailable()){
    double d0 = fabs(iMuon.muonBestTrack()->dxy(pvposition));
    double dZ = fabs(iMuon.muonBestTrack()->dz(pvposition));
    passesBestTrackID=(d0<0.2 && dZ<0.5);
    }

    int nMatchedStations = iMuon.numberOfMatchedStations();
        
    if(iMuon.track().isAvailable()){
    int nTrackerLayers = iMuon.track()->hitPattern().trackerLayersWithMeasurement();
    passesTrackID=(nTrackerLayers>5);
    }

    
    if(iMuon.innerTrack().isAvailable()){
    int nValidPiyles = iMuon.innerTrack()->hitPattern().numberOfValidPixelHits();
    passesinnerTrackID=(nValidPiyles>0);
    }
    double relIso = (iMuon.pfIsolationR04().sumChargedHadronPt + std::max( iMuon.pfIsolationR04().sumNeutralHadronEt + iMuon.pfIsolationR04().sumPhotonEt - 0.5 * iMuon.pfIsolationR04().sumPUPt,0.0)) / iMuon.pt();
    
    isPFandGlobal = (iMuon.isPFMuon() && iMuon.isGlobalMuon());
    passesKinematics = (iMuon.pt()>20 && fabs(iMuon.eta())<2.4);
    passesISO = (relIso<0.12);
    passesID = (passesinnerTrackID && passesGlobalTrackID && passesBestTrackID && passesTrackID && nMatchedStations>1);
    
    if(passesKinematics && passesISO && passesID && isPFandGlobal){
    newselectedMuons.push_back(*it);
    }
    else std::cout << " Bad Muon found. " << std::endl;
  }
  return newselectedMuons;
}


std::vector<pat::Electron> tHqAnalyzer::ElectronSelection( std::vector<pat::Electron> selectedElectrons, const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> pvposition){
  //-------------------------------------------------------------------------                                                                              
  //do electrons                              
  
  std::vector<pat::Electron> newselectedElectrons;
  
  for( std::vector<pat::Electron>::const_iterator it = selectedElectrons.begin(), ed = selectedElectrons.end(); it != ed; ++it ){
    pat::Electron iElectron = *it;
    bool passesKinematics=false;
    bool inCrack=true;
    bool passesConversion=false;
    bool passesISO=false;
    bool passessID=false;

    //check if barrel or endcap supercluster                                                                                                             
    double SCeta = (iElectron.superCluster().isAvailable()) ? iElectron.superCluster()->position().eta() : 99;
    double absSCeta = fabs(SCeta);
    bool isEB = ( absSCeta <= 1.479 );
    bool isEE = ( absSCeta > 1.479 && absSCeta < 2.5 );

    //isolation                                                                                                                                          
    double pfIsoCharged = iElectron.pfIsolationVariables().sumChargedHadronPt;
    double pfIsoNeutralHadron = iElectron.pfIsolationVariables().sumNeutralHadronEt;
    double pfIsoNeutralPhoton = iElectron.pfIsolationVariables().sumPhotonEt;
    double pfIsoSumPUPt = iElectron.pfIsolationVariables().sumPUPt;

    double relIso = (pfIsoCharged + max( pfIsoNeutralHadron + pfIsoNeutralPhoton - 0.5*pfIsoSumPUPt, 0.0 ))/iElectron.pt();


    //other stuff                                                                                                                                        
    double full5x5_sigmaIetaIeta = iElectron.full5x5_sigmaIetaIeta();
    double dEtaIn = fabs( iElectron.deltaEtaSuperClusterTrackAtVtx() );
    double dPhiIn = fabs( iElectron.deltaPhiSuperClusterTrackAtVtx() );
    double hOverE = iElectron.hcalOverEcal();

    double ooEmooP = 999;
    if( iElectron.ecalEnergy() == 0 ) ooEmooP = 1e30;
    else if( !std::isfinite(iElectron.ecalEnergy()) ) ooEmooP = 1e30;
    else ooEmooP = fabs(1.0/iElectron.ecalEnergy() - iElectron.eSuperClusterOverP()/iElectron.ecalEnergy() );

    double d0 = 999;
    double dZ = 999;
    double expectedMissingInnerHits = 999;
    if( iElectron.gsfTrack().isAvailable() ){
      d0 = fabs(iElectron.gsfTrack()->dxy(pvposition));
      dZ = fabs(iElectron.gsfTrack()->dz(pvposition));
      expectedMissingInnerHits = iElectron.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
    }

    //do the checks                                                                                                                                        
    passesConversion = ( iElectron.passConversionVeto() );
    passesKinematics = (iElectron.pt()>20 && fabs(iElectron.eta())<2.4);
    if(iElectron.superCluster().isAvailable()){
      inCrack = (fabs(iElectron.superCluster()->position().eta())>1.4442 && fabs(iElectron.superCluster()->position().eta())<1.5660);
    }

    //dZ and d0 cuts are different than in the posted recipe                                                                                               
    if(isEB){
      passessID = (full5x5_sigmaIetaIeta<0.010399 && dEtaIn<0.007641 && dPhiIn<0.032643 &&hOverE<0.060662 && d0<0.011811 && dZ<0.070775 && expectedMissingInnerHits<=1 && ooEmooP<0.153897);
      passesISO = (relIso<0.097213);
    }
    else if(isEE){
      passessID = (full5x5_sigmaIetaIeta<0.029524 && dEtaIn<0.009285 && dPhiIn<0.042447 &&hOverE<0.104263 && d0<0.051682 && dZ<0.180720 && expectedMissingInnerHits<=1 && ooEmooP<0.137468);
      passesISO = (relIso<0.116708);
    }
    else{
      passessID=false;
      passesISO=false;
    }

    //    std::cout<<relIso<<" "<<isEB<<std::endl;
    //    std::cout<<dZ<<" "<<d0<<std::endl;
    //    std::cout<<passesKinematics<<" "<<passesConversion<<" "<<passesISO<<" "<<passessID<<" "<<inCrack<<std::endl;

    if(passesKinematics && passesConversion && passesISO && passessID && !inCrack){
      newselectedElectrons.push_back(*it);
    }
    else std::cout << " Bad Electron found. " << std::endl;
  }
  return newselectedElectrons;
}


std::vector<pat::Jet> tHqAnalyzer::JetSelection( std::vector<pat::Jet> selectedJets, std::vector<pat::Electron> selectedElectrons, std::vector<pat::Muon> selectedMuons, InputCollections input){
  
  //----------------------------------------------------------------------------------  
  //do jets

  std::vector<pat::Jet> cleanEleJets;
  std::vector<pat::Jet> cleanMuonJets;
  std::vector<pat::Jet> correctedJets;
  std::vector<pat::Jet> sortedJets;
  //  std::vector<pat::Jet> selectedJets;
  std::vector<pat::Jet> taggedJets;

  //first uncorrect the jets
                                                                                                                               
  std::vector<pat::Jet> bufferJets;
  for( std::vector<pat::Jet>::const_iterator it = input.selectedJets.begin(), ed = input.selectedJets.end(); it != ed; ++it ){
    pat::Jet iJet = *it;
    //    double originalPt = iJet.pt();
    //   std::cout<<iJet.currentJECLevel()<<" "<<iJet.currentJECSet()<<std::endl;
    math::XYZTLorentzVector uncorrectedP4 = iJet.correctedP4("Uncorrected");
    //   math::XYZTLorentzVector ALTuncorrectedP4 = iJet.correctedJet(0).p4();                                                                           
    iJet.setP4(uncorrectedP4);
    //    double uncorrectedPt = iJet.pt();
    //   std::cout<<originalPt<<" "<<uncorrectedPt<<std::endl;                        
    //   iJet.setP4(ALTuncorrectedP4);
    //   uncorrectedPt = iJet.pt();
    //   std::cout<<iJet.currentJECLevel()<<" "<<iJet.currentJECSet()<<std::endl;

    //    std::cout<< "Original JetPt: " << originalPt<<"   -   Uncorrected JetPt: "<<uncorrectedPt<<std::endl;
    bufferJets.push_back(iJet);
  }
  
  bool doNewCleaning = false;
  
  if(doNewCleaning){
    //cleaning like in miniAODhelper 
    //    std::cout<<"doing new cleaning"<<std::endl;
    // std::cout<<"Number of Leptons: Muons: "<< selectedMuons.size() << "   Electrons:  " << selectedElectrons.size() << std::endl;
    for( std::vector<pat::Jet>::const_iterator it = bufferJets.begin(), ed = bufferJets.end(); it != ed; ++it ){
      //std::cout<< "BufferJet pt: " << it->pt()<<std::endl;
    }
    cleanEleJets = tHqUtils::RemoveOverlaps(selectedElectrons, bufferJets);
    for( std::vector<pat::Jet>::const_iterator it = cleanEleJets.begin(), ed = cleanEleJets.end(); it != ed; ++it ){
      //std::cout<< "After electron removal: " << it->pt()<<std::endl;       
    }
    cleanMuonJets = tHqUtils::RemoveOverlaps(selectedMuons, cleanEleJets);
    for( std::vector<pat::Jet>::const_iterator it = cleanMuonJets.begin(), ed = cleanMuonJets.end(); it != ed; ++it ){
      //std::cout<< "After muon removal: " << it->pt()<<std::endl;                                                                                                             
    }
  }
  else{
    //clean from electrons  
    for(std::vector<pat::Electron>::const_iterator iEle = selectedElectrons.begin(), ed = selectedElectrons.end(); iEle != ed; ++iEle ){
      pat::Electron Ele = *iEle;
      double maxDeltaR=0.4;
      double minDeltaR=maxDeltaR;
      double dR=999.0;
      int matchindex=-1;
      int counter=0;

      //find jet closes to electron
      for( std::vector<pat::Jet>::const_iterator it = bufferJets.begin(), ed = bufferJets.end(); it != ed; ++it ){
	pat::Jet iJet = *it;
	dR = reco::deltaR(Ele.eta(), Ele.phi(), iJet.eta(), iJet.phi());
	if(dR<minDeltaR){
	  minDeltaR = dR;
	  matchindex=counter;
        }
	counter++;
      }

      // now clean the closest jet from the electron and put it back with the rest                                                                           
      counter=0;
      for( std::vector<pat::Jet>::const_iterator it = bufferJets.begin(), ed = bufferJets.end(); it != ed; ++it ){
	pat::Jet iJet = *it;
	//	std::cout<< "BufferJet pt: " << it->pt()<<std::endl;
	if(matchindex==counter){
	  //       std::cout<<"ele jet "<<Ele.pt()<<" "<<iJet.pt()<<std::endl;
	  //       std::cout<<"ele jet "<<Ele.energy()<<" "<<iJet.energy()<<std::endl;
	  math::XYZTLorentzVector original = iJet.p4();
	  original -= Ele.p4();
	  iJet.setP4(original);
	}
	if(iJet.pt()>0.0 && iJet.energy()>0.0){
	  cleanEleJets.push_back(iJet);
	}
	//	std::cout<< "After Electron removal pt: " << iJet.pt()<<std::endl;
	counter++;


	//     std::cout<<Ele.eta()<<" "<<Ele.phi()<<" "<<iJet.eta()<<" "<<iJet.phi()<<std::endl;                                                                
      }
      bufferJets.clear();
      for( std::vector<pat::Jet>::const_iterator it = cleanEleJets.begin(), ed = cleanEleJets.end(); it != ed; ++it ){
	bufferJets.push_back(*it);
      }
      cleanEleJets.clear();
    }
    for( std::vector<pat::Jet>::const_iterator it = bufferJets.begin(), ed = bufferJets.end(); it != ed; ++it ){
      cleanEleJets.push_back(*it);
    }
    //now we should have the jets cleaned from the electrons                                                                                                 
    //clean from muos                                                                                                                                   

    bufferJets.clear();
    for( std::vector<pat::Jet>::const_iterator it = cleanEleJets.begin(), ed = cleanEleJets.end(); it != ed; ++it ){
      bufferJets.push_back(*it);
    }

    for(std::vector<pat::Muon>::const_iterator iMuon = selectedMuons.begin(), ed = selectedMuons.end(); iMuon != ed; ++iMuon ){
      pat::Muon Muon = *iMuon;
      double maxDeltaR=0.4;
      double minDeltaR=maxDeltaR;
      double dR=999.0;
      int matchindex=-1;
      int counter=0;

      //find jet closes to muon

      for( std::vector<pat::Jet>::const_iterator it = bufferJets.begin(), ed = bufferJets.end(); it != ed; ++it ){
	pat::Jet iJet = *it;
	dR = reco::deltaR(Muon.eta(), Muon.phi(), iJet.eta(), iJet.phi());
	if(dR<minDeltaR){
	  minDeltaR = dR;
	  matchindex=counter;
        }
	counter++;
      }
      // now clean the closest jet from the muon and put it back with the rest
      counter=0;
      for( std::vector<pat::Jet>::const_iterator it = bufferJets.begin(), ed = bufferJets.end(); it != ed; ++it ){
	pat::Jet iJet = *it;
	if(matchindex==counter){
	  //         std::cout<<"muon jet "<<Muon.pt()<<" "<<iJet.pt()<<std::endl;
	  //         std::cout<<"muon jet "<<Muon.energy()<<" "<<iJet.energy()<<std::endl;

	  math::XYZTLorentzVector original = iJet.p4();
	  original -= Muon.p4();
	  iJet.setP4(original);
	  //	  std::cout<< "After Muon removal pt: " << iJet.pt()<<std::endl;
	}
	if(iJet.pt()>0.0 && iJet.energy()>0.0){
	  cleanMuonJets.push_back(iJet);
	}
	counter++;
      }
      bufferJets.clear();
      for( std::vector<pat::Jet>::const_iterator it = cleanMuonJets.begin(), ed = cleanMuonJets.end(); it != ed; ++it ){
	bufferJets.push_back(*it);
      }
      cleanMuonJets.clear();
    }
    for( std::vector<pat::Jet>::const_iterator it = bufferJets.begin(), ed = bufferJets.end(); it != ed; ++it ){
      cleanMuonJets.push_back(*it);
      //   std::cout<<it->pt()<<std::endl;
    }
  }
  //now correct the jets again
  //using the PHYS14_25_V2 corrections i hope :)  
  /*  "ak4PFchsL1L2L3"*/

    //  const JetCorrector* jetCorrector = JetCorrector::getJetCorrector("ak4PFchsL1L2L3",input.setup);
  correctedJets.clear();

  //  std::cout<<"setup the corrector"<<std::endl;
  for( std::vector<pat::Jet>::const_iterator it = cleanMuonJets.begin(), ed = cleanMuonJets.end(); it != ed; ++it ){
    pat::Jet iJet = *it;
    double jec =1.0;
    //jec = jetCorrector->correction(iJet, input.edmevent, input.setup );
    //    std::cout<<"uncorrected "<<jec<<" "<<iJet.pt()<<std::endl; 
    iJet.scaleEnergy(jec);
    ///    std::cout<<"corrected "<<iJet.pt()<<std::endl;     
    //    std::cout<<iJet.currentJECLevel()<<" "<<iJet.currentJECSet()<<" "<<std::endl;
    
    correctedJets.push_back(iJet);
    //   double uncorrectedPt = iJet.pt();    
    //   math::XYZTLorentzVector correctedP4 = iJet.correctedP4("L1FastJetL2RelativeL3Absolute","","PHYS14_25_V2");                                  
    //   iJet.setP4(correctedP4);  
    //   pat::Jet newJet = iJet.correctedJet("L3Absolute","none","patJetCorrFactors"); 
    //   double correctedPt = newJet.pt();    
    //   std::cout<<uncorrectedPt<<" "<<correctedPt<<std::endl;  
    //   bufferJets.push_back(iJet);  
    //   for(unsigned int i=0; i<iJet.availableJECSets().size(); i++){   
    //     std::cout<<iJet.availableJECSets()[i]<<std::endl;   
    //   }
    //   for(unsigned int i=0; i<iJet.availableJECLevels().size(); i++){
    //     std::cout<<iJet.availableJECLevels("patJetCorrFactors")[i]<<std::endl;  
    //   } 
  }

  // std::cout<<"The Parameter Set"<<std::endl;                                                                      
  // std::cout<<Config.dump()<<std::endl;;                                                                                                
  //now sort the jets by pt                                                                                                               
  
  std::sort(correctedJets.begin(), correctedJets.end(), tHqUtils::FirstJetIsHarder);
  //  std::cout<<correctedJets.at(0).pt()<<" "<<correctedJets.at(1).pt()<<" "<<correctedJets.at(2).pt()<<std::endl; 
  

  //now select only the good jets                                                                
                                                          
  selectedJets.clear();
  for( std::vector<pat::Jet>::const_iterator it = correctedJets.begin(), ed = correctedJets.end(); it != ed; ++it ){
    pat::Jet iJet = *it;

    bool isselected=false;
    isselected=(iJet.pt()>25.0 && abs(iJet.eta())<2.4 && iJet.neutralHadronEnergyFraction()<0.99 && iJet.chargedHadronEnergyFraction()>0.0 && iJet.chargedMultiplicity()>0.0 && iJet.chargedEmEnergyFraction()<0.99 && iJet.neutralEmEnergyFraction()<0.99 && iJet.numberOfDaughters()>1 );
    if(isselected){
      selectedJets.push_back(iJet);
    }
  }

  //get the btagged jets                                                                                                                                   
  taggedJets.clear();
  for( std::vector<pat::Jet>::const_iterator it = selectedJets.begin(), ed = selectedJets.end(); it != ed; ++it ){
    pat::Jet iJet = *it;

    bool istagged=false;
    double workingPoint=0.814;
    istagged=(iJet.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags") > workingPoint);
    if(istagged){
      taggedJets.push_back(iJet);
    }
  }

  return selectedJets;

}



boosted::Event tHqAnalyzer::FillEvent(const edm::Event& iEvent, const edm::Handle<GenEventInfoProduct>& genEvtInfo, const edm::Handle<reco::BeamSpot>& beamSpot, const edm::Handle<HcalNoiseSummary>& hcalNoiseSummary, const edm::Handle< std::vector<PileupSummaryInfo> >& puSummaryInfo){
 
  boosted::Event event = boosted::Event();
  
  event.evt         = iEvent.id().event();
  event.run         = iEvent.id().run();
  event.sample      = sampleID;
  event.lumiBlock   = iEvent.luminosityBlock();
  
  
  if(genEvtInfo.isValid()){
    std::vector<double> genWeights = genEvtInfo->weights();
    for(size_t i=0;i<genWeights.size();i++){
      event.weight *= genWeights[i];
    }

    event.qScale = genEvtInfo->qScale();
    event.alphaQCD = genEvtInfo->alphaQCD();
    event.alphaQED = genEvtInfo->alphaQED();
    event.pthat = ( genEvtInfo->hasBinningValues() ? (genEvtInfo->binningValues())[0] : 0.0);
    event.scalePDF = genEvtInfo->pdf()->scalePDF;
    event.x1 = genEvtInfo->pdf()->x.first;
    event.x2 = genEvtInfo->pdf()->x.second;
    event.xPDF1 = genEvtInfo->pdf()->xPDF.first;
    event.xPDF2 = genEvtInfo->pdf()->xPDF.second;
    event.id1 = genEvtInfo->pdf()->id.first;
    event.id2 = genEvtInfo->pdf()->id.second;
  }
  
  if(beamSpot.isValid()){
    event.BSx = beamSpot->x0();
    event.BSy = beamSpot->y0();
    event.BSz = beamSpot->z0();
  }
  
 
  if( hcalNoiseSummary.isValid() ){
    event.hcalnoiseLoose = hcalNoiseSummary->passLooseNoiseFilter();
    event.hcalnoiseTight = hcalNoiseSummary->passTightNoiseFilter();
  }
  
  if( puSummaryInfo.isValid() ){
    for(std::vector<PileupSummaryInfo>::const_iterator PVI = puSummaryInfo->begin(); PVI != puSummaryInfo->end(); ++PVI) {

      int BX = PVI->getBunchCrossing();

      event.sumNVtx  += float(PVI->getPU_NumInteractions());
      event.sumTrueNVtx += float(PVI->getTrueNumInteractions());

      if( BX==0 ){
	      event.numGenPV = PVI->getPU_NumInteractions();
	      event.numTruePV = PVI->getTrueNumInteractions();
      }

      if(BX == -1) { 
	      event.nm1 = PVI->getPU_NumInteractions();
	      event.nm1_true = PVI->getTrueNumInteractions();
      }
      else if(BX == 0) { 
	      event.n0 = PVI->getPU_NumInteractions();
	      event.n0_true = PVI->getTrueNumInteractions();
      }
      else if(BX == 1) { 
	      event.np1 = PVI->getPU_NumInteractions();
	      event.np1_true = PVI->getTrueNumInteractions();
      }
    }
  }
  
  return event;
}


map<string,float> tHqAnalyzer::GetWeights(const boosted::Event& event, const reco::VertexCollection& selectedPVs, const std::vector<pat::Jet>& selectedJets, const std::vector<pat::Electron>& selectedElectrons, const std::vector<pat::Muon>& selectedMuons, const std::vector<reco::GenParticle>& genParticles){
  map<string,float> weights;
  
  if(isData){
    weights["Weight"] = 1.0;
    weights["Weight_XS"] = 1.0;
    weights["Weight_CSV"] = 1.0;
    weights["Weight_PU"] = 1.0;
    weights["Weight_TopPt"] = 1.0;
    return weights;
  }

  // not sure why the BNevent weight is !=+-1 but we dont want it that way
  float weight = event.weight;  
  float xsweight = xs*luminosity/totalMCevents;
  float csvweight = 1.;
  float puweight = 1.;
  float topptweight = 1.;
  //float csvweight = beanHelper.GetCSVweight(selectedJets,jsystype);
  //float puweight = beanHelper.GetPUweight(event[0].numTruePV);
  //float topptweight = beanHelper.GetTopPtweight(mcparticlesStatus3);
  //float q2scaleweight = beanHelper.GetQ2ScaleUp(const BNevent&);
  
  //NA, JERup, JERdown, JESup, JESdown, hfSFup, hfSFdown, lfSFdown, lfSFup, TESup, TESdown, 
  //CSVLFup, CSVLFdown, CSVHFup, CSVHFdown, CSVHFStats1up, CSVHFStats1down, CSVLFStats1up, CSVLFStats1down, CSVHFStats2up, CSVHFStats2down, CSVLFStats2up, CSVLFStats2down, CSVCErr1up, CSVCErr1down, CSVCErr2up, CSVCErr2down
  
  weight *= xsweight*csvweight*puweight*topptweight;
  weights["Weight"] = weight;
  weights["Weight_XS"] = xsweight;
  weights["Weight_CSV"] = csvweight;
  weights["Weight_PU"] = puweight;
  weights["Weight_TopPt"] = topptweight;
  
  /*
  if(doSystematics && jsystype != sysType::JESup && jsystype != sysType::JESup){
    weights["Weight_TopPtup"] = beanHelper.GetTopPtweightUp(mcparticlesStatus3)/topptweight;
    weights["Weight_TopPtdown"] = beanHelper.GetTopPtweightDown(mcparticlesStatus3)/topptweight;

    weights["Weight_CSVLFup"] = beanHelper.GetCSVweight(selectedJets,sysType::CSVLFup)/csvweight;
    weights["Weight_CSVLFdown"] = beanHelper.GetCSVweight(selectedJets,sysType::CSVLFdown)/csvweight;
    weights["Weight_CSVHFup"] = beanHelper.GetCSVweight(selectedJets,sysType::CSVHFup)/csvweight;
    weights["Weight_CSVHFdown"] = beanHelper.GetCSVweight(selectedJets,sysType::CSVHFdown)/csvweight;
    weights["Weight_CSVHFStats1up"] = beanHelper.GetCSVweight(selectedJets,sysType::CSVHFStats1up)/csvweight;
    weights["Weight_CSVHFStats1down"] = beanHelper.GetCSVweight(selectedJets,sysType::CSVHFStats1down)/csvweight;
    weights["Weight_CSVLFStats1up"] = beanHelper.GetCSVweight(selectedJets,sysType::CSVLFStats1up)/csvweight;
    weights["Weight_CSVLFStats1down"] = beanHelper.GetCSVweight(selectedJets,sysType::CSVLFStats1down)/csvweight;
    weights["Weight_CSVHFStats2up"] = beanHelper.GetCSVweight(selectedJets,sysType::CSVHFStats2up)/csvweight;
    weights["Weight_CSVHFStats2down"] = beanHelper.GetCSVweight(selectedJets,sysType::CSVHFStats2down)/csvweight;
    weights["Weight_CSVLFStats2up"] = beanHelper.GetCSVweight(selectedJets,sysType::CSVLFStats2up)/csvweight;
    weights["Weight_CSVLFStats2down"] = beanHelper.GetCSVweight(selectedJets,sysType::CSVLFStats2down)/csvweight;
    weights["Weight_CSVCErr1up"] = beanHelper.GetCSVweight(selectedJets,sysType::CSVCErr1up)/csvweight;
    weights["Weight_CSVCErr1down"] = beanHelper.GetCSVweight(selectedJets,sysType::CSVCErr1down)/csvweight;
    weights["Weight_CSVCErr2up"] = beanHelper.GetCSVweight(selectedJets,sysType::CSVCErr2up)/csvweight;
    weights["Weight_CSVCErr2down"] = beanHelper.GetCSVweight(selectedJets,sysType::CSVCErr2down)/csvweight;
    weights["Weight_Q2up"] = beanHelper.GetQ2ScaleUp(event[0]);
    weights["Weight_Q2down"] = beanHelper.GetQ2ScaleDown(event[0]);

  }
  */
  
  return weights;
}


// ------------ method called once each job just before starting event loop  ------------
void 
tHqAnalyzer::beginJob()
{
  eventcount=0;

  watch.Start();
}


// ------------ method called once each job just after ending the event loop  ------------
void 
tHqAnalyzer::endJob() 
{
  treewriter.AddSampleInformation();
  cutflow.Print();
}
// ------------ method called when starting to processes a run ------------
// needed for the hlt_config_
void
tHqAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
std::string hltTag="HLT";
bool hltchanged = true;
if (!hlt_config.init(iRun, iSetup, hltTag, hltchanged)) {
std::cout << "Warning, didn't find trigger process HLT,\t" << hltTag << std::endl;
return;
}
}
 

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
tHqAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(tHqAnalyzer);
