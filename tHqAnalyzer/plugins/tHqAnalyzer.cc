// -*- C++ -*-
//
// Package:    tHqAnalysis/tHqAnalyzer
// Class:      tHqAnalyzer
// 
/**\class tHqAnalyzer tHqAnalyzer.cc tHqAnalysis/tHqAnalyzer/plugins/tHqAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Shawn Williamson, Hannes Mildner
// Edited for thq analysis: Simon Fink
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
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "MiniAOD/MiniAODHelper/interface/MiniAODHelper.h"

#include "tHqAnalysis/tHqAnalyzer/interface/tHqUtils.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/InputCollections.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/Cutflow.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/EventInfo.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/TreeWriter.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/TriggerInfo.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/HistoReweighter.hpp"

#include "tHqAnalysis/tHqAnalyzer/interface/Selection.hpp"

#include "tHqAnalysis/tHqAnalyzer/interface/WeightProcessor.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/TopGenVarProcessor.hpp"

#include "tHqAnalysis/tHqAnalyzer/interface/MVAVarProcessor.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/BaseVarProcessor.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/RecoVarProcessor.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/GenTopEvent.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/GentHqEvent.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/tHqGenVarProcessor.hpp"


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
  map<string,float> GetWeights(const GenEventInfoProduct& genEventInfo, const EventInfo& eventInfo, const reco::VertexCollection& selectedPVs, const std::vector<pat::Jet>& selectedJets, const std::vector<pat::Electron>& selectedElectrons, const std::vector<pat::Muon>& selectedMuons, const std::vector<reco::GenParticle>& genParticles, sysType::sysType systype=sysType::NA);//, edm::Run const& iRun);
  void GetSystWeights(const LHEEventProduct& LHEEvent, vector<string> &weight_syst_id, vector<float> &weight_syst, float &Weight_orig);
  std::vector<pat::Electron> ElectronSelection( std::vector<pat::Electron> selectedElectrons,  const reco::VertexCollection& selectedPVs);
  std::vector<pat::Muon> MuonSelection( std::vector<pat::Muon> selectedMuons,  const reco::VertexCollection& selectedPVs);      
  std::vector<pat::Jet> JetSelection( std::vector<pat::Jet> selectedJets, std::vector<pat::Electron> selectedElectrons, std::vector<pat::Muon> selectedMuons, InputCollections input);
      
      
      
      // ----------member data ---------------------------
      
      /** the beanhelper is used for selections and reweighting */
      MiniAODHelper helper;
      
      /** Reweighter to match the PV distribution in data*/
  //      HistoReweighter pvWeight;


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

      /** calculate and store systematic weights? */
      bool doSystematics;
  
      /** use GenBmatching info? this is only possible if the miniAOD contains them */
      bool useGenHadronMatch;
    
      /** jet systematic that is applied (the outher systematics are done at a different place with reweighting)*/
      sysType::sysType jsystype;
      
      
      /** pu summary data access token **/
      edm::EDGetTokenT< std::vector<PileupSummaryInfo> > EDMPUInfoToken;
      
      /** pileup density data access token **/
      edm::EDGetTokenT <double> EDMRhoToken;


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

      /** gen info data access token **/
      edm::EDGetTokenT< LHEEventProduct > EDMLHEEventToken;
      edm::EDGetTokenT< LHEEventProduct >   EDMLHEEventToken_alt;      

      /** gen particles data access token **/
      edm::EDGetTokenT< std::vector<reco::GenParticle> > EDMGenParticlesToken;
      
      /** gen jets data access token **/
      edm::EDGetTokenT< std::vector<reco::GenJet> > EDMGenJetsToken;
   
      // custom genjets for tt+X categorization
      edm::EDGetTokenT< std::vector<reco::GenJet> > EDMCustomGenJetsToken;
        
      /** tt+X categorization tokens **/
      edm::EDGetTokenT<std::vector<int> > genBHadJetIndexToken;
      edm::EDGetTokenT<std::vector<int> > genBHadFlavourToken;
      edm::EDGetTokenT<std::vector<int> > genBHadFromTopWeakDecayToken;
      edm::EDGetTokenT<std::vector<reco::GenParticle> > genBHadPlusMothersToken;
      edm::EDGetTokenT<std::vector<std::vector< int > > > genBHadPlusMothersIndicesToken;
      edm::EDGetTokenT<std::vector<int> > genBHadIndexToken;
      edm::EDGetTokenT<std::vector<int> > genBHadLeptonHadronIndexToken;
      edm::EDGetTokenT<std::vector<int> > genBHadLeptonViaTauToken;
      edm::EDGetTokenT<std::vector<int> > genCHadJetIndexToken;
      edm::EDGetTokenT<std::vector<int> > genCHadFlavourToken;
      edm::EDGetTokenT<std::vector<int> > genCHadFromTopWeakDecayToken;
      edm::EDGetTokenT<std::vector<int> > genCHadBHadronIdToken;
      edm::EDGetTokenT<std::vector<int> > genCHadIndexToken;
      edm::EDGetTokenT<std::vector<reco::GenParticle> > genCHadPlusMothersToken;
      edm::EDGetTokenT<int> genTtbarIdToken;
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
tHqAnalyzer::tHqAnalyzer(const edm::ParameterSet& iConfig){
  //  pvWeight = ((tHqUtils::GetAnalyzerPath()+"/data/pvweights/PUhistos.root").c_str(),"data","mc")
  std::string era = iConfig.getParameter<std::string>("era");
  string analysisType = iConfig.getParameter<std::string>("analysisType");
  analysisType::analysisType iAnalysisType = analysisType::LJ;
  if(analysisType == "LJ") iAnalysisType = analysisType::LJ;
  else cerr << "No matching analysis type found for: " << analysisType << endl;
  luminosity = iConfig.getParameter<double>("luminostiy");
  sampleID = iConfig.getParameter<int>("sampleID");
  xs = iConfig.getParameter<double>("xs");
  totalMCevents = iConfig.getParameter<int>("nMCEvents");
  isData = iConfig.getParameter<bool>("isData");

  useGenHadronMatch = iConfig.getParameter<bool>("useGenHadronMatch");

  //  useFatJets = iConfig.getParameter<bool>("useFatJets");

  string outfileName = iConfig.getParameter<std::string>("outfileName");
  
  std::cout << "Outfile Name: " << outfileName << std::endl;
  
  // REGISTER DATA ACCESS
  EDMPUInfoToken          = consumes< std::vector<PileupSummaryInfo> >(edm::InputTag("addPileupInfo","",""));
  EDMRhoToken             = consumes <double> (edm::InputTag(std::string("fixedGridRhoFastjetAll")));
  EDMHcalNoiseToken       = consumes< HcalNoiseSummary >(edm::InputTag("hcalnoise","",""));
  EDMSelectedTriggerToken = consumes< pat::TriggerObjectStandAloneCollection > (edm::InputTag("selectedPatTrigger","",""));
  EDMTriggerResultToken   = consumes< edm::TriggerResults > (edm::InputTag("TriggerResults","","HLT"));
  EDMBeamSpotToken        = consumes< reco::BeamSpot > (edm::InputTag("offlineBeamSpot","",""));
  EDMVertexToken          = consumes< reco::VertexCollection > (edm::InputTag("offlineSlimmedPrimaryVertices","","PAT"));
  EDMMuonsToken           = consumes< std::vector<pat::Muon> >(edm::InputTag("slimmedMuons","","PAT"));
  EDMElectronsToken       = consumes< std::vector<pat::Electron> >(edm::InputTag("slimmedElectrons","","PAT"));
  EDMJetsToken            = consumes< std::vector<pat::Jet> >(edm::InputTag("slimmedJets","","PAT"));
  EDMPuppiJetsToken       = consumes< std::vector<pat::Jet> >(edm::InputTag("slimmedJetsPuppi","","PAT"));
  EDMMETsToken            = consumes< std::vector<pat::MET> >(edm::InputTag("slimmedMETs","","PAT"));
  //EDMHEPTopJetsToken      = consumes< boosted::HEPTopJetCollection >(edm::InputTag("HEPTopJetsPFMatcher","heptopjets","p"));
  // EDMSubFilterJetsToken   = consumes< boosted::SubFilterJetCollection >(edm::InputTag("CA12JetsCA3FilterjetsPFMatcher","subfilterjets","p"));
  EDMGenInfoToken         = consumes< GenEventInfoProduct >(edm::InputTag("generator","","SIM"));
  EDMLHEEventToken        = consumes< LHEEventProduct >(edm::InputTag("source","","LHEFile"));
  EDMLHEEventToken_alt  = consumes< LHEEventProduct >(edm::InputTag("externalLHEProducer","","LHE"));
  if(!isData){
    EDMGenParticlesToken    = consumes< std::vector<reco::GenParticle> >(edm::InputTag("prunedGenParticles","","PAT"));
    EDMGenJetsToken         = consumes< std::vector<reco::GenJet> >(edm::InputTag("slimmedGenJets","","PAT"));
    EDMCustomGenJetsToken   = consumes< std::vector<reco::GenJet> >(edm::InputTag("ak4GenJetsCustom","",""));

  
    // tt+X CATEGORIZATION data
    genBHadJetIndexToken           = consumes<std::vector<int> >(edm::InputTag("matchGenBHadron","genBHadJetIndex",""));
    genBHadFlavourToken            = consumes<std::vector<int> >(edm::InputTag("matchGenBHadron","genBHadFlavour",""));
    genBHadFromTopWeakDecayToken   = consumes<std::vector<int> >(edm::InputTag("matchGenBHadron","genBHadFromTopWeakDecay",""));
    genBHadPlusMothersToken        = consumes<std::vector<reco::GenParticle> >(edm::InputTag("matchGenBHadron","genBHadPlusMothers",""));
    genBHadPlusMothersIndicesToken = consumes<std::vector<std::vector<int> > >(edm::InputTag("matchGenBHadron","genBHadPlusMothersIndices",""));
    genBHadIndexToken              = consumes<std::vector<int> >(edm::InputTag("matchGenBHadron","genBHadIndex"));
    genBHadLeptonHadronIndexToken  = consumes<std::vector<int> >(edm::InputTag("matchGenBHadron","genBHadLeptonHadronIndex",""));
    genBHadLeptonViaTauToken       = consumes<std::vector<int> >(edm::InputTag("matchGenBHadron","genBHadLeptonViaTau",""));
    genCHadJetIndexToken           = consumes<std::vector<int> >(edm::InputTag("matchGenCHadron","genCHadJetIndex",""));
    genCHadFlavourToken            = consumes<std::vector<int> >(edm::InputTag("matchGenCHadron","genCHadFlavour",""));
    genCHadFromTopWeakDecayToken   = consumes<std::vector<int> >(edm::InputTag("matchGenCHadron","genCHadFromTopWeakDecay",""));
    genCHadBHadronIdToken          = consumes<std::vector<int> >(edm::InputTag("matchGenCHadron","genCHadBHadronId",""));
    genCHadIndexToken              = consumes<std::vector<int> >(edm::InputTag("matchGenCHadron","genCHadIndex"));
    genCHadPlusMothersToken        = consumes<std::vector<reco::GenParticle> >(edm::InputTag("matchGenCHadron","genCHadPlusMothers",""));
    
    genTtbarIdToken                = consumes<int>              (edm::InputTag("categorizeGenTtbar","genTtbarId",""));
  }
  // INITIALIZE MINIAOD HELPER
  helper.SetUp(era, sampleID, iAnalysisType, isData);
  helper.SetJetCorrectorUncertainty(); 

  // INITIALIZE SELECTION & CUTFLOW
  cutflow.Init((outfileName+"_Cutflow.txt").c_str());
  cutflow.AddStep("all");
  
  std::vector<std::string> selectionNames = iConfig.getParameter< std::vector<std::string> >("selectionNames");
  
  relevantTriggers = iConfig.getParameter< std::vector<std::string> >("relevantTriggers");

  // INITIALIZE TREEWRITER
  treewriter.Init(outfileName);
  std::vector<std::string> processorNames = iConfig.getParameter< std::vector<std::string> >("processorNames");
  for(vector<string>::const_iterator itPro = processorNames.begin();itPro != processorNames.end();++itPro) {
    cout << "This is Processor " << *itPro << endl;
    treewriter.FillProcessorName(*itPro);
    if(*itPro == "WeightProcessor") treewriter.AddTreeProcessor(new WeightProcessor());    
    else if(*itPro == "BaseVarProcessor") treewriter.AddTreeProcessor(new BaseVarProcessor());
    else if(*itPro == "RecoVarProcessor") treewriter.AddTreeProcessor(new RecoVarProcessor());
    else if(*itPro == "MVAVarProcessor") treewriter.AddTreeProcessor(new MVAVarProcessor());
    else if(*itPro == "TopGenVarProcessor") treewriter.AddTreeProcessor(new TopGenVarProcessor());
    else if(*itPro == "tHqGenVarProcessor") treewriter.AddTreeProcessor(new tHqGenVarProcessor());
    else cout << "No matching processor found for: " << *itPro << endl;    
    } 
  cout << "The total size of ProcessorNames is " << processorNames.size() << endl;
  treewriter.FillProcessorMap();
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
void tHqAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  cout << "########### NEW EVENT ##################" << endl;
  if(eventcount<10||eventcount%1000==0){
    cout << "Analyzing event " << eventcount << endl;
    watch.Print();
    watch.Continue();
  }
  
  eventcount++;
  
  /**** GET PILEUPSUMMARYINFO ****/
  edm::Handle< std::vector<PileupSummaryInfo> >  h_puinfosummary;
  iEvent.getByToken( EDMPUInfoToken, h_puinfosummary);
  
  /**** GET RHO ****/
  edm::Handle<double> h_rho;
  iEvent.getByToken(EDMRhoToken,h_rho);
  helper.SetRho(*h_rho);

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
  bool firstVertexIsGood=false;
  bool isFirst=true;  
  if( h_vtxs.isValid() ){
    for( reco::VertexCollection::const_iterator itvtx = vtxs.begin(); itvtx!=vtxs.end(); ++itvtx ){
      
      bool isGood = ( !(itvtx->isFake()) &&
		                  (itvtx->ndof() >= 4.0) &&
		                  (abs(itvtx->z()) <= 24.0) &&
		                  (abs(itvtx->position().Rho()) <= 2.0) 
		                );
      if( isGood ) selectedPVs.push_back(*itvtx);
      if( isFirst ) firstVertexIsGood=isGood;
      isFirst=false;
    }
  }
  if( vtxs.size()>0 ) helper.SetVertex( vtxs[0] );
  
  /**** GET LEPTONS ****/
  // MUONS
  edm::Handle< std::vector<pat::Muon> > h_muons;
  iEvent.getByToken( EDMMuonsToken,h_muons );
  std::vector<pat::Muon> const &muons = *h_muons; 
  std::vector<pat::Muon> rawMuons = muons;
  std::vector<pat::Muon> selectedMuons = helper.GetSelectedMuons( muons, 15., muonID::muonTight );
  std::vector<pat::Muon> selectedMuonsLoose = helper.GetSelectedMuons( muons, 10., muonID::muonLoose );

  // ELECTRONS`
  edm::Handle< std::vector<pat::Electron> > h_electrons;
  iEvent.getByToken( EDMElectronsToken,h_electrons );
  std::vector<pat::Electron> const &electrons = *h_electrons;
  std::vector<pat::Electron> rawElectrons = electrons;
  std::vector<pat::Electron> selectedElectrons = helper.GetSelectedElectrons( electrons, 15., electronID::electronTight );
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

  
  const JetCorrector* corrector = JetCorrector::getJetCorrector( "ak4PFchsL1L2L3", iSetup );   
  helper.SetJetCorrector(corrector);
  
  // selected jets with jet ID cuts
  std::vector<pat::Jet> idJets = helper.GetSelectedJets(pfjets, 0., 9999., jetID::jetLoose, '-' );
  // Get raw jets
  std::vector<pat::Jet> rawJets = helper.GetUncorrectedJets(idJets);
  // Clean muons from jets
  std::vector<pat::Jet> jetsNoMu = helper.RemoveOverlaps(selectedMuonsLoose, rawJets);
  // Clean electrons from jets
  std::vector<pat::Jet> jetsNoEle = helper.RemoveOverlaps(selectedElectronsLoose, jetsNoMu);
  // Apply jet corrections
  std::vector<pat::Jet> correctedJets = helper.GetCorrectedJets(jetsNoEle, iEvent, iSetup, sysType::NA);
  // Get jet Collection which pass selection
  std::vector<pat::Jet> selectedJets = helper.GetSelectedJets(correctedJets, 20., 4.7, jetID::jetLoose, '-' );
  // Get jet Collection which pass loose selection
  std::vector<pat::Jet> selectedJetsLoose = helper.GetSelectedJets(correctedJets, 15., 4.7, jetID::jetLoose, '-' );

  std::vector<pat::Jet> selectedJets_uncorrected = helper.GetSelectedJets(jetsNoEle, 15., 4.7, jetID::jetLoose, '-');



  // Get raw puppi jets
  std::vector<pat::Jet> rawPuppiJets = helper.GetUncorrectedJets(pfpuppijets);
  // Clean muons from jets
  std::vector<pat::Jet> puppiJetsNoMu = helper.RemoveOverlaps(selectedMuonsLoose, rawPuppiJets);
  // Clean electrons from jets
  std::vector<pat::Jet> puppiJetsNoEle = helper.RemoveOverlaps(selectedElectronsLoose, puppiJetsNoMu);
  // Apply jet corrections 
  std::vector<pat::Jet> correctedPuppiJets = helper.GetCorrectedJets(puppiJetsNoEle, iEvent, iSetup, sysType::NA);

  // Get puppi jet Collection which pass selection
  std::vector<pat::Jet> selectedPuppiJets = helper.GetSelectedJets(correctedPuppiJets, 20., 4.7, jetID::jetLoose, '-' );

  /**** GET MET ****/
  edm::Handle< std::vector<pat::MET> > h_pfmet;
  iEvent.getByToken( EDMMETsToken,h_pfmet );
  std::vector<pat::MET> const &pfMETs = *h_pfmet;
  // type I met corrections?
  assert(pfMETs.size()>0);

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
  if(!isData) iEvent.getByToken( EDMGenInfoToken, h_geneventinfo );

  /**** GET GENPARTICLES ****/
  
  edm::Handle< std::vector<reco::GenParticle> > h_genParticles;
  if(!isData){
    cout << "Ja was gehjt." << endl;
    iEvent.getByToken( EDMGenParticlesToken,h_genParticles );
    std::vector<reco::GenParticle> const &genParticles = *h_genParticles;
  }
  
  /**** GET GENJETS ****/
  edm::Handle< std::vector<reco::GenJet> > h_genjets;
  std::vector<reco::GenJet> selectedGenJets;
  if(!isData){
    iEvent.getByToken( EDMGenJetsToken,h_genjets );
    std::vector<reco::GenJet> const &genjets = *h_genjets;
    for(size_t i=0; i<genjets.size();i++){
      if(genjets[i].pt()>10&&fabs(genjets[i].eta())<5.)
	selectedGenJets.push_back(genjets[i]);
    }
  }
  
  // custom genjets for tt+X categorization
  edm::Handle< std::vector<reco::GenJet> > h_customgenjets;
  if(!isData){
    iEvent.getByToken( EDMCustomGenJetsToken,h_customgenjets );
  }
  
  // Fill Event Info Object
  EventInfo eventInfo(iEvent,h_beamspot,h_hcalnoisesummary,h_puinfosummary,firstVertexIsGood,*h_rho);
  
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

  //for(auto name=triggerResults.begin(); name!=triggerResults.end();name++){
    //    unsigned int TriggerID =  hlt_config.triggerIndex(*name);
  //   std::cout << "Jo: Name: " << *name << endl;
  // }

  TriggerInfo triggerInfo(triggerMap);

 
  // FIGURE OUT SAMPLE
  
  bool foundT=false;
  bool foundTbar=false;
  bool foundHiggs=false;
  HiggsDecay::HiggsDecay higgsdecay=HiggsDecay::NA;
  if(!isData){
    std::vector<reco::GenParticle> const &genParticles = *h_genParticles;
    for(size_t i=0; i<genParticles.size();i++){
      if(genParticles[i].pdgId()==6) foundT=true;
      if(genParticles[i].pdgId()==-6) foundTbar=true;
      if(genParticles[i].pdgId()==25){
	foundHiggs=true;
	if(higgsdecay==HiggsDecay::NA)higgsdecay=HiggsDecay::nonbb;
	for(uint j=0;j<genParticles[i].numberOfDaughters();j++){
	  if (abs(genParticles[i].daughter(j)->pdgId())==5){
	    higgsdecay=HiggsDecay::bb;
	  }
	}
      }
    }
  }
  GenTopEvent genTopEvt;
  GentHqEvent gentHqEvt;
  int ttid=-1;
  int ttid_full=-1;

  if(!isData&&useGenHadronMatch&&foundT&&foundTbar){

    std::cout << "Doing tt+X categorization" << endl;
    /**** tt+X categorization ****/
    // Reading gen jets from the event
    // Reading B hadrons related information
    edm::Handle<std::vector<int> > genBHadFlavour;
    edm::Handle<std::vector<int> > genBHadJetIndex;
    edm::Handle<std::vector<int> > genBHadFromTopWeakDecay;
    edm::Handle<std::vector<reco::GenParticle> > genBHadPlusMothers;
    edm::Handle<std::vector<std::vector<int> > > genBHadPlusMothersIndices;
    edm::Handle<std::vector<reco::GenParticle> > genCHadPlusMothers;
    edm::Handle<std::vector<int> > genBHadIndex;
    edm::Handle<std::vector<int> > genBHadLeptonHadronIndex;
    edm::Handle<std::vector<int> > genBHadLeptonViaTau;
    // Reading C hadrons related information
    edm::Handle<std::vector<int> > genCHadIndex;
    edm::Handle<std::vector<int> > genCHadFlavour;
    edm::Handle<std::vector<int> > genCHadJetIndex;
    edm::Handle<std::vector<int> > genCHadFromTopWeakDecay;
    edm::Handle<std::vector<int> > genCHadBHadronId;
    edm::Handle<int> genTtbarId;
    iEvent.getByToken(genCHadBHadronIdToken, genCHadBHadronId);
    iEvent.getByToken(genBHadFlavourToken, genBHadFlavour);
    iEvent.getByToken(genBHadJetIndexToken, genBHadJetIndex);  
    iEvent.getByToken(genBHadFromTopWeakDecayToken, genBHadFromTopWeakDecay);  
    iEvent.getByToken(genBHadPlusMothersToken, genBHadPlusMothers);    
    iEvent.getByToken(genBHadPlusMothersIndicesToken, genBHadPlusMothersIndices);
    iEvent.getByToken(genCHadPlusMothersToken, genCHadPlusMothers);    
    iEvent.getByToken(genBHadIndexToken, genBHadIndex);
    iEvent.getByToken(genBHadLeptonHadronIndexToken, genBHadLeptonHadronIndex);
    iEvent.getByToken(genBHadLeptonViaTauToken, genBHadLeptonViaTau);
    iEvent.getByToken(genCHadFlavourToken, genCHadFlavour);
    iEvent.getByToken(genCHadJetIndexToken, genCHadJetIndex);
    iEvent.getByToken(genCHadFromTopWeakDecayToken, genCHadFromTopWeakDecay);
    iEvent.getByToken(genCHadIndexToken, genCHadIndex);
    iEvent.getByToken(genTtbarIdToken, genTtbarId);
    ttid_full = *genTtbarId;
    ttid = ttid_full%100;
    // fill additional jet info in genTopEvt
    genTopEvt.FillTTxDetails(*h_customgenjets, 
			     *genBHadIndex, *genBHadJetIndex, *genBHadFlavour, *genBHadFromTopWeakDecay, *genBHadPlusMothers, 
			     *genCHadIndex, *genCHadJetIndex, *genCHadFlavour, *genCHadFromTopWeakDecay, *genCHadPlusMothers,
			     *genCHadBHadronId,
			     15,4.7); 
  }

  /*** FIGURE OUT SAMPLETYPE ***/
  SampleType sampleType= SampleType::nonttbkg;
  if(isData) sampleType = SampleType::data;
  else if(foundT&&foundTbar&&foundHiggs) sampleType = SampleType::tth;
  else if(foundT&&foundTbar){ 
    sampleType =SampleType::ttl;
    //if(ttid==51||ttid==52) sampleType = SampleType::ttb;
    if(ttid==51) sampleType = SampleType::ttb;
    else if(ttid==52) sampleType = SampleType::tt2b;
    else if(ttid==53||ttid==54||ttid==55) sampleType = SampleType::ttbb;
    else if(ttid==41||ttid==42) sampleType = SampleType::ttcc;
    else if(ttid==43||ttid==44||ttid==45) sampleType = SampleType::ttcc;    
  }
  else if(((foundT&&!foundTbar)||(!foundT&&foundTbar))&&foundHiggs) sampleType = SampleType::thq;

  /**** GET LHEINFO ****/
  edm::Handle<LHEEventProduct> h_lheeventinfo;
  if(!isData&&(sampleType==SampleType::thq)) iEvent.getByToken( EDMLHEEventToken, h_lheeventinfo );
  else iEvent.getByToken( EDMLHEEventToken_alt, h_lheeventinfo);

  
  /*** KICK OUT WRONG PROCESSORS ***/
  if(sampleType!=SampleType::thq) treewriter.RemoveTreeProcessor("tHqGenVarProcessor"); 
  if(sampleType==SampleType::thq || sampleType == SampleType::nonttbkg || sampleType==SampleType::data) treewriter.RemoveTreeProcessor("TopGenVarProcessor");
  if(!isData&&foundT&&foundTbar) {
    // fill genTopEvt with tt(H) information
    genTopEvt.Fill(*h_genParticles,ttid_full);
    cout << "APPARENTLY, THIS IS A SAMPLE WITH TWO TOP QUARKS" << endl;
  }
  
  if(sampleType ==  SampleType::thq){
    cout << "APPARENTLY, THIS IS A THQ SAMPLE" << endl;
    gentHqEvt.Fill(*h_genParticles);
  }

  // DO REWEIGHTING

  vector<string> syst_weights_id;
  vector<float> syst_weights;
  float Weight_orig=1;
  map<string,float> weights;
  map<string,float> weights_uncorrjets;
  if(!isData){
    weights = GetWeights(*h_geneventinfo,eventInfo,selectedPVs,selectedJets,selectedElectrons,selectedMuons,*h_genParticles,sysType::NA);
    weights_uncorrjets = GetWeights(*h_geneventinfo,eventInfo,selectedPVs,selectedJets_uncorrected,selectedElectrons,selectedMuons,*h_genParticles,sysType::NA);
  
 
    GetSystWeights(*h_lheeventinfo,syst_weights_id,syst_weights,Weight_orig);
  }
  /*
  std::cout << "Weight_orig: " << Weight_orig << std::endl;

  for (size_t i=0;i<syst_weights.size();i++){
    weights_syst[i]=syst_weights[i];
    weights_syst_id[i]=syst_weights_id[i];
  }

  std::cout << weights_syst[0] << std::endl;
  std::cout << weights_syst_id[0] << std::endl;

  */

  // Define INPUT
  InputCollections input( eventInfo,
			  //selectedTrigger,
			  //triggerResults,
                          triggerInfo,
			  selectedPVs,
			  rawMuons,
			  selectedMuons,
			  selectedMuonsLoose,
			  rawElectrons,
			  selectedElectrons,
                          selectedElectronsLoose,
			  rawJets,
                          selectedJets,
			  rawPuppiJets,
			  selectedPuppiJets,
                          selectedJetsLoose,
                          pfMETs[0],
			  genTopEvt,
			  gentHqEvt,
                          selectedGenJets,
                          sampleType,
                          weights,
			  Weight_orig,
			  syst_weights,
			  syst_weights_id
			  );
       
  

  // DO SELECTION
  cutflow.EventSurvivedStep("all");
  bool selected=true;
  for(size_t i=0; i<selections.size() && selected; i++){
    if(!selections.at(i)->IsSelected(input,cutflow)) selected=false;
  }

  if(!selected) return;    

  selectedElectrons = ElectronSelection(selectedElectrons,input.selectedPVs);
  selectedMuons = MuonSelection(selectedMuons,input.selectedPVs);
  selectedJets = JetSelection(selectedJets,selectedElectrons, selectedMuons, input);
  
  // WRITE TREE
  treewriter.Process(input);  
}


std::vector<pat::Muon> tHqAnalyzer::MuonSelection( std::vector<pat::Muon> selectedMuons,  const reco::VertexCollection& selectedPVs){
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
    
    const float minPt=5;
    const float maxEta=2.4;
    
    
    
    if(iMuon.globalTrack().isAvailable()){
      double chi2ndof = iMuon.globalTrack()->normalizedChi2();
      int nValidHits = iMuon.globalTrack()->hitPattern().numberOfValidMuonHits();
      passesGlobalTrackID=(chi2ndof<10.0 && nValidHits>0);
    } 
    
    if(iMuon.muonBestTrack().isAvailable()){
      double d0 = fabs(iMuon.muonBestTrack()->dxy(selectedPVs.at(0).position()));
      double dZ = fabs(iMuon.muonBestTrack()->dz(selectedPVs.at(0).position()));
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
    passesKinematics = (iMuon.pt()>minPt && fabs(iMuon.eta())<maxEta);
    passesISO = (relIso<0.12);
    passesID = (passesinnerTrackID && passesGlobalTrackID && passesBestTrackID && passesTrackID && nMatchedStations>1);
    
    if(passesKinematics && passesISO && passesID && isPFandGlobal){
    newselectedMuons.push_back(*it);
    }
    else std::cout << " Bad Muon found. " << std::endl;
  }
  return newselectedMuons;
}


  std::vector<pat::Electron> tHqAnalyzer::ElectronSelection( std::vector<pat::Electron> selectedElectrons,  const reco::VertexCollection& selectedPVs){
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

    const float minPt=5;
    const float maxEta=2.4;


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
      d0 = fabs(iElectron.gsfTrack()->dxy(selectedPVs.at(0).position()));
      dZ = fabs(iElectron.gsfTrack()->dz(selectedPVs.at(0).position()));
      expectedMissingInnerHits = iElectron.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
    }

    //do the checks                                                                                                                                        
    passesConversion = ( iElectron.passConversionVeto() );
    passesKinematics = (iElectron.pt()>minPt && fabs(iElectron.eta())<maxEta);
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

  const float minPt=10;
  const float maxEta=5;

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
    if(iJet.pt()>=minPt && abs(iJet.eta())<maxEta){
      if(abs(iJet.eta())<=3.0) isselected = (iJet.neutralHadronEnergyFraction()<0.99 && iJet.neutralEmEnergyFraction()<0.99 && (iJet.chargedMultiplicity()+iJet.neutralMultiplicity())>1) && ((abs(iJet.eta())<=2.4 && iJet.chargedHadronEnergyFraction()>0 && iJet.chargedMultiplicity()>0 && iJet.chargedEmEnergyFraction()<0.99) || abs(iJet.eta())>2.4) && abs(iJet.eta())<=3.0;
      //tightJetID = (iJet.neutralHadronEnergyFraction()<0.90 && iJet.neutralEmEnergyFraction()<0.90 && (iJet.chargedMultiplicity()+iJet.neutralMultiplicity())>1) && ((abs(iJet.eta())<=2.4 && iJet.chargedHadronEnergyFraction()>0 && iJet.chargedMultiplicity()>0 && iJet.chargedEmEnergyFraction()<0.99) ||  abs(iJet.eta())>2.4) && abs(iJet.eta())<=3.0;
      else isselected = (iJet.neutralEmEnergyFraction()<0.90 && iJet.neutralMultiplicity()>10 && abs(iJet.eta())>3.0 );
    }
    // isselected=(iJet.pt()>minPt && abs(iJet.eta())<maxEta && iJet.neutralHadronEnergyFraction()<0.99 && iJet.chargedHadronEnergyFraction()>0.0 && iJet.chargedMultiplicity()>0.0 && iJet.chargedEmEnergyFraction()<0.99 && iJet.neutralEmEnergyFraction()<0.99 && iJet.numberOfDaughters()>1 );
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

map<string,float> tHqAnalyzer::GetWeights(const GenEventInfoProduct&  genEventInfo,const EventInfo& eventInfo, const reco::VertexCollection& selectedPVs, const std::vector<pat::Jet>& selectedJets, const std::vector<pat::Electron>& selectedElectrons, const std::vector<pat::Muon>& selectedMuons, const std::vector<reco::GenParticle>& genParticles, sysType::sysType systype){//, edm::Run const& iRun){
  
  map<string,float> weights;
  
  if(isData){
    weights["Weight"] = 1.0;
    weights["Weight_XS"] = 1.0;
    weights["Weight_CSV"] = 1.0;
    weights["Weight_PU"] = 1.0;
    weights["Weight_PV"] = 1.0;
    weights["Weight_TopPt"] = 1.0;
    return weights;
  }

  float weight = 1.;
  assert(genEventInfo.weights().size()<=1); // before we multiply any weights we should understand what they mean
  for(size_t i=0;i<genEventInfo.weights().size();i++){
    
    weight *= (genEventInfo.weights()[i]>0 ? 1.: -1.); // overwrite intransparent MC weights, use \pm 1 instead
  }   // NOTE: IS THIS EVEN CORRECT?

  weight = genEventInfo.weight();

  //  double csvWgtHF, csvWgtLF, csvWgtCF;

  float xsweight = xs*luminosity/totalMCevents;
  float csvweight = 1.;
  float puweight = 1.;
  float topptweight = 1.;
  
  // ADD CSV WEIGHTS HERE OR SEPARATLY?


  //weight *= xsweight*csvweight*puweight*topptweight;
  weights["Weight"] = weight;
  weights["Weight_XS"] = xsweight;
  weights["Weight_CSV"] = csvweight;
  weights["Weight_PU"] = puweight;
  weights["Weight_TopPt"] = topptweight;
  //  weights["Weight_PV"] = pvWeight.GetWeight(selectedPVs.size());

  
  
  
  return weights;
}



// Systematic weights
void tHqAnalyzer::GetSystWeights(const LHEEventProduct&  LHEEvent, vector<string> &weight_syst_id, vector<float> &weight_syst, float &Weight_orig ){

  //  map<string,float> syst_weights;
  
  //  syst_weights["Weight_orig"] = LHEEvent.originalXWGTUP();

  

  Weight_orig = LHEEvent.originalXWGTUP();

  //cout << "lalala " << LHEEvent.weights().size() <<endl;
  //cout << "lululu " << LHEEvent.weights()[491].id.c_str() <<endl;
  cout << "original weight:" << LHEEvent.originalXWGTUP() << endl;
  cout << "Weight_orig2:" << Weight_orig << endl;

  for (size_t i=0; i<LHEEvent.weights().size();i++){ 
    weight_syst_id.push_back(LHEEvent.weights()[i].id.c_str());
    weight_syst.push_back(LHEEvent.weights()[i].wgt);
  }
  //  for (size_t i=0;i<LHEEvent.weights().size();i++){
  // cout << "weight id i: " << LHEEvent.weights()[i].c_str() << endl;
  // }

  //  return syst_weights;
  
  return;

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
  std::cout << "Original File had " << eventcount << " Entries." << std::endl;
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
