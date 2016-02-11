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
#include "FWCore/Framework/interface/TriggerNamesService.h"
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
#include "MiniAOD/MiniAODHelper/interface/CSVHelper.h"

#include "tHqAnalysis/tHqAnalyzer/interface/tHqUtils.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/InputCollections.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/Cutflow.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/EventInfo.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/TreeWriter.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/TriggerInfo.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/HistoReweighter.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/PUWeights.hpp"

#include "tHqAnalysis/tHqAnalyzer/interface/Selection.hpp"

#include "tHqAnalysis/tHqAnalyzer/interface/WeightProcessor.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/TopGenVarProcessor.hpp"

#include "tHqAnalysis/tHqAnalyzer/interface/MVAVarProcessor.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/BaseVarProcessor.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/RecoVarProcessor.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/GenTopEvent.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/GentHqEvent.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/tHqGenVarProcessor.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/TriggerVarProcessor.hpp"

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
  

  map<string,float> GetWeights(const GenEventInfoProduct& genEventInfo, 
			       const EventInfo& eventInfo, 
			       const reco::VertexCollection& selectedPVs, 
			       const std::vector<pat::Jet>& selectedJets, 
			       const std::vector<pat::Electron>& selectedElectrons, 
			       const std::vector<pat::Muon>& selectedMuons, 
			       const std::vector<reco::GenParticle>& genParticles, 
			       const GenTopEvent& genTopEvt,
			       const sysType::sysType systype=sysType::NA, 
			       const SampleType sampletype=SampleType::nonttbkg
			       );

  void GetSystWeights(const LHEEventProduct& LHEEvent, vector<string> &weight_syst_id, vector<float> &weight_syst, float &Weight_orig);
  float GetTopPtWeight(float toppt1, float toppt2);
      
      
  // ----------member data ---------------------------
  
  /** the beanhelper is used for selections and reweighting */
  MiniAODHelper helper;
  
  CSVHelper csvReweighter;
  
  /** Reweighter to match the PV distribution in data*/
  HistoReweighter pvWeight;
  PUWeights puWeights_;
  
  /** writes flat trees for MVA analysis */
  TreeWriter treewriter_nominal;
  
  
  /** Declare TreeWriters for JES and JER systematics */
  TreeWriter treewriter_jesup;
  TreeWriter treewriter_jesdown;
  TreeWriter treewriter_jerup;
  TreeWriter treewriter_jerdown;
  
  /** store output file names */

  std::string outfileName;
  std::string outfileName_nominal;
  std::string outfileNameJESup;
  std::string outfileNameJESdown;
  std::string outfileNameJERup;
  std::string outfileNameJERdown;

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

  /** is LHE info available */
  bool useLHE;

  /** triggers that are checked */
  vector<std::string> relevantTriggers;

  /** calculate and store systematic weights? */
  bool doSystematics = 1;

  /** calculate JER systematics **/
  bool doJERsystematic =1;
  
  /** use GenBmatching info? this is only possible if the miniAOD contains them */
  bool useGenHadronMatch;
    
  /** recorrect the jets and MET that were in MiniAOD? */
  bool recorrectMET;

  /** jet systematic that is applied (the outher systematics are done at a different place with reweighting)*/
  sysType::sysType jsystype;
      
  /** SampleType needed to see what gen variables need to be filled */
  SampleType sampleType = SampleType::nonttbkg;
      
  /** pu summary data access token **/
  edm::EDGetTokenT< std::vector<PileupSummaryInfo> > EDMPUInfoToken;
      
  /** pileup density data access token **/
  edm::EDGetTokenT <double> EDMRhoToken;


  /** hcal noise data access token **/
  edm::EDGetTokenT< HcalNoiseSummary > EDMHcalNoiseToken;
      


  edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesToken;
  /** selected trigger data access token **/
  //  edm::EDGetTokenT< pat::TriggerObjectStandAloneCollection > EDMSelectedTriggerToken;
      
  /** trigger results data access token **/
  //  edm::EDGetTokenT< edm::TriggerResults > EDMTriggerResultToken;
  // HLTConfigProvider hlt_config;

  /** beam spot data access token **/
  edm::EDGetTokenT< reco::BeamSpot > EDMBeamSpotToken;
      
  /** vertex data access token **/
  edm::EDGetTokenT< reco::VertexCollection > EDMVertexToken;
      
  /** muons data access token **/
  edm::EDGetTokenT< std::vector<pat::Muon> > EDMMuonsToken;
      
  /** electrons data access token **/
  edm::EDGetTokenT< edm::View<pat::Electron> > EDMElectronsToken;
      
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
        

  // MVA values and categories
  edm::EDGetTokenT<edm::ValueMap<float> > EDMeleMVAvaluesToken;
  edm::EDGetTokenT<edm::ValueMap<int> >   EDMeleMVAcategoriesToken;


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
tHqAnalyzer::tHqAnalyzer(const edm::ParameterSet& iConfig):csvReweighter(CSVHelper("MiniAOD/MiniAODHelper/data/csv_rwt_fit_hf_2015_11_20.root","MiniAOD/MiniAODHelper/data/csv_rwt_fit_lf_2015_11_20.root",5)),pvWeight((tHqUtils::GetAnalyzerPath()+"/data/pvweights/data.root").c_str(),"data",(tHqUtils::GetAnalyzerPath()+"/data/pvweights/mc.root").c_str(),"mc"){
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
  useLHE = iConfig.getParameter<bool>("useLHE");
  useGenHadronMatch = iConfig.getParameter<bool>("useGenHadronMatch");
  recorrectMET = iConfig.getParameter<bool>("recorrectMET");

  if(isData) doSystematics = false;

  
  //  useFatJets = iConfig.getParameter<bool>("useFatJets");
  outfileName = iConfig.getParameter<std::string>("outfileName");
  outfileName_nominal=outfileName;
  treewriter_nominal.SetTreeName("Nominal");
  
  size_t stringIndex = outfileName.find("nominal");
  if(stringIndex==std::string::npos) outfileName_nominal=outfileName+"_nominal";

  if(doSystematics){
    treewriter_jesup.SetTreeName("JESup");
    treewriter_jesdown.SetTreeName("JESdown");
    treewriter_jerup.SetTreeName("JERup");
    treewriter_jerdown.SetTreeName("JERdown");

    outfileNameJESup=outfileName_nominal;
    outfileNameJESdown=outfileName_nominal;
    outfileNameJERup=outfileName_nominal;
    outfileNameJERdown=outfileName_nominal;
  
    if(stringIndex!=std::string::npos){
      outfileNameJESup.replace(stringIndex,7,"JESUP");
      outfileNameJESdown.replace(stringIndex,7,"JESDOWN");
      outfileNameJERup.replace(stringIndex,7,"JERUP");
      outfileNameJERdown.replace(stringIndex,7,"JERDOWN");
    }
    else{
      outfileNameJESup=outfileName+"_JESUP";
      outfileNameJESdown=outfileName+"_JESDOWN";
      outfileNameJERup=outfileName+"_JERUP";
      outfileNameJERdown=outfileName+"_JERDOWN";
    }
  }
  std::cout << "Outfile Name: " << outfileName_nominal << std::endl;
  
  // REGISTER DATA ACCESS
  EDMPUInfoToken          = consumes< std::vector<PileupSummaryInfo> >(edm::InputTag("slimmedAddPileupInfo","",""));
  EDMRhoToken             = consumes <double> (edm::InputTag(std::string("fixedGridRhoFastjetAll")));
  EDMHcalNoiseToken       = consumes< HcalNoiseSummary >(edm::InputTag("hcalnoise","",""));
  //  EDMSelectedTriggerToken = consumes< pat::TriggerObjectStandAloneCollection > (edm::InputTag("selectedPatTrigger","",""));
  // EDMTriggerResultToken   = consumes< edm::TriggerResults > (edm::InputTag("TriggerResults","","HLT"));

  triggerBitsToken        = consumes< edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"));
  triggerObjectsToken     = consumes< pat::TriggerObjectStandAloneCollection>(edm::InputTag("selectedPatTrigger","",""));
  triggerPrescalesToken   = consumes< pat::PackedTriggerPrescales>(edm::InputTag("patTrigger","",""));

  EDMBeamSpotToken        = consumes< reco::BeamSpot > (edm::InputTag("offlineBeamSpot","",""));
  EDMVertexToken          = consumes< reco::VertexCollection > (edm::InputTag("offlineSlimmedPrimaryVertices"));
  EDMMuonsToken           = consumes< std::vector<pat::Muon> >(edm::InputTag("slimmedMuons"));
  EDMElectronsToken       = consumes< edm::View<pat::Electron> >(edm::InputTag("slimmedElectrons"));
  EDMJetsToken            = consumes< std::vector<pat::Jet> >(edm::InputTag("slimmedJets"));
  EDMPuppiJetsToken       = consumes< std::vector<pat::Jet> >(edm::InputTag("slimmedJetsPuppi"));
  EDMMETsToken            = consumes< std::vector<pat::MET> >(edm::InputTag("slimmedMETs"));
  EDMGenInfoToken         = consumes< GenEventInfoProduct >(edm::InputTag("generator"));

  // electron MVA info
  // TODO: these (and many of the names above) shouldn't be hard coded but set in python cfg
  EDMeleMVAvaluesToken           = consumes<edm::ValueMap<float> >(edm::InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Values",""));
  EDMeleMVAcategoriesToken       = consumes<edm::ValueMap<int> >(edm::InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories",""));

  if(useLHE){
    EDMLHEEventToken      = consumes< LHEEventProduct >(edm::InputTag("source"));
    EDMLHEEventToken_alt  = consumes< LHEEventProduct >(edm::InputTag("externalLHEProducer"));
  }
  if(!isData){
    EDMGenParticlesToken    = consumes< std::vector<reco::GenParticle> >(edm::InputTag("prunedGenParticles"));
    EDMGenJetsToken         = consumes< std::vector<reco::GenJet> >(edm::InputTag("slimmedGenJets"));
    EDMCustomGenJetsToken   = consumes< std::vector<reco::GenJet> >(edm::InputTag("ak4GenJetsCustom"));

  
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

  // INITIALIZE PU WEIGHTS
  puWeights_.init(iConfig);

  // INITIALIZE SELECTION & CUTFLOW
  cutflow.Init((outfileName_nominal+"_Cutflow.txt").c_str());
  cutflow.AddStep("all");
  
  std::vector<std::string> selectionNames = iConfig.getParameter< std::vector<std::string> >("selectionNames");
  
  relevantTriggers = iConfig.getParameter< std::vector<std::string> >("relevantTriggers");

  // INITIALIZE TREEWRITER
  treewriter_nominal.Init(outfileName_nominal);
  // in case of systematics
  if(doSystematics){
    treewriter_jesup.Init(outfileNameJESup);
    treewriter_jesdown.Init(outfileNameJESdown);
    treewriter_jerup.Init(outfileNameJERup);
    treewriter_jerdown.Init(outfileNameJERdown);
  }
  

  std::vector<std::string> processorNames = iConfig.getParameter< std::vector<std::string> >("processorNames");
  for(vector<string>::const_iterator itPro = processorNames.begin();itPro != processorNames.end();++itPro) {
    treewriter_nominal.FillProcessorName(*itPro);
    if(*itPro == "WeightProcessor") treewriter_nominal.AddTreeProcessor(new WeightProcessor());    
    else if(*itPro == "BaseVarProcessor") treewriter_nominal.AddTreeProcessor(new BaseVarProcessor());
    else if(*itPro == "RecoVarProcessor") treewriter_nominal.AddTreeProcessor(new RecoVarProcessor());
    else if(*itPro == "TopGenVarProcessor") treewriter_nominal.AddTreeProcessor(new TopGenVarProcessor());
    else if(*itPro == "tHqGenVarProcessor") treewriter_nominal.AddTreeProcessor(new tHqGenVarProcessor());
    else if(*itPro == "TriggerVarProcessor") treewriter_nominal.AddTreeProcessor(new TriggerVarProcessor(relevantTriggers));
    else cout << "No matching processor found for: " << *itPro << endl;    
  } 
  treewriter_nominal.FillProcessorMap();

  //the systematics tree writers use the same processors that are used for the nominal trees
  // it might improve the performance to turn some of them off
  if(doSystematics){
    std::vector<TreeProcessor*> tps = treewriter_nominal.GetTreeProcessors();
    std::vector<string> tpsn = treewriter_nominal.GetTreeProcessorNames();
    for(uint i=0; i<tps.size();i++){
      treewriter_jesup.AddTreeProcessor(tps[i],tpsn[i]);
      treewriter_jesdown.AddTreeProcessor(tps[i],tpsn[i]);
      treewriter_jerup.AddTreeProcessor(tps[i],tpsn[i]);
      treewriter_jerdown.AddTreeProcessor(tps[i],tpsn[i]);
    }
    treewriter_jesup.FillProcessorMap();
    treewriter_jesdown.FillProcessorMap();
    treewriter_jerup.FillProcessorMap();
    treewriter_jerdown.FillProcessorMap();
    
  }
  //}
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
  if(eventcount<10||eventcount%1000==0){
    cout << "----------------------------------------" << endl;
    cout << "-           Analyzing event " << eventcount << "          -" << endl;
    cout << "----------------------------------------" << endl;
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
  
  //  edm::Handle<edm::TriggerResults> h_triggerresults;
  //iEvent.getByToken( EDMTriggerResultToken,h_triggerresults );
  //edm::TriggerResults const &triggerResults = *h_triggerresults;
  
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

  TriggerInfo triggerInfo(iEvent,triggerBitsToken,triggerObjectsToken,triggerPrescalesToken);
  if(eventcount==1) triggerInfo.ListTriggers();
  
  /**** GET LEPTONS ****/
  // MUONS
  edm::Handle< std::vector<pat::Muon> > h_muons;
  iEvent.getByToken( EDMMuonsToken,h_muons );
  std::vector<pat::Muon> const &muons = *h_muons; 
  std::vector<pat::Muon> rawMuons = muons;
  std::vector<pat::Muon> selectedMuons = helper.GetSelectedMuons( muons, 15., muonID::muonTight );
  std::vector<pat::Muon> selectedMuonsLoose = helper.GetSelectedMuons( muons, 10., muonID::muonLoose );

  // ELECTRONS

  edm::Handle< edm::View<pat::Electron> > h_electrons;
  iEvent.getByToken( EDMElectronsToken,h_electrons );

  // add electron mva info to electrons
  edm::Handle<edm::ValueMap<float> > h_mvaValues; 
  edm::Handle<edm::ValueMap<int> > h_mvaCategories;
  iEvent.getByToken(EDMeleMVAvaluesToken,h_mvaValues);
  iEvent.getByToken(EDMeleMVAcategoriesToken,h_mvaCategories);  
  std::vector<pat::Electron> electrons = helper.GetElectronsWithMVAid(h_electrons,h_mvaValues,h_mvaCategories);

  helper.AddElectronRelIso(electrons,coneSize::R03, corrType::rhoEA,effAreaType::spring15,"relIso");
  std::vector<pat::Electron> rawElectrons = electrons;
  //  std::vector<pat::Electron> selectedElectrons = helper.GetSelectedElectrons( electrons, 15., electronID::electronEndOf15NonTrigMVA80iso0p15 );
  //  std::vector<pat::Electron> selectedElectronsLoose = helper.GetSelectedElectrons( electrons, 10., electronID::electronEndOf15NonTrigMVA80iso0p15 );

  std::vector<pat::Electron> selectedElectrons = helper.GetSelectedElectrons( electrons, 20., electronID::electronEndOf15MVA80iso0p15, 2.1 );
  std::vector<pat::Electron> selectedElectronsLoose = helper.GetSelectedElectrons( electrons, 15., electronID::electronEndOf15MVA80iso0p15, 2.4 );


  /**** GET MET ****/
  edm::Handle< std::vector<pat::MET> > h_pfmet;
  iEvent.getByToken( EDMMETsToken,h_pfmet );
  std::vector<pat::MET> const &pfMETs = *h_pfmet;

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
  
  const double jetptcut=20.0; 
  const double jetetacut=4.7; 
  const double jetptcut_loose=15.0; 
  const double jetetacut_loose=4.7;

  // selected jets with jet ID cuts
  std::vector<pat::Jet> idJets = helper.GetSelectedJets(pfjets, 0., 9999., jetID::jetLoose, '-' );
  // Get raw jets
  std::vector<pat::Jet> rawJets = helper.GetUncorrectedJets(idJets);

  // Clean muons and electrons from jets
  std::vector<pat::Jet> cleanJets = helper.GetDeltaRCleanedJets(rawJets,selectedMuonsLoose,selectedElectronsLoose,0.4);
  // Apply nominal jet corrections
  std::vector<pat::Jet> correctedJets_nominal = helper.GetCorrectedJets(cleanJets, iEvent, iSetup, sysType::NA);
  // Get jet Collection which pass selection
  std::vector<pat::Jet> selectedJets_nominal = helper.GetSelectedJets(correctedJets_nominal, jetptcut, jetetacut, jetID::none, '-' );
  // Sort jets
  selectedJets_nominal = helper.GetSortedByPt(selectedJets_nominal);

  // Get jet Collection which pass loose selection
  std::vector<pat::Jet> selectedJetsLoose_nominal = helper.GetSelectedJets(correctedJets_nominal, jetptcut_loose ,jetetacut_loose, jetID::none, '-' );
  // Sort jets
  selectedJetsLoose_nominal = helper.GetSortedByPt(selectedJetsLoose_nominal);

  //Get selected jets for MET filters
  std::vector<pat::Jet> idJetsForMET = helper.GetSelectedJets(pfjets, 0., 5.4, jetID::jetMETcorrection, '-' );
  std::vector<pat::Jet> rawJetsForMET = helper.GetUncorrectedJets(idJetsForMET);
  std::vector<pat::Jet> correctedJetsForMET_nominal = helper.GetCorrectedJets(rawJetsForMET, iEvent, iSetup, sysType::NA);
  //correct MET 
  std::vector<pat::MET> correctedMETs_nominal;
  if(recorrectMET) correctedMETs_nominal = helper.CorrectMET(idJetsForMET,correctedJetsForMET_nominal,pfMETs);
  else correctedMETs_nominal = pfMETs;

  // Get raw puppi jets
  std::vector<pat::Jet> rawPuppiJets = helper.GetUncorrectedJets(pfpuppijets);
  // DeltaR clean puppi jets
  std::vector<pat::Jet> cleanedPuppiJets = helper.GetDeltaRCleanedJets(rawPuppiJets,selectedMuonsLoose, selectedElectronsLoose, 0.4);
  // Apply jet corrections 
  std::vector<pat::Jet> correctedPuppiJets = helper.GetCorrectedJets(cleanedPuppiJets, iEvent, iSetup, sysType::NA);

  // Get puppi jet Collection which pass selection
  std::vector<pat::Jet> selectedPuppiJets = helper.GetSelectedJets(correctedPuppiJets, jetptcut, jetetacut, jetID::none, '-' );

  /** JES Correcteions for Jets */
  // Apply systematically shifted jet corrections -- these vector stay empty if you dont use makeSystematicsTrees
  std::vector<pat::Jet> correctedJets_unsorted_jesup;
  std::vector<pat::Jet> correctedJets_unsorted_jesdown;

  std::vector<pat::Jet> correctedJetsForMET_jesup;
  std::vector<pat::Jet> correctedJetsForMET_jesdown;
  std::vector<pat::MET> correctedMETs_jesup;
  std::vector<pat::MET> correctedMETs_jesdown;

  std::vector<pat::Jet> selectedJets_unsorted_jesup;
  std::vector<pat::Jet> selectedJets_unsorted_jesdown;
  std::vector<pat::Jet> selectedJets_unsorted_uncorrected;
  std::vector<pat::Jet> selectedJetsLoose_unsorted_jesup;
  std::vector<pat::Jet> selectedJetsLoose_unsorted_jesdown;
  std::vector<pat::Jet> selectedJetsLoose_unsorted_uncorrected;

  std::vector<pat::Jet> selectedJets_jesup;
  std::vector<pat::Jet> selectedJets_jesdown;
  std::vector<pat::Jet> selectedJets_uncorrected;
  std::vector<pat::Jet> selectedJetsLoose_jesup;
  std::vector<pat::Jet> selectedJetsLoose_jesdown;
  std::vector<pat::Jet> selectedJetsLoose_uncorrected;
  
  if(doSystematics){
    correctedJets_unsorted_jesup = helper.GetCorrectedJets(cleanJets, iEvent, iSetup, sysType::JESup);
    correctedJets_unsorted_jesdown = helper.GetCorrectedJets(cleanJets, iEvent, iSetup, sysType::JESdown);
 
    correctedJetsForMET_jesup = helper.GetCorrectedJets(rawJetsForMET, iEvent, iSetup, sysType::JESup);
    correctedJetsForMET_jesdown = helper.GetCorrectedJets(rawJetsForMET, iEvent, iSetup, sysType::JESdown);
    correctedMETs_jesup = recorrectMET ? helper.CorrectMET(idJetsForMET,correctedJetsForMET_jesup,pfMETs) : pfMETs;
    correctedMETs_jesdown = recorrectMET ? helper.CorrectMET(idJetsForMET,correctedJetsForMET_jesdown,pfMETs) : pfMETs;

    selectedJets_unsorted_jesup = helper.GetSelectedJets(correctedJets_unsorted_jesup, jetptcut, jetetacut, jetID::none, '-' );
    selectedJets_unsorted_jesdown = helper.GetSelectedJets(correctedJets_unsorted_jesdown, jetptcut, jetetacut, jetID::none, '-' );
    selectedJets_unsorted_uncorrected = helper.GetSelectedJets(cleanJets, jetptcut, jetetacut, jetID::none, '-' );
    selectedJetsLoose_unsorted_jesup = helper.GetSelectedJets(correctedJets_unsorted_jesup, jetptcut_loose, jetetacut_loose, jetID::none, '-' ); 
    selectedJetsLoose_unsorted_jesdown = helper.GetSelectedJets(correctedJets_unsorted_jesdown, jetptcut_loose,jetetacut_loose, jetID::none, '-' ); 
    selectedJetsLoose_unsorted_uncorrected = helper.GetSelectedJets(cleanJets, jetptcut_loose, jetetacut_loose, jetID::none, '-' ); 

    selectedJets_jesup = helper.GetSortedByPt(selectedJets_unsorted_jesup);
    selectedJets_jesdown = helper.GetSortedByPt(selectedJets_unsorted_jesdown);
    selectedJets_uncorrected = helper.GetSortedByPt(selectedJets_unsorted_uncorrected);
    selectedJetsLoose_jesup = helper.GetSortedByPt(selectedJetsLoose_unsorted_jesup);
    selectedJetsLoose_jesdown = helper.GetSortedByPt(selectedJetsLoose_unsorted_jesdown);
    selectedJetsLoose_uncorrected = helper.GetSortedByPt(selectedJetsLoose_unsorted_uncorrected);
  }


  // Apply systematically shifted jet resolution -- these vector stay empty if you dont use makeSystematicsTrees
  std::vector<pat::Jet> correctedJets_unsorted_jerup;
  std::vector<pat::Jet> correctedJets_unsorted_jerdown;

  std::vector<pat::Jet> correctedJetsForMET_jerup;
  std::vector<pat::Jet> correctedJetsForMET_jerdown;
  std::vector<pat::MET> correctedMETs_jerup;
  std::vector<pat::MET> correctedMETs_jerdown;

  std::vector<pat::Jet> selectedJets_unsorted_jerup;
  std::vector<pat::Jet> selectedJets_unsorted_jerdown;
  std::vector<pat::Jet> selectedJetsLoose_unsorted_jerup;
  std::vector<pat::Jet> selectedJetsLoose_unsorted_jerdown;

  std::vector<pat::Jet> selectedJets_jerup;
  std::vector<pat::Jet> selectedJets_jerdown;
  std::vector<pat::Jet> selectedJetsLoose_jerup;
  std::vector<pat::Jet> selectedJetsLoose_jerdown;
  
  if(doSystematics){
    correctedJets_unsorted_jerup = helper.GetCorrectedJets(cleanJets, iEvent, iSetup, sysType::JERup);
    correctedJets_unsorted_jerdown = helper.GetCorrectedJets(cleanJets, iEvent, iSetup, sysType::JERdown);

    correctedJetsForMET_jerup = helper.GetCorrectedJets(rawJetsForMET, iEvent, iSetup, sysType::JERup);
    correctedJetsForMET_jerdown = helper.GetCorrectedJets(rawJetsForMET, iEvent, iSetup, sysType::JERdown);
    correctedMETs_jerup = recorrectMET ? helper.CorrectMET(idJetsForMET,correctedJetsForMET_jerup,pfMETs) : pfMETs;
    correctedMETs_jerdown = recorrectMET ? helper.CorrectMET(idJetsForMET,correctedJetsForMET_jerdown,pfMETs) : pfMETs;

    selectedJets_unsorted_jerup = helper.GetSelectedJets(correctedJets_unsorted_jerup, jetptcut, jetetacut, jetID::none, '-' );
    selectedJets_unsorted_jerdown = helper.GetSelectedJets(correctedJets_unsorted_jerdown, jetptcut, jetetacut, jetID::none, '-' );
    selectedJetsLoose_unsorted_jerup = helper.GetSelectedJets(correctedJets_unsorted_jerup, jetptcut_loose, jetetacut_loose, jetID::none, '-' ); 
    selectedJetsLoose_unsorted_jerdown = helper.GetSelectedJets(correctedJets_unsorted_jerdown, jetptcut_loose,jetetacut_loose, jetID::none, '-' ); 

    selectedJets_jerup = helper.GetSortedByPt(selectedJets_unsorted_jerup);
    selectedJets_jerdown = helper.GetSortedByPt(selectedJets_unsorted_jerdown);
    selectedJetsLoose_jerup = helper.GetSortedByPt(selectedJetsLoose_unsorted_jerup);
    selectedJetsLoose_jerdown = helper.GetSortedByPt(selectedJetsLoose_unsorted_jerdown);
    
  }



  
  /**** GET GENEVENTINFO ****/
  edm::Handle<GenEventInfoProduct> h_geneventinfo;
  if(!isData) iEvent.getByToken( EDMGenInfoToken, h_geneventinfo );

  /**** GET GENPARTICLES ****/
  
  edm::Handle< std::vector<reco::GenParticle> > h_genParticles;
  if(!isData){
    iEvent.getByToken( EDMGenParticlesToken,h_genParticles );
    //    std::vector<reco::GenParticle> const &genParticles = *h_genParticles;
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
  if(eventcount==1){
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
    else if (((foundT&&!foundTbar)||(!foundT&&foundTbar))&&foundHiggs){
      sampleType = SampleType::thq;
      cout << "### SampleType ---> SingleTop + Higgs " << endl;
    } 
    else if (((foundT&&!foundTbar)||(!foundT&&foundTbar))&&!foundHiggs){
      sampleType = SampleType::st;
      cout << "### SampleType ---> SingleTop " << endl;
    }
  }

  /**** GET LHEINFO ****/
  edm::Handle<LHEEventProduct> h_lheeventinfo;
  if(useLHE){
    if(!isData&&(sampleType==SampleType::thq||sampleType==SampleType::st)) iEvent.getByToken( EDMLHEEventToken, h_lheeventinfo );
    else iEvent.getByToken( EDMLHEEventToken_alt, h_lheeventinfo);
  }
  
  if(eventcount==1){
  /*** KICK OUT WRONG PROCESSORS ***/
    if(sampleType!=SampleType::thq) {
      treewriter_nominal.RemoveTreeProcessor("tHqGenVarProcessor"); 
      treewriter_jesup.RemoveTreeProcessor("tHqGenVarProcessor"); 
      treewriter_jesdown.RemoveTreeProcessor("tHqGenVarProcessor"); 
      treewriter_jerup.RemoveTreeProcessor("tHqGenVarProcessor"); 
      treewriter_jerdown.RemoveTreeProcessor("tHqGenVarProcessor"); 
    }
    if(sampleType==SampleType::thq || sampleType == SampleType::nonttbkg || sampleType==SampleType::data || sampleType==SampleType::st){
      treewriter_nominal.RemoveTreeProcessor("TopGenVarProcessor");
      treewriter_jesup.RemoveTreeProcessor("TopGenVarProcessor");
      treewriter_jesdown.RemoveTreeProcessor("TopGenVarProcessor");
      treewriter_jerup.RemoveTreeProcessor("TopGenVarProcessor");
      treewriter_jerdown.RemoveTreeProcessor("TopGenVarProcessor");
    }
  }

  if(!isData&&foundT&&foundTbar) {
    // fill genTopEvt with tt(H) information
    genTopEvt.Fill(*h_genParticles,ttid_full);
    //cout << "APPARENTLY, THIS IS A SAMPLE WITH TWO TOP QUARKS" << endl;
  }
  
  if(sampleType ==  SampleType::thq){
    //cout << "APPARENTLY, THIS IS A THQ SAMPLE" << endl;
    gentHqEvt.Fill(*h_genParticles);
  }
  // DO REWEIGHTING

  vector<string> syst_weights_id;
  vector<float> syst_weights;
  float Weight_orig=1;
  map<string,float> weights_uncorrjets;
  map<string,float> weights;
  map<string,float> weights_jesup;
  map<string,float> weights_jerup;
  map<string,float> weights_jesdown;
  map<string,float> weights_jerdown;
  if(!isData){
    weights       = GetWeights(*h_geneventinfo,eventInfo,selectedPVs,selectedJets_nominal,selectedElectrons,selectedMuons,*h_genParticles,genTopEvt,sysType::NA,sampleType);
    weights_jesup = GetWeights(*h_geneventinfo,eventInfo,selectedPVs,selectedJets_jesup,selectedElectrons,selectedMuons,*h_genParticles,genTopEvt,sysType::JESup,sampleType);
    weights_jesdown = GetWeights(*h_geneventinfo,eventInfo,selectedPVs,selectedJets_jesdown,selectedElectrons,selectedMuons,*h_genParticles,genTopEvt,sysType::JESdown,sampleType);
    weights_jerup = GetWeights(*h_geneventinfo,eventInfo,selectedPVs,selectedJets_jerup,selectedElectrons,selectedMuons,*h_genParticles,genTopEvt,sysType::JERup,sampleType);
    weights_jerdown = GetWeights(*h_geneventinfo,eventInfo,selectedPVs,selectedJets_jerdown,selectedElectrons,selectedMuons,*h_genParticles,genTopEvt,sysType::JERdown,sampleType);
    weights_uncorrjets = GetWeights(*h_geneventinfo,eventInfo,selectedPVs,selectedJets_uncorrected,selectedElectrons,selectedMuons,*h_genParticles,genTopEvt,sysType::NA,sampleType);
    if(useLHE){
      GetSystWeights(*h_lheeventinfo,syst_weights_id,syst_weights,Weight_orig);
    }
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
  InputCollections input_nominal( eventInfo,
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
				  idJets,
				  rawJets,
				  selectedJets_nominal,
				  rawPuppiJets,
				  selectedPuppiJets,
				  selectedJetsLoose_nominal,
				  correctedMETs_nominal[0],
				  genTopEvt,
				  gentHqEvt,
				  selectedGenJets,
				  sampleType,
				  weights,
				  Weight_orig,
				  syst_weights,
				  syst_weights_id
				  );
       
  // define systematically shifted input (replace quantaties affected by jets)
  InputCollections input_jesup( input_nominal,selectedJets_jesup,selectedJetsLoose_jesup,correctedMETs_jesup[0],weights_jesup);
  InputCollections input_jesdown( input_nominal,selectedJets_jesdown,selectedJetsLoose_jesdown,correctedMETs_jesdown[0],weights_jesdown);
  InputCollections input_jerup( input_nominal,selectedJets_jerup,selectedJetsLoose_jerup,correctedMETs_jerup[0],weights_jerup);
  InputCollections input_jerdown( input_nominal,selectedJets_jerdown,selectedJetsLoose_jerdown,correctedMETs_jerdown[0],weights_jerdown);
  //  InputCollections input_uncorrjets( input_nominal,selectedJets_uncorrected,selectedJetsLoose_uncorrected,pfMETs[0],weights_uncorrjets);


  // DO SELECTION
  cutflow.EventSurvivedStep("all");
  bool selected=true;
  for(size_t i=0; i<selections.size() && selected; i++){
    if(!selections.at(i)->IsSelected(input_nominal,cutflow)) selected=false;
  }
  if(!selected) return;    

  // WRITE TREE
  treewriter_nominal.Process(input_nominal);  
  if (doSystematics) {
    treewriter_jesup.Process(input_jesup);
    treewriter_jesdown.Process(input_jesdown);
    treewriter_jerup.Process(input_jerup);
    treewriter_jerdown.Process(input_jerdown);
  }
}


map<string,float> tHqAnalyzer::GetWeights(const GenEventInfoProduct&  genEventInfo,
					  const EventInfo& eventInfo, 
					  const reco::VertexCollection& selectedPVs, 
					  const std::vector<pat::Jet>& selectedJets, 
					  const std::vector<pat::Electron>& selectedElectrons, 
					  const std::vector<pat::Muon>& selectedMuons, 
					  const std::vector<reco::GenParticle>& genParticles, 
					  const GenTopEvent& genTopEvt,
					  const sysType::sysType systype, 
					  const SampleType sampletype
					  ){
  
  map<string,float> weights;
  
  if(isData){
    weights["Weight"] = 1.0;
    weights["Weight_XS"] = 1.0;
    weights["Weight_CSV"] = 1.0;
    weights["Weight_PU"] = 1.0;
    weights["Weight_TopPt"] = 1.0;
    return weights;
  }

  float weight = 1.;
  assert(genEventInfo.weights().size()<=1); // before we multiply any weights we should understand what they mean
  for(size_t i=0;i<genEventInfo.weights().size();i++){
    
    weight *= (genEventInfo.weights()[i]>0 ? 1.: -1.); // overwrite intransparent MC weights, use \pm 1 instead
  }   // NOTE: IS THIS EVEN CORRECT?

  weight = genEventInfo.weight();

  double csvWgtHF, csvWgtLF, csvWgtCF;

  float xsweight = xs*luminosity/totalMCevents;
  float csvweight = 1.;
  float puweight = 1.;
  float topptweight = 1.;
  
  //Calculate TopPt Weight
  bool isTTbar;
  isTTbar = (sampletype == SampleType::ttl || sampletype == SampleType::ttb || sampletype == SampleType::ttbb ||sampletype == SampleType::ttcc || sampletype == SampleType::tt2b) ? true : false;
  topptweight = (isTTbar) ? GetTopPtWeight(genTopEvt.GetTop().pt(),genTopEvt.GetTopBar().pt()) : 1.;


  //get vectors of jet properties
  std::vector<double> jetPts;
  std::vector<double> jetEtas;
  std::vector<double> jetCSVs;
  std::vector<int> jetFlavors;
  for(std::vector<pat::Jet>::const_iterator itJet = selectedJets.begin(); itJet != selectedJets.end(); ++itJet){
    jetPts.push_back(itJet->pt());
    jetEtas.push_back(itJet->eta());
    jetCSVs.push_back(helper.GetJetCSV(*itJet,"pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    jetFlavors.push_back(itJet->hadronFlavour());
  }  


  // ADD CSV WEIGHTS HERE OR SEPARATLY?


  if(systype==sysType::JESup)csvweight= csvReweighter.getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,7, csvWgtHF, csvWgtLF, csvWgtCF);
  else if(systype==sysType::JESdown)csvweight= csvReweighter.getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,8, csvWgtHF, csvWgtLF, csvWgtCF);
  else if(systype==sysType::JERup)csvweight= csvReweighter.getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,0, csvWgtHF, csvWgtLF, csvWgtCF); //there are no SF for JER yet!!
  else if(systype==sysType::JERdown)csvweight= csvReweighter.getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,0, csvWgtHF, csvWgtLF, csvWgtCF); //there are no SF for JER yet!!
  else csvweight= csvReweighter.getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,0, csvWgtHF, csvWgtLF, csvWgtCF);

  // compute PU weights, and set nominal weight
  puWeights_.compute(eventInfo);
  puweight = puWeights_.nominalWeight();

  
  //weight *= xsweight*csvweight*puweight*topptweight;
 
  weights["Weight"] = weight;
  weights["Weight_XS"] = xsweight;
  weights["Weight_CSV"] = csvweight;
  weights["Weight_TopPt"] = topptweight;
  weights["Weight_PU"] = puweight;

  //  weights["Weight_PV"] = pvWeight.GetWeight(selectedPVs.size());

    
  cout << "PU weight :" << weights["Weight_PU"] << endl;
  cout << "CSV weight :" << weights["Weight_CSV"] << endl;
  cout << "TopPt weight :" << weights["Weight_TopPt"] << endl;

  if(systype != sysType::JESup && systype != sysType::JESup && systype != sysType::JERup && systype != sysType::JERdown) {

    weights["Weight_CSVLFup"] = csvReweighter.getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,9, csvWgtHF, csvWgtLF, csvWgtCF)/csvweight;
    weights["Weight_CSVLFdown"] = csvReweighter.getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,10, csvWgtHF, csvWgtLF, csvWgtCF)/csvweight;
    weights["Weight_CSVHFup"] = csvReweighter.getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,11, csvWgtHF, csvWgtLF, csvWgtCF)/csvweight;
    weights["Weight_CSVHFdown"] = csvReweighter.getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,12, csvWgtHF, csvWgtLF, csvWgtCF)/csvweight;
    weights["Weight_CSVHFStats1up"] = csvReweighter.getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,13, csvWgtHF, csvWgtLF, csvWgtCF)/csvweight;
    weights["Weight_CSVHFStats1down"] = csvReweighter.getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,14, csvWgtHF, csvWgtLF, csvWgtCF)/csvweight;
    weights["Weight_CSVLFStats1up"] = csvReweighter.getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,17, csvWgtHF, csvWgtLF, csvWgtCF)/csvweight;
    weights["Weight_CSVLFStats1down"] = csvReweighter.getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,18, csvWgtHF, csvWgtLF, csvWgtCF)/csvweight;
    weights["Weight_CSVHFStats2up"] = csvReweighter.getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,15, csvWgtHF, csvWgtLF, csvWgtCF)/csvweight;
    weights["Weight_CSVHFStats2down"] = csvReweighter.getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,16, csvWgtHF, csvWgtLF, csvWgtCF)/csvweight;
    weights["Weight_CSVLFStats2up"] = csvReweighter.getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,19, csvWgtHF, csvWgtLF, csvWgtCF)/csvweight;
    weights["Weight_CSVLFStats2down"] = csvReweighter.getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,20, csvWgtHF, csvWgtLF, csvWgtCF)/csvweight;
    weights["Weight_CSVCErr1up"] = csvReweighter.getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,21, csvWgtHF, csvWgtLF, csvWgtCF)/csvweight;
    weights["Weight_CSVCErr1down"] = csvReweighter.getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,22, csvWgtHF, csvWgtLF, csvWgtCF)/csvweight;
    weights["Weight_CSVCErr2up"] = csvReweighter.getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,23, csvWgtHF, csvWgtLF, csvWgtCF)/csvweight;
    weights["Weight_CSVCErr2down"] = csvReweighter.getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,24, csvWgtHF, csvWgtLF, csvWgtCF)/csvweight;
  }
  
  // set optional additional PU weights
  for(std::vector<PUWeights::Weight>::const_iterator it = puWeights_.additionalWeightsBegin();
      it != puWeights_.additionalWeightsEnd(); ++it) {
    weights[it->name()] = it->value();
  }
  if(eventInfo.puSummaryIsValid==0) weights["Weight_PU"]=1;

  return weights;
}



// Systematic weights
void tHqAnalyzer::GetSystWeights(const LHEEventProduct&  LHEEvent, vector<string> &weight_syst_id, vector<float> &weight_syst, float &Weight_orig ){

  //  map<string,float> syst_weights;
  
  //  syst_weights["Weight_orig"] = LHEEvent.originalXWGTUP();

  

  Weight_orig = LHEEvent.originalXWGTUP();

  //cout << "lalala " << LHEEvent.weights().size() <<endl;
  //cout << "lululu " << LHEEvent.weights()[491].id.c_str() <<endl;
  //  cout << "original weight:" << LHEEvent.originalXWGTUP() << endl;
  //  cout << "Weight_orig2:" << Weight_orig << endl;

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

float tHqAnalyzer::GetTopPtWeight(float toppt1, float toppt2){
  float sf1=exp(0.156-0.00137*toppt1);
  float sf2=exp(0.156-0.00137*toppt2);
  return sqrt(sf1*sf2);
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
 
  
  treewriter_nominal.AddSampleInformation();
  if(doSystematics){
    treewriter_jesup.AddSampleInformation();
    treewriter_jesdown.AddSampleInformation();
    treewriter_jerup.AddSampleInformation();
    treewriter_jerdown.AddSampleInformation();
  }

  cutflow.Print();
}
// ------------ method called when starting to processes a run ------------
// needed for the hlt_config_
void
tHqAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{

  // std::string hltTag="HLT";
  //bool hltchanged = true;
  //  if (!hlt_config.init(iRun, iSetup, hltTag, hltchanged)) {
  //  std::cout << "Warning, didn't find trigger process HLT,\t" << hltTag << std::endl;
  //  return;
  // }
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
