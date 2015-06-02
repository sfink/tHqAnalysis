import FWCore.ParameterSet.Config as cms
from tHqAnalysis.tHqAnalyzer.Selection_cff import *

tHqAnalyzer = cms.EDAnalyzer(
    'tHqAnalyzer',
    LeptonSelectionNoTrigger,
    JetTagSelection,  
    relevantTriggers = cms.vstring("HLT_IsoMu24_eta2p1_IterTrk02_v1","HLT_Ele27_eta2p1_WP85_Gsf_v1"),
    era = cms.string("2012_53x"),
    analysisType = cms.string("LJ"),
    luminostiy = cms.double(19.7),
    sampleID = cms.int32(9125),
    xs = cms.double(248),
    nMCEvents = cms.int32(25000000),
    isData = cms.bool(False),
    useFatJets = cms.bool(False),
    disableObjectSelections = cms.bool(False), # disables selection of some objects for synch exe
    outfileName = cms.string("tHqAnalyzed"),
    selectionNames = cms.vstring("LeptonSelection"),
#    processorNames = cms.vstring("MVAVarProcessor","BaseVarProcessor")
    processorNames = cms.vstring("BaseVarProcessor","RecoVarProcessor")
#    processorNames = cms.vstring("WeightProcessor","MCMatchVarProcessor","MVAVarProcessor","tHqJetVarProcessor","tHqTopHiggsVarProcessor","tHqTopVarProcessor","tHqHiggsVarProcessor")
)
