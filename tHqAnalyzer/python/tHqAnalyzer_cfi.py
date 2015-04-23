import FWCore.ParameterSet.Config as cms
from tHqAnalysis.tHqAnalyzer.Selection_cff import *

tHqAnalyzer = cms.EDAnalyzer(
    'tHqAnalyzer',
    LeptonSelectionNoTrigger,
    JetTagSelection,  
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
    processorNames = cms.vstring("MVAVarProcessor")
#    processorNames = cms.vstring("WeightProcessor","MCMatchVarProcessor","MVAVarProcessor","tHqJetVarProcessor","tHqTopHiggsVarProcessor","tHqTopVarProcessor","tHqHiggsVarProcessor")
)
