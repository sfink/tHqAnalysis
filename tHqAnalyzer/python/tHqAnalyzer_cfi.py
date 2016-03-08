import FWCore.ParameterSet.Config as cms
from tHqAnalysis.tHqAnalyzer.Selection_cff import *
from tHqAnalysis.tHqAnalyzer.Weights_cff import *

import os

 
if os.getenv('USELHE')=='1':
    print "useLHE is 1 -> Filling LHE weights"
    var_useLHE=True
else: 
    var_useLHE=False
    print "USE LHE: %s (for LHE usage : \"export USELHE=1\")" % var_useLHE

if os.getenv('ISDATA')=='1':
    print "isData is 1 -> Not using GenInfo"
    var_isData=True
else: 
    var_isData=False
    print "isData: %s (for isData usage : \"export ISDATA=1\")" % var_isData


if os.getenv('USEGENHADRONMATCH')=='1':
    print "genHadronMatch is 1 -> Doing tt+x splitting + matching"
    var_genHadronMatch=True
else: 
    var_genHadronMatch=False
    print "genHadronMatch: %s (for useGenHadronMatch usage : \"export USEGENHADRONMATCH=1\")" % var_genHadronMatch

if os.getenv('RECORRECTMET')=='1':
    print "recorrectMET is 1 -> Recorrecting MiniAOD MET"
    var_recorrectMET=True
else: 
    var_recorrectMET=False
    print "Recorrect MET: %s (To activate use : \"export RECORRECTMET=1\")" % var_recorrectMET


Trigger_data = cms.vstring("HLT_IsoMu18_v*",
                         "HLT_Ele23_WPLoose_Gsf_v*",
                         "HLT_Mu20_v*",
                         "HLT_Ele27_eta2p1_WPLoose_Gsf_v*",
                         "HLT_IsoMu20_v*",
                         "HLT_IsoTkMu20_v*",
                         "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
                         "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
                         "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*",
                         "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
                         "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*")

Trigger_mc = cms.vstring("HLT_IsoMu18_v*",
                         "HLT_Ele23_WPLoose_Gsf_v*",
                         "HLT_Mu20_v*",
                         "HLT_Ele27_eta2p1_WPLoose_Gsf_v*",
                         "HLT_IsoMu20_v*",
                         "HLT_IsoTkMu20_v*",
                         "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
                         "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
                         "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*",
                         "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
                         "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*")

if(var_isData):
    var_relevantTrigger = Trigger_data
else:
    var_relevantTrigger = Trigger_mc



tHqAnalyzer = cms.EDAnalyzer(
    'tHqAnalyzer',
    LeptonSelectionNoTrigger,
    JetTagSelection,  
    isData = cms.bool(var_isData),
    relevantTriggers = var_relevantTrigger,
    era = cms.string("2015_74x"),
    analysisType = cms.string("LJ"),
    luminostiy = cms.double(19.7),
    sampleID = cms.int32(9125),
    xs = cms.double(248),
    nMCEvents = cms.int32(25000000),
#    doSystematics = options.doSystematics,
    # PU weights, defined in Weights_cff
    nominalPUWeight = cms.PSet(NominalPUWeight),
    additionalPUWeights = cms.VPSet(AdditionalPUWeights),
    doSystematics = cms.bool(False),

    useFatJets = cms.bool(False),
    useLHE = cms.bool(var_useLHE),
    useGenHadronMatch = cms.bool(var_genHadronMatch),
    recorrectMET = cms.bool(var_recorrectMET),

    disableObjectSelections = cms.bool(False), # disables selection of some objects for synch exe
    outfileName = cms.string("tHqAnalyzed"),
    selectionNames = cms.vstring("LeptonSelection"),


    processorNames = cms.vstring("WeightProcessor","BaseVarProcessor","RecoVarProcessor","tHqGenVarProcessor","TopGenVarProcessor","TriggerVarProcessor")
#    processorNames = cms.vstring("WeightProcessor","TopGenVarProcessor","MVAVarProcessor","tHqJetVarProcessor","tHqTopHiggsVarProcessor","tHqTopVarProcessor","tHqHiggsVarProcessor")
)
