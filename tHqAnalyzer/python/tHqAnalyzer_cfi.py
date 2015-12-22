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


tHqAnalyzer = cms.EDAnalyzer(
    'tHqAnalyzer',
    LeptonSelectionNoTrigger,
    JetTagSelection,  
    relevantTriggers = cms.vstring("HLT_Ele27_eta2p1_WP85_Gsf_HT200_v1","HLT_IsoMu24_eta2p1_v1"),
    era = cms.string("2015_74x"),
    analysisType = cms.string("LJ"),
    luminostiy = cms.double(19.7),
    sampleID = cms.int32(9125),
    xs = cms.double(248),
    nMCEvents = cms.int32(25000000),

    # PU weights, defined in Weights_cff
    nominalPUWeight = cms.PSet(NominalPUWeight),
    additionalPUWeights = cms.VPSet(AdditionalPUWeights),

    isData = cms.bool(var_isData),
    useFatJets = cms.bool(False),
    useLHE = cms.bool(var_useLHE),
    useGenHadronMatch = cms.bool(var_genHadronMatch),
    recorrectMET = cms.bool(var_recorrectMET),

    disableObjectSelections = cms.bool(False), # disables selection of some objects for synch exe
    outfileName = cms.string("tHqAnalyzed"),
    selectionNames = cms.vstring("LeptonSelection"),


    processorNames = cms.vstring("WeightProcessor","BaseVarProcessor","RecoVarProcessor","tHqGenVarProcessor","TopGenVarProcessor")
#    processorNames = cms.vstring("WeightProcessor","TopGenVarProcessor","MVAVarProcessor","tHqJetVarProcessor","tHqTopHiggsVarProcessor","tHqTopVarProcessor","tHqHiggsVarProcessor")
)
