import FWCore.ParameterSet.Config as cms
import sys
import os
# To execute test, run
#  cmsRun boostedAnalysis_cfg.py isData=False outputFile=testrun maxEvents=100 inputFiles=file:/pnfs/desy.de/cms/tier2/store/user/hmildner/ttHTobb_M125_13TeV_powheg_pythia8/Boostedv2MiniAOD/151017_154254/0000/BoostedTTH_MiniAOD_1.root,file:/pnfs/desy.de/cms/tier2/store/user/hmildner/ttHTobb_M125_13TeV_powheg_pythia8/Boostedv2MiniAOD/151017_154254/0000/BoostedTTH_MiniAOD_2.root

# parse command-line arguments
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideCommandLineParsing
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
# The following variables are already defined in VarParsing class:
# maxEvents: singleton, int; default = -1
# inputFiles: (comma separated, no spaces!) list, string: default empty
options.register( "outName", "tHqAnalyzed", VarParsing.multiplicity.singleton, VarParsing.varType.string, "name and path of the output files (without extension)" )
#options.register( "inputFiles", "bla", VarParsing.multiplicity.singleton, VarParsing.varType.string, "name and path of the input file" )

options.register( "skipEvents", 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Number of events to skip" )
options.register( "isData", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "is it data or MC?" )
options.register( "globalTag", "76X_mcRun2_asymptotic_RunIIFall15DR76_v0", VarParsing.multiplicity.singleton, VarParsing.varType.string, "global tag" )
options.register( "doSystematics", True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Produce different files for JES and JER systematics?" )

options.parseArguments()

# re-set some defaults
if options.maxEvents is -1: # maxEvents is set in VarParsing class by default to -1
    options.maxEvents = 10000 # reset to 10000 for testing
if not options.inputFiles:
#    options.inputFiles=['file:/pnfs/desy.de/cms/tier2/store/user/hmildner/ttHTobb_M125_13TeV_powheg_pythia8/Boostedv3MiniAOD/151120_183808/0000/BoostedTTH_MiniAOD_10.root']
    options.inputFiles=['root://cmsxrootd.fnal.gov///store/mc/RunIIFall15MiniAODv2/THQ_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1/MI\
NIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/0AA8FC13-75C7-E511-B2D8-38EAA78E2C94.root']


# checks for correct values and consistency
if "data" in options.globalTag.lower() and not options.isData:
    print "\n\nConfig ERROR: GT contains seems to be for data but isData==False\n\n"
    sys.exit()
if "mc" in options.globalTag.lower() and options.isData:
    print "\n\nConfig ERROR: GT contains seems to be for MC but isData==True\n\n"
    sys.exit()
if not options.inputFiles:
    print "\n\nConfig ERROR: no inputFiles specified\n\n"
    sys.exit()
    
# print settings
print "\n\n***** JOB SETUP *************************"
for key in options._register:
    # do not print unused default options
    if key not in ["secondaryInputFiles","section","tag","totalSections","outputFile","secondaryOutputFile","filePrepend"]:
        print str(key)+" : "+str( options.__getattr__(key) )
print "*****************************************\n\n"


process = cms.Process("analysis")

# cmssw options
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = options.globalTag
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(int(options.maxEvents)))
process.source = cms.Source(  "PoolSource",
                              fileNames = cms.untracked.vstring(options.inputFiles),
                              skipEvents=cms.untracked.uint32(int(options.skipEvents)),
)


## Set up JetCorrections chain to be used in MiniAODHelper
## Note: name is hard-coded to ak4PFchsL1L2L3 and does not
## necessarily reflect actual corrections level
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *
from JetMETCorrections.Configuration.JetCorrectionCondDB_cff import *
process.ak4PFCHSL1Fastjet = cms.ESProducer(
  'L1FastjetCorrectionESProducer',
  level = cms.string('L1FastJet'),
  algorithm = cms.string('AK4PFchs'),
  srcRho = cms.InputTag( 'fixedGridRhoFastjetAll' )
  )
process.ak4PFchsL2Relative = ak4CaloL2Relative.clone( algorithm = 'AK4PFchs' )
process.ak4PFchsL3Absolute = ak4CaloL3Absolute.clone( algorithm = 'AK4PFchs' )
process.ak4PFchsResidual = ak4CaloResidual.clone( algorithm = 'AK4PFchs' )
process.ak4PFchsL1L2L3 = cms.ESProducer("JetCorrectionESChain",
  correctors = cms.vstring(
    'ak4PFCHSL1Fastjet',
    'ak4PFchsL2Relative',
    'ak4PFchsL3Absolute')
)
if options.isData:
  process.ak4PFchsL1L2L3.correctors.append('ak4PFchsResidual') # add residual JEC for data


#=================================== JEC from DB file for data ===============
if options.isData:
    process.load("CondCore.DBCommon.CondDBCommon_cfi")
    from CondCore.DBCommon.CondDBSetup_cfi import *
    process.jec = cms.ESSource("PoolDBESSource",
                               DBParameters = cms.PSet(
            messageLevel = cms.untracked.int32(0)
            ),
                               timetype = cms.string('runnumber'), # what does this do?
                               toGet = cms.VPSet(
            cms.PSet(
                record = cms.string('JetCorrectionsRecord'),
                tag    = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_DATA_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs')
                ),
#        ..................................................
            ## here you add as many jet types as you need
            ## note that the tag name is specific for the particular sqlite file 
            ), 
                               connect = cms.string('sqlite:///'+os.environ.get('CMSSW_BASE')+'/src/tHqAnalysis/tHqAnalyzer/data/jecs/Fall15_25nsV2_DATA.db')
#                               connect = cms.string('sqlite:../data/jecs/Fall15_25nsV2_DATA.db')
                               )
## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

#===============================================================


#process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
#process.selectedElectrons = cms.EDFilter("PATElectronSelector", src = cms.InputTag("slimmedElectrons"), cut = cms.string("pt > 5 && abs(eta)<2.5") ) 

# load and run the boosted analyzer
if options.isData:
    
    process.load("tHqAnalysis.tHqAnalyzer.tHqAnalyzer_cfi")
else:
    process.load("tHqAnalysis.tHqAnalyzer.tHqAnalyzer_cfi")
    

    # Supplies PDG ID to real name resolution of MC particles
    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
    # Needed to determine tt+x category -- is usually run when producing boosted jets in miniAOD 
    process.load("tHqAnalysis.tHqProducer.genHadronMatching_cfi")

process.tHqAnalyzer.outfileName=options.outName
process.tHqAnalyzer.doSystematics=options.doSystematics
#if not options.isData:
#    process.tHqAnalyzer.eventWeight = options.weight
#process.BoostedAnalyzer.generatorName=options.generatorName

#process.BoostedAnalyzer.doJERsystematic = False


#process.tHqAnalyzer.selectionNames = ["VertexSelection","LeptonSelection","JetTagSelection"]
#process.tHqAnalyzer.processorNames = ["WeightProcessor","MCMatchVarProcessor","BasicVarProcessor","MVAVarProcessor","BDTVarProcessor","TTbarReconstructionVarProcessor","ReconstructionMEvarProcessor","MEMProcessor","tHqTopHiggsVarProcessor","AdditionalJetProcessor"]


### electron MVA ####
# Load the producer for MVA IDs
process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
## check the event content 
process.content = cms.EDAnalyzer("EventContentAnalyzer")
#process.p = cms.Path(process.calibratedPatElectrons * process.electronMVAValueMapProducer * process.tHqAnalyzer)
process.p = cms.Path(process.electronMVAValueMapProducer * process.tHqAnalyzer)
