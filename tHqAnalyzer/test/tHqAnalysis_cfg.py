import FWCore.ParameterSet.Config as cms
import os
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *
from JetMETCorrections.Configuration.JetCorrectionCondDB_cff import *

# Set Process and Variables
#------------------------------------------------------------------------------------------------------------------------------------
process = cms.Process("boosted")

# grid-control variables
gc = {}
gc['nickname'] = '__NICK__'
gc['filenames'] = '__FILE_NAMES__'
gc['outfilename'] = '__OUT_FILE__'
gc['skip'] = '__SKIP_EVENTS__'
gc['max'] = '__MAX_EVENTS__'

gc['sampletype'] = '__SAMPLE_TYPE__'
gc['xs'] = '__XS__'
gc['mcevents'] = '__MCEVENTS__'
gc['globaltag'] = '__GLOBALTAG__'

# environment variables
env = {}
env['nickname'] = os.getenv('NICK_NAME')
env['filenames'] = os.getenv('FILE_NAMES')
env['outfilename'] = os.getenv('OUTFILE_NAME')
env['skip'] = os.getenv('SKIP_EVENTS')
env['max'] = os.getenv('MAX_EVENTS')

env['sampletype'] = os.getenv('SAMPLE_TYPE')
env['xs'] = os.getenv('XS')
env['mcevents'] = os.getenv('MCEVENTS')
env['globaltag'] = os.getenv('GLOBALTAG')

isData = os.getenv('ISDATA')

# default variables
default = {}
default['nickname'] = 'Lalala'
#default['filenames'] = 'root://cmsxrootd.fnal.gov///store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/00DF0A73-17C2-E511-B086-E41D2D08DE30.root'
default['filenames'] = 'root://cmsxrootd.fnal.gov///store/mc/RunIIFall15MiniAODv2/THQ_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/0AA8FC13-75C7-E511-B2D8-38EAA78E2C94.root'


default['outfilename'] = None
default['skip'] = '0'
default['max'] = '100'

default['sampletype'] = '9125'
default['xs'] = '248'
default['mcevents'] = '3500000'
default['globaltag'] = '76X_mcRun2_asymptotic_v12'


# fill in default values if not set by gc
values = gc.copy()
for key, value in values.iteritems():
    if value.startswith('__'):
        if env[key] is None:
            values[key] = default[key]
        else:
            values[key] = env[key]

print "The Global Tag is %s \n " % values['globaltag']
# convert strings
values['filenames'] = values['filenames'].strip(',')
values['filenames'] = map(lambda s: s.strip('" '), values['filenames'].split(","))

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = values['globaltag']

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(int(values['max'])))

process.source = cms.Source(  "PoolSource",
                              fileNames = cms.untracked.vstring(values['filenames']),
                              skipEvents = cms.untracked.uint32(int(values['skip']))                        
)



process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

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
if isData:
  process.ak4PFchsL1L2L3.correctors.append('ak4PFchsResidual') # add residual JEC for data


process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("tHqAnalysis.tHqProducer.genHadronMatching_cfi")

process.load("tHqAnalysis.tHqAnalyzer.tHqAnalyzer_cfi")
process.tHqAnalyzer.useFatJets=False
if values['outfilename'] is not None:
    process.tHqAnalyzer.outfileName=values['outfilename']
if values['sampletype'] is not None:
    process.tHqAnalyzer.sampleID=cms.int32(int(values['sampletype']))
if values['xs'] is not None:
    process.tHqAnalyzer.xs=cms.double(float(values['xs']))
if values['mcevents'] is not None:
    process.tHqAnalyzer.nMCEvents=cms.int32(int(values['mcevents']))      
    

### electron MVA ####
# Load the producer for MVA IDs
process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
## check the event content 
process.content = cms.EDAnalyzer("EventContentAnalyzer")


process.p = cms.Path(process.electronMVAValueMapProducer * process.tHqAnalyzer)
