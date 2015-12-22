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
#default['filenames'] = 'root://cmsxrootd.fnal.gov///store/user/shwillia/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/BoostedTTH_MiniAOD/150227_111650/0000/BoostedTTH_MiniAOD_15.root'
#default['filenames'] = 'root://cmsxrootd.fnal.gov///store/mc/Phys14DR/TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1E2D2522-A46A-E411-9C55-002590D0AFDC.root'
default['filenames'] = 'root://cmsxrootd.fnal.gov///store/mc/RunIISpring15MiniAODv2/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v3/60000/00181849-176A-E511-8B11-848F69FD4C94.root'
#default['filenames'] = 'root://cmsxrootd.fnal.gov///store/mc/RunIISpring15DR74/ZZ_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/0C479546-7209-E511-BA6A-3417EBE8862E.root'
#default['filenames'] = 'root://cmsxrootd.fnal.gov///store/mc/RunIISpring15DR74/QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/5084EE68-9507-E511-A8ED-0025905C3D40.root'
#default['filenames'] = 'file:/nfs/dust/cms/user/bmaier/CMSSW_7_4_6_patch6/src/tHqAnalysis/06249A8D-FE54-E511-825E-008CFA1111EC.root'
#default['filenames'] = 'root://cmsxrootd.fnal.gov///store/data/Run2015D/SingleMuon/MINIAOD/05Oct2015-v1/10000/025A01CA-8B6F-E511-B7A5-0025905A6084.root'
default['outfilename'] = None
default['skip'] = '0'
default['max'] = '100'

default['sampletype'] = '9125'
default['xs'] = '248'
default['mcevents'] = '3500000'
default['globaltag'] = '74X_mcRun2_asymptotic_v4'


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
