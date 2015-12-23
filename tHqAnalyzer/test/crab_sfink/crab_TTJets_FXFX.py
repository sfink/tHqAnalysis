from WMCore.Configuration import Configuration
import os
os.environ['GLOBALTAG'] = '74X_mcRun2_asymptotic_v4'
os.environ['ISDATA'] = "0"
os.environ['USELHE'] = "1"
os.environ['USEGENHADRONMATCH'] = "1"

config = Configuration()

config.section_("General")
config.General.requestName = 'TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/desy.de/user/f/fink/xxl/af-cms/13TeV/CMSSW_7_4_15_patch1/src/tHqAnalysis/tHqAnalyzer/test/tHqAnalysis_cfg.py'
config.JobType.outputFiles = ['tHqAnalyzed_Tree.root']

config.section_("Data")
config.Data.inputDataset = '/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v3/MINIAODSIM'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.publication = False
#config.Data.totalUnits = 5
#config.Data.publishDbsUrl = 'phys03'
config.Data.outputDatasetTag = 'THQ_MiniAOD'

config.General.transferOutputs = True

config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
config.Site.blacklist = ['T2_US_Florida']
os.environ['RECORRECTMET'] = "1"
