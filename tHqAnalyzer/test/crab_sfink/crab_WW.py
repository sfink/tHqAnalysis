from WMCore.Configuration import Configuration
import os
os.environ['GLOBALTAG'] = '74X_mcRun2_asymptotic_v4'

config = Configuration()

config.section_("General")
config.General.requestName = 'WW_TuneCUETP8M1_13TeV-pythia8'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/desy.de/user/f/fink/xxl/af-cms/13TeV/CMSSW_7_4_15_patch1/src/tHqAnalysis/tHqAnalyzer/test/tHqAnalysis_cfg.py'
config.JobType.outputFiles = ['tHqAnalyzed_Tree.root']

config.section_("Data")
config.Data.inputDataset = '/WW_TuneCUETP8M1_13TeV-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v4-v1/MINIAODSIM'
                            
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.publication = False
#config.Data.totalUnits = 5
#config.Data.publishDbsUrl = 'phys03'
config.Data.outputDatasetTag = 'tHqAnalysis_MiniAOD'

config.General.transferOutputs = True

config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
