from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'THQ_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/nfs/dust/cms/user/bmaier/CMSSW_7_4_15/src/tHqAnalysis/tHqAnalyzer/test/tHqAnalysis_cfg.py'
config.JobType.outputFiles = ['THQ-madgraph-pythia8.root']

config.section_("Data")
config.Data.inputDataset = '/THQ_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.publication = False
#config.Data.totalUnits = 5
#config.Data.publishDbsUrl = 'phys03'
config.Data.outputDatasetTag = 'Gridding_Round_One'

config.General.transferOutputs = True

config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
