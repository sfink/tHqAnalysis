from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_wtruthcontainer_weights_wmass_test5jobs2'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/desy.de/user/f/fink/xxl/af-cms/13TeV/CMSSW_7_4_6_patch6/src/tHqAnalysis/tHqAnalyzer/test/tHqAnalysis_cfg.py'
config.JobType.outputFiles = ['tHqAnalyzed_Tree.root']

config.section_("Data")
config.Data.inputDataset = '/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.publication = False
config.Data.totalUnits = 5
#config.Data.publishDbsUrl = 'phys03'
config.Data.outputDatasetTag = 'THQ_MiniAOD'

config.General.transferOutputs = True

config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
config.Site.blacklist = ['T2_US_Florida']
