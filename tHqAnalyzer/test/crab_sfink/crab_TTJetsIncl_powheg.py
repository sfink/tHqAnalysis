from WMCore.Configuration import Configuration
import os
os.environ['GLOBALTAG'] = '76X_mcRun2_asymptotic_RunIIFall15DR76_v0'
os.environ['ISDATA'] = "0"
os.environ['USELHE'] = "1"
os.environ['USEGENHADRONMATCH'] = "1"
os.environ['RECORRECTMET'] = "1"

config = Configuration()

config.section_("General")
config.General.requestName = 'TT_TuneCUETP8M1_13TeV-powheg-pythia8_ext3'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/desy.de/user/f/fink/xxl/af-cms/13TeV/CMSSW_7_6_3/src/tHqAnalysis/tHqAnalyzer/test/tHqAnalysis_cfg.py'
#config.JobType.outputFiles = ['tHqAnalyzed_Tree.root']

config.section_("Data")
config.Data.inputDataset = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/MINIAODSIM'

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
config.JobType.outputFiles = ['tHqAnalyzed_JERDOWN_Tree.root', 'tHqAnalyzed_JERUP_Tree.root', 'tHqAnalyzed_JESDOWN_Tree.root', 'tHqAnalyzed_JESUP_Tree.root', 'tHqAnalyzed_nominal_Tree.root']
