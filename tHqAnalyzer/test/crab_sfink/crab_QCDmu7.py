from WMCore.Configuration import Configuration
import os
os.environ['GLOBALTAG'] = '76X_mcRun2_asymptotic_RunIIFall15DR76_v0'
os.environ['ISDATA'] = "0"
os.environ['USELHE'] = "0"

config = Configuration()

config.section_("General")
config.General.requestName = 'QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/desy.de/user/f/fink/xxl/af-cms/13TeV/CMSSW_7_6_3/src/tHqAnalysis/tHqAnalyzer/test/tHqAnalysis_cfg.py'
#config.JobType.outputFiles = ['tHqAnalyzed_Tree.root']

config.section_("Data")
config.Data.inputDataset = '/QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM'
                            
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.publication = False
#config.Data.totalUnits = 2
#config.Data.publishDbsUrl = 'phys03'
config.Data.outputDatasetTag = 'tHqAnalysis_MiniAOD'

config.General.transferOutputs = True

config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
os.environ['RECORRECTMET'] = "1"
config.JobType.outputFiles = ['tHqAnalyzed_JERDOWN_Tree.root', 'tHqAnalyzed_JERUP_Tree.root', 'tHqAnalyzed_JESDOWN_Tree.root', 'tHqAnalyzed_JESUP_Tree.root', 'tHqAnalyzed_nominal_Tree.root']
