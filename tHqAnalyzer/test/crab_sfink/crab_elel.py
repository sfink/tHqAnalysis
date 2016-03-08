from WMCore.Configuration import Configuration
import os
os.environ['GLOBALTAG'] = '76X_dataRun2_v15'
os.environ['ISDATA'] = "1"
os.environ['USELHE'] = "1"
os.environ['RECORRECTMET'] = "1"

config = Configuration()

config.section_("General")
config.General.requestName = 'DoubleEG_json'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/desy.de/user/f/fink/xxl/af-cms/13TeV/CMSSW_7_6_3/src/tHqAnalysis/tHqAnalyzer/test/tHqAnalysis_cfg.py'
config.JobType.outputFiles = ['tHqAnalyzed_nominal_Tree.root']

config.section_("Data")
config.Data.inputDataset = '/DoubleEG/Run2015D-16Dec2015-v2/MINIAOD'

                            
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.publication = False
config.Data.totalUnits = 2
#config.Data.publishDbsUrl = 'phys03'
config.Data.outputDatasetTag = 'tHqAnalysis_MiniAOD'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_Silver_v2.txt'
config.General.transferOutputs = True

config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
os.environ['RECORRECTMET'] = "1"
