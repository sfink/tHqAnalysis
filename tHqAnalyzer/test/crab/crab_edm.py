from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'DATA_edm_all_2'
config.General.workArea = 'crab_projects'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'edm.py'
#config.JobType.psetName = 'single_top_t-chan.py'
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = '/SingleMuon/Run2015B-PromptReco-v1/MINIAOD'

config.Data.splitting = 'LumiBased'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/DCSOnly/json_DCSONLY_Run2015B.txt'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-251252_13TeV_PromptReco_Collisions15_JSON.txt'
config.Data.inputDBS = 'global'
config.Data.unitsPerJob = 50
config.Data.outLFNDirBase = '/store/user/nfalterm/single_top/' # or '/store/group/<subdir>'
config.Data.publication = True
config.Data.publishDBS = 'phys03'
config.Data.publishDataName = 'DATA_new_all_EDM_2'

config.Site.storageSite = 'T2_DE_DESY'
config.Data.ignoreLocality =  True
#config.Site.whitelist = ['T2_CH_*']
config.Site.blacklist = ['T2_US_Vanderbilt']

config.User.voGroup = 'dcms'





