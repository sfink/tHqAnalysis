tHqAnalysis
=======

Framework of the tHq group of Karlsruhe, much credit to the ttH group

## Installation
Follow These Steps:

    export SCRAM_ARCH=slc6_amd64_gcc491
    scram project CMSSW_7_4_6_patch6
    cd CMSSW_7_4_6_patch6/src
    cmsenv   
    git cms-addpkg PhysicsTools/JetMCAlgos/
    cd PhysicsTools/JetMCAlgos/plugins/
    rm GenHFHadronMatcher.cc
    wget https://twiki.cern.ch/twiki/pub/CMSPublic/GenHFHadronMatcher/GenHFHadronMatcher.cc
    wget https://twiki.cern.ch/twiki/pub/CMSPublic/GenHFHadronMatcher/GenTtbarCategorizer.cc
    cd -
    cd PhysicsTools/JetMCAlgos/python/
    wget https://twiki.cern.ch/twiki/pub/CMSPublic/GenHFHadronMatcher/GenTtbarCategorizer_cfi.py.txt
    mv GenTtbarCategorizer_cfi.py.txt GenTtbarCategorizer_cfi.py
    cd -
    git cms-merge-topic gkasieczka:htt-v2-74X

    git clone https://github.com/cms-ttH/MiniAOD.git
    git clone https://github.com/sfink/tHqAnalysis.git
    
    scram b -j 12


