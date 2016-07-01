tHqAnalysis
=======

Framework of the tHq group of Karlsruhe, much credit to the ttH group

## Installation
Follow These Steps:

    scram project CMSSW_7_4_15_patch1
    cd CMSSW_7_4_15_patch1/src/
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
    git cms-merge-topic gkasieczka:htt-v2-76X
    
    git clone git@github.com:sfink/tHqAnalysis.git -b CMSSW_76x
    git clone git@github.com:sfink/MiniAOD.git -b CMSSW_76x
    
    scram b -j 12

    


