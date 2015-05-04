#ifndef THQANALYZER_THQEVENT_HPP
#define THQANALYZER_THQEVENT_HPP

#include "tHqAnalysis/tHqObjects/interface/Event.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "tHqAnalysis/tHqAnalyzer/interface/InputCollections.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/tHqUtils.hpp"

class tHqEvent{
  
  public:
    
    // Constructor & Destructor
    tHqEvent(const InputCollections& input);
    ~tHqEvent();
    
    // Manage Event
    void ResetEvent();
    
    // Reconstruction Functions
    
    // Charged Lepton Reconstruction
    void LeptonRec();
    
    // Neutrino Reconstruction
    void NeutrinoRec();
    
    // Anti-kt 5 Jets Reconstruction
    void ak5JetsRec();
    void ak5JetsClean(bool cleanHiggsCand = false, bool cleanTopHadCand = false, bool cleanTopLepCand = false);    
    
    // Get Functions
    // Input Collection
    const InputCollections& GetInput();
    
    // Charged Lepton Reconstruction
    math::XYZTLorentzVector  GetLeptonVec();
    
    // Neutrino Reconstruction
    math::XYZTLorentzVector  GetNeutrinoVec();

    // W Reconstruction
    math::XYZTLorentzVector  GetWVec();

    
    // Anti-kt 5 Jets
    pat::JetCollection Getak5JetsAll();
    int             GetNJets();
    int             GetNBTagL();
    int             GetNBTagM();
    int             GetNBTagT();
    float           GetAverageCSV();
    
    pat::JetCollection Getak5JetsCleaned();
    int             GetNCleanedak5Jets();
    int             GetNCleanedBTagL();
    int             GetNCleanedBTagM();
    int             GetNCleanedBTagT();
    float           GetAverageCSVClean();
    
  private:
    
    // Inhput Collection
    const InputCollections& input;
    
    // Analysis Settings
    
    // Charged Lepton
    math::XYZTLorentzVector lepVecCand;
    
    // Neutrino
    math::XYZTLorentzVector nuVecCand;

    // W Boson
    math::XYZTLorentzVector wVecCand;

    
    // Anti-kt 5 Jets
    pat::JetCollection selectedJets;
    std::vector<bool> BTagL;
    std::vector<bool> BTagM;
    std::vector<bool> BTagT;
    int nJets;
    int nBTagL;
    int nBTagM;
    int nBTagT;
    
    pat::JetCollection cleanedak5Jets;
    int nCleanedak5Jets;
    int nCleanedBTagL;
    int nCleanedBTagM;
    int nCleanedBTagT;
    
    bool verbose;
    const char* btagger;
};

#endif
