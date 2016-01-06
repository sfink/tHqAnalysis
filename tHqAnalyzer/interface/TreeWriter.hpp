#ifndef TREEWRITER_HPP
#define TREEWRITER_HPP

#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"

#include "tHqAnalysis/tHqAnalyzer/interface/VariableContainer.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/InputCollections.hpp"
#include "tHqAnalysis/tHqAnalyzer/interface/TreeProcessor.hpp"

/*
  The TreeWriter is used to write all the variables that are used in the MVA analysis in flat TTrees. Different Processors can be loaded that write a certain class of variables. Calculating BDT outputs and histograms from these variables will be supported soon.
 */
class TreeWriter{
  public:
    
    /**
       Creates a TreeWriter and the associated TTree/ 
    */
    TreeWriter();
    ~TreeWriter();
    
    /**
       Process a single event.
       @param outfileName filename of the created TTree
       @param input input collections (cannot be changed, this is left to the BEANRunner)
       @param sampleType type of the sample, decides wheter data or mc is analyzed, which mc-matching is possible and which ttbar subsample the event belongs to
       @param weights the nominal and systematics weights of the event
    */
    void Init(std::string fileName);
    bool Process(const InputCollections& input);
    void AddTreeProcessor(TreeProcessor* processor);
    void AddTreeProcessor(TreeProcessor* processor,string name);
    void AddSampleInformation();
    void FillProcessorMap();
    void RemoveTreeProcessor(string name);
    void FillProcessorName(string name);

    std::vector<TreeProcessor*> GetTreeProcessors() const;
    std::vector<std::string> GetTreeProcessorNames() const;


  private:
  
    void Init();
    void FillTree();
    bool initialized;
    TTree* tree;
    TDirectory* dir;
    TFile* outFile;
    VariableContainer vars;
    std::vector<TreeProcessor*> processors;
    std::vector<TStopwatch> stopwatches;
    std::vector<std::string> processorNames;
    std::map<std::string,TreeProcessor*>  ProcessorMap;
};
#endif
