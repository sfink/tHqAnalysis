#include "DataFormats/Common/interface/Wrapper.h"

//Add includes for your classes here
#include "tHqAnalysis/tHqObjects/interface/Event.h"

#include <vector>

namespace {

  struct BEAN_Collections {

    //add 'dummy' Wrapper variable for each class type you put into the Event
    boosted::Event eventdummy0;
    edm::Wrapper<boosted::Event> eventdummy1;
    std::vector<boosted::Event> eventdummy2;
    edm::Wrapper<std::vector<boosted::Event> > eventdummy3;
  };
}
