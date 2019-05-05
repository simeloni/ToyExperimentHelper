#include <iostream>
#include "TFile.h"
#include "ToyExperiment.hpp"
#include "ToyModule.hpp"

int main(int, char**) {
    
    ToyExperiment* experiment = new ToyExperiment();
    experiment->setNRepetitions(10);
    
    ToyModule* module = new ToyModule();

    experiment->addModule(module);
    experiment->run();

    return 0;
}
