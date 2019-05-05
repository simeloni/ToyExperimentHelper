#include <iostream>
#include "TFile.h"
#include "ToyExperiment.hpp"
#include "ToyModule.hpp"

int main(int, char**) {
    
    ToyExperiment* experiment = new ToyExperiment();
    experiment->setNRepetitions(10);
    experiment->setNData(1000.);
    experiment->setNMC(10000.);
    
    ToyModule* module = new ToyModule();

    experiment->addModule(module);
    experiment->run();

    return 0;
}
