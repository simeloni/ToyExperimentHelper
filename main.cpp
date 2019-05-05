#include <iostream>
#include "TFile.h"
#include "ToyExperiment.hpp"
#include "ToyModule.hpp"

int main(int, char**) {
    
    //If you want to parallelize the job execution you can do it here
    ToyExperiment* experiment = new ToyExperiment();
    experiment->setNRepetitions(10);
    experiment->setNData(1000.);
    experiment->setNMC(10000.);
    experiment->setOutputFileName("ResultsFile_0.root");
    
    ToyModule* module = new ToyModule();

    experiment->addModule(module);
    experiment->run();

    return 0;
}
