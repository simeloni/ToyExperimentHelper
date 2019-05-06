#include <iostream>
#include "TFile.h"
#include "ToyExperiment.hpp"
#include "ToyModule.hpp"
#include "ExampleParameters.hpp"

int main(int, char**) {
    
    //If you want to parallelize the job execution you can do it here
    ToyExperiment* experiment = new ToyExperiment();
    ToyExperiment* experiment2 = new ToyExperiment();
    experiment->setNRepetitions(10);
    experiment->setNData(1000.);
    experiment->setNMC(10000.);
    experiment->setOutputFileName("ResultsFile_0.root");
    
    ToyModule* module = new ToyModule();
    module->setName("module");
    ToyModule* module1 = new ToyModule();
    module1->setName("module1");

    ExampleParameters* pars = new ExampleParameters();

    experiment->addModule(module);
    experiment2->addModule(module);
    experiment->addModule(module1);
    experiment->run();
    experiment2->run();

    return 0;
}
