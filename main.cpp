#include <iostream>
#include "TFile.h"
#include "ToyExperiment.hpp"
#include "ToyModule.hpp"
#include "ExampleParameters.hpp"
#include "ExampleExperiment.hpp"
#include "ExampleModule.hpp"

int main(int, char**) {
    
    //If you want to parallelize the job execution you can do it here
    ExampleExperiment* experiment = new ExampleExperiment();
    
    ExampleParameters* pars = new ExampleParameters();
    experiment->addParameters(*pars);
    experiment->setNRepetitions(10000);

    ExampleModule* module = new ExampleModule();
    module->setName("ExampleModule");

    experiment->addModule(*module);
    experiment->run();

    return 0;
}
