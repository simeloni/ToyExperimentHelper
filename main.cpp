#include <iostream>
#include "TFile.h"
#include "ToyExperiment.hpp"
#include "ToyModule.hpp"

int main(int, char**) {
    std::cout << "Hello, world!\n";
    
    ToyExperiment* experiment = new ToyExperiment();
    
    ToyModule* module = new ToyModule();

    experiment->addModule(module);
    experiment->run();

    return 0;
}
