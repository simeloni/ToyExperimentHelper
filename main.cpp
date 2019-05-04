#include <iostream>
#include "TFile.h"
#include "ToyExperiment.hpp"

int main(int, char**) {
    std::cout << "Hello, world!\n";
    
    ToyExperiment* experiment = new ToyExperiment();
    
    return 0;
}
