#include "ToyExperiment.hpp"
#include "TFile.h"

ToyExperiment::ToyExperiment() {
    _nData = 10;
    _nRepetitions = 1;
    _outputFile = new TFile("outputFile.root", "RECREATE");
}

ToyExperiment::~ToyExperiment() {

}

void ToyExperiment::run() {

    std::vector<ToyModule*>::iterator it;

    for(it=_modules.begin(); it != _modules.end(); ++it)
        (*it)->beforeInitialize();

    initialize();

    for(it=_modules.begin(); it != _modules.end(); ++it)
        (*it)->beforeGenerate();
    
    generate();
    generateMC();
    
    for(it=_modules.begin(); it != _modules.end(); ++it)
        (*it)->beforeFit();

    fit();
    
    for(it=_modules.begin(); it != _modules.end(); ++it)
        (*it)->beforeSave();
    
    save();

    for(it=_modules.begin(); it != _modules.end(); ++it)
        (*it)->afterSave();
}

void ToyExperiment::addModule(ToyModule* module) {

    if (module != NULL) {
        _modules.push_back(module);
    }
    else {
        std::cout << "WARNING: Empty module passed to ToyExperiment" << std::endl;
    }
    
}