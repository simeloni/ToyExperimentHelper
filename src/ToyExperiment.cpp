#include "ToyExperiment.hpp"
#include "TFile.h"

ToyExperiment::ToyExperiment() {
    _nData = 10;
    _nRepetitions = 1;
    _outputFile = new TFile("outputFile.root", "RECREATE");
}

ToyExperiment::~ToyExperiment() {

}

void ToyExperiment::setNRepetitions(int nRepetitions) {
    
    if (nRepetitions < 0){
        std::cout << "ERROR: Number of repetitions in Toy Experiment is negative. Setting it to 0" << std::endl;
        _nRepetitions = 0;
    }
    else {
        _nRepetitions = nRepetitions;
    }
}

void ToyExperiment::setNData(int nData) {

    if (nData < 0){
        std::cout << "ERROR: Number of data events in Toy Experiment is negative. Setting it to 0" << std::endl;
        _nData = 0;
    }
    else {
        _nData = nData;
    }
}

void ToyExperiment::setNMC(int nMC) {

    if (nMC < 0){
        std::cout << "ERROR: Number of MC events in Toy Experiment is negative. Setting it to 0" << std::endl;
        _nMC = 0;
    }
    else {
        _nMC = nMC;
    }
}

void ToyExperiment::run() {

    std::vector<ToyModule*>::iterator it;

    for (int i=0; i<_nRepetitions; ++i) {

        for(it=_modules.begin(); it != _modules.end(); ++it)
            (*it)->beforeInitialize();

        initialize();

        for(it=_modules.begin(); it != _modules.end(); ++it)
            (*it)->beforeGenerate();
        
        generate(_nData);
        generateMC(_nMC);
        
        for(it=_modules.begin(); it != _modules.end(); ++it)
            (*it)->beforeFit();

        fit();
        
        for(it=_modules.begin(); it != _modules.end(); ++it)
            (*it)->beforeSave();
        
        save();

        for(it=_modules.begin(); it != _modules.end(); ++it)
            (*it)->afterSave();

    }

}

void ToyExperiment::addModule(ToyModule* module) {

    if (module != NULL) {
        _modules.push_back(module);
    }
    else {
        std::cout << "WARNING: Empty module passed to ToyExperiment" << std::endl;
    }
    
}