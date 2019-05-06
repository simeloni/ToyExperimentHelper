#include "ToyExperiment.hpp"
#include "TFile.h"

ToyExperiment::ToyExperiment() {
    _nData = 10;
    _nMC   = 10;
    _nRepetitions = 0;
    _outputFileName = TString("outputFile.root");
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

void ToyExperiment::setOutputFileName(TString outputFileName) {

    if(outputFileName == TString("")){
        std::cout << "WARNING: OutputFileName set to empty String! Using standard one" << std::endl;
    }
    else {
        _outputFileName = outputFileName;
    }
}

void ToyExperiment::createOutputFile() {

    _outputFile = new TFile(_outputFileName, "RECREATE");

}

void ToyExperiment::run() {

    createOutputFile();

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

    ToyModule* module_ = new ToyModule(*module);
    module_->setReferenceToExperiment(this); 

    if (module != NULL) {
        _modules.push_back(module_);
    }
    else {
        std::cout << "WARNING: Empty module passed to ToyExperiment" << std::endl;
    }
    
}
/*
void ToyExperiment::addParameters(parameters* parameters) {

    if (parameter != NULL) {
        _parameters.push_back(parameters);
    }
    else {
        std::cout << "WARNING: Empty parameters passed to ToyExperiment" << std::endl;
    }
}
*/