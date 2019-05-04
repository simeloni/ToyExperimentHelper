#include "ToyExperiment.hpp"
#include "TFile.h"

ToyExperiment::ToyExperiment() {
    _nData = 10;
    _nRepetitions = 1;
    _outputFile = new TFile("outputFile.root");
}

ToyExperiment::~ToyExperiment() {

}

void ToyExperiment::run() {
    initialize();
    generate();
    generateMC();
    fit();
    save();
}