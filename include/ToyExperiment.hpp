#ifndef TOYEXPERIMENT_HH
#define TOYEXPERIMENT_HH

#include "TFile.h"

class ToyExperiment {

public:
    ToyExperiment();
    ~ToyExperiment();

private:
    int _nRepetitions; //Number of times you want to repeat the experiment
    int _nData;        //Number of events you want to generate
    int _nMC;          //Number of MC events you want to generate. This is used if the fit is a templated fit

    TFile* _ouptutFile; //This is the file into which you want to save the output of your toy experiment
};

#endif