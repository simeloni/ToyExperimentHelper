#ifndef TOYEXPERIMENT_HH
#define TOYEXPERIMENT_HH

#include "TFile.h"
#include "TString.h"
#include <iostream>
#include <vector>
#include "ToyModule.hpp"

class ToyModule;

class ToyExperiment {

public:
    ToyExperiment();
    ~ToyExperiment();

    void addModule(ToyModule* module);
    //void addParameters(parameters* parameters);

    void setNRepetitions(int nRepetitions);
    void setNData(int nData);
    void setNMC(int nMC);

    void setOutputFileName(TString outputFileName);

    void run();

private: 
    void createOutputFile();

private:
    //TODO: Maybe you want to return some code status and catch the error codes in the running
    virtual void initialize(){std::cout << "Initializing" << std::endl;};
    virtual void generate(int nData){std::cout << "Generating data" << std::endl;};
    virtual void generateMC(int nMC){std::cout << "Generating MC" << std::endl;};
    virtual void fit(){std::cout << "Fitting" << std::endl;};
    virtual void save(){std::cout << "Saving the results" << std::endl;};

private:
    int _nRepetitions; //Number of times you want to repeat the experiment. Use it in the generateMC function
    int _nData;        //Number of events you want to generate. Use it in the generate function
    int _nMC;          //Number of MC events you want to generate. This is used if the fit is a templated fit
    TString _outputFileName;

    TFile* _outputFile; //This is the file into which you want to save the output of your toy experiment

    std::vector<ToyModule*> _modules;
    //std::vector<Parameters*> _parameters;
};

#endif