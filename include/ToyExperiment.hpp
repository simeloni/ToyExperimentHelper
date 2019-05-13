#ifndef TOYEXPERIMENT_HH
#define TOYEXPERIMENT_HH

#include "TFile.h"
#include "TString.h"
#include <iostream>
#include <vector>
#include "ToyModule.hpp"
#include "FitParameters.hpp"
#include "Parameter.hpp"
#include "TTree.h"

class ToyModule;

class ToyExperiment {

public:
    ToyExperiment();
    ~ToyExperiment();

    void addModule(ToyModule& module);
    void addParameters(FitParameters& parameters);

    void setNRepetitions(int nRepetitions);
    void setNData(int nData);
    void setNMC(int nMC);

    void setOutputFileName(TString outputFileName);

    void run();

    int getNParameters();
    std::vector<FitParameters*> getParameters();

private: 
    void createOutputFile();

private:
    //TODO: Maybe you want to return some code status and catch the error codes in the running
    void initialize();
    virtual void generate(int nData){std::cout << "Generating data" << std::endl;};
    virtual void generateMC(int nMC){std::cout << "Generating MC" << std::endl;};
    virtual void fit(){std::cout << "Fitting" << std::endl;};
    void save();

private:
    int _nRepetitions; //Number of times you want to repeat the experiment. Use it in the generateMC function
    int _nData;        //Number of events you want to generate. Use it in the generate function
    int _nMC;          //Number of MC events you want to generate. This is used if the fit is a templated fit
    TString _outputFileName;

    int _nParams; //Keep track of the number of parameters you have added

    TFile* _outputFile; //This is the file into which you want to save the output of your toy experiment
    TTree* _outputTree; //This is the outputTree into which you want to save the output of your toy experiment

    std::vector<ToyModule*> _modules;

    std::vector<double*> _initial_values;
    std::vector<double*> _fitted_values;
    std::vector<double*> _errorHigh_values;
    std::vector<double*> _errorLow_values;
    std::vector<double*> _error_values;
    std::vector<double*> _pull_values;
    std::vector<double*> _residual_values;    

protected: 
    std::vector<FitParameters*> _parameters;    
};

#endif