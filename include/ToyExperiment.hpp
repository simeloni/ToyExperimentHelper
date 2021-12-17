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

    int getNRepetitions();
    int getNData();
    int getNMC();

    void setOutputDirectory(const char* path);
    char* getOutputDirectory() const;

    void setOutputFileName(TString outputFileName);
    bool buildDirectoryTree();

    void run();

    int getNParameters();
    std::vector<FitParameters*> getParameters();


    void setParameterValue(const char* name, double v);
    void setParameterStartValue(const char* name, double v);
    void setParameterTrueValue(const char* name, double v);
    void setParameterErrorHigh(const char* name, double v);
    void setParameterErrorLow(const char* name, double v);

    double getParameterValue(const char* name);
    double getParameterStartValue(const char* name);
    double getParameterTrueValue(const char* name);
    double getParameterErrorHigh(const char* name);
    double getParameterErrorLow(const char* name);

    int getRepetition(){ return _currentRepetition; };

private: 
    void createOutputFile();

    parameter* getParameter(const char* name);

private:
    //TODO: Maybe you want to return some code status and catch the error codes in the running
    void initialize();
    virtual void setup(){std::cout << "INFO: Setup not implemented" << std::endl;};
    virtual void closing(){std::cout << "INFO: Closing not implemented" << std::endl;}; 
    virtual void generate(int nData){std::cout << "INFO: Generating data not implemented" << std::endl;};
    virtual void generateMC(int nMC){std::cout << "INFO: Generating MC not implemented" << std::endl;};
    virtual void fit(){std::cout << "INFO: Fitting not implemented" << std::endl;};
    void save();

private:
    int _nRepetitions; //Number of times you want to repeat the experiment. Use it in the generateMC function
    int _nData;        //Number of events you want to generate. Use it in the generate function
    int _nMC;          //Number of MC events you want to generate. This is used if the fit is a templated fit
    
    int _currentRepetition; //Current repetition number

    TString _outputFileName ;
    TString _dirPath        ;
    TString _xmlPath        ;
    TString _xmlPathTrue    ;
    TString _xmlPathInit    ;
    TString _xmlPathFinal   ;

    int _nParams; //Keep track of the number of parameters you have added

    TFile* _outputFile; //This is the file into which you want to save the output of your toy experiment
    TTree* _outputTree; //This is the outputTree into which you want to save the output of your toy experiment

    std::vector<ToyModule*> _modules;

    std::vector<double*> _initial_values;
    std::vector<double*> _fitted_values;
    std::vector<double*> _true_values;
    std::vector<double*> _errorHigh_values;
    std::vector<double*> _errorLow_values;
    std::vector<double*> _error_values;
    std::vector<double*> _pull_values;
    std::vector<double*> _residual_values;    

protected: 
    std::vector<FitParameters*> _parameters;    
};

#endif