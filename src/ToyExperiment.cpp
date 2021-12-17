#include "ToyExperiment.hpp"
#include "TFile.h"
#include "TBranch.h"
#include "TTree.h"
#include "TString.h"
#include <sys/stat.h>
#include <cmath>

ToyExperiment::ToyExperiment() {
    _nData = 10;
    _nMC   = 10;
    _nRepetitions = 0;
    _currentRepetition = 0;
    _outputFileName = TString("outputFile.root");

    _dirPath      = TString("./outputExperiment")          ;
    _xmlPath      = TString(_dirPath + "/Xml/" )  ;
    _xmlPathTrue  = TString(_xmlPath  + "/true/" );
    _xmlPathInit  = TString(_xmlPath  + "/init/" );
    _xmlPathFinal = TString(_xmlPath  + "/final/");
    
    _nParams = 0;
}

ToyExperiment::~ToyExperiment() {

}

void ToyExperiment::setOutputDirectory(const char* path){
    _dirPath      = TString(path)                 ;
    _xmlPath      = TString(_dirPath + "/Xml/" )  ;
    _xmlPathTrue  = TString(_xmlPath  + "/true/" );
    _xmlPathInit  = TString(_xmlPath  + "/init/" );
    _xmlPathFinal = TString(_xmlPath  + "/final/");
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

    _outputFile = new TFile(_dirPath + "/" + _outputFileName, "RECREATE");
    _outputTree = new TTree("ExperimentTree", "ExperimentTree");

    double val;

    double par_counter = 0;

    for (int i = 0; i<_parameters.size(); i++) {
        FitParameters* pars;
        pars = _parameters.at(i);
    
        for (int j=0; j<pars->nparameters(); j++){

            _fitted_values.push_back(new double);
            _true_values.push_back(new double);
            _initial_values.push_back(new double);
            _errorLow_values.push_back(new double);
            _errorHigh_values.push_back(new double);
            _error_values.push_back(new double);
            _pull_values.push_back(new double);
            _residual_values.push_back(new double);

            parameter* par;
            par = pars->get_parameter(j);

            TString par_name(par->get_name());

            _outputTree->Branch(par_name+TString("_fitted"), (_fitted_values.at(par_counter)), par_name+TString("_fitted")+TString("/D"));
            _outputTree->Branch(par_name+TString("_true"), (_true_values.at(par_counter)), par_name+TString("_true")+TString("/D"));
            _outputTree->Branch(par_name+TString("_initial"), (_initial_values.at(par_counter)), par_name+TString("_initial")+TString("/D"));
            _outputTree->Branch(par_name+TString("_errorLow"), (_errorLow_values.at(par_counter)), par_name+TString("errorLow")+TString("/D"));
            _outputTree->Branch(par_name+TString("_errorHigh"), (_errorHigh_values.at(par_counter)), par_name+TString("errorHigh")+TString("/D"));
            _outputTree->Branch(par_name+TString("_error"), (_error_values.at(par_counter)), par_name+TString("error")+TString("/D"));
            _outputTree->Branch(par_name+TString("_pull"), (_pull_values.at(par_counter)), par_name+TString("pull")+TString("/D"));
            _outputTree->Branch(par_name+TString("_residual"), (_residual_values.at(par_counter)), par_name+TString("residual")+TString("/D"));

            par_counter += 1;
        }
    }
}

bool ToyExperiment::buildDirectoryTree(){
    int idx = 1;
    struct stat info;

    setOutputDirectory(_dirPath) ;

    TString aux = _dirPath;
    while(stat(aux.Data(), &info) == 0){
        aux = _dirPath + std::to_string(idx++);
        if (idx > 100){
            std::cout << "ERROR: more than 100 directory with the same name exists" << std::endl;
            return false;
        }
    }
    setOutputDirectory(_dirPath = aux);
    
    if (mkdir(_dirPath.Data(), S_IRWXU) != 0){
        return false;
    }
    
    mkdir((_xmlPath     ).Data(), S_IRWXU);
    mkdir((_xmlPathTrue ).Data(), S_IRWXU);
    mkdir((_xmlPathInit ).Data(), S_IRWXU);
    mkdir((_xmlPathFinal).Data(), S_IRWXU);

    createOutputFile();

    return true;
}

void ToyExperiment::run() {

    if (! buildDirectoryTree()) return ;

    std::vector<ToyModule*>::iterator it;

    setup();

    for (int i=0; i<_nRepetitions; ++i) {

        _currentRepetition = i;
        
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

    closing();
    
    _outputFile->cd();
    _outputTree->Write();
    _outputFile->Close();

}

void ToyExperiment::addModule(ToyModule& module) {

    module.setReferenceToExperiment(this); 
    _modules.push_back(&module);    
}

void ToyExperiment::addParameters(FitParameters& parameters) {

    _parameters.push_back(&parameters);
    _nParams += parameters.nparameters();

}

void ToyExperiment::initialize() {


    double par_counter = 0;
    for (int i = 0; i<_parameters.size(); i++){
        FitParameters* pars;
        pars = _parameters.at(i);

        for (int j=0; j<pars->nparameters(); j++){

            parameter* par;
            par = pars->get_parameter(j);

            *(_initial_values.at(par_counter)) = par->get_start_value();
            par_counter += 1;
        }
    }

    return;
}

void ToyExperiment::save() {

    double par_counter = 0;
    for (int i = 0; i<_parameters.size(); i++){
        FitParameters* pars;
        pars = _parameters.at(i);

        for (int j=0; j<pars->nparameters(); j++){
            parameter* par;
            par = pars->get_parameter(j);

            *(_fitted_values.at(par_counter))    = par->get_value();
            *(_true_values.at(par_counter))    = par->get_true_value();
            *(_error_values.at(par_counter))     = std::abs(par->get_error());
            *(_errorHigh_values.at(par_counter)) = std::abs(par->get_error_up());
            *(_errorLow_values.at(par_counter))  = std::abs(par->get_error_down());
            *(_residual_values.at(par_counter))  = par->get_value() - par->get_true_value();
            double pull;
            if (par->get_value() > par->get_true_value()) pull = (par->get_value() - par->get_true_value())/par->get_error_down();
            else                                           pull = (par->get_value() - par->get_true_value())/par->get_error_up();
            *(_pull_values.at(par_counter))      = pull;
        
            par_counter += 1;
        }
    }

    _outputFile->cd();
    _outputTree->Fill();

    return;
}

int ToyExperiment::getNParameters() {
    return _nParams;
}

std::vector<FitParameters*> ToyExperiment::getParameters() {
    return _parameters;
}

char* ToyExperiment::getOutputDirectory() const { return const_cast<char*> (_dirPath.Data()) ; }

int ToyExperiment::getNRepetitions(){

    return _nRepetitions;
}

int ToyExperiment::getNData(){

    return _nMC;
}


int ToyExperiment::getNMC(){
    
    return _nData;
}

parameter* ToyExperiment::getParameter(const char* name){
    parameter* par = NULL; 
    
    for (int i=0; i<_parameters.size(); ++i){
        FitParameters* pars = _parameters.at(i);
        par = pars->get_parameter(name);
        if (par != NULL){
            return par;
        }
    }
    return par;
}

void ToyExperiment::setParameterValue(const char* name, double v){
    parameter* par = this->getParameter(name);
    if (par != NULL){
        par->set_value(v);
        return;
    }
    std::cout << "ERROR: Parameter " << name << " not found. Not setting its value!" << std::endl;
}
void ToyExperiment::setParameterStartValue(const char* name, double v){
    parameter* par = this->getParameter(name);
    if (par != NULL){
        par->set_start_value(v);
        return;
    }
    std::cout << "ERROR: Parameter " << name << " not found. Not setting its start value!" << std::endl;
}
void ToyExperiment::setParameterTrueValue(const char* name, double v){
    parameter* par = this->getParameter(name);
    if (par != NULL){
        par->set_true_value(v);
        return;
    }
    std::cout << "ERROR: Parameter " << name << " not found. Not setting its true value!" << std::endl;
}
void ToyExperiment::setParameterErrorHigh(const char* name, double v){
    parameter* par = this->getParameter(name);
    if (par != NULL){
        par->set_error_up(v);
        return;
    }
    std::cout << "ERROR: Parameter " << name << " not found. Not setting its error up value!" << std::endl;
}
void ToyExperiment::setParameterErrorLow(const char* name, double v){
    parameter* par = this->getParameter(name);
    if (par != NULL){
        par->set_error_down(v);
        return;
    }
    std::cout << "ERROR: Parameter " << name << " not found. Not setting its error down value!" << std::endl;
}

double ToyExperiment::getParameterValue(const char* name){
    parameter* par = this->getParameter(name);
    if (par != NULL){
        return par->get_value();
    }
    std::cout << "ERROR: Parameter " << name << " not found. Returning 0 as its value!" << std::endl; 
    return 0.;
}
double ToyExperiment::getParameterStartValue(const char* name){
    parameter* par = this->getParameter(name);
    if (par != NULL){
        return par->get_start_value();
    }
    std::cout << "ERROR: Parameter " << name << " not found. Returning 0 as its start value!" << std::endl; 
    return 0.;
}
double ToyExperiment::getParameterTrueValue(const char* name){
    parameter* par = this->getParameter(name);
    if (par != NULL){
        return par->get_true_value();
    }
    std::cout << "ERROR: Parameter " << name << " not found. Returning 0 as its true value!" << std::endl; 
    return 0.;
}
double ToyExperiment::getParameterErrorHigh(const char* name){
    parameter* par = this->getParameter(name);
    if (par != NULL){
        return par->get_error_up();
    }
    std::cout << "ERROR: Parameter " << name << " not found. Returning 0 as its error up value!" << std::endl; 
    return 0.;
}
double ToyExperiment::getParameterErrorLow(const char* name){
    parameter* par = this->getParameter(name);
    if (par != NULL){
        return par->get_error_down();
    }
    std::cout << "ERROR: Parameter " << name << " not found. Returning 0 as its error down value!" << std::endl; 
    return 0.;
}
