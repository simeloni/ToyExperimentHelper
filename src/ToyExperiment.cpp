#include "ToyExperiment.hpp"
#include "TFile.h"
#include "TBranch.h"
#include "TTree.h"
#include "TString.h"

ToyExperiment::ToyExperiment() {
    _nData = 10;
    _nMC   = 10;
    _nRepetitions = 0;
    _outputFileName = TString("outputFile.root");
    _nParams = 0;
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
    _outputTree = new TTree("ExperimentTree", "ExperimentTree");

    double val;

    double par_counter = 0;

    for (int i = 0; i<_parameters.size(); i++) {
        FitParameters* pars;
        pars = _parameters.at(i);
    
        for (int j=0; j<pars->nparameters(); j++){

            _fitted_values.push_back(new double);
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
    
    _outputFile->cd();
    _outputTree->Write();
    _outputFile->Close();

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

void ToyExperiment::addParameters(FitParameters* parameters) {

    if (parameters != NULL) {
        _parameters.push_back(parameters);
        _nParams += parameters->nparameters();
    }
    else {
        std::cout << "WARNING: Empty parameters passed to ToyExperiment" << std::endl;
    }
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
            *(_error_values.at(par_counter))     = par->get_error();
            *(_errorHigh_values.at(par_counter)) = par->get_error_up();
            *(_errorLow_values.at(par_counter))  = par->get_error_down();
            *(_residual_values.at(par_counter))  = par->get_value() - par->get_start_value();
            double pull;
            if (par->get_value() > par->get_start_value()) pull = (par->get_value() - par->get_start_value())/par->get_error_up();
            else                                           pull = (par->get_value() - par->get_start_value())/par->get_error_down();
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