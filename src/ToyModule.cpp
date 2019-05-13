#include "ToyModule.hpp"
#include <string>

ToyModule::ToyModule() {
    _name = std::string("ToyModule");
    _MCExperiment = NULL;
}

ToyModule::ToyModule(const ToyModule& anotherModule) {
    _name         = anotherModule._name;
    _MCExperiment = anotherModule._MCExperiment;
}

ToyModule::~ToyModule() {

}

void ToyModule::setName(std::string name) {
    _name = name;
}

void ToyModule::setReferenceToExperiment(ToyExperiment* MCExp) {
    if (MCExp == NULL) {
        std::cout << "ERROR: in setReferenceToExperiment, ToyModule.hpp. Null pointer to ToyMCExperiment" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (getExperiment()) { 
        std::cout << "ERROR: in setReferenceToExperiment, ToyModule.hpp. Module already used from another experiment. Please use another instance!" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    _MCExperiment = MCExp;
}

std::string ToyModule::getName() {
    return _name;
}

ToyExperiment* ToyModule::getExperiment() {
    return _MCExperiment;
}

std::vector<FitParameters*>  ToyModule::getExperimentParameters() {
    return _MCExperiment->getParameters();
}