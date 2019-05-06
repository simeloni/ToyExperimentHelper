#include "ToyModule.hpp"

ToyModule::ToyModule() {
    _name = TString("ToyModule");
    _MCExperiment = NULL;
}

ToyModule::ToyModule(const ToyModule& anotherModule) {

    _name         = anotherModule._name;
    _MCExperiment = anotherModule._MCExperiment;

}

ToyModule::~ToyModule() {

}

void ToyModule::setName(TString name) {
    _name = name;
}

void ToyModule::setReferenceToExperiment(ToyExperiment* MCExp) {
        if (MCExp == NULL) {
            std::cout << "ERROR: in setReferenceToExperiment, ToyModule.hpp. Null pointer to ToyMCExperiment" << std::endl;
            exit(EXIT_FAILURE);
        }
        _MCExperiment = MCExp;
}