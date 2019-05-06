#include "FitParameters.hpp"
#include <iostream>

FitParameters::FitParameters() {

}

FitParameters::~FitParameters() {

}

void FitParameters::addParameter(parameter* param) {
    if (param == NULL) {
        std::cout << "ERROR: in addParameter, FitParameters class. The pointer is pointing to NULL. Not adding any parameter " << std::endl;
    }

    else {
        param->set_index(_parameters.size());
        _parameters.push_back(param);
    }
}

bool FitParameters::is_blind() {

    for (unsigned int i = 0; i<_parameters.size(); i++) {
        if (_parameters.at(i)->is_blind()) return true;
    }
    
    return false;
}
