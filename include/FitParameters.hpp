#ifndef FITPARAMETERS_HH
#define FITPARAMETERS_HH

#include <vector>
#include "Parameter.hpp"

class FitParameters {

public:
    FitParameters();
    ~FitParameters();

private: 
    std::vector<parameter*> _parameters;
};


#endif 