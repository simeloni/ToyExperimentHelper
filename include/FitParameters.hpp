#ifndef FITPARAMETERS_HH
#define FITPARAMETERS_HH

#include <vector>
#include "Parameter.hpp"

class FitParameters {

public:
    FitParameters();
    ~FitParameters();

    bool is_blind(); //return true if at least one of the parameters is blind

protected:
    //This method is protected since I want the specializations of these classed
    // to use it only internally!
    void addParameter(parameter* param);


private: 
    std::vector<parameter*> _parameters;
};


#endif 