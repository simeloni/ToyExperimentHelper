#ifndef FITPARAMETERS_HH
#define FITPARAMETERS_HH

#include <vector>
#include "Parameter.hpp"

class FitParameters {

public:
    FitParameters();
    ~FitParameters();

    bool is_blind() const; //return true if at least one of the parameters is blind
    int nparameters() const; 
    parameter* get_parameter(unsigned int i) const;
    parameter* get_parameter(std::string name) const;
    void reset_parameters();
    void print_parameters(bool latex_output=true) const;
    void fix_parameters(); //fix all the parameters;
    void fix_param(std::string tobefixed);
    void take_current_as_start(); //set all parameters to current values.
    
//protected:
    //This method is protected since I want the specializations of these classed
    // to use it only internally!
    void addParameter(parameter* param);

private: 
    std::vector<parameter*> _parameters;
};


#endif 