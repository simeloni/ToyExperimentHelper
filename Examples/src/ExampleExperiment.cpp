#include "ExampleExperiment.hpp"
#include "FitParameters.hpp"
#include "Parameter.hpp"

ExampleExperiment::ExampleExperiment() {
}

ExampleExperiment::~ExampleExperiment() {

}

void ExampleExperiment::fit() {

    //Setting the parameters and the parameters errors. You would like to do this at the end of your fitter
    for (int i=0; i<_parameters.size(); ++i) {
        FitParameters* pars = _parameters.at(i); 

        parameter* par;
        par = pars->get_parameter("delta_RhoSq");
        if (par != NULL) {
            par->set_value(0.1);
            par->set_error_up(0.1*0.10);
            par->set_error_down(0.1*0.10);
        }
        par = pars->get_parameter("delta_R1");
        if (par != NULL) {
            par->set_value(0.2);
            par->set_error_up(0.2*0.10);
            par->set_error_down(0.2*0.10);
        }
        par = pars->get_parameter("delta_R2");
        if (par != NULL) {
            par->set_value(0.3);
            par->set_error_up(0.3*0.10);
            par->set_error_down(0.3*0.10);
        }
        par = pars->get_parameter("delta_R0");
        if (par != NULL) {
            par->set_value(0.4);
            par->set_error_up(0.4*0.10);
            par->set_error_down(0.4*0.10);
        }
    }

    return;
}

int ExampleExperiment::getAdditional() {
    return _additional;
}

void ExampleExperiment::setAdditional(double value) {
    _additional = value;
}