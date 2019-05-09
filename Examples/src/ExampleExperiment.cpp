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
            par->set_value(0.2);
            par->set_error_up(0.2);
            par->set_error_down(0.2);
            par->set_error(0.2);
        }
        par = pars->get_parameter("delta_R1");
        if (par != NULL) {
            par->set_value(0.4);
            par->set_error_up(0.4);
            par->set_error_down(0.4);
            par->set_error(0.4);
        }
        par = pars->get_parameter("delta_R2");
        if (par != NULL) {
            par->set_value(0.6);
            par->set_error_up(0.6);
            par->set_error_down(0.6);
            par->set_error(0.6);
        }
        par = pars->get_parameter("delta_R0");
        if (par != NULL) {
            par->set_value(0.8);
            par->set_error_up(0.8);
            par->set_error_down(0.8);
            par->set_error(0.8);
        }
    }

    return;
}