#include "ExampleParameters.hpp"
#include "Parameter.hpp"

ExampleParameters::ExampleParameters() {
    
    //Initialize the parameters
    parameter* par_of_interest = new parameter();
    par_of_interest->set_name("par_of_interest");
    par_of_interest->set_min(0.);
    par_of_interest->set_max(10.);
    par_of_interest->set_start_value(5.);
    par_of_interest->set_step_size(0.005);

    parameter* par_nuisance = new parameter();
    par_nuisance->set_name("par_nuisance");
    par_nuisance->set_min(0.);
    par_nuisance->set_max(1.);
    par_nuisance->set_start_value(0.5);
    par_nuisance->set_step_size(0.0005);

    // add the parameters to the vector
    addParameter(par_of_interest);
    addParameter(par_nuisance);
}

ExampleParameters::~ExampleParameters() {

}