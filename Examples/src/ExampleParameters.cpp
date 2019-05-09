#include "ExampleParameters.hpp"
#include "Parameter.hpp"

ExampleParameters::ExampleParameters() {
    
    parameter* delta_RhoSq = new parameter();
    delta_RhoSq->set_name("delta_RhoSq");
    delta_RhoSq->set_start_value(0.1);
    delta_RhoSq->set_min(-0.5);
    delta_RhoSq->set_max(0.5);
    delta_RhoSq->set_value(0.);

    parameter* delta_R1 = new parameter();
    delta_R1->set_name("delta_R1");
    delta_R1->set_start_value(0.2);
    delta_R1->set_min(-0.5);
    delta_R1->set_max(0.5);
    delta_R1->set_value(0.);

    parameter* delta_R2  = new parameter();
    delta_R2->set_name("delta_R2");
    delta_R2->set_start_value(0.3);
    delta_R2->set_min(-0.5);
    delta_R2->set_max(0.5);
    delta_R2->set_value(0.);

    parameter* delta_R0 = new parameter();
    delta_R0->set_name("delta_R0");
    delta_R0->set_start_value(0.4);
    delta_R0->set_min(-0.5);
    delta_R0->set_max(0.5);
    delta_R0->set_value(0.);

    addParameter(delta_RhoSq);
    addParameter(delta_R1);
    addParameter(delta_R2);
    addParameter(delta_R0);

}

ExampleParameters::~ExampleParameters() {

}