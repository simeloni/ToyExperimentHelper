#include "ExampleModule.hpp"
#include "ExampleExperiment.hpp"
#include <vector>
#include "ExampleParameters.hpp"


ExampleModule::ExampleModule() {
    generator = new TRandom();
    generator->SetSeed(42);
}

ExampleModule::~ExampleModule() {

}

void ExampleModule::beforeInitialize() {
    std::cout << "Running Derived module " << getName() << ": Before initialize" << std::endl;
    
    //This is a good place to float the parameters randomly for example 
    //Access the data of the experiment
    
    std::vector<FitParameters*> parameters = getExperimentParameters();
    for (unsigned int i=0; i<parameters.size(); i++){
        
        //If you are sure the parameters are of the type you think, cast them
        ExampleParameters* params = (ExampleParameters*) parameters.at(i);
        std::cout << params << std::endl;  

        //Loop over the parameters and add a 10% fluctuation  
        for (unsigned j=0; j<params->nparameters(); j++) {
            parameter* par = params->get_parameter(j);
            double val = par->get_value();
            par->set_start_value(val + generator->Gaus(0, 0.10*val));
        
        }
    }
}


void ExampleModule::beforeGenerate(){
    std::cout << "Running Derived module " << getName() <<": Before Generate" << std::endl; 

    std::cout << "---Reading the protected members of experiment from the example module---" << std::endl;

    std::cout << "---Accessing the Experiment private members using the standard getters of the Module---" << std::endl;
    ((ExampleExperiment*)getExperiment())->setAdditional(42.);
    std::cout << ((ExampleExperiment*)getExperiment())->getAdditional() << std::endl;

    std::cout << "---Access the data of the experiment" << std::endl;
    std::cout << "The experiment has " << getExperimentParameters().size() << " sets of parameters" << std::endl;
}

ExampleModule::ExampleModule(const ExampleModule& another){
    
}