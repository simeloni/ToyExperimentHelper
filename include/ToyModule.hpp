#ifndef TOYMODULE_HH
#define TOYMODULE_HH

#include <iostream>
#include <stdlib.h>
#include <string>
#include "ToyExperiment.hpp"
#include <vector>
#include "FitParameters.hpp" 

class ToyExperiment;
class ToyModule {

friend ToyExperiment;

public:
    ToyModule();
    ToyModule(const ToyModule&);
    ~ToyModule();

    virtual void setName(std::string name);

    virtual void beforeInitialize(){;};
    virtual void beforeGenerate(){;};
    virtual void beforeFit(){;};
    virtual void beforeSave(){;};
    virtual void afterSave(){;};

    std::string getName();

private:
    std::string          _name;
    ToyExperiment* _MCExperiment;

protected:
    void setReferenceToExperiment(ToyExperiment* MCExp);
    ToyExperiment* getExperiment(); //Return the experiment I am belonging to. 
    std::vector<FitParameters*> getExperimentParameters(); //Return the list of parameters added to the experiment you belong to 
};

#endif