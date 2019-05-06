#ifndef TOYMODULE_HH
#define TOYMODULE_HH

#include <iostream>
#include <stdlib.h>
#include "TString.h"

class ToyExperiment;
class ToyModule {

friend ToyExperiment;

public:
    ToyModule();
    ToyModule(const ToyModule &anotherModule); //Copy constructor
    ~ToyModule();

    void setName(TString name);

    virtual void beforeInitialize(){std::cout << "Running Module " << _name << ": before Initialize" << std::endl;};
    virtual void beforeGenerate(){std::cout << "Running Module " << _name << ": before Generate" << std::endl;};
    virtual void beforeFit(){std::cout << "Running Module " << _name << ": before Fit" << std::endl;};
    virtual void beforeSave(){std::cout << "Running Module " << _name << ": before Save" << std::endl;};
    virtual void afterSave(){std::cout << "Running Module " << _name << ": after Save" << std::endl;};

private:
    void setReferenceToExperiment(ToyExperiment* MCExp);

private:
    TString        _name;
    ToyExperiment* _MCExperiment;
};

#endif