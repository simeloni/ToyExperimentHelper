#ifndef TOYMODULE_HH
#define TOYMODULE_HH

#include <iostream>

class ToyMCExperiment;
class ToyModule {

friend ToyMCExperiment;

public:
    ToyModule();
    ~ToyModule();

    virtual void beforeInitialize(){std::cout << "Running Module: before Initialize" << std::endl;};
    virtual void beforeGenerate(){std::cout << "Running Module: before Generate" << std::endl;};
    virtual void beforeFit(){std::cout << "Running Module: before Fit" << std::endl;};
    virtual void beforeSave(){std::cout << "Running Module: before Save" << std::endl;};
    virtual void afterSave(){std::cout << "Running Module: after Save" << std::endl;};

private:
    void setReferenceToExperiment(*ToyMCExperiment MCExp){
        //TODO: Check everything is okay
        _MCExperiment = MCExp;
    };

private:
    ToyMCExperiment* _MCExperiment;


    };

};

#endif