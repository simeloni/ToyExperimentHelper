#ifndef EXAMPLEMODULE_HH
#define EXAMPLEMODULE_HH

#include "ToyModule.hpp"
#include "TRandom.h"

class ExampleModule : public ToyModule{

public:

    ExampleModule(const ExampleModule&);
    ExampleModule();
    ~ExampleModule();

    void beforeGenerate();
    void beforeInitialize();

private:
    TRandom* generator;

};

#endif