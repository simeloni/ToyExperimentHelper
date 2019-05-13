#ifndef EXAMPLEEXPERIMENT_HH
#define EXAMPLEEXPERIMENT_HH

#include "ToyExperiment.hpp"

class ExampleExperiment : public ToyExperiment {

public:

    void fit();

    ExampleExperiment();
    ~ExampleExperiment();

    int getAdditional(); 
    void setAdditional(double value);

private: 
    //Additional private parameters
    int _additional;
};

#endif