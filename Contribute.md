#How to contribute
For one feature, open a branch and when you are happy, open a merge request.

# TODO list
##General coding style
 * [ ] Use CXX 11
 * [ ] One default constructor with standard values
 * [ ] One constructor with non standard parameters 'a la python'
 
##Parameters class
 * [ ] It has to be a vector of parameters
 * [ ] It has to initialize the vector of parameters wither with a addParameter(*parameter) or with a common I/O method (e.g. xml)
 * [ ] It has to save itself in the same format used by the I/O method

##ToyMCExperiment
 * [ ] It has to have a addParameters() method
 * [ ] Init must NOT be virtual
 * [ ] Init must initialize the I/O file and save to it the initial values of the parameters
 * [ ] (The structure of the I/O file should be:
  A tree for each type: True, Initial, Final.
  For each tree you save all the parameters properties.
  The properties are: value, error, is blinded, error up, error down, step size, is constrained, number of the experiment)
 * [ ] I/O should be a common Root file
 * [ ] Save must NOT be virtual. If someone wants to to something more to the I/O file, he should do it with another Module
 * [ ] The functions should have a signature which returns an error object and catches the errors in the run function
 * [ ] It must have some getParameters object, through which the modules can access to the parameters vector
 * [ ] The init() function should initialize the output folder structure
 * [ ] Everything should be executable even though it has not been filled
 * [ ] Write a DryRun, which executes everything without saving anything to disk

##ToyModule
* [ ] One of the private members of the module class should be the reference to the experiment to which it has been assigned
* [ ] The refence should be initialized when the module is added in the addModule function of ToyExperiment

##Example modules
* [ ] Plot Module
* [ ] Smearing module: defines some smearing factors for some of the parameters

