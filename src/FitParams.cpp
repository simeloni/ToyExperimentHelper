#ifndef FITPARAMS_H
#define FITPARAMS_H

#include <string>
#include <iomanip>
#include <iostream>
#include <vector>
#include "FitParams.hpp"
//#include "funcs.hh"

class parameters;

/// The parameter class stores the different attributes of a parameter used in the analysis, e.g. the name, start value and the parameter range.
class parameter {
  friend class parameters;

private:
  ///blind this parameter
  bool blind;
  ///blinding string     
  std::string blinding_string;
  ///blinding delta, never look at this
  double blinding_delta;
  //blinding range
  double blinding_scale;
  //This name is supplied as the parameter name to minuit. Please give every parameter a different name.                                                                                            
  std::string name;
  ///The parameter description is used by the plotter for plot captions and axis titles.                                                                                                             
  std::string description;
  ///The start value of the parameter. After every step in the toy study, the parameter is reset to this value.                                                                                      
  double start_value;
  ///The current value of the parameter                                                                                                                                                              
  double value;
  ///The minimum possible value for the parameter                                                                                                                                                    
  double min;
  ///The maximum possible value for the parameter                                                                                                                                                    
  double max;
  ///The initial step size used by minuit                                                                                                                                                            
  double step_size;
  ///The error of the fitted result                                                                                                                                                                  
  double error;
  ///asymmetric errors if a minos error analysis has been performed                                                                                                                                  
  double error_up;
  ///asymmetric errors if a minos error analysis has been performed                                                                                                                                  
  double error_down;
  ///the index of this parameter in the params array of the parameter set                                                                                                                            
  int index;
  ///should the parameter float without limits?                                                                                                                                                      
  bool unlimited;
  ///this has been the error of a previous measurement                                                                                                                                               
  double previous_error;
  ///previous measurement, this can be different from the start value, e.g. when running a toy study                                                                                                 
  double previous_measurement;
  ///do we want to vary the parameter so that it is compatible with an earlier measurement                                                                                                           
  bool gaussian_constraint;

public:
  //correlations                                                                                                                                                                                     
  std::vector<double> correlations;
  ///public set routines                 
  void set_name(std::string n) {name = n;};
  void set_min(double m) {min = m;};
  void set_max(double m) {max = m;};
  void set_step_size(double s) {step_size = s;};
  void set_value(double v) {value = v;};
  void set_start_value(double v) {start_value = v;};
  void set_description(std::string d) {description = d;};
  void set_error(double e) {error = e;};
  void set_error_down(double e) {error_down = e;};
  void set_error_up(double e) {error_up = e;};
  void set_unlimited(bool u) {unlimited = u;};
  void set_index(int idx) {index = idx;};
  void set_previous_measurement(double m) {previous_measurement = m;};
  ///set blinding string                                                                                                                                                                             
  void set_blinding_string(std::string s) {
    blinding_string = s;
  };
    void set_blinding(bool b, double b_scale, bool is_angle=false, std::string b_string="DefaultBlinder") {
    blind = b;
    blinding_scale = b_scale;
    blinding_string = b_string;
    //this is new, apparently the shi(f)t needs to be subtracted, who knew?                                                                                                                          
    if (!is_angle)
      blinding_delta = -evaluate_unblind_uniform(0.0, blinding_string.c_str(), blinding_scale);
    else
      blinding_delta = -evaluate_unblind_uniform_angle(0.0, blinding_string.c_str(), blinding_scale);
  };
  ///set all values which make sense on initialisation                                                                                                                                               
  void init(std::string n, std::string d, double v, double min, double max)
  {
    set_name(n);
    set_description(d);
    set_value(v);
    start_value = v;
    set_min(min);
    set_max(max);
    set_error(0.0);
    set_error_up(0.0);
    set_error_down(0.0);
    set_unlimited(false);
    gaussian_constraint = false;
    previous_measurement = v;
    correlations.clear();
  };
  void init(std::string n, std::string d, double v, double min, double max, double stepsize)
  {
    set_name(n);
    set_description(d);
    set_value(v);
    start_value = v;
    set_min(min);
    set_max(max);
    set_step_size(stepsize);
    set_error(0.0);
    set_error_up(0.0);
    set_error_down(0.0);
    set_unlimited(false);
    gaussian_constraint = false;
    previous_measurement = v;
    correlations.clear();
  };
  void init(std::string n, std::string d, double v, double min, double max, double stepsize, double previous_measurement_error)
  {
    init(n, d, v, min, max, stepsize);
    previous_error = previous_measurement_error;
    gaussian_constraint = true;
  }
  void init(std::string n, std::string d, double v, double min, double max, double stepsize, bool nolimits)
  {
    init(n, d, v, min, max, stepsize);
    unlimited = nolimits;
  };
  ///public get routines
    std::string get_name() const {return name;};
  const char* get_mn_name() const {return get_name().c_str();};
  double get_min() const {return min;};
  double get_max() const {return max;};
  double get_blinding_delta() const {return blinding_delta;};
  double get_step_size() const {return step_size;};
  double get_value() const {return value;};
  double get_start_value() const {return start_value;};
  double get_previous_error() const {return previous_error;};
  double get_previous_measurement() const {return previous_measurement;};
  std::string get_description() const {return description;};
  bool get_unlimited() const {return unlimited;};
  bool get_gaussian_constraint() const {return gaussian_constraint;};
  std::string get_root_name() const {
    std::string rootstring = get_description();
    size_t j;
    while ((j = rootstring.find("\\mathrm ")) != std::string::npos)
      {
        rootstring.replace(j,8,"");
      }
    while ((j = rootstring.find("\\")) != std::string::npos)
      {
        rootstring.replace(j,1,"#");
      }
    return rootstring;
  };
  double operator() () const {return value;};
  double get_error() const {return error;};
  double get_error_up() const {return error_up;};
  double get_error_down() const {return error_down;};
  int get_index() const {return index;};
  void reset_start() {value = start_value; error = 0.0;};
  ///default constructor                                                                                                                                                                             
  parameter():
    blind(false),
    blinding_string("DefaultBlinder"),
    blinding_delta(0.0),
    name("default"),
    description("default"),
    start_value(0.0),
    value(0.0),
    min(-1.0),
    max(+1.0),
    step_size(0.001),
    error(0.0),
    error_up(0.0),
    error_down(0.0),
    index(-1),
    unlimited(false),
    previous_error(0.0),
    previous_measurement(0.0),
    gaussian_constraint(false)
  {};
  ///normal constructor
    parameter(std::string n, std::string d, double v,
            double minimum, double maximum, double stepsize = 0.0, bool nolimits=false):
    blind(false),
    blinding_string("DefaultBlinder"),
    blinding_delta(0.0),
    name(n),
    description(d),
    start_value(v),
    value(v),
    min(minimum),
    max(maximum),
    step_size(stepsize),
    error(0.0),
    error_up(0.0),
    error_down(0.0),
    index(-1),
    unlimited(nolimits),
    gaussian_constraint(false)
  {};
  bool is_blind() {
    return blind;
  };
};

///The parameters class stores a collection of parameters used in the analysis. Every analysis should implement its own set e.g. bs2jpsiphi_parameters                                               
class parameters
{
private:
  ///this vector holds the different parameters                                                                                                                                                      
  std::vector<parameter*> params;
protected:
  ///Add a certain parameter to the params vector.                                                                                                                                                   
  void add_parameter(parameter* param) {
    param->set_index(params.size());
    params.push_back(param);
  };
public:
  bool is_blind() const {
    for (unsigned int i = 0; i<params.size(); i++)
      if (params.at(i)->is_blind())
        return true;
    return false;
  };
    virtual ~parameters() {};
  ///Method gives back he number of parameters                                                                                                                                                       
  unsigned int nparameters() const
  {
    return params.size();
  };
  ///Methos returns a pointer to parameter with index i                                                                                                                                              
  parameter* get_parameter(unsigned int i) const
  {
    return params.at(i);
  };
  ///returns parameter with specific name                                                                                                                                                            
  parameter* get_parameter(std::string name) const
  {
    parameter* res=0;
    for (unsigned int i =0; i<nparameters();i++)
      if (get_parameter(i)->get_name() == name)
      {
        res = get_parameter(i);
        break;
      }
    return res;
  };
  ///Resets all parameters to their start values                                                                                                                                                     
  void reset_parameters()
  {
    for (unsigned int i = 0; i < nparameters(); i++)
      params.at(i)->reset_start();
  };
  ///Method that prints the parameter set, the current values, errors and the deviation from the start value in sigma
    void print_parameters(bool latex_output=true) const
  {
    if (!latex_output)
    {
      std::cout << std::endl << "Parameters" << std::endl;

      std::cout << "Blinded Parameters: ";
      for (unsigned int j = 0; j < nparameters(); j++)
        if (params[j]->is_blind())
          std::cout << params[j]->get_name() << " ";
      std::cout << std::endl;

      for (unsigned int j = 0; j < nparameters(); j++)
      {
        std::string parname(params[j]->get_name());//the_minimizer->fCpnam[j]);                                                                                                                      
        std::cout.setf(std::ios::left);
        double nsigma = 0.0;
        if (params[j]->get_error() > 0.0)
          nsigma = (params[j]->get_value() - params[j]->get_start_value())/params[j]->get_error();
        double value = params[j]->get_value();
        if (params[j]->is_blind())
          value += params[j]->blinding_delta;//todo                                                                                                                                                  
        double error = params[j]->get_error();
        double error_up = params[j]->get_error_up();
        double error_down = params[j]->get_error_down();
        if (error_up == error_down)
          std::cout << std::setw(22) << parname
            //<< std::setw(12) << (params[j]->is_blind() ? "(*)" : format_value(value, error))                                                                                                       
                    << std::setw(12) << (format_value(value, error))
                    << std::setw(12) << "+-" << format_error(value, error)
                    << (params[j]->is_blind() ? " (*)": "")
                    << " " << (params[j]->is_blind() ? "(*)" : format_double(nsigma)) << " sigma" << std::endl;
        else
          std::cout << std::setw(22) << parname
            //<< std::setw(12) << (params[j]->is_blind() ? "(*)" : format_value(value, error))                                                                                                       
                    << std::setw(12) << (format_value(value, error))
                    << std::setw(12) << " " << format_error(value, error_down)
                    << " " << format_error(value, error_up)
                    << (params[j]->is_blind() ? "(*)": "")
                    << " " << (params[j]->is_blind() ? "(*)" : format_double(nsigma)) << " sigma" << std::endl;
      }
    }
    else
        {
                std::cout << "Parameters" << std::endl;
                std::cout << "Blinded Parameters: ";
                for (unsigned int j = 0; j < nparameters(); j++)
                        if (params[j]->is_blind())
                                std::cout << params[j]->get_name() << " ";
                std::cout << std::endl;

                if (params[0]->get_error_up() == params[0]->get_error_down())
                {
                        std::cout << "\\begin{tabular}{|c|c|c|} \\hline" << std::endl;
                        std::cout << "parameter & result & $\\sigma$ from nominal \\\\ \\hline \\hline" << std::endl;
                }
                else
                {
                        std::cout << "\\begin{tabular}{|c|c|c|c|c|} \\hline" << std::endl;
                        std::cout << "parameter & result & $\\sigma_\\text{down}$ & $\\sigma_\\text{up}$ & $\\sigma$ from nominal \\\\ \\hline \\hline" << std::endl;
                }
                for (unsigned int j = 0; j < nparameters(); j++)
                {
                        std::string pardesc(params[j]->get_description());//the_minimizer->fCpnam[j]);                                                                                               
                        std::string partex("");
                        for (std::string::iterator it = pardesc.begin(); it != pardesc.end(); it++)
                                if (*it != '#')
                                        partex += *it;
                                else
                                        partex += "\\";
                        double value = params[j]->get_value();
                        if (params[j]->is_blind())
                                value += params[j]->blinding_delta;//todo                                                                                                                            
                        double error = params[j]->get_error();
                        double error_up = params[j]->get_error_up();
                        double error_down = params[j]->get_error_down();
                        double nsigma = 0.0;
                        if (error > 0.0)
                                //nsigma = params[j]->get_start_value();                                                                                                                             
                                nsigma = (value - params[j]->get_start_value())/error;
                        std::cout.setf(std::ios::left);

                        if (error_up == error_down)
                                std::cout << "$" << std::setw(30) << partex << "$ & $"
                                        //<< std::setw(12) << (params[j]->is_blind() ? "(*)" : format_value(value, error))                                                                           
                                        << std::setw(12) << ((params[j]->get_step_size() == 0.0 && params[j]->is_blind()) ? "" : format_value(value, error))
                                        << "\\pm " << std::setw(12) << format_error(value, error) << (params[j]->is_blind() ? "(*)": "") << "$ & "
                                        << (params[j]->is_blind() ? "(*)" : format_double(nsigma)) << "\\\\" << std::endl;
                        else
                                std::cout << "$" << std::setw(30) << partex << "$ & $"
                                        << std::setw(12)
                                        //<< (params[j]->is_blind() ? "(*)" : format_value(value, error))                                                                                            
                                        << (format_value(value, error))
                                        << "\\pm "
                                        << std::setw(12) << format_error(value, error) << (params[j]->is_blind() ? "(*)": "") << "$ & "
                                        << std::setw(12) << format_error(value, error_down) << " & "
                                        << std::setw(12) << "+" << format_error(value, error_up) << " & "
                                        << (params[j]->is_blind() ? "(*)" : format_double(nsigma)) << "\\\\" << std::endl;
               }
                std::cout << "\\hline" << std::endl;
                std::cout << "\\end{tabular}" << std::endl;
        }
  };
  ///fix all parameters                                                                                                                                                                              
  void fix_parameters() {
    for (unsigned int i = 0; i < params.size(); i++)
        params.at(i)->set_step_size(0.0);
  };
  ///fix desired parameter
    void fix_param(std::string tobefixed){
    for(unsigned int i = 0; i < params.size(); i++){
      std::string parname =params.at(i)->get_name();
        if(tobefixed==parname)
          params.at(i)->set_step_size(0.0);
    }
  };
  ///set all parameters to current values! Careful with this, please.                                                                                                                                
  void take_current_as_start() {
    for (unsigned int i = 0; i < params.size(); i++)
      {
        parameter* p = params.at(i);
        p->set_start_value(p->get_value());
      }
  };

};


#endif
