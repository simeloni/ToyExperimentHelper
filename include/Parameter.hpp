#ifndef PARAMETER_HH
#define PARAMETER_HH

#include <string>
#include <iomanip>
#include <iostream>
#include <vector>

class parameters;
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
  void set_blinding(bool b, double b_scale, bool is_angle=false, std::string b_string="DefaultBlinder");
  ///set all values which make sense on initialisation                                                                                                                                               
  void init(std::string n, std::string d, double v, double min, double max);
  void init(std::string n, std::string d, double v, double min, double max, double stepsize);
  void init(std::string n, std::string d, double v, double min, double max, double stepsize, double previous_measurement_error);
  void init(std::string n, std::string d, double v, double min, double max, double stepsize, bool nolimits);
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
  std::string get_root_name() const;
  double operator() () const {return value;};
  double get_error() const {return error;};
  double get_error_up() const {return error_up;};
  double get_error_down() const {return error_down;};
  int get_index() const {return index;};
  void reset_start() {value = start_value; error = 0.0;};
  bool is_blind() {return blind;};
};

#endif