#include "Parameter.hpp"
#include "funcs.hh"

void parameter::set_blinding(bool b, double b_scale, bool is_angle, std::string b_string) {
    blind = b;
    blinding_scale = b_scale;
    blinding_string = b_string;
    //this is new, apparently the shi(f)t needs to be subtracted, who knew?                                                                                                                          
    if (!is_angle)
      blinding_delta = -evaluate_unblind_uniform(0.0, blinding_string.c_str(), blinding_scale);
    else
      blinding_delta = -evaluate_unblind_uniform_angle(0.0, blinding_string.c_str(), blinding_scale);
}

//set all values which make sense on initialisation                                                                                                                                               
void parameter::init(std::string n, std::string d, double v, double min, double max)
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
    blind = false;
    constant = false;
    previous_measurement = v;
    correlations.clear();
}

void parameter::init(std::string n, std::string d, double v, double min, double max, double stepsize)
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
    blind = false;
    constant = false;
    previous_measurement = v;
    correlations.clear();
}

void parameter::init(std::string n, std::string d, double v, double min, double max, double stepsize, double previous_measurement_error)
{
    init(n, d, v, min, max, stepsize);
    previous_error = previous_measurement_error;
    gaussian_constraint = false;
    constant = false;
    blind = false;
}

void parameter::init(std::string n, std::string d, double v, double min, double max, double stepsize, bool nolimits)
{
    init(n, d, v, min, max, stepsize);
    unlimited = nolimits;
}

std::string parameter::get_root_name() const {
  
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
}