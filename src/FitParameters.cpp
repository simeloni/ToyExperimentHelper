#include "FitParameters.hpp"
#include <iostream>
#include "funcs.hh"

FitParameters::FitParameters() {

}

FitParameters::~FitParameters() {

}

void FitParameters::addParameter(parameter* param) {
    if (param == NULL) {
        std::cout << "ERROR: in addParameter, FitParameters class. The pointer is pointing to NULL. Not adding any parameter " << std::endl;
    }

    else {
        param->set_index(_parameters.size());
        _parameters.push_back(param);
    }
}

bool FitParameters::is_blind() const {

    for (unsigned int i = 0; i<_parameters.size(); i++) {
        if (_parameters.at(i)->is_blind()) return true;
    }
    
    return false;
}

int FitParameters::nparameters() const {
    return _parameters.size();
}

parameter* FitParameters::get_parameter(unsigned int i) const {
    return _parameters.at(i);
}
    
parameter* FitParameters::get_parameter(std::string name) const {
    parameter* res=NULL;
    for (unsigned int i=0; i<nparameters(); i++){
        if (get_parameter(i)->get_name() == name) {
            res = get_parameter(i);
            break;
        }
    }
    return res;
}

void FitParameters::reset_parameters() {
    for (unsigned int i = 0; i < nparameters(); i++) {
        _parameters.at(i)->reset_start();
    }
}

void FitParameters::print_parameters(bool latex_output) const {
    if (!latex_output){
      std::cout << std::endl << "Parameters" << std::endl;

      std::cout << "Blinded Parameters: ";
      for (unsigned int j = 0; j < nparameters(); j++)
        if (_parameters[j]->is_blind())
          std::cout << _parameters[j]->get_name() << " ";
      std::cout << std::endl;

      for (unsigned int j = 0; j < nparameters(); j++)
      {
        std::string parname(_parameters[j]->get_name());                                                                                                                   
        std::cout.setf(std::ios::left);
        double nsigma = 0.0;
        if (_parameters[j]->get_error() > 0.0)
          nsigma = (_parameters[j]->get_value() - _parameters[j]->get_start_value())/_parameters[j]->get_error();
        double value = _parameters[j]->get_value();
        if (_parameters[j]->is_blind())
          value += _parameters[j]->blinding_delta;//todo                                                                                                                                                  
        double error = _parameters[j]->get_error();
        double error_up = _parameters[j]->get_error_up();
        double error_down = _parameters[j]->get_error_down();
        if (error_up == error_down)
          std::cout << std::setw(22) << parname
            //<< std::setw(12) << (_parameters[j]->is_blind() ? "(*)" : format_value(value, error))                                                                                                       
                    << std::setw(12) << (format_value(value, error))
                    << std::setw(12) << "+-" << format_error(value, error)
                    << (_parameters[j]->is_blind() ? " (*)": "")
                    << " " << (_parameters[j]->is_blind() ? "(*)" : format_double(nsigma)) << " sigma" << std::endl;
        else
          std::cout << std::setw(22) << parname
            //<< std::setw(12) << (_parameters[j]->is_blind() ? "(*)" : format_value(value, error))                                                                                                       
                    << std::setw(12) << (format_value(value, error))
                    << std::setw(12) << " " << format_error(value, error_down)
                    << " " << format_error(value, error_up)
                    << (_parameters[j]->is_blind() ? "(*)": "")
                    << " " << (_parameters[j]->is_blind() ? "(*)" : format_double(nsigma)) << " sigma" << std::endl;
      }
    }
    else
        {
                std::cout << "Parameters" << std::endl;
                std::cout << "Blinded Parameters: ";
                for (unsigned int j = 0; j < nparameters(); j++)
                        if (_parameters[j]->is_blind())
                                std::cout << _parameters[j]->get_name() << " ";
                std::cout << std::endl;

                if (_parameters[0]->get_error_up() == _parameters[0]->get_error_down())
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
                        std::string pardesc(_parameters[j]->get_description());//the_minimizer->fCpnam[j]);                                                                                               
                        std::string partex("");
                        for (std::string::iterator it = pardesc.begin(); it != pardesc.end(); it++)
                                if (*it != '#')
                                        partex += *it;
                                else
                                        partex += "\\";
                        double value = _parameters[j]->get_value();
                        if (_parameters[j]->is_blind())
                                value += _parameters[j]->blinding_delta;//todo                                                                                                                            
                        double error = _parameters[j]->get_error();
                        double error_up = _parameters[j]->get_error_up();
                        double error_down = _parameters[j]->get_error_down();
                        double nsigma = 0.0;
                        if (error > 0.0)
                                //nsigma = _parameters[j]->get_start_value();                                                                                                                             
                                nsigma = (value - _parameters[j]->get_start_value())/error;
                        std::cout.setf(std::ios::left);

                        if (error_up == error_down)
                                std::cout << "$" << std::setw(30) << partex << "$ & $"
                                        //<< std::setw(12) << (_parameters[j]->is_blind() ? "(*)" : format_value(value, error))                                                                           
                                        << std::setw(12) << ((_parameters[j]->get_step_size() == 0.0 && _parameters[j]->is_blind()) ? "" : format_value(value, error))
                                        << "\\pm " << std::setw(12) << format_error(value, error) << (_parameters[j]->is_blind() ? "(*)": "") << "$ & "
                                        << (_parameters[j]->is_blind() ? "(*)" : format_double(nsigma)) << "\\\\" << std::endl;
                        else
                                std::cout << "$" << std::setw(30) << partex << "$ & $"
                                        << std::setw(12)
                                        //<< (_parameters[j]->is_blind() ? "(*)" : format_value(value, error))                                                                                            
                                        << (format_value(value, error))
                                        << "\\pm "
                                        << std::setw(12) << format_error(value, error) << (_parameters[j]->is_blind() ? "(*)": "") << "$ & "
                                        << std::setw(12) << format_error(value, error_down) << " & "
                                        << std::setw(12) << "+" << format_error(value, error_up) << " & "
                                        << (_parameters[j]->is_blind() ? "(*)" : format_double(nsigma)) << "\\\\" << std::endl;
               }
                std::cout << "\\hline" << std::endl;
                std::cout << "\\end{tabular}" << std::endl;
        }
}
    
void FitParameters::fix_parameters() {
    for (unsigned int i = 0; i<nparameters(); i++) {
        _parameters.at(i)->set_step_size(0.0);
    }
}

void FitParameters::fix_param(std::string tobefixed) {
    for (unsigned int i = 0; i < nparameters(); i++) {
        std::string parname = _parameters.at(i)->get_name();
        if (tobefixed == parname) {
            _parameters.at(i)->set_step_size(0.0);
        }
    }
}

void FitParameters::take_current_as_start() {
    for (unsigned int i = 0; i<nparameters(); i++) {
        parameter* p = _parameters.at(i);
        p->set_start_value(p->get_value());
    }
} 

