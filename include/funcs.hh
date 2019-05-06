/**
 * Borrowed from HD fitter
 */

#ifndef FUNCS_H
#define FUNCS_H

//#include "hists.hh"
//#include "event.hh"
#include <complex>
#include <Math/Vector4D.h>
#include <Math/Boost.h>


  typedef ROOT::Math::PxPyPzEVector LorentzVector;
  typedef ROOT::Math::Boost LorentzBoost;
  typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> > Vector3;

  void test_diff();

  //blinding
  double evaluate_unblind_uniform(double valueSeenByAnalist, const char* blindingString, double scale);
  double bring_to_2pi_domain(double angle);
  double evaluate_unblind_uniform_angle(double valueSeenByAnalist, const char* blindingString, double scale);

  std::string format_double(double value, unsigned int digits=3);
  std::string format_value(double value, double error, unsigned int digits=3);
  std::string format_error(double value, double error, unsigned int digits=3);
  std::string format_result(double value, double error, unsigned int digits=3, std::string pm=" \\pm ");

  double sqr(double x);

  double diff_term(unsigned int n, double lambda, double tau, double sigma);
  double diff_term_b(unsigned int n, double lambda, double tau, double sigma);

  double integrate_t_to_the_n_convoluted_exp_old(unsigned int n, double tau, double sigma);

  void test_integrations();
  //new stuff
  double integrate_t_to_the_n_convoluted_exp(unsigned int n, double tau, double sigma);
  double integrate_t_to_the_n_convoluted_exp(double from, double to, int n, double tau, double sigma);
  double diff_rest(double from, double to, int n, double tau, double sigma); 
  double diff_erf_lambda(double x, int n, double tau, double sigma); 
  double diff_exp_one_minus_lambda_squared(int n, double tau, double sigma); 
  //new stuff ends here

  double diff_tau_over_lambda(int n, double tau);

  double integrate_convoluted_exp(double a, double b, double tau, double sigma);
  double integrate_convoluted_exp_cos(double a, double b, double sigma, double tau, double deltam);

  std::vector<double> integrate_legendre(unsigned int order, double a, double b);
  //double integrate_legendre(int order, double a, double b);
  //double integrate_legendre(std::vector<double> coeffs, double a, double b);
  double integrate_legendre_4d(std::vector<double> coeffs, 
			       double xa, double xb, double ya, double yb, double za, double zb, double wa, double wb, 
			       unsigned int orderx, unsigned int ordery,
			       unsigned int orderz, unsigned int orderw);

  //void save_to_txt(std::vector<event> events, std::string filename);
  //std::vector<event> load_from_txt(std::string filename);

  //void save_to_root(std::vector<event> events, std::string filename);
  //std::vector<event> load_from_root(std::string filename);

  //entry_vector calc_legendre_coefficients_4d(const entry_vector& data, 
  //					     unsigned int nx, unsigned int ny, 
  //					     unsigned int nz, unsigned int nw, 
  //					     unsigned int orderx, unsigned int ordery,
  //					     unsigned int orderz, unsigned int orderw);
  //entry_vector calc_legendre_coefficients_4d_viachi2(const hist_4d& data, const hist_4d& prediction, 
  //						     unsigned int orderx, unsigned int ordery,
  //						     unsigned int orderz, unsigned int orderw);
  //entry_vector calc_legendre_coefficients_3d_viachi2(const hist_3d& data,  
  //						     unsigned int orderx, unsigned int ordery, unsigned int orderz);
  //entry_vector calc_legendre_coefficients_1d(const entry_vector& data, 
  //					     unsigned int n, unsigned int order);
  //entry_vector calc_legendre_coefficients_1d_viachi2(const entry_vector& data, 
  //						     unsigned int n, unsigned int order);
  //
  ////std::vector<double> constants_from_legendre(const std::vector<double>& coeffs, bool orthonormal=false);
  //entry_vector constants_from_legendre_1d(const entry_vector& coeffs);
  //entry_vector constants_from_legendre_4d(const entry_vector& coeffs, 
  //					  unsigned int order_x, unsigned int order_y, 
  //					  unsigned int order_z, unsigned int order_w);
  //entry_vector constants_from_legendre_3d(const entry_vector& coeffs, 
  //					  unsigned int order_x, unsigned int order_y, unsigned int order_z);
  //
  //
  //entry_vector correct_constants_1d(const entry_vector& constants, 
  //				    unsigned int order_t, double t_min, double t_max);
  //entry_vector correct_constants_4d(const entry_vector& constants, 
  //				    unsigned int order_t, unsigned order_cos_theta, 
  //				    unsigned int order_phi, unsigned int order_cos_psi,
  //				    double t_min, double t_max, double phi_min, double phi_max);
  //entry_vector correct_constants_3d(const entry_vector& constants, 
  //				    unsigned order_cos_theta, unsigned int order_phi, unsigned int order_cos_psi,
  //				    double phi_min, double phi_max);


  double integral_x_to_n_times_exp_minus_x(double x, double tau, int n);

  double integral_x_to_n_times_sin_x(double x, int n);
  double integral_x_to_n_times_cos_x(double x, int n);

  double integral_x_to_n_times_sin_2x(double x, int n);
  double integral_x_to_n_times_cos_2x(double x, int n);
  double integral_x_to_n_times_cos_x_2(double x, int n);
  double integral_x_to_n_times_sin_x_2(double x, int n);

  double integral_x_to_n_times_asin_x(double x, int n);
  double integral_x_to_n_times_sqrt_1_minus_x2(double x, int n);

  void test_legendre();
  void test_legendre_2pi();
  void test_legendre_2d();
  void test_chebychev();

  double chebychev(double x, int n);
  void chebychev(double x, int n, std::vector<double>& results);
  void legendre(double x, int n, std::vector<double>& results);
  void orthonormal_legendre(double x, int n, std::vector<double>& results);

  //void calc_data_eff(const bs2jpsiphi_pdf* prob, const std::vector<event>& events, unsigned int n);

  void test_integrals();

  /*
    void integral_ijk_fx(double thetaa, double thetab, 
    double phia, double phib, 
    double psia, double psib, 
    int i, int j, int k,
    double& f1, double& f2, double& f3, double& f4, double& f5, double& f6);

    double integral_ijk_f1(double thetaa, double thetab, double phia, double phib, double psia, double psib, 
    int i, int j, int k);
    double integral_ijk_f2(double thetaa, double thetab, double phia, double phib, double psia, double psib, 
    int i, int j, int k);
    double integral_ijk_f3(double thetaa, double thetab, double phia, double phib, double psia, double psib, 
    int i, int j, int k);
    double integral_ijk_f4(double thetaa, double thetab, double phia, double phib, double psia, double psib, 
    int i, int j, int k);
    double integral_ijk_f5(double thetaa, double thetab, double phia, double phib, double psia, double psib, 
    int i, int j, int k);
    double integral_ijk_f6(double thetaa, double thetab, double phia, double phib, double psia, double psib, 
    int i, int j, int k);

    void integral_fx(double thetaa, double thetab, 
    double phia, double phib, 
    double psia, double psib, 
    double& f1, double& f2, double& f3, double& f4, double& f5, double& f6);
    double integral_f1(double thetaa, double thetab, double phia, double phib, double psia, double psib);
    double integral_f2(double thetaa, double thetab, double phia, double phib, double psia, double psib);
    double integral_f3(double thetaa, double thetab, double phia, double phib, double psia, double psib);
    double integral_f4(double thetaa, double thetab, double phia, double phib, double psia, double psib);
    double integral_f5(double thetaa, double thetab, double phia, double phib, double psia, double psib);
    double integral_f6(double thetaa, double thetab, double phia, double phib, double psia, double psib);

    void int_test();
  */

  //simple progress bar
  void progressbar(unsigned int percent, unsigned int nchar=50, std::string prefix="", unsigned int bar_index=0);
  ///compares the different implementations of the complex error function
  void test_errf();
  void test_convolutions();
  void test_int_convolutions();
  ///calculates the legendre polynomial of lth order
  double legendre_poly(unsigned int l, double x);
  ///calculates the associated legendre polynomial for l,m order
  double assoc_legendre(unsigned int l, unsigned int m, double x);
  ///calculates the real part of the spherical harmonics for l, m 
  double spherical_harmonic_re(unsigned int l, unsigned int m, double theta, double phi);
  ///calculates the imaginary part of the spherical harmonics for l, m 
  double spherical_harmonic_im(unsigned int l, unsigned int m, double theta, double phi);
  ///function calculates the angles in the transversity base using the momenta of the final state particles 
  /*void calculate_transversity_angles(const LorentzVector& mu_minus, 
				     const LorentzVector& mu_plus, 
				     const LorentzVector& k_minus, 
				     const LorentzVector& k_plus, 
				     double & cos_theta, double & angle_phi, double & cos_psi);*/
  ///this function calculates the convolution sin(angle) from 0 to pi
  double convoluted_sin(double angle, double sigma);
  double convoluted_cos_sin(double angle, double sigma);
  double convoluted_cos(double angle, double sigma);//not used
  double convoluted_cos_2(double angle, double sigma);//not used
  double convoluted_sin_2(double angle, double sigma);//used by s wave
  double convoluted_2_sin(double angle, double sigma);//used by s wave
  double convoluted_2_cos(double angle, double sigma);//not used
  double convoluted_cos_2_sin(double angle, double sigma);
  double convoluted_2_sin_sin(double angle, double sigma);
  double convoluted_sin_3(double angle, double sigma);
  //this integrates the above expressions (vonvolution before has been from 0 to pi
  double int_convoluted_sin(double angle, double sigma);
  double int_convoluted_cos_sin(double angle, double sigma);
  double int_convoluted_cos_2_sin(double angle, double sigma);
  double int_convoluted_2_sin_sin(double angle, double sigma);
  double int_convoluted_sin_3(double angle, double sigma);
  double int_convoluted_sin_2(double angle, double sigma);//used by s wave
  double int_convoluted_2_sin(double angle, double sigma);//used by s wave
  //convolutions for phi. this is easier due to the wraparound at +-pi -> convolution from -int to inf
  double inf_convoluted_sin(double angle, double sigma);
  double inf_convoluted_cos(double angle, double sigma);
  double inf_convoluted_cos_2(double angle, double sigma);
  double inf_convoluted_sin_2(double angle, double sigma);
  double inf_convoluted_2_sin(double angle, double sigma);
  //integral of above
  double int_inf_convoluted_const(double angle, double sigma);
  double int_inf_convoluted_sin(double angle, double sigma);
  double int_inf_convoluted_cos(double angle, double sigma);
  double int_inf_convoluted_cos_2(double angle, double sigma);
  double int_inf_convoluted_sin_2(double angle, double sigma);
  double int_inf_convoluted_2_sin(double angle, double sigma);
  //integrals from -inf to inf over the convolution from 0 to pi with acceptance
  double int_convoluted_sin_cos_n(int n, double sigma);
  double int_convoluted_sin_3_cos_n(int n, double sigma);
  double int_convoluted_cos_2_sin_cos_n(int n, double sigma);
  double int_convoluted_2_sin_sin_cos_n(int n, double sigma);
  double int_convoluted_cos_sin_cos_n(int n, double sigma);
  double int_convoluted_sin_2_cos_n(int n, double sigma);
  //integrals from 0 to pi over the convolution from -inf to inf with acceptance
  double int_inf_convoluted_const_x_n(int n, double sigma);
  double int_inf_convoluted_cos_2_x_n(int n, double sigma);
  double int_inf_convoluted_sin_2_x_n(int n, double sigma);
  double int_inf_convoluted_sin_x_n(int n, double sigma);
  double int_inf_convoluted_2_sin_x_n(int n, double sigma);
  double int_inf_convoluted_cos_x_n(int n, double sigma);

  ///implementation of Crystal Ball function
  double CB_function(double m, double mean, double alpha, double n, double sigma);
  double norm_CB_function(double mean, double alpha, double n, double sigma, double m_min, double m_max);
  ///function gives the result of exp(-ct/tau)*sin(delta_m*ct) and exp(-ct/tau)*cos(delta_m*ct) convoluted with a gaussian with width sigma 
  void convoluted_exp_sincos(double ct, double sigma, double tau, double deltam, double& one, double& two);
  ///function gives the norm of the expression exp(-ct/tau)*sin(delta_m*ct) and exp(-ct/tau)*cos(delta_m*ct) convoluted with a gaussian with width sigma 
  void norm_convoluted_exp_sincos(double tau, double deltam, double& one, double& two);
  //function gives the result of exp(-ct/tau)*sin(delta_m*ct) and exp(-ct/tau)*cos(delta_m*ct) convoluted with a gaussian with width sigma and mean m and k-factor k
  void convoluted_exp_sincos_sl_k(double ct, double sigma, double tau, double deltam, double k, double m, double& one, double& two);
  ///function calculates exp(-ct/tau) convoluted with a gaussian with width sigma 
  double convoluted_exp(double ct, double sigma_ct, double S);
  ///function calculates exp(-ct/tau) convoluted with a gaussian with width sigma and mean m_reso, for semileptonics decays, with K-FACTORS
  double convoluted_exp_sl_k(double ct, double sigma_ct, double S, double m_reso, double k);
  //integrated the function exp(-ct/tau)* cos(dm t)convoluted with a gaussian with width sigma and mean mean, for semileptonics decays, with K-FACTORS 
  double integrate_convoluted_exp_cos_sl_k(double a, double b, double sigma, double tau, double deltam, double mean, double k);
  //integrated the function t *exp(-ct/tau)* cos(dm t)convoluted with a gaussian with width sigma and mean mean, for semileptonics decays, with K-FACTORS 
  double integrate_convoluted_t_exp_cos_sl_k(double a, double b, double sigma, double tau, double deltam, double mean, double k);
  //integrated full time pdf, mixed, unmixed separately, including acceptance

  double integrate_convoluted_exp_cos_sl_k_generic2(double t, double  sigma, double tau, double deltam, double mean, double k, double alpha, double beta, double shift);
  //integrated full time pdf, mixed, unmixed separately, including acceptance
  double integrate_convoluted_exp_cos_sl_k_ab(double a, double b, double  sigma, double tau, double deltam, double mean, double k, double alpha, double beta, double shift);
  //integrated time exp pdf, mixed, unmixed separately, including acceptance
  double integrate_convoluted_exp_sl_k_generic(double t, double  sigma, double tau, double deltam, double mean, double k, double alpha, double beta, double gamma, double shift);
  //integrated time exp pdf, mixed, unmixed separately, including acceptance
  double integrate_convoluted_exp_sl_k_ab(double a, double b, double  sigma, double tau, double deltam, double mean, double k, double alpha, double beta, double gamma, double shift);

  double integrate_convoluted_exp_cos_sl_k_generic2_TurboAcc(double t, double  sigma, double tau, double deltam, double mean, double k, double alpha, double beta, double gamma, double shift, double alpha2, double frac);
  //integrated full time pdf, mixed, unmixed separately, including acceptance with 2 exponentials
  double integrate_convoluted_exp_cos_sl_k_ab_TurboAcc(double a, double b, double  sigma, double tau, double deltam, double mean, double k, double alpha, double beta, double gamma, double shift, double alpha2, double frac);
  double integrate_convoluted_exp_sin_sl_k_generic2_TurboAcc(double t, double  sigma, double tau, double deltam, double mean, double k, double alpha, double beta, double gamma, double shift, double alpha2, double frac);
  //integrated full time pdf, mixed, unmixed separately, including acceptance with 2 exponentials
  double integrate_convoluted_exp_sin_sl_k_ab_TurboAcc(double a, double b, double  sigma, double tau, double deltam, double mean, double k, double alpha, double beta, double gamma, double shift, double alpha2, double frac);
  //integrated time exp pdf, mixed, unmixed separately, including acceptance with 2 exponentials
  double integrate_convoluted_exp_sl_k_generic_TurboAcc(double t, double  sigma, double tau, double deltam, double mean, double k, double alpha, double beta, double gamma, double shift, double alpha2, double frac);
  //integrated time exp pdf, mixed, unmixed separately, including acceptance
  double integrate_convoluted_exp_sl_k_ab_TurboAcc(double a, double b, double  sigma, double tau, double deltam, double mean, double k, double alpha, double beta, double gamma, double shift, double alpha2, double frac);

  double integrate_convoluted_exp_cos_sl_k_NoAcc(double t, double  sigma, double tau, double deltam, double mean, double k);
  //integrated full time pdf, mixed, unmixed separately, including acceptance
  double integrate_convoluted_exp_cos_sl_k_ab_NoAcc(double a, double b, double  sigma, double tau, double deltam, double mean, double k);
  //integrated time exp pdf, mixed, unmixed separately, including acceptance
  double integrate_convoluted_exp_sl_k_NoAcc(double t, double  sigma, double tau, double deltam, double mean, double k);
  //integrated time exp pdf, mixed, unmixed separately, including acceptance
  double integrate_convoluted_exp_sl_k_ab_NoAcc(double a, double b, double  sigma, double tau, double deltam, double mean, double k);

  ///function calculates the norm of exp(-ct/tau) convoluted with a gaussian with width sigma 
  double norm_convoluted_exp(double ct, double sigma_ct, double S);
  ///johnson function (for delta Mass sig?) 
  double johnson(double x, double mu, double sigma, double gamma, double delta, double x_min, double x_max);
  ///DstD0BG - roofit function - function (for delta Mass bkg) 
  double DstD0BG(double x, double dm0, double C, double A, double B, double x_min, double x_max);
  ///function calculates exp(-ct/tau) depending on the parameters normalized or not 
  double expdecay(double t, double tau, bool normalized);
  ///function calculates 1/(sqrt(2*pi)*exp(-(x-mean)^2/(2*sigma^2)))
  double gauss(double x, double sigma, double mean);
  ///function calculates a normed gaussian between the values min and max used in the mass pdfs 
  double normed_gauss(double x, double sigma, double mean, double min, double max);

  ///function calculates the normalization of a gaussian used in the mass pdfs 
  double norm_gauss(double sigma, double mean, double min, double max);
  ///function calculates a NOTnormed gaussian between the values min and max used in the mass pdfs 
  double NOTnormed_gauss(double x, double sigma, double mean, double min, double max);

  double integrate_normed_gauss(double sigma, double mean, double min, double max, double a, double b);

  ///linear probability, not used in bs2jpsiphi atm. 
  double linear(double x, double min, double max, double slope);
  ///very important helper function, the error function with complex argument. This is probably one of the performance bottlenecks. Probably much faster with a lookup table and interpolation
  std::complex<double> cErrF(const std::complex<double>& x);
  ///implementation 2 of the complex error function
  std::complex<double> cErrF_2(const std::complex<double>& x);
  std::complex<double> ErrF_2(const std::complex<double>& x);

  ///implementation 3 of the complex error function
  std::complex<double> cErrF_3(const std::complex<double>& x);
  ///the faddeeva function function
  std::complex<double> Faddeeva_2 (const std::complex<double>& z);
  ///used by cErrF, does the actual calculations
  std::complex<double> wErrF(const std::complex<double>& arg);
  ///This is another version of the complex error function, used in the newer cdf code
  std::complex<double> nwwerf(const std::complex<double> z);



#endif
