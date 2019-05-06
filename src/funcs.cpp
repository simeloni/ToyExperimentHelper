/**
 * Colletion of functions borrowed  from asld fitter
 */

#include <complex>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <TMath.h>
#include <TLegend.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <Math/SpecFuncMathCore.h>
#include <Math/SpecFuncMathMore.h>
#include <Math/Vector3D.h>
#include <Math/Boost.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>
//#include <TDecompBK.h>
#include <TDecompChol.h>
#include "funcs.hh"
//#include "event.hh"
#include <RooFit.h>
#include <RooRealVar.h>
#include <RooUnblindUniform.h>
//#include "symbolic.hh"
//#include "bs2jpsiphi_pdf.hh"

using namespace std;

const double M_4_PI = 1.2732395447351626861510701069801;   // 4/Pi
const double M_SQRTPI = 1.7724538509055160272981674833411;   // sqrt(Pi)
const double M_SQRT2PI = 2.5066282746310005024157652848111;   // sqrt(2Pi)
const double M_SQRT_PI_2 = 1.2533141373155002512078826424055;   // sqrt(Pi/2)
const double M_SQRT_2_PI = 0.79788456080286535587989211986876;  // sqrt(2/Pi)
const double M_1_SQRTPI = 0.56418958354775628694807945156077;  // 1/sqrt(Pi)
//const double M_2_SQRTPI = 1.1283791670955125738961589031216;   // 2/sqrt(Pi)
const double m_2_sqrtpi = 1.1283791670955125738961589031216;   // 2/sqrt(Pi)
const double M_SQRT_2PI = 0.39894228040143267793994605993438;  // 1/sqrt(2Pi)
const double M_SQRT_2 = 0.70710678118654752440084436210485;  // sqrt(1/2)
const double M_2PI = 6.2831853071795864769252867665590;   // 2Pi
const double M_4PI = 12.5663706143591729538505735331180;//4Pi

double  evaluate_unblind_uniform(double valueSeenByAnalist, const char* blindingString, double scale)
{
  // note: this is really slow because RooBlindTools generates _many_ random numbers
  // for all of the blinders for every instance. Need to discuss with Wouter V.
  // how we can solve this.
  double rc = valueSeenByAnalist;
  bool gBlind = true;
  if(gBlind) {
    RooRealVar blindedx("bx","bx",valueSeenByAnalist);
    RooUnblindUniform unblindedx("ubx","ubx",blindingString,scale,blindedx);
    RooAbsHiddenReal* ar = &unblindedx;
    rc = ar->getHiddenVal();
  }
  return rc;
}

double bring_to_2pi_domain(double angle)
{
  const double pi = TMath::Pi();
  double multipleoftwopi = (angle+pi)/(2*pi);
  if( multipleoftwopi < 0 ) {
    multipleoftwopi += 1-int( multipleoftwopi );
  } else {
    multipleoftwopi -= int( multipleoftwopi );
  }
  return pi*(2*multipleoftwopi-1);
}

double evaluate_unblind_uniform_angle(double valueSeenByAnalist, const char* blindingString, double scale)
{
  double rc = valueSeenByAnalist;
  bool gBlind = true;
  if(gBlind) {
    RooRealVar blindedx("bx","bx",valueSeenByAnalist);
    RooUnblindUniform unblindedx("ubx","ubx",blindingString,scale,blindedx);
    RooAbsReal& ar = unblindedx;
    rc = bring_to_2pi_domain(ar.getVal());
  }
  return rc;
}


std::string format_double(double value, unsigned int digits)
{
  stringstream out;
  out << fixed << setprecision(digits) << value;
  return out.str();
};

std::string format_value(double value, double error, unsigned int digits)
{
  //format value such that the relevant number of digits is given back
  //determine first "digits" significant places of error
  stringstream out;
  if (error == 0.0)
    {
      out << value;
    }
  else
    {
      int first_significant = floor(log10(error));
      if (first_significant > 0)
	first_significant = 0;
      //ln(e^x)=x log10(c*10^x)=x+log10(c) with log10(c) < 1
      out << fixed << setprecision(fabs(first_significant)+digits-1)
	  << value;
    }
  return out.str();
};

std::string format_error(double value, double error, unsigned int digits)
{
  //format value such that the relevant number of digits is given back
  //determine first "digits" significant places of error
  stringstream out;
  if (error == 0.0)
    {
      out << error;
    }
  else
    {
      int first_significant = floor(log10(error));
      if (first_significant > 0)
	first_significant = 0;
      //ln(e^x)=x log10(c*10^x)=x+log10(c) with log10(c) < 1
      out << fixed << setprecision(fabs(first_significant)+digits-1)
	  << error;
    }
  return out.str();
};

std::string format_result(double value, double error, unsigned int digits, std::string pm)
{
  //format value such that the relevant number of digits is given back
  //determine first "digits" significant places of error
  stringstream out;
  if (error == 0.0)
    {
      out << value << pm << error;
    }
  else
    {
      int first_significant = floor(log10(error));
      if (first_significant > 0)
	first_significant = 0;
      //ln(e^x)=x log10(c*10^x)=x+log10(c) with log10(c) < 1
      out << fixed << setprecision(fabs(first_significant)+digits-1)
	  << value
	  << pm
	  << error;
    }
  return out.str();
};


void test_diff() 
{
  double sigma = 0.03;
  double tau = 1.5;
  double lambda = 1.0;
  unsigned int degree = 20;

  std::cout << std::endl;
  std::cout << "tau=" << tau << " sigma=" << sigma << std::endl;
  for (unsigned int i = 0; i<= degree; i++)
    {
      std::cout << "via recursion: " << "(d/dlambda)^" << i << " f(lambda)|lambda=1 = " << diff_term(i, lambda, tau, sigma) << std::endl;

      std::cout << "int t^" << i << " exp(-t/tau) x gauss(sigma) " << integrate_t_to_the_n_convoluted_exp(i, tau, sigma) << std::endl;
    }
}

//calculates int_a^b int_0^inf  1/sqrt(2*pi)/sigma exp(-(t-d)^2/2/sigma^2) exp(-d/tau) dd dt
double integrate_convoluted_exp(double a, double b, double tau, double sigma)
{
  /*
    -%e^(-b/tau-a/tau)*(tau*%e^(a/tau+sigma^2/(2*tau^2))*erf((sqrt(2)*b*tau-\
    sqrt(2)*sigma^2)/(2*sigma*tau))-tau*%e^(b/tau+sigma^2/(2*tau^2))*erf((sqrt(2)*\
    a*tau-sqrt(2)*sigma^2)/(2*sigma*tau))+tau*%e^(a/tau+sigma^2/(2*tau^2))+((erf(a\
    /(sqrt(2)*sigma))-erf(b/(sqrt(2)*sigma)))*tau*%e^(a/tau)-tau*%e^(sigma^2/(2*ta\
    u^2)))*%e^(b/tau))/2

    -exp(-b/tau-a/tau)*
    (tau
    *exp(a/tau+sigma^2/(2*tau^2))
    *erf((sqrt(2)*b*tau-sqrt(2)*sigma^2)/(2*sigma*tau))
    -tau
    *exp(b/tau+sigma^2/(2*tau^2))
    *erf((sqrt(2)*a*tau-sqrt(2)*sigma^2)/(2*sigma*tau))
    
    +tau
    *exp(a/tau+sigma^2/(2*tau^2))
    
    +((erf(a/(sqrt(2)*sigma))-erf(b/(sqrt(2)*sigma)))
    *tau
    *exp(a/tau)
    -tau
    *exp(sigma^2/(2*tau^2)))
    *exp(b/tau))/2

  */
  double result =
    -exp(-b/tau-a/tau)*
    (tau 
     * exp(a/tau+sigma*sigma/(2.0*tau*tau))
     * erf((sqrt(2.0)*b*tau-sqrt(2.0)*sigma*sigma)/(2.0*sigma*tau))
     - tau
     * exp(b/tau+sigma*sigma/(2.0*tau*tau))
     * erf((sqrt(2.0)*a*tau-sqrt(2.0)*sigma*sigma)/(2.0*sigma*tau))
     + tau
     * exp(a/tau+sigma*sigma/(2.0*tau*tau))
     + ((erf(a/(sqrt(2.0)*sigma))-erf(b/(sqrt(2.0)*sigma)))
	* tau
	* exp(a/tau)
	- tau
	* exp(sigma*sigma/(2.0*tau*tau)))
     * exp(b/tau))/2.0;
  return result;
}

//this is to differentiate 1/lambda * exp(a*(1-lambda)^2) recursively
double diff_term(unsigned int n, double lambda, double tau, double sigma)
{
  assert(n >= 0);
  //double a = 0.5*sigma*sigma/(tau*tau);
  if (n == 0)
    //return 1.0/lambda*exp(a*(1.0-lambda)*(1.0-lambda));
    //assume lambda=1
    return 1.0;
  else
    {
      double result = 0.0;
      for (unsigned int k = 0; k<n; k++)
	{
	  result += diff_term_b(k, lambda, tau, sigma) 
	    * diff_term(n - k - 1, lambda, tau, sigma) 
	    * TMath::Binomial(n-1, k);
	}
      return result;
    }
}



//is this for exp(a*(1-lambda)^2)?
double diff_term_b(unsigned int n, double lambda, double tau, double sigma)
{
  assert(n >= 0);
  double a = 0.5*sigma*sigma/(tau*tau);
  //assume lambda = 1.0;
  if (n == 0)
    return -1.0;
  else if (n == 1)
    return 1.0 + 2.0*a;
  else 
    return (n % 2 == 0 ? -1.0 : 1.0) * TMath::Factorial(n);
}

void test_integrations() 
{
  double tau = 1.5;
  double sigma = 0.03;
  double min=-2.0, max=200.0;
  double w = max - min;
  unsigned int divisions = 100;
  double dw = w/divisions;
  double a, b, value = 0.0, value2 = 0.0, current, current2;
  for (unsigned int i = 0; i<15; i++)
    {
      value = 0.0;
      value2 = 0.0;
      for (unsigned int j = 0; j < divisions; j++)
	{
	  a = min + j*dw;
	  b = min + (j+1)*dw;
	  current = integrate_t_to_the_n_convoluted_exp(a, b, i, tau, sigma);
	  //cout << "Int_" << a << "^" << b << " t^" << i << " smeared exp = " << current << endl;
	  if (i == 0)
	    {
	      current2 = integrate_convoluted_exp(a, b, tau, sigma);
	      //cout << "Int_" << a << "^" << b << " t^0 smeared exp = " << current2  << endl;
	      value2 += current2;
	    }
	  value += current;
	}
      cout << "0. Int_-inf^+inf t^" << i << " smeared exp = " << value << endl;
      cout << "1. Int_-inf^+inf t^" << i << " smeared exp = " << value2 << endl;
      cout << "2. Int_-inf^+inf t^" << i << " smeared exp = " << integrate_t_to_the_n_convoluted_exp(i, tau, sigma) << endl;
      cout << "3. Int_-inf^+inf t^" << i << " smeared exp = " << integrate_t_to_the_n_convoluted_exp_old(i, tau, sigma) << endl;
      cout << endl;
    }
}

//maxima understands something like this: 
//for n:0 thru 10 do print(diff(erf(a+(1-lambda)*b),lambda,n));
//for n:0 thru 10 do print(string(diff(erf(a+(1-lambda)*b),lambda,n)));

//actually this works, too
//for n:0 thru 3 do print(string(subst(1,lambda,diff(tau/lambda*exp(c*(1-lambda)^2)*erf(a+(1-lambda)*b),lambda,n))));

//for n:0 thru 14 do print(string(integrate(sqrt(alpha/%pi)*exp(-alpha*(t-d)^2)*cos(t)^n, t, -inf, inf)));
//for n:0 thru 8 do print("n=",n,": ",expand(integrate(sqrt(alpha/%pi)*exp(-alpha*(t-d)^2)*cos(t)^n, t, -inf, inf)));

//for k:0 thru 3 do print("k=",k,":  ",expand(exp(-(2*k+1)^2/4/alpha)*cos((2*k+1)*t)/2^7*binomial(3,k)));
//n=7
//for k:0 thru 3 do print(expand(exp(-(2*k+1)^2/4/alpha)*cos((7-2*k+1)*t)/2^6*binomial(7,k)));
//n=10
//for k:0 thru 5 do print(expand(exp(-(2*k+1)^2/4/alpha)*cos((10-2*k)*t)/2^9*binomial(10,k)));
//for k:0 thru n/2 do print(expand(1/2^(n-1)*binomial(n,k)*exp(-1/4/alpha*(n-2*k)^2)*cos((n-2*k)*t)));
//for k:0 thru 5 do print(expand(1/2^(10-1)*binomial(10,k)*exp(-1/4/alpha*(10-2*k)^2)*cos((10-2*k)*t)));
//for k:0 thru 4 do print(expand(1/2^(9-1)*binomial(9,k)*exp(-1/4/alpha*(9-2*k)^2)*cos((9-2*k)*t)));

//differentiates tau/lambda n times
double diff_tau_over_lambda(int n, double tau) 
{
  if (n % 2 == 0)
    return tau*TMath::Factorial(n);
  else
    return -tau*TMath::Factorial(n);
}

//for n:0 thru 10 do print(string(subst(1,lambda,diff(exp(a*(1-lambda)^2),lambda,n))));
//this differentiates exp(a*(1-lambda)^2) n times and then sets lambda = 1
double diff_exp_one_minus_lambda_squared(int n, double tau, double sigma) 
{
  if (n == 0)
    return 1.0;
  else if (n % 2 == 0)
    {
      double a = sigma*sigma/(2.0*tau*tau);
      return TMath::Factorial(n)/TMath::Factorial(n/2)*pow(a, n/2);
    }
  else
    return 0.0;
}

//this differentiates erf(a+(1-lambda)*b) n times to lambda and then evaluates for lambda = 1
double diff_erf_lambda(double x, int n, double tau, double sigma) 
{
  /*
  //for n:0 thru 15 do print(string(factor(subst(1,lambda,diff(erf(a+(1-lambda)*b),lambda,n)))));
  erf(a) 
  -2*%e^-a^2*b/sqrt(%pi) 
  -4*a*%e^-a^2*b^2/sqrt(%pi) 
  -4*(2*a^2-1)*%e^-a^2*b^3/sqrt(%pi) 
  -8*a*(2*a^2-3)*%e^-a^2*b^4/sqrt(%pi) 
  -8*(4*a^4-12*a^2+3)*%e^-a^2*b^5/sqrt(%pi) 
  -16*a*(4*a^4-20*a^2+15)*%e^-a^2*b^6/sqrt(%pi) 
  -16*(8*a^6-60*a^4+90*a^2-15)*%e^-a^2*b^7/sqrt(%pi) 
  -32*a*(8*a^6-84*a^4+210*a^2-105)*%e^-a^2*b^8/sqrt(%pi) 
  -32*(16*a^8-224*a^6+840*a^4-840*a^2+105)*%e^-a^2*b^9/sqrt(%pi) 
  -64*a*(16*a^8-288*a^6+1512*a^4-2520*a^2+945)*%e^-a^2*b^10/sqrt(%pi) 
  -64*(32*a^10-720*a^8+5040*a^6-12600*a^4+9450*a^2-945)*%e^-a^2*b^11/sqrt(%pi) 
  -128*a*(32*a^10-880*a^8+7920*a^6-27720*a^4+34650*a^2-10395)*%e^-a^2*b^12/sqrt(%pi) 
  -128*(64*a^12-2112*a^10+23760*a^8-110880*a^6+207900*a^4-124740*a^2+10395)*%e^-a^2*b^13/sqrt(%pi) 
  -256*a*(64*a^12-2496*a^10+34320*a^8-205920*a^6+540540*a^4-540540*a^2+135135)*%e^-a^2*b^14/sqrt(%pi) 
  -256*(128*a^14-5824*a^12+96096*a^10-720720*a^8+2522520*a^6-3783780*a^4+1891890*a^2-135135)*%e^-a^2*b^15/sqrt(%pi) 
  */
  const double sqrt_2 = sqrt(2.0);
  const double pi = TMath::Pi();
  double a = x/(sqrt_2*sigma);
  double b = sigma/(sqrt_2*tau);
  switch(n) {
  case 0: 
    return erf(a);
    break;
  case 1: 
    return -2.0*exp(-a*a)*b/sqrt(pi); 
    break;
  case 2: 
    return -4.0*a*exp(-a*a)*pow(b,2)/sqrt(pi); 
    break;
  case 3: 
    return -4.0*(2*pow(a,2)-1.0)*exp(-a*a)*pow(b,3)/sqrt(pi); 
    break;
  case 4: 
    return -8.0*a*(2.0*pow(a,2)-3.0)*exp(-a*a)*pow(b,4)/sqrt(pi); 
    break;
  case 5: 
    return -8*(4*pow(a,4)-12*pow(a,2)+3)*exp(-a*a)*pow(b,5)/sqrt(pi); 
    break;
  case 6: 
    return -16*a*(4*pow(a,4)-20*pow(a,2)+15)*exp(-a*a)*pow(b,6)/sqrt(pi); 
    break;
  case 7: 
    return -16*(8*pow(a,6)-60*pow(a,4)+90*pow(a,2)-15)*exp(-a*a)*pow(b,7)/sqrt(pi); 
    break;
  case 8: 
    return -32*a*(8*pow(a,6)-84*pow(a,4)+210*pow(a,2)-105)*exp(-a*a)*pow(b,8)/sqrt(pi); 
    break;
  case 9: 
    return -32*(16*pow(a,8)-224*pow(a,6)+840*pow(a,4)-840*pow(a,2)+105)*exp(-a*a)*pow(b,9)/sqrt(pi); 
    break;
  case 10: 
    return -64*a*(16*pow(a,8)-288*pow(a,6)+1512*pow(a,4)-2520*pow(a,2)+945)*exp(-a*a)*pow(b,10)/sqrt(pi); 
    break;
  case 11: 
    return -64*(32*pow(a,10)-720*pow(a,8)+5040*pow(a,6)-12600*pow(a,4)+9450*pow(a,2)-945)*exp(-a*a)*pow(b,11)/sqrt(pi); 
    break;
  case 12: 
    return -128*a*(32*pow(a,10)-880*pow(a,8)+7920*pow(a,6)-27720*pow(a,4)+34650*pow(a,2)-10395)*exp(-a*a)*pow(b,12)/sqrt(pi); 
    break;
  case 13: 
    return -128*(64*pow(a,12)-2112*pow(a,10)+23760*pow(a,8)-110880*pow(a,6)+207900*pow(a,4)-124740*pow(a,2)+10395)*exp(-a*a)*pow(b,13)/sqrt(pi); 
    break;
  case 14: 
    return -256*a*(64*pow(a,12)-2496*pow(a,10)+34320*pow(a,8)-205920*pow(a,6)+540540*pow(a,4)-540540*pow(a,2)+135135)*exp(-a*a)*pow(b,14)/sqrt(pi); 
    break;
  case 15: 
    return -256*(128*pow(a,14)-5824*pow(a,12)+96096*pow(a,10)-720720*pow(a,8)+2522520*pow(a,6)-3783780*pow(a,4)+1891890*pow(a,2)-135135)*exp(-a*a)*pow(b,15)/sqrt(pi); 
    break;
  default: 
    cout << "Degree " << n << " to large" << endl;
    return 0;
  }
}

//differentiates the large sum 
double diff_rest(double from, double to, int n, double tau, double sigma) 
{
  double result=0.0, factor=0.0;
  const double sqrt_2 = sqrt(2.0);
  double p1 = pow(-from/tau,n)*exp(-from/tau)*(TMath::Erf(from/(sqrt_2*sigma)-sigma/(sqrt_2*tau))+1.0);
  double p2 = pow(-to/tau,n)*exp(-to/tau)*(TMath::Erf(to/(sqrt_2*sigma)-sigma/(sqrt_2*tau))+1.0);
  //cout << "p1=" << p1 << " p2=" << p2 << " for from=" << from << " to=" << to << " order=" << n << endl; 
  result += p1;
  result -= p2;
  factor = exp(-sigma*sigma/(2.0*tau*tau));
  for (int k = 0; k <= n; k++)
    {
      result += 
	factor 
	* diff_exp_one_minus_lambda_squared(k, tau, sigma) 
	//* (diff_erf_lambda(to, n-k-1, tau, sigma) - diff_erf_lambda(from, n-k-1, tau, sigma))
	//* TMath::Binomial(n-1, k);
	* (diff_erf_lambda(to, n-k, tau, sigma) - diff_erf_lambda(from, n-k, tau, sigma))
	* TMath::Binomial(n, k);
    }
  //return factor*result;
  return result;
}

//calculates int_from^to int_0^inf  1/sqrt(2*pi)/sigma * t^n exp(-(t-d)^2/2/sigma^2) exp(-d/tau) dd dt
//this formula is correct
double integrate_t_to_the_n_convoluted_exp(double from, double to, int n, double tau, double sigma)
{
  double result=0.0, factor=0.0;//, factor2=0.0;
  //factor2 = exp(-sigma*sigma/(2.0*tau*tau));
  for (int k = 0; k <= n; k++)
    {
      result += 
	//diff_tau_over_lambda(from, to, k, tau, sigma) 
	diff_tau_over_lambda(k, tau) 
	//(k % 2 == 0 ? tau : -tau)
	//* diff_rest(from, to, n-k-1, tau, sigma) 
	//* TMath::Binomial(n-1, k);
	* diff_rest(from, to, n-k, tau, sigma) 
	* TMath::Binomial(n, k);
    }
  factor = 0.5*exp(+sigma*sigma/(2.0*tau*tau))*pow(-tau, n);
  return factor*result;
}

//calculates int_-inf^+inf int_0^inf  1/sqrt(2*pi)/sigma * t^n exp(-(t-d)^2/2/sigma^2) exp(-d/tau) dd dt
double integrate_t_to_the_n_convoluted_exp(unsigned int n, double tau, double sigma)
{
  double result=0.0, factor=0.0;
  for (unsigned int k = 0; k <= n; k++)
    {
      result += 
	//diff_tau_over_lambda(from, to, k, tau, sigma) 
	//(k % 2 == 0 ? tau : -tau)
	diff_tau_over_lambda(k, tau) 
	//* diff_rest(from, to, n-k-1, tau, sigma) 
	//* TMath::Binomial(n-1, k);
	* diff_exp_one_minus_lambda_squared(n-k, tau, sigma) 
	//* diff_rest(from, to, n-k, tau, sigma) 
	* TMath::Binomial(n, k);
    }
  factor = pow(-tau, int(n));
  return factor*result;
}

//calculates int_-inf^+inf int_0^inf  1/sqrt(2*pi)/sigma * t^n exp(-(t-d)^2/2/sigma^2) exp(-d/tau) dd dt
double integrate_t_to_the_n_convoluted_exp_old(unsigned int n, double tau, double sigma)
{
  
  //simple, eh?
  double lambda = 1.0;
  double result = diff_term(n, lambda, tau, sigma);

  result *= ((n%2 == 0 ? 1.0 : -1.0) * pow(tau, int(n+1)));
  return result;
}

std::vector<double> integrate_legendre(unsigned int order, double a, double b)
{
  std::vector<double> results(order, 0.0);
  std::vector<double> leg_a, leg_b;
  legendre(a, order+1, leg_a);
  legendre(b, order+1, leg_b);

  results.at(0) = (b-a);
  for (unsigned int i = 1; i<order; i++)
    results.at(i) =  (leg_b.at(i+1.0) - leg_b.at(i-1.0) - leg_a.at(i+1.0) + leg_a.at(i-1.0))/(2.0*i+1.0);
  return results;
}

double integrate_legendre_4d(std::vector<double> coeffs, 
			     double xa, double xb, double ya, double yb, double za, double zb, double wa, double wb, 
			     unsigned int orderx, unsigned int ordery,
			     unsigned int orderz, unsigned int orderw)
{
  assert(orderx*ordery*orderz*orderw == coeffs.size());

  std::vector<double> legendre_x = integrate_legendre(orderx, xa, xb);
  std::vector<double> legendre_y = integrate_legendre(ordery, ya, yb);
  std::vector<double> legendre_z = integrate_legendre(orderz, za, zb);
  std::vector<double> legendre_w = integrate_legendre(orderw, wa, wb);

  double coeff, result=0.0;
  unsigned int coeff_bin;
  for (unsigned int m = 0; m<orderx; m++)
    for (unsigned int n = 0; n<ordery; n++)
      for (unsigned int o = 0; o<orderz; o++)
	for (unsigned int p = 0; p<orderw; p++)
	  {
	    coeff_bin = m + orderx*n + orderx*ordery*o + orderx*ordery*orderz*p;
	    coeff = coeffs.at(coeff_bin);
	    result += coeff 
	      * legendre_x.at(m)
	      * legendre_y.at(n)
	      * legendre_z.at(o)
	      * legendre_w.at(p);
	  }

  return result;
}

double sqr(double x)
{
  return x*x;
}

/*
  chi^2 = sum_i (data_i - fitted(x_i))^2 / sigma(x_i)^2
  = sum_i (data_i - sum_j (c_j*P_j(x_i)) )^2 / sigma(x_i)^2

  dchi^2/ dc_k = sum_i  2*P_k(x_i)*(data_i - sum_j (c_j*P_j(x_i)) ) / sigma(x_i)^2 = 0

  sum_i P_k(x_i)*(data_i) / sigma(x_i)^2 = sum_i P_k(x_i)*(sum_j (c_j*P_j(x_i)) ) / sigma(x_i)^2
*/
// entry_vector calc_legendre_coefficients_1d_viachi2(const entry_vector& data, 
// 						   unsigned int n, unsigned int order)
// {
//   assert(data.size() == n);//not strictly necessary

//   std::vector<double> legendre_res;
//   std::vector<std::vector<double> > legendre_x;

//   std::vector<double> ajk(order*order, 0.0);
//   std::vector<double> bk(order, 0.0);

//   double x, value, error, sqr_error;

//   for (unsigned int i=0; i<n; i++)//for every x
//     {
//       x = -1.0 + (i+0.5)*2.0/n;
//       //cout << "funcs.cc: " << x << endl;
//       legendre(x, order, legendre_res);
//       legendre_x.push_back(legendre_res);
//     }

//   unsigned int data_bin;
//   for (unsigned int i=0; i<n; i++)
//     {
//       data_bin = i;
//       value = data.at(data_bin).value;
//       error = data.at(data_bin).error;      
//       //calculate b_k and a_jk
//       for (unsigned int k = 0; k < order; k++)
// 	{
// 	  if (error == 0.0) 
// 	    error = 1.0;
// 	  sqr_error = error*error;
// 	  bk[k] += 
// 	    legendre_x.at(i).at(k)
// 	    * value 
// 	    / sqr_error;
// 	  for (unsigned int j = 0; j < order; j++)
// 	    {
// 	      ajk[j+k*order] += 
// 		legendre_x.at(i).at(k)
// 		*legendre_x.at(i).at(j)
// 		/ sqr_error;  
// 	    }
// 	}//end calculations of b_k and a_jk
//     }
//   cout << "Matrix a_jk" << endl;
//   double mean_diagonal=0.0;
//   unsigned int nlarger=0, nsmaller =0;
//   for (unsigned int i=0; i<order; i++)
//     mean_diagonal += ajk.at(i*order)/order;
//   cout << "mean diagonal element: " << mean_diagonal << endl;
//   for (unsigned int i=0; i<order; i++)
//     {
//       for (unsigned int j=0; j<order; j++)
// 	{
// 	  if (fabs(ajk.at(j+i*order)) > 0.01*mean_diagonal)
// 	    {
// 	      cout << scientific << setprecision(2) << setw(12) 
// 		   << ajk.at(j+i*order)
// 		   << " ";
// 	      nlarger++;
// 	    }
// 	  else
// 	    {
// 	      cout << scientific << setprecision(2) << setw(12) 
// 		   << "-"
// 		   << " ";
// 	      nsmaller++;
// 	    }
// 	}
//       cout << endl;
//     }
//   cout << "Vector b_k" << endl;
//   for (unsigned int i=0; i<order; i++)
//     cout << scientific << setprecision(2) << setw(12) 
// 	 << bk.at(i) << endl;
    
//   cout << "Matrix elements which can possibly be cut: " << nsmaller << " of " << nsmaller+nlarger << endl;
//   cout << "Using ROOT Matrices to solve a_ik c_i = b_k" << endl;
//   TMatrixDSym a(order);
//   TVectorD b(order);
//   for (unsigned int i=0; i<order; i++)
//     {
//       for (unsigned int j=i; j<order; j++)
// 	a(i,j) = ajk.at(i+order*j);
//       b[i] = bk.at(i);
//     }
//   bool fine = false;
//   TDecompChol decomp(a);
//   TVectorD result  = decomp.Solve(b, fine);
//   cout << "Decomposition successfull? " << (fine ? "YES" : "NO") << endl;
//   entry_vector result_coeffs(order, 0.0);
//   for (unsigned int i=0; i<order; i++)
//     {
//       cout << "c(" << setw(3) << i << ") = " << result[i] << endl;
//       result_coeffs.at(i) = result[i];
//     }  
  
//   return result_coeffs;
// }

// entry_vector calc_legendre_coefficients_1d(const entry_vector& data, 
// 					   unsigned int n, unsigned int order)
// {
//   assert(data.size() == n);
//   std::vector<double> legendre_x;
//   entry_vector coeffs(order);
//   double dv = 2.0/n;
//   double x, value, error, current_error;
//   for (unsigned int i=0; i<n; i++)
//     {
//       x = -1.0 + (i+0.5)*2.0/n;
//       //legendre(x, n, legendre_x);
//       legendre(x, order, legendre_x);
//       value = data.at(i).value;
//       error = data.at(i).error;
	    
//       for (unsigned int m = 0; m<order; m++)
// 	{
// 	  coeffs.at(m).value += 
// 	    value * dv * (0.5*(2.0*m+1.0)) * legendre_x.at(m);
// 	  current_error = coeffs.at(m).error;
// 	  coeffs.at(m).error =
// 	    sqrt(sqr(current_error) + sqr(error * dv * (0.5*(2.0*m+1.0)) * legendre_x.at(m)));
// 	}
//     }
//   return coeffs;

// }

// entry_vector calc_legendre_coefficients_4d(const entry_vector& data, 
// 					   unsigned int nx, unsigned int ny, 
// 					   unsigned int nz, unsigned int nw, 
// 					   unsigned int orderx, unsigned int ordery,
// 					   unsigned int orderz, unsigned int orderw)
// {
//   unsigned int n = orderx*ordery*orderz*orderw;//nx*ny*nz*nw;
//   //this was wrong
//   assert(data.size() == nx*ny*nz*nw);

//   std::vector<double> legendre_x, legendre_y, legendre_z, legendre_w;
//   entry_vector coeffs(n);
//   double dv = 2.0*2.0*2.0*2.0/(nx*ny*nz*nw);
//   double x, y, z, w, value, error;
//   unsigned int data_bin, legendre_bin;
//   for (unsigned int i=0; i<nx; i++)
//     {
//       progressbar(double(i)/nx*100, 50, "Determining Legendre coefficients ");
//       x = -1.0 + (i+0.5)*2.0/nx;
//       for (unsigned int j=0; j<ny; j++)
// 	{
// 	  y = -1.0 + (j+0.5)*2.0/ny;
// 	  for (unsigned int k=0; k<nz; k++)
// 	    {
// 	      z = -1.0 + (k+0.5)*2.0/nz;
// 	      for (unsigned int l=0; l<nw; l++)
// 		{
// 		  w = -1.0 + (l+0.5)*2.0/nw;
// 		  legendre(x, orderx, legendre_x);
// 		  legendre(y, ordery, legendre_y);
// 		  legendre(z, orderz, legendre_z);
// 		  legendre(w, orderw, legendre_w);
		  
// 		  data_bin = i+ nx*j + nx*ny*k + nx*ny*nz*l;
// 		  value = data.at(data_bin).value;
// 		  error = data.at(data_bin).error;
// 		  if (!std::isnan(value) && !std::isinf(value))
// 		    {
// 		      for (unsigned int m = 0; m<orderx; m++)
// 			for (unsigned int n = 0; n<ordery; n++)
// 			  for (unsigned int o = 0; o<orderz; o++)
// 			    for (unsigned int p = 0; p<orderw; p++)
// 			      {
// 				legendre_bin = m 
// 				  + orderx*n 
// 				  + orderx*ordery*o 
// 				  + orderx*ordery*orderz*p;
			    
// 				coeffs.at(legendre_bin) += 
// 				  value * dv 
// 				  * (0.5*(2.0*m+1.0)) * (0.5*(2.0*n+1.0)) 
// 				  * (0.5*(2.0*o+1.0)) * (0.5*(2.0*p+1.0))
// 				  * legendre_x.at(m) * legendre_y.at(n)
// 				  * legendre_z.at(o) * legendre_w.at(p);
			    
// 			      }
// 		    }
// 		}
// 	    }
// 	}
//     }
//   progressbar(100, 50, "Determining Legendre coefficients ");
//   std::cout << std::endl;
//   return coeffs;
// }

// entry_vector calc_legendre_coefficients_3d_viachi2(const hist_3d& data,  
// 						   unsigned int orderx, unsigned int ordery, unsigned int orderz)
// {
//   unsigned int orderxy = orderx*ordery;
//   unsigned int orderxyz = orderx*ordery*orderz;
//   unsigned int order = orderxyz;
//   unsigned int nx = data.get_bins_x(), 
//     ny = data.get_bins_y(), 
//     nz = data.get_bins_z();
//   unsigned int nxy = nx*ny;
//   //unsigned int nxyz = nx*ny*nz;
//   assert(data.size() == nx*ny*nz);

//   std::vector<double> legendre_res;
//   std::vector<std::vector<double> > legendre_x, legendre_y, legendre_z;
//   std::vector<double> aik(order*order, 0.0);
//   std::vector<double> bk(order, 0.0);

//   double x, y, z, value, error, legendre_product, sqr_error;

//   for (unsigned int i=0; i<nx; i++)
//     {
//       x = -1.0 + (i+0.5)*2.0/nx;
//       legendre(x, orderx, legendre_res);
//       legendre_x.push_back(legendre_res);
//     }
//   for (unsigned int j=0; j<ny; j++)
//     {
//       y = -1.0 + (j+0.5)*2.0/ny;
//       legendre(y, ordery, legendre_res);
//       legendre_y.push_back(legendre_res);
//     }
//   for (unsigned int k=0; k<nz; k++)
//     {
//       z = -1.0 + (k+0.5)*2.0/nz;
//       legendre(z, orderz, legendre_res);
//       legendre_z.push_back(legendre_res);
//     }

//   unsigned int data_bin;
//   std::vector<double> legendre_products(order, 0.0);
//   double curr_x, curr_y, curr_z;
//   for (unsigned int k=0; k<nz; k++)
//     {
//       for (unsigned int j=0; j<ny; j++)
// 	{
// 	  for (unsigned int i=0; i<nx; i++)
// 	    {
// 	      data_bin = i + nx*j + nxy*k;
// 	      value = data.at(data_bin).value;
// 	      //pred = prediction.at(data_bin).value;
// 	      error = data.at(data_bin).error;
// 	      //need to add more bins if value is < 5!
// 	      //need to fit eps*theory to data instead!
// 	      if (!std::isnan(value) && !std::isinf(value))// && error != 0.0 && value != 0.0)
// 		{
// 		  //calculate coefficients
// 		  for (unsigned int o = 0; o<orderz; o++)
// 		    {
// 		      curr_z = legendre_z[k][o];
// 		      for (unsigned int n = 0; n<ordery; n++)
// 			{
// 			  curr_y = legendre_y[j][n];
// 			  for (unsigned int m = 0; m<orderx; m++)
// 			    {
// 			      curr_x = legendre_x[i][m];
// 			      legendre_products[m + orderx*n + orderxy*o] 
// 				= curr_x*curr_y*curr_z;
// 			    }
// 			}
// 		    }//end calculations of legendre products
		    
// 		  //calculate b_k and a_ik
// 		  for (unsigned int m = 0; m < order; m++)
// 		    {
// 		      if (error == 0.0)
// 			//error = pred;
// 			error = 1.0;
// 		      sqr_error = error*error;
		      
// 		      legendre_product = legendre_products[m];

// 		      bk[m] += 
// 			legendre_product 
// 			* value / sqr_error;
// 		      //the prediction already has the dv factor!
// 		      for (unsigned int n = 0; n < order; n++)
// 			{
// 			  /*
// 			    aik[n+m*order] += 
// 			    legendre_product
// 			    * pred * pred
// 			    * legendre_products[n]
// 			    / sqr_error;
// 			  */
// 			  if (n > m)
// 			    {
// 			      double temp =
// 				legendre_product
// 				* legendre_products[n]
// 				/ sqr_error;

// 			      aik[n+m*order] += temp;
// 			      aik[m+n*order] += temp;
// 			    }
// 			  else if (m == n)
// 			    aik[n+m*order] += 
// 			      legendre_product
// 			      * legendre_products[n]
// 			      / sqr_error;
// 			}
// 		    }//end calculations of b_k and a_ik

// 		}//std::isnan(value) end

// 	    }
// 	}
//     }//end loop over bins
//   progressbar(100, 50, "Determining Legendre coefficients ");
//   std::cout << std::endl;
  
//   cout << "Matrix a_ik" << endl;
//   double mean_diagonal=0.0;
//   unsigned int nlarger=0, nsmaller =0;
//   for (unsigned int i=0; i<order; i++)
//     mean_diagonal += aik.at(i*order)/order;
//   //cout << "mean diagonal element: " << mean_diagonal << endl;
//   for (unsigned int i=0; i<order; i++)
//     {
//       for (unsigned int j=0; j<order; j++)
// 	{
// 	  if (fabs(aik.at(j+i*order)) > 0.01*mean_diagonal)
// 	    {
// 	      cout << scientific << setprecision(2) << setw(12) 
// 		   << aik.at(j+i*order)
// 		   << " ";
// 	      nlarger++;
// 	    }
// 	  else
// 	    {
// 	      cout << scientific << setprecision(2) << setw(12) 
// 		   << "-"
// 		   << " ";
// 	      nsmaller++;
// 	    }
// 	}
//       cout << endl;
//     }
  
//   cout << "Vector b_k" << endl;
//   for (unsigned int i=0; i<order; i++)
//     cout << scientific << setprecision(2) << setw(12) 
// 	 << bk.at(i) << endl;
    
//   //cout << "Matrix elements which can possibly be cut: " << nsmaller << " of " << nsmaller+nlarger << endl;
//   cout << "Using ROOT Matrices to solve a_ik c_i = b_k" << endl;
//   TMatrixDSym a(order);
//   TVectorD b(order);
//   for (unsigned int i=0; i<order; i++)
//     {
//       for (unsigned int j=i; j<order; j++)
// 	a(i,j) = aik.at(i+order*j);
//       b[i] = bk.at(i);
//     }

//   bool fine = false;
//   TDecompChol decomp(a);
//   TVectorD result  = decomp.Solve(b, fine);
//   cout << "Decomposition successfull? " << (fine ? "YES" : "NO") << endl;
//   entry_vector result_coeffs(order, 0.0);
//   for (unsigned int i=0; i<order; i++)
//     {
//       cout << "c(" << setw(3) << i << ") = " << result[i] << endl;
//       result_coeffs.at(i) = result[i];
//     }  
//   return result_coeffs;
// };

// entry_vector calc_legendre_coefficients_4d_viachi2(const hist_4d& data, const hist_4d& prediction, 
// 							       unsigned int orderx, unsigned int ordery,
// 							       unsigned int orderz, unsigned int orderw)
// {
//   unsigned int orderxy = orderx*ordery;
//   unsigned int orderxyz = orderx*ordery*orderz;
//   unsigned int order = orderx*ordery*orderz*orderw;
//   unsigned int nx = data.get_bins_w(), 
//     ny = data.get_bins_x(), 
//     nz = data.get_bins_y(),
//     nw = data.get_bins_z();
//   unsigned int nxy = nx*ny;
//   unsigned int nxyz = nx*ny*nz;
//   assert(data.size() == nx*ny*nz*nw);

//   std::vector<double> legendre_res;
//   std::vector<std::vector<double> > legendre_x, legendre_y, legendre_z, legendre_w;
//   std::vector<double> aik(order*order, 0.0);
//   std::vector<double> bk(order, 0.0);

//   double x, y, z, w, value, error, pred, legendre_product, sqr_error;

//   for (unsigned int l=0; l<nw; l++)
//     {
//       w = -1.0 + (l+0.5)*2.0/nw;
//       legendre(w, orderw, legendre_res);
//       legendre_w.push_back(legendre_res);
//     }
//   for (unsigned int k=0; k<nz; k++)
//     {
//       z = -1.0 + (k+0.5)*2.0/nz;
//       legendre(z, orderz, legendre_res);
//       legendre_z.push_back(legendre_res);
//     }
//   for (unsigned int j=0; j<ny; j++)
//     {
//       y = -1.0 + (j+0.5)*2.0/ny;
//       legendre(y, ordery, legendre_res);
//       legendre_y.push_back(legendre_res);
//     }
//   for (unsigned int i=0; i<nx; i++)
//     {
//       x = -1.0 + (i+0.5)*2.0/nx;
//       legendre(x, orderx, legendre_res);
//       legendre_x.push_back(legendre_res);
//     }

//   double dv = 1.0;
//   unsigned int data_bin;
//   std::vector<double> legendre_products(order, 0.0);
//   double curr_x, curr_y, curr_z, curr_w;
//   for (unsigned int l=0; l<nw; l++)
//     {
//       progressbar(double(l)/nw*100, 50, "Determining Legendre coefficients ");
//       for (unsigned int k=0; k<nz; k++)
// 	{
// 	  for (unsigned int j=0; j<ny; j++)
// 	    {
// 	      for (unsigned int i=0; i<nx; i++)
// 		//x is proper time (i is first index, this is w in hists, this means w in hists.hh
// 		{
// 		  //problem is that the t bin should be the last (need to check this in hist_4d!!)
// 		  //then we can add up bins until value > 5.0
// 		  // if (data.irregular_binning())
// 		  //   dv = 2.0*2.0*2.0/(nw*ny*nz)*
// 		  //     (data.get_bin_hi_w(i)-data.get_bin_lo_w(i));
// 		  data_bin = i + nx*j + nxy*k + nxyz*l;
// 		  value = data.at(data_bin).value;
// 		  pred = prediction.at(data_bin).value;
// 		  error = data.at(data_bin).error;
// 		  //need to add more bins if value is < 5!
// 		  //need to fit eps*theory to data instead!
// 		  if (!std::isnan(value) && !std::isinf(value) && value > 5.0)// && error != 0.0 && value != 0.0)
// 		    {
// 		      //calculate coefficients
// 		      for (unsigned int p = 0; p<orderw; p++)
// 			{
// 			  curr_w = legendre_w[l][p];
// 			  for (unsigned int o = 0; o<orderz; o++)
// 			    {
// 			      curr_z = legendre_z[k][o];
// 			      for (unsigned int n = 0; n<ordery; n++)
// 				{
// 				  curr_y = legendre_y[j][n];
// 				  for (unsigned int m = 0; m<orderx; m++)
// 				    {
				      
// 				      curr_x = legendre_x[i][m];
// 				      legendre_products[m + orderx*n + orderxy*o + orderxyz*p] 
// 					= curr_x*curr_y*curr_z*curr_w;
// 				    }
// 				}
// 			    }
// 			}//end calculations of legendre products
		    
// 		      //calculate b_k and a_ik
// 		      for (unsigned int m = 0; m < order; m++)
// 			{
// 			  if (error == 0.0)
// 			    //error = pred;
// 			    error = 1.0;
// 			  sqr_error = error*error;
		      
// 			  legendre_product = legendre_products[m];

// 			  bk[m] += 
// 			    legendre_product 
// 			    * value * pred / sqr_error;
// 			  //the prediction already has the dv factor!
// 			  for (unsigned int n = 0; n < order; n++)
// 			    {/*
// 			       aik[n+m*order] += 
// 			       legendre_product
// 			       * pred * pred
// 			       * legendre_products[n]
// 			       / sqr_error;
// 			     */
// 			      if (n > m)
// 				{
// 				  double temp =
// 				    legendre_product
// 				    * pred * pred
// 				    * legendre_products[n]
// 				    / sqr_error;

// 				  aik[n+m*order] += temp;
// 				  aik[m+n*order] += temp;
// 				}
// 			      else if (m == n)
// 				aik[n+m*order] += 
// 				  legendre_product
// 				  * pred * pred * dv * dv 
// 				  * legendre_products[n]
// 				  / sqr_error;
// 			    }
// 			}//end calculations of b_k and a_ik

// 		    }//std::isnan(value) end

// 		}
// 	    }
// 	}
//     }//end loop over bins
//   progressbar(100, 50, "Determining Legendre coefficients ");
//   std::cout << std::endl;
  
//   cout << "Matrix a_ik" << endl;
//   double mean_diagonal=0.0;
//   unsigned int nlarger=0, nsmaller =0;
//   for (unsigned int i=0; i<order; i++)
//     mean_diagonal += aik.at(i*order)/order;
//   cout << "mean diagonal element: " << mean_diagonal << endl;
//   for (unsigned int i=0; i<order; i++)
//     {
//       for (unsigned int j=0; j<order; j++)
// 	{
// 	  if (fabs(aik.at(j+i*order)) > 0.01*mean_diagonal)
// 	    {
// 	      cout << scientific << setprecision(2) << setw(12) 
// 		   << aik.at(j+i*order)
// 		   << " ";
// 	      nlarger++;
// 	    }
// 	  else
// 	    {
// 	      cout << scientific << setprecision(2) << setw(12) 
// 		   << "-"
// 		   << " ";
// 	      nsmaller++;
// 	    }
// 	}
//       cout << endl;
//     }
//   cout << "Vector b_k" << endl;
//   for (unsigned int i=0; i<order; i++)
//     cout << scientific << setprecision(2) << setw(12) 
// 	 << bk.at(i) << endl;
    
//   cout << "Matrix elements which can possibly be cut: " << nsmaller << " of " << nsmaller+nlarger << endl;
//   //for (unsigned int i=0; i<order*order; i++)
//   //  cout << fixed << setprecision(2) << test_vector.at(i).value << " " << endl;
//   cout << "Using ROOT Matrices to solve a_ik c_i = b_k" << endl;
//   TMatrixDSym a(order);
//   TVectorD b(order);
//   for (unsigned int i=0; i<order; i++)
//     {
//       for (unsigned int j=i; j<order; j++)
// 	a(i,j) = aik.at(i+order*j);
//       b[i] = bk.at(i);
//     }
//   bool fine = false;
//   TDecompChol decomp(a);
//   TVectorD result  = decomp.Solve(b, fine);
//   cout << "Decomposition successfull? " << (fine ? "YES" : "NO") << endl;
//   entry_vector result_coeffs(order, 0.0);
//   for (unsigned int i=0; i<order; i++)
//     {
//       cout << "c(" << setw(3) << i << ") = " << result[i] << endl;
//       result_coeffs.at(i) = result[i];
//     }  
//   //results = coeffs/coeffs2;
//   //return results;
  
  
//   return result_coeffs;
// }

// //this is probably unsafe for large orders!
// entry_vector constants_from_legendre_1d(const entry_vector& coeffs)
// {
//   unsigned int n = coeffs.size();
//   entry_vector poly_coeffs(n, 0.0);
//   //this is the general formula to extract the coefficients for the monomials, can be found, for example, here:
//   //http://functions.wolfram.com/Polynomials/LegendreP/06/01/02/0006/
//   for (unsigned int i=0; i<n; i++)//for every P_i(x)
//     for (unsigned int k = 0; k <= i/2; k++)//go from 0 to i/2
//       {
// 	poly_coeffs.at(i-2*k) += 
// 	  (k%2==0 ? 1.0 : -1.0) //pow(-1.0,k) 
// 	  * TMath::Binomial(i, k)
// 	  * TMath::Binomial(2*i-2*k, i)
// 	  //* (orthonormal ? sqrt(0.5*(2.0*i + 1.0)) : 1.0)//correct for orthonormal polynomials
// 	  * coeffs.at(i).value//apply legendre coefficient
// 	  / pow(2.0, int(i));
//       }
//   return poly_coeffs;
// }

// entry_vector constants_from_legendre_4d(const entry_vector& coeffs, 
// 					unsigned int order_x, unsigned int order_y, 
// 					unsigned int order_z, unsigned int order_w)
// {
//   unsigned int n = coeffs.size();
//   assert(n == order_x*order_y*order_z*order_w);
//   //this is the general formula to extract the coefficients for the monomials, can be found, for example, here:
//   //http://functions.wolfram.com/Polynomials/LegendreP/06/01/02/0006/

//   unsigned int from_bin, to_bin;
//   entry_vector corr_x(n, 0.0);
//   for (unsigned int i=0; i<order_x; i++)
//     for (unsigned int j=0; j<order_y; j++)
//       for (unsigned int k=0; k<order_z; k++)
// 	for (unsigned int l=0; l<order_w; l++)
// 	  {
// 	    //the current coefficient we are working on
// 	    from_bin = i+order_x*j + order_x*order_y*k + order_x*order_y*order_z*l;
// 	    for (unsigned int m = 0; m <= i/2; m++)
// 	      {
// 		to_bin = (i-2*m)+order_x*j + order_x*order_y*k + order_x*order_y*order_z*l;
// 		corr_x.at(to_bin) += 
// 		  (m%2==0 ? 1.0 : -1.0)
// 		  * TMath::Binomial(i, m)
// 		  * TMath::Binomial(2*i-2*m, i)
// 		  * coeffs.at(from_bin).value
// 		  / pow(2.0, int(i));
// 	      }
// 	  }
//   entry_vector corr_y(n, 0.0);
//   for (unsigned int i=0; i<order_x; i++)
//     for (unsigned int j=0; j<order_y; j++)
//       for (unsigned int k=0; k<order_z; k++)
// 	for (unsigned int l=0; l<order_w; l++)
// 	  {
// 	    //the current coefficient we are working on
// 	    from_bin = i+order_x*j + order_x*order_y*k + order_x*order_y*order_z*l;
// 	    for (unsigned int m = 0; m <= j/2; m++)
// 	      {
// 		to_bin = i+order_x*(j-2*m) + order_x*order_y*k + order_x*order_y*order_z*l;
// 		corr_y.at(to_bin) += 
// 		  (m%2==0 ? 1.0 : -1.0)
// 		  * TMath::Binomial(j, m)
// 		  * TMath::Binomial(2*j-2*m, j)
// 		  * corr_x.at(from_bin).value
// 		  / pow(2.0, int(j));
// 	      }
// 	  }
//   entry_vector corr_z(n, 0.0);
//   for (unsigned int i=0; i<order_x; i++)
//     for (unsigned int j=0; j<order_y; j++)
//       for (unsigned int k=0; k<order_z; k++)
// 	for (unsigned int l=0; l<order_w; l++)
// 	  {
// 	    //the current coefficient we are working on
// 	    from_bin = i+order_x*j + order_x*order_y*k + order_x*order_y*order_z*l;
// 	    for (unsigned int m = 0; m <= k/2; m++)
// 	      {
// 		to_bin = i+order_x*j + order_x*order_y*(k-2*m) + order_x*order_y*order_z*l;
// 		corr_z.at(to_bin) += 
// 		  (m%2==0 ? 1.0 : -1.0)
// 		  * TMath::Binomial(k, m)
// 		  * TMath::Binomial(2*k-2*m, k)
// 		  * corr_y.at(from_bin).value
// 		  / pow(2.0, int(k));
// 	      }
// 	  }
//   entry_vector poly_coeffs(n, 0.0);
//   for (unsigned int i=0; i<order_x; i++)
//     for (unsigned int j=0; j<order_y; j++)
//       for (unsigned int k=0; k<order_z; k++)
// 	for (unsigned int l=0; l<order_w; l++)
// 	  {
// 	    //the current coefficient we are working on
// 	    from_bin = i+order_x*j + order_x*order_y*k + order_x*order_y*order_z*l;
// 	    for (unsigned int m = 0; m <= l/2; m++)
// 	      {
// 		to_bin = i+order_x*j + order_x*order_y*k + order_x*order_y*order_z*(l-2*m);
// 		poly_coeffs.at(to_bin) += 
// 		  (m%2==0 ? 1.0 : -1.0)
// 		  * TMath::Binomial(l, m)
// 		  * TMath::Binomial(2*l-2*m, l)
// 		  * corr_z.at(from_bin).value
// 		  / pow(2.0, int(l));
// 	      }
// 	  }

//   return poly_coeffs;
// }

// entry_vector constants_from_legendre_3d(const entry_vector& coeffs, 
// 					unsigned int order_x, unsigned int order_y, unsigned int order_z)
// {
//   unsigned int n = coeffs.size();
//   assert(n == order_x*order_y*order_z);
//   //this is the general formula to extract the coefficients for the monomials, can be found, for example, here:
//   //http://functions.wolfram.com/Polynomials/LegendreP/06/01/02/0006/
//   unsigned int from_bin, to_bin;
//   entry_vector corr_x(n, 0.0);
//   for (unsigned int i=0; i<order_x; i++)
//     for (unsigned int j=0; j<order_y; j++)
//       for (unsigned int k=0; k<order_z; k++)
// 	{
// 	  //the current coefficient we are working on
// 	  from_bin = i+order_x*j + order_x*order_y*k;
// 	  for (unsigned int m = 0; m <= i/2; m++)
// 	    {
// 	      to_bin = (i-2*m)+order_x*j + order_x*order_y*k;
// 	      corr_x.at(to_bin) += 
// 		(m%2==0 ? 1.0 : -1.0)
// 		* TMath::Binomial(i, m)
// 		* TMath::Binomial(2*i-2*m, i)
// 		* coeffs.at(from_bin).value
// 		/ pow(2.0, int(i));
// 	    }
// 	}
//   entry_vector corr_y(n, 0.0);
//   for (unsigned int i=0; i<order_x; i++)
//     for (unsigned int j=0; j<order_y; j++)
//       for (unsigned int k=0; k<order_z; k++)
// 	{
// 	  //the current coefficient we are working on
// 	  from_bin = i+order_x*j + order_x*order_y*k;
// 	  for (unsigned int m = 0; m <= j/2; m++)
// 	    {
// 	      to_bin = i+order_x*(j-2*m) + order_x*order_y*k;
// 	      corr_y.at(to_bin) += 
// 		(m%2==0 ? 1.0 : -1.0)
// 		* TMath::Binomial(j, m)
// 		* TMath::Binomial(2*j-2*m, j)
// 		* corr_x.at(from_bin).value
// 		/ pow(2.0, int(j));
// 	    }
// 	}
//   entry_vector corr_z(n, 0.0);
//   for (unsigned int i=0; i<order_x; i++)
//     for (unsigned int j=0; j<order_y; j++)
//       for (unsigned int k=0; k<order_z; k++)
// 	{
// 	  //the current coefficient we are working on
// 	  from_bin = i+order_x*j + order_x*order_y*k;
// 	  for (unsigned int m = 0; m <= k/2; m++)
// 	    {
// 	      to_bin = i+order_x*j + order_x*order_y*(k-2*m);
// 	      corr_z.at(to_bin) += 
// 		(m%2==0 ? 1.0 : -1.0)
// 		* TMath::Binomial(k, m)
// 		* TMath::Binomial(2*k-2*m, k)
// 		* corr_y.at(from_bin).value
// 		/ pow(2.0, int(k));
// 	    }
// 	}
//   return corr_z;
// };

// entry_vector correct_constants_3d(const entry_vector& constants, 
// 				  unsigned order_cos_theta, unsigned int order_phi, unsigned int order_cos_psi,
// 				  double phi_min, double phi_max)
// {
//   cout << "Correcting polynomial constants for phi" << endl;
//   entry_vector phi_result(constants.size(), 0.0);
//   assert(constants.size() == order_cos_theta*order_phi*order_cos_psi);
//   unsigned int the_bin;
//   for (unsigned int j = 0; j<order_cos_theta; j++)
//     for (unsigned int k = 0; k<order_phi; k++)
//       for (unsigned int l = 0; l<order_cos_psi; l++)
// 	{
// 	  the_bin = j + order_cos_theta*k + order_cos_theta*order_phi*l;
// 	  phi_result.at(the_bin) = constants.at(the_bin)/pow(TMath::Pi(),int(k));//correct for -pi..pi instead of -1..+1	      
// 	}
//   return phi_result;
// };

// entry_vector correct_constants_1d(const entry_vector& constants, 
// 				  unsigned int order_t, double t_min, double t_max)
// {
//   assert(constants.size() == order_t);
//   cout << "Correcting polynomial constants for t" << endl;
//   double c = 2.0/(t_max-t_min);
//   double d = -2.0*t_min/(t_max-t_min)-1.0;
//   double e = d/c;
//   entry_vector t_result(constants.size(), 0.0);
//   unsigned int the_bin, new_bin;
//   for (unsigned int i = 0; i<order_t; i++)
//     {
//       the_bin = i;
//       //correct for t_min..t_max instead of -1..+1
//       entry constant_dash = constants.at(the_bin);
//       for (unsigned int m=0; m<=i; m++)
// 	{
// 	  new_bin = m;
// 	  t_result.at(new_bin) += 
// 	    constant_dash
// 	    *pow(c,int(i))
// 	    *TMath::Binomial(i,m)
// 	    *pow(e,int(i-m));
// 	}
//     }
//   return t_result;
// }

// entry_vector correct_constants_4d(const entry_vector& constants, 
// 				  unsigned int order_t, unsigned order_cos_theta, 
// 				  unsigned int order_phi, unsigned int order_cos_psi,
// 				  double t_min, double t_max, double phi_min, double phi_max)
// {
//   //entry_vector result(constants);
//   cout << "Correcting polynomial constants for phi" << endl;
//   entry_vector phi_result(constants.size(), 0.0);
//   assert(constants.size() == order_t*order_cos_theta*order_phi*order_cos_psi);
//   unsigned int the_bin;
//   for (unsigned int i = 0; i<order_t; i++)
//     for (unsigned int j = 0; j<order_cos_theta; j++)
//       for (unsigned int k = 0; k<order_phi; k++)
// 	for (unsigned int l = 0; l<order_cos_psi; l++)
// 	  {
// 	    the_bin = i + order_t*j + order_t*order_cos_theta*k + order_t*order_cos_theta*order_phi*l;
// 	    phi_result.at(the_bin) = constants.at(the_bin)/pow(TMath::Pi(),int(k));//correct for -pi..pi instead of -1..+1	      
// 	  }
//   cout << "Correcting polynomial constants for t" << endl;
//   double c = 2.0/(t_max-t_min);
//   double d = -2.0*t_min/(t_max-t_min)-1.0;
//   double e = d/c;
//   unsigned int new_bin;
//   entry_vector t_result(constants.size(), 0.0);
//   for (unsigned int i = 0; i<order_t; i++)
//     for (unsigned int j = 0; j<order_cos_theta; j++)
//       for (unsigned int k = 0; k<order_phi; k++)
// 	for (unsigned int l = 0; l<order_cos_psi; l++)
// 	  {
// 	    the_bin = i + order_t*j + order_t*order_cos_theta*k + order_t*order_cos_theta*order_phi*l;
// 	    //correct for t_min..t_max instead of -1..+1
// 	    entry constant_dash = phi_result.at(the_bin);
// 	    for (unsigned int m=0; m<=i; m++)
// 	      {
// 		new_bin = m + order_t*j + order_t*order_cos_theta*k + order_t*order_cos_theta*order_phi*l;
// 		t_result.at(new_bin) += 
// 		  constant_dash
// 		  *pow(c,int(i))
// 		  *TMath::Binomial(i,m)
// 		  *pow(e,int(i-m));
// 	      }
// 	  }
//   return t_result;
// }


//calculates all chebychev polynomials up to including degree n
void chebychev(double x, int n, std::vector<double>& results)
{
  assert(n >= 0);
  results.clear();
  results.reserve(n+1);
  results.push_back(1.0);
  results.push_back(x);
  for (int i=2; i<n+1; i++)
    results.push_back(2.0*x*results.at(i-1) - results.at(i-2));
}

//calculates all legendre polynomials up to including degree n
void legendre(double x, int n, std::vector<double>& results)
{
  assert(n >= 0);
  results.clear();
  //
  //results = std::vector<double> ();
  results.reserve(n+1);
  results.push_back(1.0);
  results.push_back(x);//ok
  for (int i=2; i<n+1; i++)
    results.push_back( ((2.0*i-1.0) * x * results.at(i-1) - (i-1.0) * results.at(i-2)) /double(i)  );
}

void orthonormal_legendre(double x, int n, std::vector<double>& results)
{
  legendre(x, n, results);
  for (unsigned int i=0; i<results.size(); i++)
    results.at(i) *= sqrt(0.5*(2.0*i+1.0));
}


//calculates chebychev polynomial of order n
double chebychev(double x, int n)
{
  assert(n >= 0);
  if (n == 0)
    return 1;
  else if (n == 1)
    return x;
  else
    return 2.0*x*chebychev(x, n-1) - chebychev(x, n-2);
}

double integral_x_to_n_times_exp_minus_x(double x, double tau, int n)
{
  assert(n>=0);
  if (n == 0)
    return -tau*exp(-x/tau);
  else
    return pow(x, n)*(-tau*exp(-x/tau)) + n * tau * integral_x_to_n_times_exp_minus_x(x, tau, n-1);
}

//this methods are helpful for the phi integration
//calculates int x^n * sin(x) dx
double integral_x_to_n_times_sin_x(double x, int n)
{
  assert(n >= 0);
  if (n == 0)
    return -cos(x);
  else
    return -pow(x, n)*cos(x) + n*integral_x_to_n_times_cos_x(x, n-1);
}

//calculates int x^n * cos(x) dx
double integral_x_to_n_times_cos_x(double x, int n)
{
  assert(n >= 0);
  if (n == 0)
    return sin(x);
  else
    return pow(x,n)*sin(x) - n*integral_x_to_n_times_sin_x(x, n-1);
}

//calculates int x^n * sin(2x) dx
double integral_x_to_n_times_sin_2x(double x, int n)
{
  assert(n >= 0);
  if (n == 0)
    return -0.5*cos(2.0*x);
  else
    return -pow(x,n)*0.5*cos(2.0*x)
      +0.5*n*integral_x_to_n_times_cos_2x(x,n-1);
}

//calculates int x^n * cos(2x) dx
double integral_x_to_n_times_cos_2x(double x, int n)
{
  assert(n >= 0);
  if (n == 0)
    return 0.5*sin(2.0*x);
  else
    return +0.5*pow(x,n)*sin(2.0*x)
      -0.5*n*integral_x_to_n_times_sin_2x(x,n-1);
}

//calculates int x^n * cos(x)^2 dx
double integral_x_to_n_times_cos_x_2(double x, int n)
{
  assert(n >= 0);
  return +1.0/(1.0 + n)*pow(x,n+1)*cos(x)*cos(x)
    +1.0/(1.0+n)*integral_x_to_n_times_sin_2x(x, n+1);
}

//calculates int x^n * sin(x)^2 dx
double integral_x_to_n_times_sin_x_2(double x, int n)
{
  assert(n >= 0);
  return +1.0/(1.0 + n)*pow(x,n+1)*sin(x)*sin(x)
    -1.0/(1.0+n)*integral_x_to_n_times_sin_2x(x, n+1);
}

//calculates int x^n * asin(x) dx
double integral_x_to_n_times_asin_x(double x, int n)
{
  assert(n >= 0);
  if (n == 0)
    return x*asin(x)+sqrt(1-x*x);
  else
    return 1.0/(n+1.0)*pow(x,n)*(x*asin(x)+sqrt(1-x*x))
      -n*integral_x_to_n_times_sqrt_1_minus_x2(x, n-1);
}

//calculates int x^n * sqrt(1-x^2) dx
double integral_x_to_n_times_sqrt_1_minus_x2(double x, int n)
{
  assert(n >= 0);
  if (n == 0)
    return 0.5*asin(x)+0.5*x*sqrt(1-x*x);
  else
    return 2.0/(n+2.0)*pow(x, n)*(0.5*asin(x)+0.5*x*sqrt(1-x*x))
      -n/(n+2.0)*integral_x_to_n_times_asin_x(x, n-1);
}


void test_integrals()
{
  cout << "Starting Test" << endl;
  //sin test
  unsigned int nbins=100;
  for (unsigned int i = 0; i < nbins; i++)
    {
      double x = -1.0 + 2.0*i/double(nbins); 
      double a = integral_x_to_n_times_sin_x_2(x, 0);
      double b = 0.5*(x-0.5*sin(2.0*x));
      if (fabs(a-b) > 1.0e-8)
	{
	  cout << a;
	  cout << " != ";
	  cout << b;
	  cout << endl;
	}
    }
  for (unsigned int i = 0; i < nbins; i++)
    {
      double x = -1.0 + 2.0*i/double(nbins); 
      double a = integral_x_to_n_times_sin_x_2(x, 1);
      double b = -(2*x*sin(2*x)+cos(2*x)-2*x*x)/8;
      if (fabs(a-b) > 1.0e-8)
	{
	  cout << a;
	  cout << " != ";
	  cout << b;
	  cout << endl;
	}
    }
  for (unsigned int i = 0; i < nbins; i++)
    {
      double x = -1.0 + 2.0*i/double(nbins); 
      double a = integral_x_to_n_times_sin_x_2(x, 2);
      double b = -((6*x*x-3)*sin(2*x)+6*x*cos(2*x)-4*x*x*x)/24;
      if (fabs(a-b) > 1.0e-8)
	{
	  cout << a;
	  cout << " != ";
	  cout << b;
	  cout << endl;
	}
    }
  for (unsigned int i = 0; i < nbins; i++)
    {
      double x = -1.0 + 2.0*i/double(nbins); 
      double a = integral_x_to_n_times_sin_x_2(x, 3);
      double b = -((4*x*x*x-6*x)*sin(2*x)+(6*x*x-3)*cos(2*x)-2*x*x*x*x)/16;
      if (fabs(a-b) > 1.0e-8)
	{
	  cout << a;
	  cout << " != ";
	  cout << b;
	  cout << endl;
	}
    }
  for (unsigned int i = 0; i < nbins; i++)
    {
      double x = -1.0 + 2.0*i/double(nbins); 
      double a = integral_x_to_n_times_sin_x_2(x, 1);
      double b = -(2*x*sin(2*x)+cos(2*x)-2*x*x)/8;
      if (fabs(a-b) > 1.0e-8)
	{
	  cout << a;
	  cout << " != ";
	  cout << b;
	  cout << endl;
	}
    }
  //cos test
  for (unsigned int i = 0; i < nbins; i++)
    {
      double x = -1.0 + 2.0*i/double(nbins); 
      double a = integral_x_to_n_times_cos_x_2(x, 0);
      double b = (sin(2.0*x)*0.5+x)*0.5;
      if (fabs(a-b) > 1.0e-8)
	{
	  cout << a;
	  cout << " != ";
	  cout << b;
	  cout << endl;
	}
    }
  std::vector<double> res;
  chebychev(0.43, 13, res);
  for (unsigned int i = 0; i < 14; i++)
    {
      if (fabs(res.at(i) - chebychev(0.43,i)) > 1.0e-8)
	cout << "CHEBYCHEV: " << res.at(i) << " != " <<  chebychev(0.43,i)<< endl;
    }
  legendre(0.43, 13, res);
  for (unsigned int i = 0; i < 14; i++)
    {
      if (fabs(res.at(i) - ROOT::Math::legendre(i, 0.43)) > 1.0e-8)
	cout << "LEGENDREPOLY " << i << ": " << res.at(i) << " != " <<  ROOT::Math::legendre(i, 0.43)<< endl;
    }

  cout << "Finished Test" << endl;
}

void test_legendre()
{
  TRandom3 rnd;
  TCanvas *canvas = new TCanvas("canvas", "canvas",1600,1200);
  canvas->cd();
  unsigned int nbins = 200;
  TH1D * data = new TH1D("data", "data", nbins, -1.0, 1.0);
  for (unsigned int i =0; i<200000; i++)
    {
      data->Fill(rnd.Gaus(0.0,0.5));
      bool finished = false;
      double x;
      while (!finished)
	{
	  x = rnd.Rndm();
	  if (rnd.Rndm() < x)
	    {
	      x = x*2.0-1.0;
	      finished = true;
	    }
	}
      data->Fill(x);
    }
  //data->FillRandom("gaus", 200000);
  data->Draw("e");

  unsigned int n = 20;

  std::vector<double> l;
  double dv = 2.0/nbins;
  std::vector<double> coeffs(n+1, 0.0);
  for (unsigned int j=0; j<nbins; j++)
    {
      double x = -1.0 + (j+0.5)*2.0/nbins;
      //orthonormal_legendre(x, n, l);
      legendre(x, n, l);
      double y = data->GetBinContent(j+1);
      for (unsigned int i = 0; i<n+1; i++)
	coeffs.at(i) += y*l.at(i) * dv * (0.5*(2.0*i+1.0));
    }
  for (unsigned int i = 0; i<n+1; i++)
    cout << "coefficient " << i << "=" << coeffs.at(i) << endl;
  TH1D * interpol = new TH1D("inter", "inter", nbins, -1.0, 1.0);
  for (unsigned int j=0; j<nbins; j++)
    {
      double x = -1.0 + (j+0.5)*2.0/nbins;
      //orthonormal_legendre(x, n, l);
      legendre(x, n, l);
      double y = 0.0;
      for (unsigned int i = 0; i<n+1; i++)
	y += coeffs.at(i) * l.at(i);
      interpol->SetBinContent(j+1, y);
    }  
  interpol->SetLineColor(2);
  interpol->Draw("same l");

  TH1D * pol1d = new TH1D("pol1d", "pol1d", nbins, -1.0, 1.0);
  std::vector<double> poly_coeffs;
  for (unsigned int i=0; i<n+1; i++)
    cout << "Poly Coefficient " << i << " = " << poly_coeffs.at(i) << endl;
  for (unsigned int j=0; j<nbins; j++)
    {
      double x = -1.0 + (j+0.5)*2.0/nbins;
      double y = 0.0;
      for (unsigned int i = 0; i<n+1; i++)
	y += poly_coeffs.at(i) * pow(x, int(i));
      pol1d->SetBinContent(j+1, y);
    }  
  pol1d->SetLineColor(3);
  pol1d->SetLineStyle(2);
  pol1d->Draw("same l");

  canvas->Update();
  canvas->Print("legendre.eps", "eps");
  delete data;
  delete interpol;
  delete pol1d;
  delete canvas;
}

void test_legendre_2pi()
{
  TRandom3 rnd;
  TCanvas *canvas = new TCanvas("canvas", "canvas",1600,1200) ;
  canvas->cd();
  unsigned int nbins = 200;
  TH1D * data = new TH1D("data", "data", nbins, -TMath::Pi(), TMath::Pi());
  for (unsigned int i =0; i<200000; i++)
    {
      data->Fill(rnd.Gaus(0.0,0.5));
      bool finished = false;
      double x;
      while (!finished)
	{
	  x = rnd.Rndm();
	  if (rnd.Rndm() < x)
	    {
	      x = x*2.0*TMath::Pi()-TMath::Pi();
	      finished = true;
	    }
	}
      data->Fill(x);
    }
  //data->FillRandom("gaus", 200000);
  data->Draw("e");

  unsigned int n = 20;

  std::vector<double> l;
  double dv = 2.0/nbins;
  std::vector<double> coeffs(n+1, 0.0);
  for (unsigned int j=0; j<nbins; j++)
    {
      double x = -1.0 + (j+0.5)*2.0/nbins; //this transforms correctly to -1...+1
      //orthonormal_legendre(x, n, l);
      legendre(x, n, l);
      double y = data->GetBinContent(j+1);
      for (unsigned int i = 0; i<n+1; i++)
	coeffs.at(i) += y*l.at(i) * dv * (0.5*(2.0*i+1.0));
    }
  for (unsigned int i = 0; i<n+1; i++)
    cout << "coefficient " << i << "=" << coeffs.at(i) << endl;
  TH1D * interpol = new TH1D("inter", "inter", nbins, -TMath::Pi(), +TMath::Pi());
  for (unsigned int j=0; j<nbins; j++)
    {
      double x = -1.0 + (j+0.5)*2.0/nbins;
      //double x = (-1.0 + (j+0.5)*2.0/nbins) / TMath::Pi();
      
      //double x = -TMath::Pi() + (j+0.5)*2.0*TMath::Pi()/nbins;
      //orthonormal_legendre(x, n, l);
      legendre(x, n, l);
      double y = 0.0;
      for (unsigned int i = 0; i<n+1; i++)
	y += coeffs.at(i) * l.at(i);
      interpol->SetBinContent(j+1, y);
    }  
  interpol->SetLineColor(2);
  interpol->Draw("same l");

  TH1D * pol1d = new TH1D("pol1d", "pol1d", nbins, -TMath::Pi(), +TMath::Pi());
  std::vector<double> poly_coeffs;
  for (unsigned int i=0; i<n+1; i++)
    cout << "Poly Coefficient " << i << " = " << poly_coeffs.at(i) << endl;
  for (unsigned int j=0; j<nbins; j++)
    {
      double x = -1.0 + (j+0.5)*2.0/nbins;
      //double x = (-1.0 + (j+0.5)*2.0/nbins)/TMath::Pi();
      //double x = xc/TMath::Pi();
      double y = 0.0;
      for (unsigned int i = 0; i<n+1; i++)
	y += poly_coeffs.at(i) * pow(x, int(i));
      pol1d->SetBinContent(j+1, y);
    }  
  pol1d->SetLineColor(3);
  pol1d->SetLineStyle(2);
  pol1d->Draw("same l");

  canvas->Update();
  canvas->Print("legendre2pi.eps", "eps");
  delete data;
  delete interpol;
  delete pol1d;
  delete canvas;
}

void test_chebychev()
{
  TRandom3 rnd;
  TCanvas *canvas = new TCanvas("canvas", "canvas",1600,1200) ;
  canvas->cd();
  unsigned int nbins = 200;
  TH1D * data = new TH1D("data", "data", nbins, -1.0, 1.0);
  for (unsigned int i =0; i<200000; i++)
    {
      data->Fill(rnd.Gaus(0.0,0.5));
      bool finished = false;
      double x;
      while (!finished)
	{
	  x = rnd.Rndm();
	  if (rnd.Rndm() < x)
	    {
	      x = x*2.0-1.0;
	      finished = true;
	    }
	}
      data->Fill(x);
    }
  //data->FillRandom("gaus", 200000);
  data->Draw("e");

  unsigned int n = 20;

  std::vector<double> l;
  double dv = 2.0/nbins;
  std::vector<double> coeffs(n+1, 0.0);
  for (unsigned int j=0; j<nbins; j++)
    {
      double x = -1.0 + (j+0.5)*2.0/nbins;
      chebychev(x, n, l);
      double y = data->GetBinContent(j+1);
      for (unsigned int i = 0; i<n+1; i++)
	coeffs.at(i) += y*l.at(i) * dv / sqrt(1-x*x) / (i==0 ? TMath::Pi() : 0.5*TMath::Pi() );
    }
  for (unsigned int i = 0; i<n+1; i++)
    cout << "coefficient " << i << "=" << coeffs.at(i) << endl;
  TH1D * interpol = new TH1D("inter", "inter", nbins, -1.0, 1.0);
  for (unsigned int j=0; j<nbins; j++)
    {
      double x = -1.0 + (j+0.5)*2.0/nbins;
      chebychev(x, n, l);
      double y = 0.0;
      for (unsigned int i = 0; i<n+1; i++)
	y += coeffs.at(i) * l.at(i);
      interpol->SetBinContent(j+1, y);
    }  
  interpol->SetLineColor(2);
  interpol->Draw("same l");

  canvas->Update();
  canvas->Print("chebychev.eps", "eps");
  delete data;
  delete interpol;
  //delete pol1d;
  delete canvas;
}

void test_legendre_2d()
{
  //gStyle->SetPalette(1);
  TRandom3 rnd;
  TCanvas *canvas = new TCanvas("canvas", "canvas",1600,1200) ;
  canvas->cd();
  unsigned int nbins = 200;
  TH2D * data = new TH2D("data", "data", nbins, -1.0, 1.0, nbins, -1.0, 1.0);
  for (unsigned int i =0; i<1000000; i++)
    {
      //double r = rnd.Gaus(0.7,0.2);
      //double phi = rnd.Rndm() * 2.0 * TMath::Pi();
      //data->Fill(fabs(r)*cos(phi),fabs(r)*sin(phi));
      data->Fill(rnd.Gaus(-0.1,0.3), rnd.Gaus(-0.1,0.3));
      data->Fill(rnd.Gaus(0.5,0.3), rnd.Gaus(0.5,0.3));
    }
  //data->FillRandom("gaus", 200000);
  data->Draw("colz");
  canvas->Update();
  canvas->Print("2d.eps", "eps");
  unsigned int n = 10;
  std::vector<double> legendre_x, legendre_y;
  double dv = 2.0/nbins*2.0/nbins;
  std::vector< std::vector<double > > coeffs(n+1, std::vector<double >(n+1, 0.0));
  for (unsigned int i=0; i<nbins; i++)
    {
      double x = -1.0 + (i+0.5)*2.0/nbins;
      for (unsigned int j=0; j<nbins; j++)
	{
	  double y = -1.0 + (j+0.5)*2.0/nbins;	  
	  double z = data->GetBinContent(i+1,j+1);
	  //orthonormal_legendre(x, n, legendre_x);
	  //orthonormal_legendre(y, n, legendre_y);
	  legendre(x, n, legendre_x);
	  legendre(y, n, legendre_y);	  
	  for (unsigned int l = 0; l<n+1; l++)
	    for (unsigned int m = 0; m<n+1; m++)
	      coeffs.at(l).at(m) += z*legendre_x.at(l)*legendre_y.at(m)*dv* (0.5*(2.0*l+1.0)) * (0.5*(2.0*m+1.0));
	}
    }
  for (unsigned int l = 0; l<n+1; l++)
    for (unsigned int m = 0; m<n+1; m++)
      cout << "coefficient (" << l << ", " << m << ") =" << coeffs.at(l).at(m) << endl;


  TH2D * interpol = new TH2D("inter", "inter", nbins, -1.0, 1.0, nbins, -1.0, 1.0);
  for (unsigned int i=0; i<nbins; i++)
    {
      double x = -1.0 + (i+0.5)*2.0/nbins;
      for (unsigned int j=0; j<nbins; j++)
	{
	  double y = -1.0 + (j+0.5)*2.0/nbins;
	  //orthonormal_legendre(x, n, legendre_x);
	  //orthonormal_legendre(y, n, legendre_y);
	  legendre(x, n, legendre_x);
	  legendre(y, n, legendre_y);

	  double z = 0.0;
	  for (unsigned int l = 0; l<n+1; l++)
	    for (unsigned int m = 0; m<n+1; m++)
	      z += coeffs.at(l).at(m)*legendre_x.at(l)*legendre_y.at(m);
	  interpol->SetBinContent(i+1, j+1, z);
	}
    }
  //calculate chi^2
  double chi2 = 0.0;
  for (unsigned int i=0; i<nbins; i++)
    {
      for (unsigned int j=0; j<nbins; j++)
	{
	  double d = data->GetBinContent(i+1, j+1);
	  double dd = data->GetBinError(i+1, j+1);
	  double f = interpol->GetBinContent(i+1, j+1);
	  if (dd != 0.0)
	    chi2 += (d-f)*(d-f)/dd/dd;
	}
    }  
  cout << "Found a chi^2 of " << chi2 << " with " << nbins*nbins << " bins." << endl;
  //cleanup
  interpol->Draw("colz");
  canvas->Update();
  canvas->Print("2dinter.eps", "eps");
  delete data;
  delete interpol;
  delete canvas;
}

void progressbar(unsigned int percent, unsigned int nchar, std::string prefix, unsigned int bar_index)
{
  if (bar_index)
    cout << "\033[" << bar_index << "A";//cursor up
  cout << "\033[" << 255 << "D";//delete line
  //cout << '\r';  
  cout << prefix << setw(3) << percent << "% [";
  cout << "\033[0;32m";
  for (unsigned int i=0; i<nchar; i++)
    if (round(double(i)/nchar*100) <= percent)
      cout << "#";
    else
      cout << " ";
  cout << "\033[0m";
  cout << "]";
  if (bar_index)
    cout << "\033[" << bar_index << "B";//cursor down
  flush(cout);
}

double legendre_poly(unsigned int l, double x)
{
  return ROOT::Math::legendre(l, x);
}

double assoc_legendre(unsigned int l, unsigned int m, double x)
{
  return ROOT::Math::assoc_legendre(l, m, x);
  //return ROOT::Math::sph_legendre(l, m, x);
}

double spherical_harmonic_re(unsigned int l, unsigned int m, double cos_theta, double phi)
{
  return sqrt((2*l+1.0)/(4.0*TMath::Pi()))*sqrt(TMath::Factorial(l-m)/TMath::Factorial(l+m))
    *assoc_legendre(l, m, cos_theta)
    *cos(m*phi);
}

double spherical_harmonic_im(unsigned int l, unsigned int m, double cos_theta, double phi)
{
  return sqrt((2*l+1.0)/(4.0*TMath::Pi()))*sqrt(TMath::Factorial(l-m)/double(TMath::Factorial(l+m)))
    *assoc_legendre(l, m, cos_theta)
    *sin(m*phi);
}

/*void calculate_transversity_angles(const LorentzVector& mu_minus, 
				   const LorentzVector& mu_plus, 
				   const LorentzVector& k_minus, 
				   const LorentzVector& k_plus, 
				   double & cos_theta, double & angle_phi, double & cos_psi)
{
  LorentzVector bs = mu_plus + mu_minus + k_plus + k_minus;
  LorentzVector jpsi = mu_plus + mu_minus;
  LorentzVector phi = k_plus + k_minus;
  LorentzBoost jpsiboost(jpsi.BoostToCM());//boosttocm actually gives beta vector
  LorentzBoost phiboost(phi.BoostToCM());//pxyzm has boost to cm pxyzm not?
  //boosted quantities in jpsi system
  LorentzVector mu_minus_d = jpsiboost(mu_minus);
  LorentzVector mu_plus_d = jpsiboost(mu_plus);
  LorentzVector k_minus_d = jpsiboost(k_minus);
  LorentzVector k_plus_d = jpsiboost(k_plus);
  LorentzVector phi_d = jpsiboost(phi);
  //boosted quantities in phi system
  LorentzVector k_plus_dd = phiboost(k_plus);
  LorentzVector jpsi_dd = phiboost(jpsi);

  //cos theta is calculated in the cms of the jpsi
  //cos theta
  Vector3 normalxy = k_minus_d.Vect().Cross(k_plus_d.Vect());
  normalxy = normalxy.Unit();//normals should be normalized
  cos_theta = (normalxy.Dot(mu_plus_d.Vect()) 
	       / sqrt(normalxy.Mag2()) 
	       / sqrt(mu_plus_d.Vect().Mag2()));
  //phi is calculated in the cms of the jpsi
  //cos phi 
  Vector3 mu_perp = normalxy * mu_plus_d.Vect().Dot(normalxy);
  Vector3 mu_parallel = mu_plus_d.Vect() - mu_perp;
  angle_phi = acos(mu_parallel.Dot(phi_d.Vect()) 
		   / sqrt(mu_parallel.Mag2()) 
		   / sqrt(phi_d.Vect().Mag2()));
  //now get orientation
  Vector3 checknormal = phi_d.Vect().Cross(mu_parallel);
  bool samedirection = (checknormal.Dot(normalxy) > 0.0);
  if (!samedirection)
    angle_phi = -angle_phi;

  //cos psi is calculated in the cms of the phi
  //cos psi
  cos_psi = k_plus_dd.Vect().Dot(-jpsi_dd.Vect())
    /sqrt(k_plus_dd.Vect().Mag2())
    /sqrt(jpsi_dd.Vect().Mag2());

}*/

//this is not normalized!
double convoluted_exp(double ct, double sigma_ct, double ctau)
{
  //use error function instead of numeric integration, should be a lot faster
  if (ct > -6.0*sigma_ct)//this safety is important for the plotting
    {
      const double sqrt2 = TMath::Sqrt(2.0);
      //  const double norm_factor = 2.0/TMath::Sqrt(TMath::Pi());
      // Erf(x) = (2/sqrt(pi)) Integral(exp(-t^2))dt between 0 and x, already taken into account?
      double exp_part = exp(-ct/ctau + sigma_ct*sigma_ct/(2.0*ctau*ctau));
      double erfc = TMath::Erfc(sigma_ct/ctau/sqrt2 - ct/sigma_ct/sqrt2);
      // / norm_factor; the normalization here is already done?
      //the expression below was normalized, but it shouldnt have been
      //double prob = 0.5/ctau * exp * erfc;//is the factor 2 correct? yes.
      double prob = 0.5 * exp_part * erfc;//is the factor 2 correct? yes.
      if (std::isnan(prob) || prob < 0.0) 
	{
	  cout << "CONVEXP gives prob: " << prob 
	       << "; ct/sigma_ct/ctau/exp/erfc/exparg/erfarg: " << ct 
	       << " " << sigma_ct << endl;
	  cout << "ctau/exp" << ctau << " " << exp_part << endl;
	  cout << "erfc/exparg/erfarg" << erfc << " " 
	       << -ct/ctau + sigma_ct*sigma_ct/(2.0*ctau*ctau) << " " 
	       << sigma_ct/ctau/sqrt2 - ct/sigma_ct/sqrt2 << endl;
	  prob = 0.0;
	}
      return prob;
    }
  else return 0.0;//safety added
}

//this is not normalized! 
double convoluted_exp_sl_k(double ct, double sigma_ct, double ctau, double m_reso, double k)
{

  //use error function instead of numeric integration, should be a lot faster
  if (ct > -6.0*sigma_ct)//this safety is important for the plotting
    {
      const double sqrt2 = TMath::Sqrt(2.0);
      //  const double norm_factor = 2.0/TMath::Sqrt(TMath::Pi());
      // Erf(x) = (2/sqrt(pi)) Integral(exp(-t^2))dt between 0 and x, already taken into account?
      double exp_part = exp(-ct*k/ctau + sigma_ct*sigma_ct*k*k/(2.0*ctau*ctau) - m_reso*k/ctau);
      double erfc = TMath::Erfc(sigma_ct*k/ctau/sqrt2 - ct/sigma_ct/sqrt2 - m_reso/(sqrt2*sigma_ct));
      // / norm_factor; the normalization here is already done?
      //the expression below was normalized, but it shouldnt have been
      //double prob = 0.5/ctau * exp * erfc;//is the factor 2 correct? yes.
      double prob = 0.5 * exp_part * erfc;//is the factor 2 correct? yes.
      if (std::isnan(prob) || prob < 0.0)
        {
          cout << "CONVEXP gives prob: " << prob
               << "; ct/sigma_ct/ctau/exp/erfc/exparg/erfarg: " << ct
               << " " << sigma_ct << endl;
          cout << "ctau/exp" << ctau << " " << exp_part << endl;
          cout << "erfc/exparg/erfarg" << erfc << " "
               << -ct/ctau + sigma_ct*sigma_ct/(2.0*ctau*ctau) << " "
               << sigma_ct/ctau/sqrt2 - ct/sigma_ct/sqrt2 << endl;
          prob = 0.0;
        }
      return prob;
    }
  else return 0.0;//safety added
}

//void fitter::ct_tagging_part(measurement* meas, double& one, double& two)
//this function calculates the convolution of the product exp(-t/tau) * sin(deltam*t) and exp(-t/tau) * cos(deltam*t)
void convoluted_exp_sincos_sl_k(double ct, double sigma, double tau, double deltam, double k, double m, double& one, double& two)
{
  double sigma_2 = sigma*sigma;
  double tau_2 = tau*tau;
  double k_2 = k*k;
  double deltam_2 = deltam*deltam;
  const double sqrt2 = TMath::Sqrt(2.0);

  const bool cdf_way = false;
  if (!cdf_way)
    {
      //safety for large negative decay times 
      if (ct < - 6.0*sigma)
        {
          one = 0.0;
          two = 0.0;
          return;
        }

      double c1 = 0.5;
      //      double exp1arg = 0.5*sigma_2/tau_2 - 0.5*deltam_2*sigma_2 - ct/tau;
      //      double exp1arg = 0.5*sigma_2*(1.0/tau_2 - deltam_2) - ct/tau;
      //      double exp1arg = 0.5*sigma_2*(k_2/tau_2 - deltam_2) - ct*k/tau -k*m/tau;
      double exp1arg = 0.5*sigma_2*k_2*(1.0/tau_2 - deltam_2) - ct*k/tau -k*m/tau;

      double exp1 = exp(exp1arg);

      //      double exp2arg = -deltam*(ct - sigma_2/tau);
      //      double exp2arg = -deltam*(ct - sigma_2*k/tau +m);
      double exp2arg = -deltam*k*(ct - sigma_2*k/tau +m);
      complex<double> exp2(cos(exp2arg), sin(exp2arg));

      //      complex<double> cerfarg(sigma/(tau*sqrt2) - ct/(sigma*sqrt2), +deltam*sigma/sqrt2);
      //      complex<double> cerfarg(sigma*k/(tau*sqrt2) - ct/(sigma*sqrt2) -m/(sigma*sqrt2), +deltam*sigma/sqrt2);

      complex<double> cerfarg(sigma*k/(tau*sqrt2) - ct/(sigma*sqrt2) -m/(sigma*sqrt2), +deltam*k*sigma/sqrt2);
      complex<double> cerf;
      if  (cerfarg.real() < -20.0)
        cerf = complex<double>(2.0,0.0);
      else
        cerf = cErrF_2(cerfarg);//best complex error function
      complex<double> c2(exp2*cerf);
      double im = -c2.imag();//exp*sin
      double re = +c2.real();//exp*cos

      one = c1*exp1*im ;
      two = c1*exp1*re ;
      //      std::cout<< "k = "<< k <<std::endl;
   }
  else //all three methods are analytically the same. There might be numerical differences, though
    {
      if (ct < -6.0*sigma)
        {
          one = 0.0;
          two = 0.0;
          return;
        }

      //current cdf
      std::complex<double> z(deltam*sigma/sqrt2, (sigma/tau-ct/sigma)/sqrt2);
      if (ct<0) {//i do not quite get this part. the calculation should also be corerct for ct < 0, right? 
        one= 2.0*nwwerf(z).real()/4.0*exp(-ct*ct/2.0/sigma/sigma);
        two= 2.0*nwwerf(z).imag()/4.0*exp(-ct*ct/2.0/sigma/sigma);
      }
      else {
        one= -2.0*nwwerf(std::conj(z)).real()/tau/4*exp(-ct*ct/2.0/sigma/sigma) +
          exp(sigma*sigma/2 *(1/tau/tau - deltam*deltam) - ct/tau)*cos(deltam*ct - deltam/tau*sigma*sigma);
        two= +2.0*nwwerf(std::conj(z)).imag()/tau/4*exp(-ct*ct/2.0/sigma/sigma) +
          exp(sigma*sigma/2 *(1/tau/tau - deltam*deltam) - ct/tau)*sin(deltam*ct - deltam/tau*sigma*sigma);
      }
    }

  if (std::isnan(one) || std::isnan(two))
    {
      //void convoluted_exp_sincos(double ct, double sigma, double tau, double deltam, double& one, double& two)
      cout << "ct= " << ct << " sigma=" << sigma << " tau=" << tau
           << " deltam=" << deltam << " one=" << one << " two=" << two << endl;
      cout << "Error: NaN" << endl;
      assert(0);

    }
  return;
};

double integrate_convoluted_exp_cos_sl_k(double a, double b, double 
						     sigma, double tau, double deltam, double mean, double k)
{
double sigma_2 = sigma*sigma;
double tau_2 = tau*tau;
double m = mean;

complex<double> i(0.0, 1.0);
//std::cout<<setprecision(8)<<i.imag()<<std::endl;
//double prefexp = exp(sigma*(2.0+sigma*(-1.0+tau_2*deltam*deltam))/(2.0*tau_2));
 complex<double> pref(tau/((-1.0-i*tau*deltam)*k));
 complex<double> erfarga ((a+m)/(TMath::Sqrt(2)*sigma), 0.0);
 complex<double> erfcarga ((-tau*(a+m) + sigma_2*k*(1.0 + 
i*tau*deltam))/(TMath::Sqrt(2)*sigma*tau));
complex<double> exparga(-(tau*deltam - i)*k*(sigma_2*k*(tau*deltam - i) + 2.0*i*tau*(a+m))/(2.0*tau_2));

//std::cout<<setprecision(8)<<prefa.imag()<<std::endl;

 complex<double> erfargb  ((b+m)/(TMath::Sqrt(2)*sigma), 0.0);
 complex<double> erfcargb ((-tau*(b+m) + sigma_2*k*(1.0 + 
i*tau*deltam))/(TMath::Sqrt(2)*sigma*tau));
 complex<double> expargb (-(tau*deltam - i)*k*(sigma_2*k*(tau*deltam - i) + 2.0*i*tau*(b+m))/(2.0*tau_2));

//std::cout<<setprecision(8)<<prefb.real()<<std::endl;

complex<double> erfa = ErrF_2(erfarga);
complex<double> erfca = cErrF_2(erfcarga);
complex<double> expa = exp(exparga);

complex<double> erfb = ErrF_2(erfargb);
complex<double> erfcb = cErrF_2(erfcargb);
complex<double> expb = exp(expargb);


complex <double> za = pref * (- erfa + erfca * expa);
complex <double> zb = pref* (- erfb + erfcb * expb);

double two = 0.5*(zb.real() - za.real());

 if (std::isnan(two))
    {
      //void convoluted_exp_sincos(double ct, double sigma, double tau, double deltam, double& one, double& two)
      cout << " sigma=" << sigma << " tau=" << tau
           << " deltam=" << deltam << " two=" << two << endl;
      cout << "Error: NaN" << endl;
      assert(0);

    }
  return two;

};

double integrate_convoluted_exp_cos_sl_k_generic2(double t, double  sigma, double tau, double deltam, double mean, double kfactor, double alpha, double beta, double shift)
{
  const double sqrt2 = TMath::Sqrt(2.0);
  complex<double> i(0.0, 1.0);

  complex<double> k = kfactor*(1.0 + i*deltam*tau);
  complex<double> a (1.0/(sigma*M_SQRT2), 0.0);
  double adouble = a.real();
  complex<double> b = k/tau;
  complex<double> b_2 = b*b;
  complex<double> ctau (tau, 0.0);

  complex<double> m_reso (mean,0.0);
  complex<double> ct (t,0.0);
  complex<double> m_shift (shift, 0.0);

  complex<double> bprime=(k/ctau + 1.0/alpha);
  complex<double> c = (sigma*k/(ctau*M_SQRT2) - m_reso/(sigma*M_SQRT2));

 complex<double> erfcarg (-a*ct + c);

  complex<double> norm = (0.0,0.0);

  complex<double> I0A = 1.0/b*(std::exp(b_2/(4.0*a*a) - b*c/a)*ErrF_2(b/(2.0*a) -c +a*ct) - exp(-b*ct)*cErrF_2(erfcarg));

  complex<double> I0B = -std::exp(m_shift/alpha)*1.0/bprime*(exp(bprime*bprime/(4.0*a*a) - bprime*c/a)*ErrF_2(bprime/(2.0*a) -c +a*ct) - exp(-bprime*ct)*cErrF_2(erfcarg));

  complex<double> factor0 = 0.5 * exp(sigma*sigma*k*k/(2.0*ctau*ctau))* exp(-m_reso*k/ctau);
  complex<double> I1A =  (1.0/(b_2) - 1.0/(2.0*a*a) + c/(a*b))*exp(b_2/(4.0*a*a) - b*c/a)*ErrF_2(a*ct -c +b/(2.0*a)) - 1.0/(a*b*std::sqrt(M_PI))*exp(-(-a*ct+c)*(-a*ct+c) - b*ct)-(b*ct + 1.0)/(b_2)*exp(-b*ct)*cErrF_2(erfcarg);

  complex<double> I1B = - std::exp(m_shift/alpha)*((1.0/bprime/bprime - 1.0/(2.0*a*a) + c/(a*bprime))*exp(bprime*bprime/(4.0*a*a) - bprime*c/a)*ErrF_2(a*ct -c +bprime/(2.0*a))
						   - 1.0/(a*bprime*std::sqrt(M_PI))*std::exp(-(-a*ct+c)*(-a*ct+c) - bprime*ct) 
 						   - (bprime*ct + 1.0)/bprime/bprime*std::exp(-bprime*ct)
 						   *cErrF_2(erfcarg));


  norm = factor0 *((I0A + I0B) + beta * (I1A + I1B));

  double result = norm.real();
  //norm = factor0*I0A;

  return result;

}


double integrate_convoluted_exp_cos_sl_k_ab(double a, double b, double  sigma, double tau, double deltam, double mean, double k, double alpha, double beta, double shift)
{
    double result = integrate_convoluted_exp_cos_sl_k_generic2(b, sigma, tau, deltam, mean, k, alpha, beta, shift) -
      integrate_convoluted_exp_cos_sl_k_generic2(a, sigma, tau, deltam, mean, k, alpha, beta, shift);
  // double result = integrate_convoluted_exp_cos_sl_k(a, b, sigma, tau, deltam, mean, k);
  //double result = integrate_convoluted_exp_cos(a, b, sigma, tau, deltam);

    return result;
}

double integrate_convoluted_exp_cos_sl_k_generic2_TurboAcc(double t, double  sigma, double tau, double deltam, double mean, double kfactor, double alpha, double beta, double gamma, double shift, double alpha2, double frac)
{
  const double sqrt2 = TMath::Sqrt(2.0);
  complex<double> i(0.0, 1.0);

  complex<double> k = kfactor*(1.0 + i*deltam*tau);
  complex<double> a (1.0/(sigma*M_SQRT2), 0.0);
  double adouble = a.real();
  complex<double> a_4 = adouble*adouble*adouble*adouble;

  complex<double> b = k/tau;
  complex<double> b_2 = b*b;
  complex<double> b_4 = b*b*b*b;

  complex<double> ctau (tau, 0.0);

  complex<double> m_reso (mean,0.0);
  complex<double> ct (t,0.0);
  complex<double> m_shift (shift, 0.0);

  complex<double> bprime=(k/ctau + 1.0/alpha);
  complex<double> bprime_2 = bprime*bprime;
  complex<double> bprime_4 = bprime*bprime*bprime*bprime;


  complex<double> bprime2=(k/ctau + 1.0/alpha2);

  complex<double> c = (sigma*k/(ctau*M_SQRT2) - m_reso/(sigma*M_SQRT2));

 complex<double> erfcarg (-a*ct + c);

  complex<double> norm = (0.0,0.0);

  complex<double> I0A = 1.0/b*(std::exp(b_2/(4.0*a*a) - b*c/a)*ErrF_2(b/(2.0*a) -c +a*ct) - exp(-b*ct)*cErrF_2(erfcarg));

  complex<double> I0B = -std::exp(m_shift/alpha)*1.0/bprime*(exp(bprime*bprime/(4.0*a*a) - bprime*c/a)*ErrF_2(bprime/(2.0*a) -c +a*ct) - exp(-bprime*ct)*cErrF_2(erfcarg));

  complex<double> I0B2 = -std::exp(m_shift/alpha2)*1.0/bprime2*(exp(bprime2*bprime2/(4.0*a*a) - bprime2*c/a)*ErrF_2(bprime2/(2.0*a) -c +a*ct) - exp(-bprime2*ct)*cErrF_2(erfcarg));

  complex<double> factor0 = 0.5 * exp(sigma*sigma*k*k/(2.0*ctau*ctau))* exp(-m_reso*k/ctau);
  complex<double> I1A =  (1.0/(b_2) - 1.0/(2.0*a*a) + c/(a*b))*exp(b_2/(4.0*a*a) - b*c/a)*ErrF_2(a*ct -c +b/(2.0*a)) - 1.0/(a*b*std::sqrt(M_PI))*exp(-(-a*ct+c)*(-a*ct+c) - b*ct)-(b*ct + 1.0)/(b_2)*exp(-b*ct)*cErrF_2(erfcarg);

  complex<double> I1B = - std::exp(m_shift/alpha)*((1.0/bprime/bprime - 1.0/(2.0*a*a) + c/(a*bprime))*exp(bprime*bprime/(4.0*a*a) - bprime*c/a)*ErrF_2(a*ct -c +bprime/(2.0*a))
						   - 1.0/(a*bprime*std::sqrt(M_PI))*std::exp(-(-a*ct+c)*(-a*ct+c) - bprime*ct) 
 						   - (bprime*ct + 1.0)/bprime/bprime*std::exp(-bprime*ct)
 						   *cErrF_2(erfcarg));

  complex<double> I1B2 = - std::exp(m_shift/alpha2)*((1.0/bprime2/bprime2 - 1.0/(2.0*a*a) + c/(a*bprime2))*exp(bprime2*bprime2/(4.0*a*a) - bprime2*c/a)*ErrF_2(a*ct -c +bprime2/(2.0*a))
						   - 1.0/(a*bprime2*std::sqrt(M_PI))*std::exp(-(-a*ct+c)*(-a*ct+c) - bprime2*ct) 
 						   - (bprime2*ct + 1.0)/bprime2/bprime2*std::exp(-bprime2*ct)
 						   *cErrF_2(erfcarg));

  complex<double> I2A =1.0/(4.0*a_4*b*b*b)*std::exp(-b*(c/a + ct))*((8.0*a_4+b_4+8.0*a*a*a*b*c-4.0*a*b*b*b*c+2.0*a*a*b_2*(-1.0+2.0*c*c))*std::exp(b_2/(4.0*a*a)+b*ct)*ErrF_2(b/(2.0*a)-c+a*ct)+1.0/std::sqrt(M_PI)*2.0*a*std::exp(b*c/a)*(b*std::exp(-((c - a*ct)*(c - a*ct)))*(b_2 - 2.0*a*(2.0*a+b*c+a*b*ct))-2.0*a*a*a*std::sqrt(M_PI)*(2.0+b*ct*(2.0+b*ct))*cErrF_2(erfcarg)));

  complex<double> I2B = - std::exp(m_shift/alpha)*1.0/(4.0*a_4*bprime*bprime*bprime)*std::exp(-bprime*(c/a + ct))*((8.0*a_4+bprime_4+8.0*a*a*a*bprime*c-4.0*a*bprime*bprime*bprime*c+2.0*a*a*bprime_2*(-1.0+2.0*c*c))*std::exp(bprime_2/(4.0*a*a)+bprime*ct)*ErrF_2(bprime/(2.0*a)-c+a*ct)+1.0/std::sqrt(M_PI)*2*a*std::exp(bprime*c/a)*(bprime*std::exp(-((c - a*ct)*(c - a*ct)))*(bprime_2 - 2.0*a*(2.0*a+bprime*c+a*bprime*ct))-2.0*a*a*a*std::sqrt(M_PI)*(2.0+bprime*ct*(2.0+bprime*ct))*cErrF_2(erfcarg)));


  norm = factor0 *((I0A + frac*I0B + (1.0-frac)*I0B2) + beta * (I1A + frac*I1B+(1.0-frac)*I1B2) + gamma * (I2A + I2B));


  //  norm = factor0 *((I0A + I0B) + beta * (I1A + I1B));
  norm = factor0 *((I0A + frac*I0B + (1.0-frac)*I0B2) + beta * (I1A + frac*I1B+(1.0-frac)*I1B2));

  double result = norm.real();
  //norm = factor0*I0A;

  return result;

}


double integrate_convoluted_exp_cos_sl_k_ab_TurboAcc(double a, double b, double  sigma, double tau, double deltam, double mean, double k, double alpha, double beta, double gamma, double shift, double alpha2, double frac)
{
  double result = integrate_convoluted_exp_cos_sl_k_generic2_TurboAcc(b, sigma, tau, deltam, mean, k, alpha, beta, gamma, shift, alpha2, frac) -
    integrate_convoluted_exp_cos_sl_k_generic2_TurboAcc(a, sigma, tau, deltam, mean, k, alpha, beta, gamma, shift, alpha2, frac);
  // double result = integrate_convoluted_exp_cos_sl_k(a, b, sigma, tau, deltam, mean, k);
  //double result = integrate_convoluted_exp_cos(a, b, sigma, tau, deltam);

    return result;
}


double integrate_convoluted_exp_sin_sl_k_generic2_TurboAcc(double t, double  sigma, double tau, double deltam, double mean, double kfactor, double alpha, double beta, double gamma, double shift, double alpha2, double frac)
{
  const double sqrt2 = TMath::Sqrt(2.0);
  complex<double> i(0.0, 1.0);

  complex<double> k = kfactor*(1.0 + i*deltam*tau);
  complex<double> a (1.0/(sigma*M_SQRT2), 0.0);
  double adouble = a.real();
  complex<double> a_4 = adouble*adouble*adouble*adouble;

  complex<double> b = k/tau;
  complex<double> b_2 = b*b;
  complex<double> b_4 = b*b*b*b;

  complex<double> ctau (tau, 0.0);

  complex<double> m_reso (mean,0.0);
  complex<double> ct (t,0.0);
  complex<double> m_shift (shift, 0.0);

  complex<double> bprime=(k/ctau + 1.0/alpha);
  complex<double> bprime_2 = bprime*bprime;
  complex<double> bprime_4 = bprime*bprime*bprime*bprime;


  complex<double> bprime2=(k/ctau + 1.0/alpha2);

  complex<double> c = (sigma*k/(ctau*M_SQRT2) - m_reso/(sigma*M_SQRT2));

 complex<double> erfcarg (-a*ct + c);

  complex<double> norm = (0.0,0.0);

  complex<double> I0A = 1.0/b*(std::exp(b_2/(4.0*a*a) - b*c/a)*ErrF_2(b/(2.0*a) -c +a*ct) - exp(-b*ct)*cErrF_2(erfcarg));

  complex<double> I0B = -std::exp(m_shift/alpha)*1.0/bprime*(exp(bprime*bprime/(4.0*a*a) - bprime*c/a)*ErrF_2(bprime/(2.0*a) -c +a*ct) - exp(-bprime*ct)*cErrF_2(erfcarg));

  complex<double> I0B2 = -std::exp(m_shift/alpha2)*1.0/bprime2*(exp(bprime2*bprime2/(4.0*a*a) - bprime2*c/a)*ErrF_2(bprime2/(2.0*a) -c +a*ct) - exp(-bprime2*ct)*cErrF_2(erfcarg));

  complex<double> factor0 = 0.5 * exp(sigma*sigma*k*k/(2.0*ctau*ctau))* exp(-m_reso*k/ctau);
  complex<double> I1A =  (1.0/(b_2) - 1.0/(2.0*a*a) + c/(a*b))*exp(b_2/(4.0*a*a) - b*c/a)*ErrF_2(a*ct -c +b/(2.0*a)) - 1.0/(a*b*std::sqrt(M_PI))*exp(-(-a*ct+c)*(-a*ct+c) - b*ct)-(b*ct + 1.0)/(b_2)*exp(-b*ct)*cErrF_2(erfcarg);

  complex<double> I1B = - std::exp(m_shift/alpha)*((1.0/bprime/bprime - 1.0/(2.0*a*a) + c/(a*bprime))*exp(bprime*bprime/(4.0*a*a) - bprime*c/a)*ErrF_2(a*ct -c +bprime/(2.0*a))
						   - 1.0/(a*bprime*std::sqrt(M_PI))*std::exp(-(-a*ct+c)*(-a*ct+c) - bprime*ct) 
 						   - (bprime*ct + 1.0)/bprime/bprime*std::exp(-bprime*ct)
 						   *cErrF_2(erfcarg));

  complex<double> I1B2 = - std::exp(m_shift/alpha2)*((1.0/bprime2/bprime2 - 1.0/(2.0*a*a) + c/(a*bprime2))*exp(bprime2*bprime2/(4.0*a*a) - bprime2*c/a)*ErrF_2(a*ct -c +bprime2/(2.0*a))
						   - 1.0/(a*bprime2*std::sqrt(M_PI))*std::exp(-(-a*ct+c)*(-a*ct+c) - bprime2*ct) 
 						   - (bprime2*ct + 1.0)/bprime2/bprime2*std::exp(-bprime2*ct)
 						   *cErrF_2(erfcarg));

  complex<double> I2A =1.0/(4.0*a_4*b*b*b)*std::exp(-b*(c/a + ct))*((8.0*a_4+b_4+8.0*a*a*a*b*c-4.0*a*b*b*b*c+2.0*a*a*b_2*(-1.0+2.0*c*c))*std::exp(b_2/(4.0*a*a)+b*ct)*ErrF_2(b/(2.0*a)-c+a*ct)+1.0/std::sqrt(M_PI)*2.0*a*std::exp(b*c/a)*(b*std::exp(-((c - a*ct)*(c - a*ct)))*(b_2 - 2.0*a*(2.0*a+b*c+a*b*ct))-2.0*a*a*a*std::sqrt(M_PI)*(2.0+b*ct*(2.0+b*ct))*cErrF_2(erfcarg)));

  complex<double> I2B = - std::exp(m_shift/alpha)*1.0/(4.0*a_4*bprime*bprime*bprime)*std::exp(-bprime*(c/a + ct))*((8.0*a_4+bprime_4+8.0*a*a*a*bprime*c-4.0*a*bprime*bprime*bprime*c+2.0*a*a*bprime_2*(-1.0+2.0*c*c))*std::exp(bprime_2/(4.0*a*a)+bprime*ct)*ErrF_2(bprime/(2.0*a)-c+a*ct)+1.0/std::sqrt(M_PI)*2*a*std::exp(bprime*c/a)*(bprime*std::exp(-((c - a*ct)*(c - a*ct)))*(bprime_2 - 2.0*a*(2.0*a+bprime*c+a*bprime*ct))-2.0*a*a*a*std::sqrt(M_PI)*(2.0+bprime*ct*(2.0+bprime*ct))*cErrF_2(erfcarg)));


  norm = factor0 *((I0A + frac*I0B + (1.0-frac)*I0B2) + beta * (I1A + frac*I1B+(1.0-frac)*I1B2) + gamma * (I2A + I2B));


  //  norm = factor0 *((I0A + I0B) + beta * (I1A + I1B));
  norm = factor0 *((I0A + frac*I0B + (1.0-frac)*I0B2) + beta * (I1A + frac*I1B+(1.0-frac)*I1B2));

  double result = norm.imag();
  //norm = factor0*I0A;

  return result;

}


double integrate_convoluted_exp_sin_sl_k_ab_TurboAcc(double a, double b, double  sigma, double tau, double deltam, double mean, double k, double alpha, double beta, double gamma, double shift, double alpha2, double frac)
{
  double result = integrate_convoluted_exp_sin_sl_k_generic2_TurboAcc(b, sigma, tau, deltam, mean, k, alpha, beta, gamma, shift, alpha2, frac) -
    integrate_convoluted_exp_sin_sl_k_generic2_TurboAcc(a, sigma, tau, deltam, mean, k, alpha, beta, gamma, shift, alpha2, frac);
  // double result = integrate_convoluted_exp_cos_sl_k(a, b, sigma, tau, deltam, mean, k);
  //double result = integrate_convoluted_exp_cos(a, b, sigma, tau, deltam);

    return result;
}

double integrate_convoluted_exp_cos_sl_k_NoAcc(double t, double  sigma, double tau, double deltam, double mean, double kfactor)
{
  const double sqrt2 = TMath::Sqrt(2.0);
  complex<double> i(0.0, 1.0);

  complex<double> k = kfactor*(1.0 + i*deltam*tau);
  complex<double> a (1.0/(sigma*M_SQRT2), 0.0);
  double adouble = a.real();
  complex<double> b = k/tau;
  complex<double> b_2 = b*b;
  complex<double> ctau (tau, 0.0);

  complex<double> m_reso (mean,0.0);
  complex<double> ct (t,0.0);

  complex<double> c = (sigma*k/(ctau*M_SQRT2) - m_reso/(sigma*M_SQRT2));

 complex<double> erfcarg (-a*ct + c);

  complex<double> norm = (0.0,0.0);

  complex<double> I0A = 1.0/b*(std::exp(b_2/(4.0*a*a) - b*c/a)*ErrF_2(b/(2.0*a) -c +a*ct) - exp(-b*ct)*cErrF_2(erfcarg));


  complex<double> factor0 = 0.5 * exp(sigma*sigma*k*k/(2.0*ctau*ctau))* exp(-m_reso*k/ctau);


  //norm = factor0 *((I0A + I0B) + beta * (I1A + I1B));
  norm = factor0*I0A;

  double result = norm.real();

  return result;

}

double integrate_convoluted_exp_cos_sl_k_ab_NoAcc(double a, double b, double  sigma, double tau, double deltam, double mean, double k)
{
    double result = integrate_convoluted_exp_cos_sl_k_NoAcc(b, sigma, tau, deltam, mean, k) -
      integrate_convoluted_exp_cos_sl_k_NoAcc(a, sigma, tau, deltam, mean, k);

    return result;
}


double integrate_convoluted_exp_sl_k_generic(double t, double  sigma, double tau, double deltam, double mean, double k, double alpha, double beta, double gamma, double shift)
{
  const double sqrt2 = TMath::Sqrt(2.0);

  double a = 1.0/(sigma*M_SQRT2); 
  double a_4 = a*a*a*a;
  double b = k/tau;
  double b_2 = b*b;
  double b_4 = b*b*b*b;

  double ctau = tau;

  double m_reso = mean;
  double ct = t;
  double m_shift = shift;

  double bprime=(k/ctau + 1.0/alpha);
  double bprime_2 = bprime*bprime;
  double bprime_4 = bprime*bprime*bprime*bprime;

  double c = (sigma*k/(ctau*M_SQRT2) - m_reso/(sigma*M_SQRT2));

  double norm = 0.0;

  double I0A = 1.0/b*(std::exp(b_2/(4*a*a) - b*c/a)*TMath::Erf(b/(2*a) -c +a*ct) - std::exp(-b*ct)*TMath::Erfc(-a*ct + c));

  double I0B = -std::exp(m_shift/alpha)*1.0/bprime*(std::exp(bprime*bprime/(4*a*a) - bprime*c/a)*TMath::Erf(bprime/(2*a) -c +a*ct) - std::exp(-bprime*ct)*TMath::Erfc(-a*ct + c));

  double factor0 = 0.5 * std::exp(sigma*sigma*k*k/(2*ctau*ctau))* std::exp(-m_reso*k/ctau);
  double I1A =  (1.0/(b_2) - 1.0/(2.0*a*a) + c/(a*b))*std::exp(b_2/(4*a*a) - b*c/a)*TMath::Erf(a*ct -c +b/(2*a)) - 1.0/(a*b*std::sqrt(M_PI))*std::exp(-(-a*ct+c)*(-a*ct+c) - b*ct)-(b*ct + 1.0)/(b_2)*std::exp(-b*ct)*TMath::Erfc(-a*ct + c);

  double I1B = - std::exp(m_shift/alpha)*((1.0/bprime/bprime - 1.0/(2.0*a*a) + c/(a*bprime))*std::exp(bprime*bprime/(4*a*a) - bprime*c/a)*TMath::Erf(a*ct -c +bprime/(2*a)) - 1.0/(a*bprime*std::sqrt(M_PI))*std::exp(-(-a*ct+c)*(-a*ct+c) - bprime*ct) - (bprime*ct + 1.0)/bprime/bprime*std::exp(-bprime*ct)*TMath::Erfc(-a*ct + c));

  double I2A =1.0/(4.0*a_4*b*b*b)*std::exp(-b*(c/a + ct))*((8.0*a_4+b_4+8.0*a*a*a*b*c-4.0*a*b*b*b*c+2.0*a*a*b_2*(-1.0+2.0*c*c))*std::exp(b_2/(4*a*a)+b*ct)*TMath::Erf(b/(2.0*a)-c+a*ct)+1.0/std::sqrt(M_PI)*2.0*a*std::exp(b*c/a)*(b*std::exp(-((c - a*ct)*(c - a*ct)))*(b_2 - 2.0*a*(2.0*a+b*c+a*b*ct))-2.0*a*a*a*std::sqrt(M_PI)*(2.0+b*ct*(2.0+b*ct))*TMath::Erfc(c - a*ct)));

  double I2B = - std::exp(m_shift/alpha)*1.0/(4.0*a_4*bprime*bprime*bprime)*std::exp(-bprime*(c/a + ct))*((8.0*a_4+bprime_4+8.0*a*a*a*bprime*c-4.0*a*bprime*bprime*bprime*c+2.0*a*a*bprime_2*(-1.0+2.0*c*c))*std::exp(bprime_2/(4.0*a*a)+bprime*ct)*TMath::Erf(bprime/(2.0*a)-c+a*ct)+1.0/std::sqrt(M_PI)*2*a*std::exp(bprime*c/a)*(bprime*std::exp(-((c - a*ct)*(c - a*ct)))*(bprime_2 - 2.0*a*(2.0*a+bprime*c+a*bprime*ct))-2.0*a*a*a*std::sqrt(M_PI)*(2.0+bprime*ct*(2.0+bprime*ct))*TMath::Erfc(c - a*ct)));

  norm = factor0 *((I0A + I0B) + beta * (I1A + I1B) + gamma * (I2A + I2B));

  // norm = factor0 *((I0A + frac*I0B + (1.0-frac)*I0B2) + beta * (I1A + frac*I1B+(1.0-frac)*I1B2) + gamma * (I2A + I2B));

  // norm = factor0 *((I0A + I0B) + beta * (I1A + I1B));

  //norm = factor0*I0A;

  return norm;

}

double integrate_convoluted_exp_sl_k_ab(double a, double b, double  sigma, double tau, double deltam, double mean, double k, double alpha, double beta, double gamma, double shift)
{
  double result = 0; 
 
  result = integrate_convoluted_exp_sl_k_generic(b, sigma, tau, deltam, mean, k, alpha, beta, gamma,shift) - integrate_convoluted_exp_sl_k_generic(a, sigma, tau, deltam, mean, k, alpha, beta, gamma,shift);

    return result;
}

double integrate_convoluted_exp_sl_k_generic_TurboAcc(double t, double  sigma, double tau, double deltam, double mean, double k, double alpha, double beta, double gamma, double shift, double alpha2, double frac)
{
  const double sqrt2 = TMath::Sqrt(2.0);

  double a = 1.0/(sigma*M_SQRT2);
  double a_4 = a*a*a*a;
  double b = k/tau;
  double b_2 = b*b; 
  double b_4 = b*b*b*b;
  double ctau = tau;

  double m_reso = mean;
  double ct = t;
  double m_shift = shift;

  double bprime=(k/ctau + 1.0/alpha);
  double bprime_2 = bprime*bprime;
  double bprime_4 = bprime*bprime*bprime*bprime;

  double bprime2=(k/ctau + 1.0/alpha2);

  double c = (sigma*k/(ctau*M_SQRT2) - m_reso/(sigma*M_SQRT2));

  double norm = 0.0;

  double I0A = 1.0/b*(std::exp(b_2/(4*a*a) - b*c/a)*TMath::Erf(b/(2*a) -c +a*ct) - std::exp(-b*ct)*TMath::Erfc(-a*ct + c));

  double I0B = -std::exp(m_shift/alpha)*1.0/bprime*(std::exp(bprime*bprime/(4*a*a) - bprime*c/a)*TMath::Erf(bprime/(2*a) -c +a*ct) - std::exp(-bprime*ct)*TMath::Erfc(-a*ct + c));

  double I0B2 = -std::exp(m_shift/alpha2)*1.0/bprime2*(std::exp(bprime2*bprime2/(4*a*a) - bprime2*c/a)*TMath::Erf(bprime2/(2*a) -c +a*ct) - std::exp(-bprime2*ct)*TMath::Erfc(-a*ct + c));

  double factor0 = 0.5 * std::exp(sigma*sigma*k*k/(2*ctau*ctau))* std::exp(-m_reso*k/ctau);
  double I1A =  (1.0/(b_2) - 1.0/(2.0*a*a) + c/(a*b))*std::exp(b_2/(4*a*a) - b*c/a)*TMath::Erf(a*ct -c +b/(2*a)) - 1.0/(a*b*std::sqrt(M_PI))*std::exp(-(-a*ct+c)*(-a*ct+c) - b*ct)-(b*ct + 1.0)/(b_2)*std::exp(-b*ct)*TMath::Erfc(-a*ct + c);

  double I1B = - std::exp(m_shift/alpha)*((1.0/bprime/bprime - 1.0/(2.0*a*a) + c/(a*bprime))*std::exp(bprime*bprime/(4*a*a) - bprime*c/a)*TMath::Erf(a*ct -c +bprime/(2*a)) - 1.0/(a*bprime*std::sqrt(M_PI))*std::exp(-(-a*ct+c)*(-a*ct+c) - bprime*ct) - (bprime*ct + 1.0)/bprime/bprime*std::exp(-bprime*ct)*TMath::Erfc(-a*ct + c));

  double I1B2 = - std::exp(m_shift/alpha2)*((1.0/bprime2/bprime2 - 1.0/(2.0*a*a) + c/(a*bprime2))*std::exp(bprime2*bprime2/(4*a*a) - bprime2*c/a)*TMath::Erf(a*ct -c +bprime2/(2*a)) - 1.0/(a*bprime2*std::sqrt(M_PI))*std::exp(-(-a*ct+c)*(-a*ct+c) - bprime2*ct) - (bprime2*ct + 1.0)/bprime2/bprime2*std::exp(-bprime2*ct)*TMath::Erfc(-a*ct + c));

  double I2A =1.0/(4.0*a_4*b*b*b)*std::exp(-b*(c/a + ct))*((8.0*a_4+b_4+8.0*a*a*a*b*c-4.0*a*b*b*b*c+2.0*a*a*b_2*(-1.0+2.0*c*c))*std::exp(b_2/(4*a*a)+b*ct)*TMath::Erf(b/(2.0*a)-c+a*ct)+1.0/std::sqrt(M_PI)*2.0*a*std::exp(b*c/a)*(b*std::exp(-((c - a*ct)*(c - a*ct)))*(b_2 - 2.0*a*(2.0*a+b*c+a*b*ct))-2.0*a*a*a*std::sqrt(M_PI)*(2.0+b*ct*(2.0+b*ct))*TMath::Erfc(c - a*ct)));

  double I2B = - std::exp(m_shift/alpha)*1.0/(4.0*a_4*bprime*bprime*bprime)*std::exp(-bprime*(c/a + ct))*((8.0*a_4+bprime_4+8.0*a*a*a*bprime*c-4.0*a*bprime*bprime*bprime*c+2.0*a*a*bprime_2*(-1.0+2.0*c*c))*std::exp(bprime_2/(4.0*a*a)+bprime*ct)*TMath::Erf(bprime/(2.0*a)-c+a*ct)+1.0/std::sqrt(M_PI)*2*a*std::exp(bprime*c/a)*(bprime*std::exp(-((c - a*ct)*(c - a*ct)))*(bprime_2 - 2.0*a*(2.0*a+bprime*c+a*bprime*ct))-2.0*a*a*a*std::sqrt(M_PI)*(2.0+bprime*ct*(2.0+bprime*ct))*TMath::Erfc(c - a*ct)));


  norm = factor0 *((I0A + frac*I0B + (1.0-frac)*I0B2) + beta * (I1A + frac*I1B+(1.0-frac)*I1B2) + gamma * (I2A + I2B));

  //norm = factor0 *((I0A + I0B) + beta * (I1A + I1B) + gamma * (I2A + I2B));

  // norm = factor0 *((I0A + frac*I0B + (1.0-frac)*I0B2) + beta * (I1A + frac*I1B+(1.0-frac)*I1B2));

  //norm = factor0 *((I0A + I0B) + beta * (I1A + I1B));

  //norm = factor0*I0A;

  return norm;

}

double integrate_convoluted_exp_sl_k_ab_TurboAcc(double a, double b, double  sigma, double tau, double deltam, double mean, double k, double alpha, double beta, double gamma,  double shift, double alpha2, double frac)
{
  double result = 0; 
 
  result = integrate_convoluted_exp_sl_k_generic_TurboAcc(b, sigma, tau, deltam, mean, k, alpha, beta, gamma, shift, alpha2, frac) - integrate_convoluted_exp_sl_k_generic_TurboAcc(a, sigma, tau, deltam, mean, k, alpha, beta, gamma, shift, alpha2, frac);

    return result;
}

double integrate_convoluted_exp_sl_k_NoAcc(double t, double  sigma, double tau, double deltam, double mean, double k)
{
  const double sqrt2 = TMath::Sqrt(2.0);

  double a = 1.0/(sigma*M_SQRT2);
  double b = k/tau;
  double b_2 = b*b;
  double ctau = tau;

  double m_reso = mean;
  double ct = t;
  double c = (sigma*k/(ctau*M_SQRT2) - m_reso/(sigma*M_SQRT2));

  double norm = 0.0;

  double I0A = 1.0/b*(std::exp(b_2/(4*a*a) - b*c/a)*TMath::Erf(b/(2*a) -c +a*ct) - std::exp(-b*ct)*TMath::Erfc(-a*ct + c));

  double factor0 = 0.5 * std::exp(sigma*sigma*k*k/(2*ctau*ctau))* std::exp(-m_reso*k/ctau);

  //norm = factor0 *((I0A + I0B) + beta * (I1A + I1B));

  norm = factor0*I0A;

  return norm;

}

double integrate_convoluted_exp_sl_k_ab_NoAcc(double a, double b, double  sigma, double tau, double deltam, double mean, double k)
{
  double result = 0; 
 
  result = integrate_convoluted_exp_sl_k_NoAcc(b, sigma, tau, deltam, mean, k) - integrate_convoluted_exp_sl_k_NoAcc(a, sigma, tau, deltam, mean, k);

    return result;
}

double integrate_convoluted_t_exp_cos_sl_k(double a, double b, double sigma, double tau, double deltam, double mean, double k)
{
  double sigma_2 = sigma*sigma;
  double tau_2 = tau*tau;
  double m = mean;
  double k_2 = k*k;

  complex<double> i(0.0, 1.0);
  //std::cout<<setprecision(8)<<i.imag()<<std::endl;
  //double prefexp = exp(sigma*(2.0+sigma*(-1.0+tau_2*deltam*deltam))/(2.0*tau_2));
  complex<double> pref(tau/((-1.0+i*tau*deltam)*k_2));
  
  complex<double> erfarga ((a+m)/(TMath::Sqrt(2)*sigma), 0.0);
  complex<double> erfcarga ((-tau*(a+m) + sigma_2*k*(1.0 + i*tau*deltam))/(TMath::Sqrt(2)*sigma*tau));
  complex<double> exparga(0.5*(1.0/tau - i*deltam)*k*((1.0/tau -i*deltam)*k*sigma_2 - 2.0*(a+m)));

  //std::cout<<setprecision(8)<<prefa.imag()<<std::endl;
  
  complex<double> erfargb  ((b+m)/(TMath::Sqrt(2)*sigma), 0.0);
  complex<double> erfcargb ((-tau*(b+m) + sigma_2*k*(1.0 + i*tau*deltam))/(TMath::Sqrt(2)*sigma*tau));
  complex<double> expargb(0.5*(1.0/tau - i*deltam)*k*((1.0/tau -i*deltam)*k*sigma_2 - 2.0*(b+m)));
  
  //std::cout<<setprecision(8)<<prefb.real()<<std::endl;
  
  complex<double> erfa = ErrF_2(erfarga);
  complex<double> erfca = cErrF_2(erfcarga);
  complex<double> expa = exp(exparga);
  
  complex<double> erfb = ErrF_2(erfargb);
  complex<double> erfcb = cErrF_2(erfcargb);
  complex<double> expb = exp(expargb);
  
  
  complex <double> za = pref * (- erfa*(tau/(1.0-i*deltam*tau)-k*a) + erfca * expa *(tau/(1.0-i*deltam*tau)+k*a) 
				+TMath::Sqrt(2.0/TMath::Pi())*sigma*k*exp(-(m+a)*(m+a)/(2.0*sigma_2)) );
  complex <double> zb = pref * (- erfb*(tau/(1.0-i*deltam*tau)-k*b) + erfcb * expb *(tau/(1.0-i*deltam*tau)+k*b) 
				+TMath::Sqrt(2.0/TMath::Pi())*sigma*k*exp(-(m+b)*(m+b)/(2.0*sigma_2)) );
  
  double two = 0.5*(zb.real() - za.real());

 if (std::isnan(two))
    {
      //void convoluted_exp_sincos(double ct, double sigma, double tau, double deltam, double& one, double& two)
      cout << " sigma=" << sigma << " tau=" << tau
           << " deltam=" << deltam << " two=" << two << endl;
      cout << "Error: NaN" << endl;
      assert(0);

    }
  return two;

};


double johnson(double x, double mu, double sigma, double gamma, double delta, double x_min, double x_max)
{
  if(x<x_min || x>x_max)
    return 0.;
  if(sigma==0. || delta==0.)
    return 0.;
  double z       = (x-mu)/sigma;
  double arg     = gamma+delta*TMath::ASinH(z);
  double arg_max = gamma+delta*TMath::ASinH((x_max-mu)/sigma);
  double arg_min = gamma+delta*TMath::ASinH((x_min-mu)/sigma);
  double norm    = TMath::Sqrt(TMath::Pi()/2.)*sigma/delta*(TMath::Erf(arg_max/TMath::Sqrt2())-TMath::Erf(arg_min/TMath::Sqrt2()));
  if(norm==0)
    return 0.;
  return TMath::Exp(-arg*arg/2.)/(norm*TMath::Sqrt(1.+z*z));
}

double DstD0BG(double x, double dm0, double C, double A, double B, double x_min, double x_max)
{

  double norm = 0.0;
  
  double arg= x-dm0;
  if (arg <= 0 ) return 0;
  double ratio= x/dm0;
  double val = (1- TMath::Exp(-arg/C))* TMath::Power(ratio, A) + B*(ratio-1);
  
  if (x_max <= dm0 ) return 0;
  else if (x_min < dm0) x_min = dm0;
  
  
  if( (A!=0) || ((B<0) && ( 1.0- TMath::Exp(-(x_max-dm0)/C) + B*(x_max/dm0 -1) < 0)))
    std::cout<< "numerical integration needed for DstD0BG: first check values of parameters!"<<std::endl;
  else
    norm = (x_max-x_min)+ C * TMath::Exp(dm0/C)* (TMath::Exp(-x_max/C)- TMath::Exp(-x_min/C)) + B * (0.5 * (x_max*x_max - x_min*x_min)/dm0 - (x_max- x_min));
  
  return val/norm;
  
}

double CB_function(double m, double mean, double alpha, double n, double sigma)
{
  double t = (m-mean)/sigma;
  if (alpha < 0) t = -t;

  double absAlpha = fabs(alpha);

  if (t >= -absAlpha) {
    return std::exp(-0.5*t*t);
  }
  else {
    double a = TMath::Power(n/absAlpha,n)*std::exp(-0.5*absAlpha*absAlpha);
    double b= n/absAlpha - absAlpha;

    return a/TMath::Power(b-t,n);
  }
}

double norm_CB_function(double mean, double alpha, double n, double sigma, double m_min, double m_max)
{
  static const double sqrtPiOver2 = 1.2533141373;
  static const double sqrt2 = 1.4142135624;

  double result = 0.0;
  bool useLog = false;

  if( fabs(n-1.0) < 1.0e-05 )
    useLog = true;
  
  double sig = fabs(sigma);
  
  double tmin = (m_min-mean)/sig;
  double tmax = (m_max-mean)/sig;
  
  if(alpha < 0) {
    double tmp = tmin;
    tmin = -tmax;
    tmax = -tmp;
  }

  double absAlpha = fabs(alpha);
  
  if( tmin >= -absAlpha ) {
    result += sig*sqrtPiOver2*(   TMath::Erf(tmax/sqrt2)
                                - TMath::Erf(tmin/sqrt2) );
  }
  else if( tmax <= -absAlpha ) {
    double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    double b = n/absAlpha - absAlpha;
    
    if(useLog) {
      result += a*sig*( std::log(b-tmin) - std::log(b-tmax) );
    }
    else {
      result += a*sig/(1.0-n)*(   1.0/(TMath::Power(b-tmin,n-1.0))
                                - 1.0/(TMath::Power(b-tmax,n-1.0)) );
    }
  }
  else {
    double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    double b = n/absAlpha - absAlpha;
    
    double term1 = 0.0;
    if(useLog) {
      term1 = a*sig*(  log(b-tmin) - log(n/absAlpha));
    }
    else {
      term1 = a*sig/(1.0-n)*(   1.0/(TMath::Power(b-tmin,n-1.0))
                              - 1.0/(TMath::Power(n/absAlpha,n-1.0)) );
    }
    
    double term2 = sig*sqrtPiOver2*(   TMath::Erf(tmax/sqrt2)
                                     - TMath::Erf(-absAlpha/sqrt2) );
    
    
    result += term1 + term2;
  }
  
  return result;
}



double norm_convoluted_exp(double ct, double sigma_ct, double ctau)
{
  return ctau;
}

//may be needed later...
double expdecay(double t, double tau, bool normalized)
{
  if (t >= 0.0)
    {
      if (normalized)
	return 1.0/tau*exp(-t/tau);
      else
	return exp(-t/tau);
    }
  else 
    return 0.0;
}

double gauss(double x, double sigma, double mean)
{
  //  what do we have libraries for...
  //this one is normalised
  return TMath::Gaus(x, mean, sigma, kTRUE);
}

double norm_gauss(double sigma, double mean, double min, double max)
{
  //corrected error
  double norm = 1.0/2.0
    *(TMath::Erf((max-mean)/(TMath::Sqrt(2.0)*sigma))
      -TMath::Erf((min-mean)/(TMath::Sqrt(2.0)*sigma)));
  return norm;
}

double NOTnormed_gauss(double x, double sigma, double mean, double min, double max)
{
  //corrected error
  return 1.0/(TMath::Sqrt(2.0*TMath::Pi())*sigma)*exp(-(mean-x)*(mean-x)/(2*sigma*sigma));
}

double normed_gauss(double x, double sigma, double mean, double min, double max)
{
  //corrected error
  double norm = 1.0/2.0
    *(TMath::Erf((max-mean)/(TMath::Sqrt(2.0)*sigma))
      -TMath::Erf((min-mean)/(TMath::Sqrt(2.0)*sigma)));
  return 1.0/(TMath::Sqrt(2.0*TMath::Pi())*sigma)*exp(-(mean-x)*(mean-x)/(2*sigma*sigma))/norm;
}

double integrate_normed_gauss(double sigma, double mean, double min, double max, double a, double b)
{
  //corrected error                                                                                                                                                                                                                                                            
  double norm = 1.0/2.0
    *(TMath::Erf((max-mean)/(TMath::Sqrt(2.0)*sigma))
      -TMath::Erf((min-mean)/(TMath::Sqrt(2.0)*sigma)));
  return  1.0/2.0*(TMath::Erf((b-mean)/(TMath::Sqrt(2.0)*sigma))-TMath::Erf((a-mean)/(TMath::Sqrt(2.0)*sigma)))/norm;
}


double linear(double x, double min, double max, double slope)
{
  //have corrected the slope (parameter needs to be big enough for minuit -> divide by width = 500^2)
  double width = max - min;
  double ymin = 1.0/width - 0.5*slope/(width*width)*width;//corrected the 0.5 factor
  return ymin + slope/(width*width) * (x - min);
}

/** This function tests the different versions of the complex error function.
 * Unfortunately the complex error function is not yet in the standard c library.
 **/
void test_errf()
{ 
  //gStyle->SetPalette(1);
  unsigned int divisions=1000;
  double min=-10.0;
  double max=+10.0;
  TH2D* real_1 = new TH2D("real_1", "real part", divisions, min, max, divisions, min, max);
  real_1->SetMinimum(-20.0);
  real_1->SetMaximum(+20.0);
  TH2D* imaginary_1 = new TH2D("imaginary_1", "imaginary part", divisions, min, max, divisions, min, max);
  imaginary_1->SetMinimum(-20.0);
  imaginary_1->SetMaximum(+20.0);
  TH2D* real_2 = new TH2D("real_2", "real part", divisions, min, max, divisions, min, max);
  TH2D* imaginary_2 = new TH2D("imaginary_2", "imaginary part", divisions, min, max, divisions, min, max);
  TH2D* real_3 = new TH2D("real_3", "real part", divisions, min, max, divisions, min, max);
  TH2D* imaginary_3 = new TH2D("imaginary_3", "imaginary part", divisions, min, max, divisions, min, max);
  double r, i;
  complex<double> result, arg;
  for (unsigned int re = 0; re < divisions; re++)
    for (unsigned int im = 0; im < divisions; im++)
      {
	r = min + (max - min)/divisions*re;
	i = min + (max - min)/divisions*im;
	arg = complex<double>(r, i);
	result = 1.0-cErrF(arg);
	//if (!std::isnan(result.real()))
	real_1->SetBinContent(re, im, result.real());
	imaginary_1->SetBinContent(re, im, result.imag());
	result = 1.0-cErrF_2(arg);
	//if (!std::isnan(result.real()))
	real_2->SetBinContent(re, im, result.real());
	imaginary_2->SetBinContent(re, im, result.imag());
	result = 1.0-cErrF_3(arg);
	//if (!std::isnan(result.real()))
	real_3->SetBinContent(re, im, result.real());
	imaginary_3->SetBinContent(re, im, result.imag());
      }
  
  TCanvas *canvas_1 = new TCanvas("canvas_1", "canvas 1",1600,1200) ;
  canvas_1->Divide(2,1);
  canvas_1->cd(1);
  real_1->Draw("COLZ");
  canvas_1->cd(2);
  imaginary_1->Draw("COLZ");
  canvas_1->Update();
  canvas_1->Print("error_1.eps", "eps");

  TCanvas *canvas_2 = new TCanvas("canvas_2", "canvas 2",1600,1200) ;
  canvas_2->Divide(2,1);
  canvas_2->cd(1);
  real_2->Draw("COLZ");
  canvas_2->cd(2);
  imaginary_2->Draw("COLZ");
  canvas_2->Update();
  canvas_2->Print("error_2.eps", "eps");

  TCanvas *canvas_3 = new TCanvas("canvas_3", "canvas 3",1600,1200) ;
  canvas_3->Divide(2,1);
  canvas_3->cd(1);
  real_3->Draw("COLZ");
  canvas_3->cd(2);
  imaginary_3->Draw("COLZ");
  canvas_3->Update();
  canvas_3->Print("error_3.eps", "eps");

  TH2D* real_12 = new TH2D((*real_1) / (*real_2));
  TH2D* imaginary_12 = new TH2D((*imaginary_1) / (*imaginary_2));
  TCanvas *canvas_12 = new TCanvas("canvas_12", "canvas 1 - 2",1600,1200) ;
  canvas_12->Divide(2,1);
  canvas_12->cd(1);
  real_12->Draw("COLZ");
  canvas_12->cd(2);
  imaginary_12->Draw("COLZ");
  canvas_12->Update();
  canvas_12->Print("error_1_over_2.eps", "eps");

  TH2D* real_23 = new TH2D((*real_2) / (*real_3));
  TH2D* imaginary_23 = new TH2D((*imaginary_2) / (*imaginary_3));
  TCanvas *canvas_23 = new TCanvas("canvas_23", "canvas 1 - 2",1600,1200) ;
  canvas_23->Divide(2,1);
  canvas_23->cd(1);
  real_23->Draw("COLZ");
  canvas_23->cd(2);
  imaginary_23->Draw("COLZ");
  canvas_23->Update();
  canvas_23->Print("error_2_over_3.eps", "eps");

  TH2D* real_13 = new TH2D((*real_1) / (*real_3));
  TH2D* imaginary_13 = new TH2D((*imaginary_1) / (*imaginary_3));
  TCanvas *canvas_13 = new TCanvas("canvas_13", "canvas 1 - 2",1600,1200) ;
  canvas_13->Divide(2,1);
  canvas_13->cd(1);
  real_13->Draw("COLZ");
  canvas_13->cd(2);
  imaginary_13->Draw("COLZ");
  canvas_13->Update();
  canvas_13->Print("error_1_over_3.eps", "eps");

}


void test_convolutions()
{ 
  //gStyle->SetPalette(1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  unsigned int divisions=1000;
  double min=-TMath::Pi();
  double max=+2.0*TMath::Pi();
  TH1D* conv_sin = new TH1D("conv_sin", "", divisions, min, max);
  TH1D* conv_cos = new TH1D("conv_cos", "", divisions, min, max);
  TH1D* conv_sin_2 = new TH1D("conv_sin_2", "", divisions, min, max);
  TH1D* conv_cos_2 = new TH1D("conv_cos_2", "", divisions, min, max);
  TH1D* conv_2_sin = new TH1D("conv_2_sin", "", divisions, min, max);
  TH1D* conv_2_cos = new TH1D("conv_2_cos", "", divisions, min, max);

  TH1D* conv_cos_2_sin = new TH1D("conv_cos_2_sin", "", divisions, min, max);
  TH1D* conv_sin_3 = new TH1D("conv_sin_3", "", divisions, min, max);
  TH1D* conv_2_sin_sin = new TH1D("conv_2_sin_sin", "", divisions, min, max);
  
  double sigma = 0.03;
  for (unsigned int i = 0; i < divisions; i++)
    {
      double angle = min+(max - min)/double(divisions)*i;
      conv_sin->SetBinContent(i+1, convoluted_sin(angle, sigma));
      conv_cos->SetBinContent(i+1, convoluted_cos(angle, sigma));
      conv_sin_2->SetBinContent(i+1, convoluted_sin_2(angle, sigma));
      conv_cos_2->SetBinContent(i+1, convoluted_cos_2(angle, sigma));
      conv_2_sin->SetBinContent(i+1, convoluted_2_sin(angle, sigma));
      conv_2_cos->SetBinContent(i+1, convoluted_2_cos(angle, sigma));
      
      conv_cos_2_sin->SetBinContent(i+1, convoluted_cos_2_sin(angle, sigma));
      conv_sin_3->SetBinContent(i+1, convoluted_sin_3(angle, sigma));
      //cout << convoluted_sin_3(angle, sigma) << endl;
      conv_2_sin_sin->SetBinContent(i+1, convoluted_2_sin_sin(angle, sigma));
  
    }

  conv_sin->SetLineColor(2);
  conv_cos->SetLineColor(3);
  conv_sin_2->SetLineColor(4);
  conv_cos_2->SetLineColor(5);
  conv_2_sin->SetLineColor(6);
  conv_2_cos->SetLineColor(7);
  conv_cos_2_sin->SetLineColor(8);
  conv_sin_3->SetLineColor(9);
  conv_2_sin_sin->SetLineColor(15);

  TCanvas *canvas = new TCanvas("canvas", "canvas",1600,1200) ;  
  canvas->cd();
  conv_sin->SetMinimum(-1.0);
  conv_sin->SetMaximum(1.0);
  conv_sin->Draw();
  conv_cos->Draw("same");
  conv_sin_2->Draw("same");
  conv_cos_2->Draw("same");
  conv_2_sin->Draw("same");
  conv_2_cos->Draw("same");
  conv_cos_2_sin->Draw("same");
  conv_sin_3->Draw("same");
  conv_2_sin_sin->Draw("same");
  TLegend* conv_legend = new TLegend(0.7,0.6,0.9,0.9);
  conv_legend->AddEntry(conv_sin, "sin #Theta", "L");
  conv_legend->AddEntry(conv_cos, "cos #Theta", "L");
  conv_legend->AddEntry(conv_sin_2, "sin #Theta^{2}", "L");
  conv_legend->AddEntry(conv_cos_2, "cos #Theta^{2}", "L");
  conv_legend->AddEntry(conv_2_sin, "sin 2#Theta", "L");
  conv_legend->AddEntry(conv_2_cos, "cos 2#Theta", "L");
  conv_legend->AddEntry(conv_cos_2_sin, "cos #Theta^{2} sin #Theta", "L");
  conv_legend->AddEntry(conv_sin_3, "sin #Theta ^{3}", "L");
  conv_legend->AddEntry(conv_2_sin_sin, "sin 2#Theta sin #Theta", "L");
  conv_legend->Draw();
  canvas->Print("convolutions.eps", "eps");

  delete conv_sin;
  delete conv_cos;
  delete conv_sin_2;
  delete conv_cos_2;
  delete conv_2_sin;
  delete conv_2_cos;
  delete conv_cos_2_sin;
  delete conv_sin_3;
  delete conv_2_sin_sin;

  delete canvas;

  TH1D* conv_cos_2_a = new TH1D("conv_cos_2_a", "", divisions, -1.0, +1.0);
  TH1D* conv_cos_2_b = new TH1D("conv_cos_2_b", "", divisions, -1.0, +1.0);
  TH1D* conv_cos_2_c = new TH1D("conv_cos_2_c", "", divisions, -1.0, +1.0);
  TH1D* conv_cos_2_d = new TH1D("conv_cos_2_d", "", divisions, -1.0, +1.0);
  TH1D* conv_cos_2_e = new TH1D("conv_cos_2_e", "", divisions, -1.0, +1.0);
  conv_cos_2_a->SetMinimum(-1.0);
  conv_cos_2_a->SetMaximum(1.0);
  for (unsigned int i = 0; i < divisions; i++)
    {
      double cos_theta = -1.0 + 2.0/divisions*(i+0.5);
      double angle = acos(cos_theta);
      conv_cos_2_a->SetBinContent(i+1, cos_theta*cos_theta);
      conv_cos_2_b->SetBinContent(i+1, convoluted_cos_2(angle, sigma));
      conv_cos_2_c->SetBinContent(i+1, convoluted_cos_2(angle, sigma)
				  +convoluted_cos_2(-angle, sigma)
				  +convoluted_cos_2(2.0*TMath::Pi()-angle, sigma));
      conv_cos_2_d->SetBinContent(i+1,1.0/sin(angle)*convoluted_cos_2_sin(angle, sigma));
      conv_cos_2_e->SetBinContent(i+1,1.0/sin(angle)*(convoluted_cos_2_sin(angle, sigma)
						      +convoluted_cos_2_sin(-angle, sigma)
						      +convoluted_cos_2_sin(2.0*TMath::Pi()-angle, sigma)
						      ));
    }

  TCanvas *canvas2 = new TCanvas("canvas2", "canvas",1600,1200);  
  canvas2->cd();
  conv_cos_2_b->SetLineColor(2);
  conv_cos_2_c->SetLineColor(3);
  conv_cos_2_d->SetLineColor(4);
  conv_cos_2_e->SetLineColor(5);
  conv_cos_2_a->Draw();
  conv_cos_2_b->Draw("same");
  conv_cos_2_c->Draw("same");
  conv_cos_2_d->Draw("same");
  conv_cos_2_e->Draw("same");
  canvas2->Print("overcos.eps", "eps");
  delete conv_cos_2_a;
  delete conv_cos_2_b;
  delete conv_cos_2_c;
  delete conv_cos_2_d;
  delete conv_cos_2_e;
  delete canvas2;
}

void test_int_convolutions()
{ 
  //gStyle->SetPalette(1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  unsigned int divisions=1000;
  double min=-TMath::Pi();
  double max=+2.0*TMath::Pi();
  TH1D* conv_sin = new TH1D("conv_sin", "", divisions, min, max);
  TH1D* conv_cos = new TH1D("conv_cos", "", divisions, min, max);
  TH1D* conv_sin_2 = new TH1D("conv_sin_2", "", divisions, min, max);
  TH1D* conv_cos_2 = new TH1D("conv_cos_2", "", divisions, min, max);
  TH1D* conv_2_sin = new TH1D("conv_2_sin", "", divisions, min, max);
  TH1D* conv_2_cos = new TH1D("conv_2_cos", "", divisions, min, max);

  TH1D* conv_cos_2_sin = new TH1D("conv_cos_2_sin", "", divisions, min, max);
  TH1D* conv_sin_3 = new TH1D("conv_sin_3", "", divisions, min, max);
  TH1D* conv_2_sin_sin = new TH1D("conv_2_sin_sin", "", divisions, min, max);
  
  double sigma = 0.03;
  for (unsigned int i = 0; i < divisions; i++)
    {
      double angle = min+(max - min)/double(divisions)*i;
      conv_sin->SetBinContent(i+1, int_convoluted_sin(angle, sigma)-int_convoluted_sin(min, sigma));
      
      // conv_cos->SetBinContent(i+1, convoluted_cos(angle, sigma));
      conv_sin_2->SetBinContent(i+1, int_convoluted_sin_2(angle, sigma) - int_convoluted_sin_2(min, sigma));
      // conv_cos_2->SetBinContent(i+1, convoluted_cos_2(angle, sigma));
      conv_2_sin->SetBinContent(i+1, int_convoluted_2_sin(angle, sigma)-int_convoluted_2_sin(min, sigma));
      // conv_2_cos->SetBinContent(i+1, convoluted_2_cos(angle, sigma));
      
      conv_cos_2_sin->SetBinContent(i+1, int_convoluted_cos_2_sin(angle, sigma) - int_convoluted_cos_2_sin(min, sigma));
      conv_sin_3->SetBinContent(i+1, int_convoluted_sin_3(angle, sigma)-int_convoluted_sin_3(min, sigma));
      conv_2_sin_sin->SetBinContent(i+1, int_convoluted_2_sin_sin(angle, sigma) - int_convoluted_2_sin_sin(min, sigma));
    }

  conv_sin->SetLineColor(2);
  conv_cos->SetLineColor(3);
  conv_sin_2->SetLineColor(4);
  conv_cos_2->SetLineColor(5);
  conv_2_sin->SetLineColor(6);
  conv_2_cos->SetLineColor(7);
  conv_cos_2_sin->SetLineColor(8);
  conv_sin_3->SetLineColor(9);
  conv_2_sin_sin->SetLineColor(15);

  TCanvas *canvas = new TCanvas("canvas", "canvas",1600,1200) ;  
  canvas->cd();
  conv_sin->SetMinimum(-0.5);
  conv_sin->SetMaximum(3.0);
  conv_sin->Draw();
  
  // conv_cos->Draw("same");
  conv_sin_2->Draw("same");
  // conv_cos_2->Draw("same");
  conv_2_sin->Draw("same");
  // conv_2_cos->Draw("same");
  conv_cos_2_sin->Draw("same");
  conv_sin_3->Draw("same");
  conv_2_sin_sin->Draw("same");
  TLegend* conv_legend = new TLegend(0.7,0.6,0.9,0.9);
  conv_legend->AddEntry(conv_sin, "sin #Theta", "L");
  conv_legend->AddEntry(conv_cos, "cos #Theta", "L");
  conv_legend->AddEntry(conv_sin_2, "sin #Theta^{2}", "L");
  conv_legend->AddEntry(conv_cos_2, "cos #Theta^{2}", "L");
  conv_legend->AddEntry(conv_2_sin, "sin 2#Theta", "L");
  conv_legend->AddEntry(conv_2_cos, "cos 2#Theta", "L");
  conv_legend->AddEntry(conv_cos_2_sin, "cos #Theta^{2} sin #Theta", "L");
  conv_legend->AddEntry(conv_sin_3, "sin #Theta ^{3}", "L");
  conv_legend->AddEntry(conv_2_sin_sin, "sin 2#Theta sin #Theta", "L");
  conv_legend->Draw();
  canvas->Print("int_convolutions.eps", "eps");

  delete conv_sin;
  delete conv_cos;
  delete conv_sin_2;
  delete conv_cos_2;
  delete conv_2_sin;
  delete conv_2_cos;
  delete conv_cos_2_sin;
  delete conv_sin_3;
  delete conv_2_sin_sin;
  delete canvas;
}

//computes erfc, because erfc(x)=exp(-x^2)*w(i*x)
complex<double> cErrF(const complex<double>& x)
{
  const complex<double> i(0.0,1.0);
  complex<double> z(i*x);
  //complex<double> z(-x.imag(), x.real());
  complex<double> result = exp(-x*x)*wErrF(z);
  return result;
}

//computes erfc, because erfc(x)=exp(-x^2)*w(i*x)
complex<double> cErrF_2(const complex<double>& x)
{
  const complex<double> i(0.0,1.0);
  complex<double> z(i*x);
  complex<double> result = exp(-x*x)*Faddeeva_2(z);
  
  if (x.real() > 20.0)// && fabs(x.imag()) < 20.0)
    result = 0.0;
  if (x.real() < -20.0)// && fabs(x.imag()) < 20.0)
    result = 2.0;
  

  return result;
}

complex<double> ErrF_2(const complex<double>& x)
{
  return 1.0 - cErrF_2(x);
}

//computes erfc, because erfc(x)=exp(-x^2)*w(i*x)
complex<double> cErrF_3(const complex<double>& x)
{
  const complex<double> i(0.0,1.0);
  complex<double> z(i*x);
  complex<double> result = exp(-x*x)*nwwerf(z);
  return result;
  //careful! the fortran version is apparently not thread safe
}

//from matpack
//precision 1e-14
complex<double> Faddeeva_2 (const complex<double>& z)
{
  // table 1: coefficients for h = 0.5  
  static double n1[12] =
    { 0.25, 1.0, 2.25, 4.0, 6.25, 9.0, 12.25, 16.0,
      20.25, 25.0, 30.25, 36.0 };
  static double e1[12] =
    { 0.7788007830714049,    0.3678794411714423,
      1.053992245618643e-1,  1.831563888873418e-2,
      1.930454136227709e-3,  1.234098040866795e-4,
      4.785117392129009e-6,  1.125351747192591e-7,
      1.605228055185612e-9,  1.388794386496402e-11,
      7.287724095819692e-14, 2.319522830243569e-16 };

  // table 2: coefficients for h = 0.53 
  static double n2[12] =
    { 0.2809, 1.1236, 2.5281, 4.4944, 7.0225, 10.1124,
      13.7641, 17.9776, 22.7529, 28.09, 33.9889, 40.4496 };
  static double e2[12] =
    { 0.7551038420890235,    0.3251072991205958, 
      7.981051630007964e-2,  1.117138143353082e-2,
      0.891593719995219e-3,  4.057331392320188e-5,
      1.052755021528803e-6,  1.557498087816203e-8,
      1.313835773243312e-10, 6.319285885175346e-13,
      1.733038792213266e-15, 2.709954036083074e-18 };
    
  // tables for Pade approximation 
  static double C[7] =
    { 65536.0, -2885792.0, 69973904.0, -791494704.0,
      8962513560.0, -32794651890.0, 175685635125.0 };
  static double D[7] =
    { 192192.0, 8648640.0, 183783600.0, 2329725600.0,
      18332414100.0, 84329104860.0, 175685635125.0 };


  double *n,*e,t,u,r,s,d,f,g,h;
  complex<double> c,d2,v,w,zz;
  int i;
    
  // use Pade approximation 
  s = norm(z);
  if (s < 1e-7) {
    zz = z*z;
    v  = exp(zz);
    c  = C[0];
    d2 = D[0];
    for (i = 1; i <= 6; i++) {
      c  = c  * zz + C[i];
      d2 = d2 * zz + D[i];
    }
    w = 1.0 / v + complex<double>(0.0,m_2_sqrtpi) * c/d2 * z * v;
    return w;

    // use trapezoid rule 
  } else {

    // select default table 1 
    n = n1;
    e = e1;
    r = M_1_PI * 0.5;
 
    // if z is too close to a pole select table 2 
    if (fabs(imag(z)) < 0.01 && fabs(real(z)) < 6.01) {
      h = modf(2*fabs(real(z)),&g);
      if (h < 0.02 || h > 0.98) {
	n = n2;
	e = e2;
	r = M_1_PI * 0.53;
      }
    }
        
    d = (imag(z) - real(z)) * (imag(z) + real(z));
    f = 4 * real(z) * real(z) * imag(z) * imag(z);

    g = h = 0.0;
    for (i = 0; i < 12; i++) {
      t = d + n[i];
      u = e[i] / (t * t + f);
      g += (s + n[i]) * u;
      h += (s - n[i]) * u;
    }
    u = 1 / s;
    c = r * complex<double>(imag(z) * (u + 2.0 * g),
			    real(z) * (u + 2.0 * h) );
        
    if (imag(z) < M_2PI) {
      s = 2.0 / r;
      t = s * real(z);
      u = s * imag(z);
      s = sin(t);
      h = cos(t);
      f = exp(- u) - h;
      g = 2.0 * exp(d-u) / (s * s + f * f);
      u = 2.0 * real(z) * imag(z);
      h = cos(u);
      t = sin(u);
      c += g * complex<double>( (h * f - t * s), -(h * s + t * f));
    }
    return c;
  }
}


//this is the complex error function w(z)
complex<double> wErrF(const complex<double>& arg)
{
  const double c1 = 7.4;
  const double c2 = 8.3;
  const double c3 = 0.3125;
  const double c4 = 1.6;
  
  const double c = 1.12837916709551257;
  const double p = 46768052394588893.38251791464692105662;
  
  //----------------------------------------------------------------------------
  complex<double>  r[37];
  complex<double>  zh, s, t, v;
  complex<double>  hh, den;
  
  double xl = p;
  
  double x  = arg.real();
  double y  = arg.imag();
  double xa = fabs(x);
  double ya = fabs(y);
  
  if (ya < c1 && xa < c2) {
    zh = complex<double>(ya+c4, xa);
    
    for (int n=35; n>-1; n--) {
      t    = zh + std::conj(double(n+1)*r[n+1]);
      r[n] = 0.5 / (t.real()*t.real() + t.imag()*t.imag()) * t;
    }
    
    for (int n=32; n>-1; n--) {
      xl =  xl*c3;
      s  = (xl+s)*r[n];
    }
    
    v = c * s;
  }
  else {
    zh = complex<double>(ya, xa);
    
    for (int n=8; n>-1; n--) {
      t    = zh + std::conj(double(n+1)*r[0]);
      r[0] = 0.5 / (t.real()*t.real() + t.imag()*t.imag()) * t;
    }
    
    v = c * r[0];
  }
  
  if (ya == 0)
    v = complex<double>(exp(-xa*xa), v.imag());
  
  if (y < 0) {
    hh = complex<double>(xa, ya);
    complex<double> hprod = -hh*hh;
    complex<double> anexp(exp(hprod.real())*cos(hprod.imag()),exp(hprod.real())*sin(hprod.imag()));
    v = 2.0 * anexp - v;
  }
  else if (x < 0)
    v = conj(v);
  
  if(std::isnan(v.real()))
    cout << "wErrF: WARNING! the result is NaN, do not trust anything!"
	 << endl;

  return v;

}

std::complex<double> nwwerf(const std::complex<double> z) {
  std::complex<double>  zh,r[38],s,t,v;

  const double z1 = 1;
  const double hf = z1/2;
  const double z10 = 10;
  const double c1 = 74/z10;
  const double c2 = 83/z10;
  const double c3 = z10/32;
  const double c4 = 16/z10;
  const double c = 1.12837916709551257;
  const double p = pow(2.0*c4,33);

  double x=z.real();
  double y=z.imag();
  double xa=(x >= 0) ? x : -x;
  double ya=(y >= 0) ? y : -y;
  if(ya < c1 && xa < c2){
    zh = std::complex<double>(ya+c4,xa);
    r[37]= std::complex<double>(0,0);
    //       do 1 n = 36,1,-1
    for(int n = 36; n>0; n--){   
      t=zh+double(n)*std::conj(r[n+1]);
      r[n]=hf*t/std::norm(t);
    }
    double xl=p;
    s=std::complex<double>(0,0);
    //    do 2 n = 33,1,-1
    for(int k=33; k>0; k--){
      xl=c3*xl;
      s=r[k]*(s+xl);
    }
    v=c*s;
  }
  else{
    zh=std::complex<double>(ya,xa);
    r[1]=std::complex<double>(0,0);
    //       do 3 n = 9,1,-1
    for(int n=9;n>0;n--){
      t=zh+double(n)*std::conj(r[1]);
      r[1]=hf*t/std::norm(t);
    }
    v=c*r[1];
  }
  if(ya == 0) v= std::complex<double>(exp(-xa*xa),v.imag());
  if(y < 0) {
    v=2.0*std::exp(std::complex<double>(-xa,-ya)*std::complex<double>(xa,ya))-v;
    if(x > 0) v=std::conj(v);
  }
  else{
    if(x < 0) v=std::conj(v);
  }
  return v;
}

void norm_convoluted_exp_sincos(double tau, double deltam, double& one, double& two)
{
  //double tau = 1.0/gamma();//0.5*(ctau_l() + ctau_h());//is this correct? 1/(0.5(gammal + gammah))
  double tau_2 = tau*tau;
  //double deltam = delta_m();
  double deltam_2 = deltam*deltam;
  one = tau_2 * deltam / (tau_2 * deltam_2 + 1.0);//sin
  two = tau / (tau_2 * deltam_2 + 1.0);//cos
}

//calculates the gaussian convolution sin(angle) from 0 to pi
double convoluted_sin(double angle, double sigma)
{
  double b = -angle/sqrt(2.0)/sigma + TMath::Pi()/sqrt(2.0)/sigma;
  double c = angle/sqrt(2.0)/sigma;
  double a = sigma/sqrt(2.0);
  double sigma_2 = sigma*sigma;
  double result = 0.0;
  complex<double> arg1(b, a);
  complex<double> arg2(c, a);
  complex<double> erf1 = 1.0-cErrF_2(arg1);
  complex<double> erf2 = 1.0-cErrF_2(arg2);
  result = 0.5*sin(angle)*exp(-0.5*sigma_2)
    *(erf1.real()+erf2.real())
    +0.5*cos(angle)*exp(-0.5*sigma_2)
    *(-erf1.imag()+erf2.imag());
  //cout << "arg1: " << b << " arg2: " << c << " argim: " << a << "*i angle: " << angle << " sigma: " << sigma << " result: " << result << endl;
  return result;
}

double convoluted_cos_sin(double angle, double sigma)
{
  return 0.5*convoluted_2_sin(angle, sigma);
}

//calculates the gaussian convolution sin(angle) from 0 to pi
double convoluted_cos(double angle, double sigma)
{
  double b = -angle/sqrt(2.0)/sigma + TMath::Pi()/sqrt(2.0)/sigma;
  double c = -angle/sqrt(2.0)/sigma;//also changed sign here
  double a = sigma/sqrt(2.0);
  double sigma_2 = sigma*sigma;
  double result = 0.0;
  complex<double> arg1(b, a);
  complex<double> arg2(c, a);
  complex<double> erf1 = 1.0-cErrF_2(arg1);
  complex<double> erf2 = 1.0-cErrF_2(arg2);
  result = 0.5*sin(angle)*exp(-0.5*sigma_2)
    *(erf1.imag()-erf2.imag())
    +0.5*cos(angle)*exp(-0.5*sigma_2)
    *(erf1.real()-erf2.real());
  return result;
}

//calculates the gaussian convolution cos(angle)^2 from 0 to pi
double convoluted_cos_2(double angle, double sigma)
{
  double b = -angle/sqrt(2.0)/sigma + TMath::Pi()/sqrt(2.0)/sigma;
  double c = -angle/sqrt(2.0)/sigma;//changed the sign here, more consistent
  double a = sigma*sqrt(2.0);//!!!
  double sigma_2 = sigma*sigma;
  double result = 0.0;
  complex<double> arg1(b, a);
  complex<double> arg2(c, a);
  complex<double> erf1 = 1.0-cErrF_2(arg1);
  complex<double> erf2 = 1.0-cErrF_2(arg2);
  result = 0.25*sin(2.0*angle)*exp(-2.0*sigma_2)
    *(erf1.imag()-erf2.imag())
    +0.25*cos(2.0*angle)*exp(-2.0*sigma_2)
    *(erf1.real()-erf2.real())
    -0.25*erf(-b)
    +0.25*erf(-c);
  return result;
}

//calculates the gaussian convolution cos(angle)^2*sin(angle) from 0 to pi
double convoluted_cos_2_sin(double angle, double sigma)
{
  double b = -angle/sqrt(2.0)/sigma + TMath::Pi()/sqrt(2.0)/sigma;
  double c = -angle/sqrt(2.0)/sigma;//changed the sign here, more consistent
  
  double a = 3.0*sigma/sqrt(2.0);
  double sigma_2 = sigma*sigma;
  double result = 0.0;
  complex<double> arg1(b, a);
  complex<double> arg2(c, a);
  complex<double> erf1 = 1.0-cErrF_2(arg1);
  complex<double> erf2 = 1.0-cErrF_2(arg2);

  result = 1.0/8.0*sin(3.0*angle)*exp(-4.5*sigma_2)
    *(erf1.real()-erf2.real())
    +1.0/8.0*cos(3.0*angle)*exp(-4.5*sigma_2)
    *(-erf1.imag()+erf2.imag());

  a = sigma/sqrt(2.0);
  arg1 = complex<double>(b, a);
  arg2 = complex<double>(c, a);
  
  erf1 = 1.0-cErrF_2(arg1);
  erf2 = 1.0-cErrF_2(arg2);

  result += 1.0/8.0*sin(angle)*exp(-0.5*sigma_2)
    *(erf1.real()-erf2.real())
    +1.0/8.0*cos(angle)*exp(-0.5*sigma_2)
    *(-erf1.imag()+erf2.imag());
  return result;
}

//gaussian convolution sin(angle)*sin(2*angle) from 0 to pi
double convoluted_2_sin_sin(double angle, double sigma)
{
  double b = -angle/sqrt(2.0)/sigma + TMath::Pi()/sqrt(2.0)/sigma;
  double c = -angle/sqrt(2.0)/sigma;//changed the sign here, more consistent
  double a = 3.0*sigma/sqrt(2.0);
  double sigma_2 = sigma*sigma;
  double result = 0.0;
  complex<double> arg1(b, a);
  complex<double> arg2(c, a);
  complex<double> erf1 = 1.0-cErrF_2(arg1);
  complex<double> erf2 = 1.0-cErrF_2(arg2);
  result = 1.0/4.0*sin(3.0*angle)*exp(-4.5*sigma_2)
    *(-erf1.imag()+erf2.imag())
    +1.0/4.0*cos(3.0*angle)*exp(-4.5*sigma_2)
    *(-erf1.real()+erf2.real());

  a = sigma/sqrt(2.0);
  arg1 = complex<double>(b, a);
  arg2 = complex<double>(c, a);
  
  erf1 = 1.0-cErrF_2(arg1);
  erf2 = 1.0-cErrF_2(arg2);

  result += 1.0/4.0*sin(angle)*exp(-0.5*sigma_2)
    *(+erf1.imag()-erf2.imag())
    +1.0/4.0*cos(angle)*exp(-0.5*sigma_2)
    *(+erf1.real()-erf2.real());
  return result;//please recheck
}

//gaussian convolution sin(angle)^3 from 0 to pi
double convoluted_sin_3(double angle, double sigma)
{
  double b = -angle/sqrt(2.0)/sigma + TMath::Pi()/sqrt(2.0)/sigma;
  double c = -angle/sqrt(2.0)/sigma;//changed the sign here, more consistent
  double a = 3.0*sigma/sqrt(2.0);
  double sigma_2 = sigma*sigma;
  double result = 0.0;
  complex<double> arg1(b, a);
  complex<double> arg2(c, a);

  complex<double> erf1 = 1.0-cErrF_2(arg1);
  complex<double> erf2 = 1.0-cErrF_2(arg2);

  result = 1.0/8.0*sin(3.0*angle)*exp(-4.5*sigma_2)
    *(-erf1.real()+erf2.real())
    +1.0/8.0*cos(3.0*angle)*exp(-4.5*sigma_2)
    *(erf1.imag()-erf2.imag());

  a = sigma/sqrt(2.0);
  arg1 = complex<double>(b, a);
  arg2 = complex<double>(c, a);
  
  erf1 = 1.0-cErrF_2(arg1);
  erf2 = 1.0-cErrF_2(arg2);

  result += 3.0/8.0*sin(angle)*exp(-0.5*sigma_2)
    *(erf1.real()-erf2.real())
    +3.0/8.0*cos(angle)*exp(-0.5*sigma_2)
    *(-erf1.imag()+erf2.imag());
  return result;//this should be correct now->recheck
}

//calculates the gaussian convolution sin(angle)^2 from 0 to pi
double convoluted_sin_2(double angle, double sigma)
{
  double b = -angle/sqrt(2.0)/sigma + TMath::Pi()/sqrt(2.0)/sigma;
  double c = -angle/sqrt(2.0)/sigma;//changed the sign here, more consistent
  double a = sigma*sqrt(2.0);//!!!
  double sigma_2 = sigma*sigma;
  double result = 0.0;
  complex<double> arg1(b, a);
  complex<double> arg2(c, a);
  complex<double> erf1 = 1.0-cErrF_2(arg1);
  complex<double> erf2 = 1.0-cErrF_2(arg2);
  result = 0.25*sin(2.0*angle)*exp(-2.0*sigma_2)
    *(-erf1.imag()+erf2.imag())
    +0.25*cos(2.0*angle)*exp(-2.0*sigma_2)
    *(-erf1.real()+erf2.real())
    -0.25*erf(-b)
    +0.25*erf(-c);
  return result;
}

//calculates the gaussian convolution cos(2*angle) from 0 to pi
double convoluted_2_cos(double angle, double sigma)
{
  double b = -angle/sqrt(2.0)/sigma + TMath::Pi()/sqrt(2.0)/sigma;
  double c = -angle/sqrt(2.0)/sigma;//changed the sign here, more consistent
  double a = sigma*sqrt(2.0);//!!!
  double sigma_2 = sigma*sigma;
  double result = 0.0;
  complex<double> arg1(b, a);
  complex<double> arg2(c, a);
  complex<double> erf1 = 1.0-cErrF_2(arg1);
  complex<double> erf2 = 1.0-cErrF_2(arg2);
  result = 
    0.5*sin(2.0*angle)*exp(-2.0*sigma_2)
    *(erf1.imag()-erf2.imag())
    +0.5*cos(2.0*angle)*exp(-2.0*sigma_2)
    *(erf1.real()-erf2.real());
  return result;
}

//calculates the gaussian convolution sin(2*angle) from 0 to pi
double convoluted_2_sin(double angle, double sigma)
{
  double b = -angle/sqrt(2.0)/sigma + TMath::Pi()/sqrt(2.0)/sigma;
  double c = -angle/sqrt(2.0)/sigma;//changed the sign here, more consistent
  double a = sigma*sqrt(2.0);//rechcek this for all other convs!
  double sigma_2 = sigma*sigma;
  double result = 0.0;
  complex<double> arg1(b, a);
  complex<double> arg2(c, a);
  complex<double> erf1 = 1.0-cErrF_2(arg1);
  complex<double> erf2 = 1.0-cErrF_2(arg2);
  result = 0.5*sin(2.0*angle)*exp(-2.0*sigma_2)
    *(erf1.real()-erf2.real())//ok
    +0.5*cos(2.0*angle)*exp(-2.0*sigma_2)
    *(-erf1.imag()+erf2.imag());

  return result;
}

//integrate(sqrt(alpha/%pi)*sin(d)*exp(-alpha*(t-d)^2),d,-inf,inf);       
//integrate(1/sqrt(2*%pi)/sigma*sin(d)*exp(-(t-d)^2/2/sigma^2),d,-inf,inf); why does this not work?
double inf_convoluted_sin(double angle, double sigma)
{
  //%e^-(1/(4*alpha))*sin(t)
  return exp(-0.5*sigma*sigma)*sin(angle);
}

//integrate(sqrt(alpha/%pi)*cos(d)*exp(-alpha*(t-d)^2),d,-inf,inf);       
double inf_convoluted_cos(double angle, double sigma)
{
  // %e^-(1/(4*alpha))*cos(t)
  return exp(-0.5*sigma*sigma)*cos(angle);
}

//integrate(sqrt(alpha/%pi)*cos(d)^2*exp(-alpha*(t-d)^2),d,-inf,inf);       
double inf_convoluted_cos_2(double angle, double sigma)
{
  //%e^-(1/alpha)*cos(2*t)/2+1/2
  return 0.5*exp(-2.0*sigma*sigma)*cos(2.0*angle)+0.5;
}

//integrate(sqrt(alpha/%pi)*sin(d)^2*exp(-alpha*(t-d)^2),d,-inf,inf);       
double inf_convoluted_sin_2(double angle, double sigma)
{
  // 1/2-%e^-(1/alpha)*cos(2*t)/2
  return 0.5-0.5*exp(-2.0*sigma*sigma)*cos(2.0*angle);
}

//integrate(sqrt(alpha/%pi)*sin(2*d)*exp(-alpha*(t-d)^2),d,-inf,inf);       
double inf_convoluted_2_sin(double angle, double sigma)
{
  //%e^-(1/alpha)*sin(2*t)
  return exp(-2.0*sigma*sigma)*sin(2.0*angle);
}

//string(expand(integrate(integrate(sqrt(alpha/%pi)*exp(-alpha*(t-d)^2),d,-inf,inf),t)));       
double int_inf_convoluted_const(double angle, double sigma)
{
  return angle;
}

//string(expand(integrate(integrate(sqrt(alpha/%pi)*sin(d)*exp(-alpha*(t-d)^2),d,-inf,inf),t)));       
double int_inf_convoluted_sin(double angle, double sigma)
{
  //-%e^-(1/(4*alpha))*cos(t)
  return -exp(-0.5*sigma*sigma)*cos(angle);
}

//string(expand(integrate(integrate(sqrt(alpha/%pi)*cos(d)*exp(-alpha*(t-d)^2),d,-inf,inf),t)));       
double int_inf_convoluted_cos(double angle, double sigma)
{
  //%e^-(1/(4*alpha))*sin(t)
  return exp(-0.5*sigma*sigma)*sin(angle);
}

//string(expand(integrate(integrate(sqrt(alpha/%pi)*cos(d)^2*exp(-alpha*(t-d)^2),d,-inf,inf),t)));       
double int_inf_convoluted_cos_2(double angle, double sigma)
{
  //%e^-(1/alpha)*sin(2*t)/4+t/2
  return 0.25*exp(-2.0*sigma*sigma)*sin(2.0*angle)+0.5*angle;
}

//string(expand(integrate(integrate(sqrt(alpha/%pi)*sin(d)^2*exp(-alpha*(t-d)^2),d,-inf,inf),t)));       
double int_inf_convoluted_sin_2(double angle, double sigma)
{
  //t/2-%e^-(1/alpha)*sin(2*t)/4
  return 0.5*angle-0.25*exp(-2.0*sigma*sigma)*sin(2.0*angle);
}

//string(expand(integrate(integrate(sqrt(alpha/%pi)*sin(2*d)*exp(-alpha*(t-d)^2),d,-inf,inf),t)));       
double int_inf_convoluted_2_sin(double angle, double sigma)
{
  // -%e^-(1/alpha)*cos(2*t)/2
  return -0.5*exp(-2.0*sigma*sigma)*cos(2.0*angle);
}

//expand(integrate(integrate(sqrt(alpha/%pi)*sin(d)*exp(-alpha*(t-d)^2),d,0,%pi), t));       
double int_convoluted_sin(double angle, double sigma)
{
  double alpha = 1.0/(2.0*sigma*sigma);
  double sa = sqrt(alpha);
  double sat = sa*angle;
  double sapi = sa*TMath::Pi();
  double overtwosa = 1.0/(2.0*sa);
  double expa = exp(-1.0/(4.0*alpha));
  double cost = cos(angle);
  double sint = sin(angle);
  double result = 0.0;
  result += 0.5*expa*cost*ErrF_2(complex<double>(sat - sapi, overtwosa)).real();
  result += -0.5*expa*sint*ErrF_2(complex<double>(sat - sapi, overtwosa)).imag();
  result += -0.5*expa*cost*ErrF_2(complex<double>(sat, overtwosa)).real();
  result += 0.5*expa*sint*ErrF_2(complex<double>(sat, overtwosa)).imag();
  result += 0.5*erf(sat - sapi);
  result += 0.5*erf(sat);
  return result;
}

double int_convoluted_cos_sin(double angle, double sigma)
{
  return 0.5*int_convoluted_2_sin(angle, sigma);
}

//expand(integrate(integrate(sqrt(alpha/%pi)*cos(d)^2*sin(d)*exp(-alpha*(t-d)^2),d,0,%pi), t));     
double int_convoluted_cos_2_sin(double angle, double sigma)
{
  double alpha = 1.0/(2.0*sigma*sigma);
  double sa = sqrt(alpha);
  double sat = sa*angle;
  double sapi = sa*TMath::Pi();
  double overtwosa = 1.0/(2.0*sa);
  double threeovertwosa = 3.0/(2.0*sa);
  double expa = exp(-1.0/(4.0*alpha));
  double expninea = exp(-9.0/(4.0*alpha));
  double cost = cos(angle);
  double sint = sin(angle);
  double costhreet = cos(3.0*angle);
  double sinthreet = sin(3.0*angle);
  double result = 0.0;
  
  result += -(1.0/24.0)*expninea*sinthreet*ErrF_2(complex<double>(sat - sapi, threeovertwosa)).imag();
  result += (1.0/24.0)*expninea*costhreet*ErrF_2(complex<double>(sat - sapi, threeovertwosa)).real();

  result += -(1.0/8.0)*expa*sint*ErrF_2(complex<double>(sat - sapi, overtwosa)).imag();
  result += (1.0/8.0)*expa*cost*ErrF_2(complex<double>(sat - sapi, overtwosa)).real();

  result += (1.0/24.0)*expninea*sinthreet*ErrF_2(complex<double>(sat, threeovertwosa)).imag();
  result += -(1.0/24.0)*expninea*costhreet*ErrF_2(complex<double>(sat, threeovertwosa)).real();

  result += (1.0/8.0)*expa*sint*ErrF_2(complex<double>(sat, overtwosa)).imag();
  result += -(1.0/8.0)*expa*cost*ErrF_2(complex<double>(sat, overtwosa)).real();

  result += (1.0/6.0)*erf(sat-sapi);
  result += (1.0/6.0)*erf(sat);

  return result;

}

//expand(integrate(integrate(sqrt(alpha/%pi)*sin(d)*sin(2*d)*exp(-alpha*(t-d)^2),d,0,%pi), t));
double int_convoluted_2_sin_sin(double angle, double sigma)
{
  double alpha = 1.0/(2.0*sigma*sigma);
  double sa = sqrt(alpha);
  double sat = sa*angle;
  double sapi = sa*TMath::Pi();
  double overtwosa = 1.0/(2.0*sa);
  double threeovertwosa = 3.0/(2.0*sa);
  double expa = exp(-1.0/(4.0*alpha));
  double expninea = exp(-9.0/(4.0*alpha));
  double cost = cos(angle);
  double sint = sin(angle);
  double costhreet = cos(3.0*angle);
  double sinthreet = sin(3.0*angle);
  double result = 0.0;
  
  result += (1.0/12.0)*expninea*sinthreet*ErrF_2(complex<double>(sat - sapi, threeovertwosa)).real();
  result += (1.0/12.0)*expninea*costhreet*ErrF_2(complex<double>(sat - sapi, threeovertwosa)).imag();

  result += -(1.0/4.0)*expa*sint*ErrF_2(complex<double>(sat - sapi, overtwosa)).real();
  result += -(1.0/4.0)*expa*cost*ErrF_2(complex<double>(sat - sapi, overtwosa)).imag();

  result += -(1.0/12.0)*expninea*sinthreet*ErrF_2(complex<double>(sat, threeovertwosa)).real();
  result += -(1.0/12.0)*expninea*costhreet*ErrF_2(complex<double>(sat, threeovertwosa)).imag();

  result += (1.0/4.0)*expa*sint*ErrF_2(complex<double>(sat, overtwosa)).real();
  result += (1.0/4.0)*expa*cost*ErrF_2(complex<double>(sat, overtwosa)).imag();

  return result;
}

//expand(integrate(integrate(sqrt(alpha/%pi)*sin(d)^3*exp(-alpha*(t-d)^2),d,0,%pi), t));  
double int_convoluted_sin_3(double angle, double sigma)
{
  double alpha = 1.0/(2.0*sigma*sigma);
  double sa = sqrt(alpha);
  double sat = sa*angle;
  double sapi = sa*TMath::Pi();
  double overtwosa = 1.0/(2.0*sa);
  double threeovertwosa = 3.0/(2.0*sa);
  double expa = exp(-1.0/(4.0*alpha));
  double expninea = exp(-9.0/(4.0*alpha));
  double cost = cos(angle);
  double sint = sin(angle);
  double costhreet = cos(3.0*angle);
  double sinthreet = sin(3.0*angle);
  double result = 0.0;
  result += (1.0/24.0)*expninea*sinthreet*ErrF_2(complex<double>(sat - sapi, threeovertwosa)).imag();
  result += -(1.0/24.0)*expninea*costhreet*ErrF_2(complex<double>(sat - sapi, threeovertwosa)).real();
  result += -(3.0/8.0)*expa*sint*ErrF_2(complex<double>(sat - sapi, overtwosa)).imag();
  result += (3.0/8.0)*expa*cost*ErrF_2(complex<double>(sat - sapi, overtwosa)).real();

  result += -(1.0/24.0)*expninea*sinthreet*ErrF_2(complex<double>(sat, threeovertwosa)).imag();
  result += (1.0/24.0)*expninea*costhreet*ErrF_2(complex<double>(sat, threeovertwosa)).real();
  result += +(3.0/8.0)*expa*sint*ErrF_2(complex<double>(sat, overtwosa)).imag();
  result += -(3.0/8.0)*expa*cost*ErrF_2(complex<double>(sat, overtwosa)).real();
  
  result += 1.0/3.0*erf(sat-sapi);
  result += 1.0/3.0*erf(sat);

  return result;  
}

//expand(integrate(integrate(sqrt(alpha/%pi)*sin(d)^2*exp(-alpha*(t-d)^2),d,0,%pi), t));
double int_convoluted_sin_2(double angle, double sigma)//used by s wave
{
  double t = angle;
  double alpha = 1.0/(2.0*sigma*sigma);
  double sa = sqrt(alpha);
  double sat = sa*angle;
  double sapi = sa*TMath::Pi();
  double oversa = 1.0/(sa);
  double expa = exp(-1.0/(alpha));
  double costwot = cos(2.0*angle);
  double sintwot = sin(2.0*angle);
  double pi = TMath::Pi();
  double spi = sqrt(pi);
  double result = 0.0;
  result += (1.0/8.0)*expa*sintwot*ErrF_2(complex<double>(sat - sapi, oversa)).real();
  result += (1.0/8.0)*expa*costwot*ErrF_2(complex<double>(sat - sapi, oversa)).imag();

  result += -(1.0/8.0)*expa*sintwot*ErrF_2(complex<double>(sat, oversa)).real();
  result += -(1.0/8.0)*expa*costwot*ErrF_2(complex<double>(sat, oversa)).imag();  
  
  result += -(t/4.0)*erf(sat-sapi);
  result += (pi/4.0)*erf(sat-sapi);
  result += (t/4.0)*erf(sat);
  
  result += -(1.0/(4.0*sa*spi))*exp(-alpha*t*t+2.0*pi*alpha*t-pi*pi*alpha);
  result += (1.0/(4.0*sa*spi))*exp(-alpha*t*t);

  return result;
  
}

//expand(integrate(integrate(sqrt(alpha/%pi)*sin(2*d)*exp(-alpha*(t-d)^2),d,0,%pi), t));
double int_convoluted_2_sin(double angle, double sigma)//used by s wave
{
  //double t = angle;
  double alpha = 1.0/(2.0*sigma*sigma);
  double sa = sqrt(alpha);
  double sat = sa*angle;
  double sapi = sa*TMath::Pi();
  double oversa = 1.0/(sa);
  double expa = exp(-1.0/(alpha));
  double costwot = cos(2.0*angle);
  double sintwot = sin(2.0*angle);
  //double pi = TMath::Pi();
  //double spi = sqrt(pi);
  double result = 0.0;
  result += -(1.0/4.0)*expa*sintwot*ErrF_2(complex<double>(sat - sapi, oversa)).imag();
  result += (1.0/4.0)*expa*costwot*ErrF_2(complex<double>(sat - sapi, oversa)).real();

  result += (1.0/4.0)*expa*sintwot*ErrF_2(complex<double>(sat, oversa)).imag();
  result += -(1.0/4.0)*expa*costwot*ErrF_2(complex<double>(sat, oversa)).real();

  result += -(1.0/4.0)*erf(sat-sapi);
  result += (1.0/4.0)*erf(sat);

  return result;

}

//for n:0 thru 8 do print("n=",n,": ",expand(integrate(sqrt(alpha/%pi)*exp(-alpha*(t-d)^2)*cos(t)^n, t, -inf, inf)));
//for k:0 thru n/2 do print(expand(1/2^(n-1)*binomial(n,k)*exp(-1/4/alpha*(n-2*k)^2)*cos((n-2*k)*t)));
//for k:0 thru 8/2 do print(expand(1/(2^(8-1))*binomial(8,k)*exp(-1/4/alpha*(8-2*k)^2)*cos((8-2*k)*t)));
//for k:0 thru 8/2 do print(expand(1/(2^(8))*binomial(8,k)*exp(-1/4/alpha*(8-2*k)^2)*cos((8-2*k)*t)));

//for k:0 thru 8/2 do print(expand(1/2^(8)*binomial(8-1,k)*exp(-1/4/alpha*(8-2*k)^2)*cos((8-2*k)*t)));

//for n:0 thru 20 do print("n=",n,": ",string(expand(integrate(integrate(sqrt(alpha/%pi)*exp(-alpha*(t)^2)*cos(t+d)^n*sin(d), t, -inf, inf), d, 0, %pi))));
double int_convoluted_sin_cos_n(int n, double sigma)
{
  assert(n >= 0);
  assert(n <= 20);
  double alpha = 1.0/(2.0*sigma*sigma);
  if (n % 2 == 1)
    return 0.0;
  switch (n) {
  case 0: return 2;
    break;
  case 2: return 1-exp(-1/alpha)/3;
    break;
  case 4: return -exp(-1/alpha)/3-exp(-4/alpha)/60+3.0/4.0;
    break;
  case 6: return -5*exp(-1/alpha)/16-exp(-4/alpha)/40-exp(-9/alpha)/560+5.0/8.0;
    break;
  case 8: return -7*exp(-1/alpha)/24-7*exp(-4/alpha)/240-exp(-9/alpha)/280-exp(-16/alpha)/4032+35.0/64.0;
    break;
  case 10: return -35*exp(-1/alpha)/128-exp(-4/alpha)/32-9*exp(-9/alpha)/1792-5*exp(-16/alpha)/8064-exp(-25/alpha)/25344+63.0/128.0;
    break;
  case 12: return -33*exp(-1/alpha)/128-33*exp(-4/alpha)/1024-11*exp(-9/alpha)/1792-11*exp(-16/alpha)/10752-exp(-25/alpha)/8448-exp(-36/alpha)/146432+231.0/512.0;
    break;
  case 14: return -1001*exp(-1/alpha)/4096-1001*exp(-4/alpha)/30720-143*exp(-9/alpha)/20480-13*exp(-16/alpha)/9216-91*exp(-25/alpha)/405504-7*exp(-36/alpha)/292864-exp(-49/alpha)/798720+429.0/1024.0;
    break;
  case 16: return -715*exp(-1/alpha)/3072-1001*exp(-4/alpha)/30720-39*exp(-9/alpha)/5120-65*exp(-16/alpha)/36864-35*exp(-25/alpha)/101376-15*exp(-36/alpha)/292864-exp(-49/alpha)/199680-exp(-64/alpha)/4177920+6435.0/16384.0;
    break;
  case 18: return -7293*exp(-1/alpha)/32768-663*exp(-4/alpha)/20480-663*exp(-9/alpha)/81920-17*exp(-16/alpha)/8192-85*exp(-25/alpha)/180224-51*exp(-36/alpha)/585728-51*exp(-49/alpha)/4259840-3*exp(-64/alpha)/2785280-exp(-81/alpha)/21168128+12155.0/32768.0;
    break;
  case 20: return -20995*exp(-1/alpha)/98304-4199*exp(-4/alpha)/131072-969*exp(-9/alpha)/114688-1615*exp(-16/alpha)/688128-323*exp(-25/alpha)/540672-4845*exp(-36/alpha)/37486592-19*exp(-49/alpha)/851968-19*exp(-64/alpha)/6684672-5*exp(-81/alpha)/21168128-exp(-100/alpha)/104595456+46189.0/131072.0;
    break;
  }
  return 0.0;
}

//for n:0 thru 20 do print("n=",n,": ",string(expand(integrate(integrate(sqrt(alpha/%pi)*exp(-alpha*(t)^2)*cos(t+d)^n*sin(d)^3, t, -inf, inf), d, 0, %pi))));
double int_convoluted_sin_3_cos_n(int n, double sigma)
{
  assert(n >= 0);
  assert(n <= 20);
  double alpha = 1.0/(2.0*sigma*sigma);
  if (n % 2 == 1)
    return 0.0;
  switch (n) {  
  case 0: return 4.0/3.0 ;
    break;
  case 2: return 2.0/3.0-2*exp(-1/alpha)/5 ;
    break;
  case 4: return -2*exp(-1/alpha)/5+exp(-4/alpha)/70+1.0/2.0 ;
    break;
  case 6: return -3*exp(-1/alpha)/8+3*exp(-4/alpha)/140+exp(-9/alpha)/2520+5.0/12.0 ;
    break;
  case 8: return -7*exp(-1/alpha)/20+exp(-4/alpha)/40+exp(-9/alpha)/1260+exp(-16/alpha)/36960+35.0/96.0 ;
    break;
  case 10: return -21*exp(-1/alpha)/64+3*exp(-4/alpha)/112+exp(-9/alpha)/896+exp(-16/alpha)/14784+exp(-25/alpha)/384384+21.0/64.0 ;
    break;
  case 12: return -99*exp(-1/alpha)/320+99*exp(-4/alpha)/3584+11*exp(-9/alpha)/8064+exp(-16/alpha)/8960+exp(-25/alpha)/128128+exp(-36/alpha)/3294720+77.0/256.0 ;
    break;
  case 14: return -3003*exp(-1/alpha)/10240+143*exp(-4/alpha)/5120+143*exp(-9/alpha)/92160+13*exp(-16/alpha)/84480+exp(-25/alpha)/67584+7*exp(-36/alpha)/6589440+exp(-49/alpha)/24893440+143.0/512.0 ;
    break;
  case 16: return -143*exp(-1/alpha)/512+143*exp(-4/alpha)/5120+13*exp(-9/alpha)/7680+13*exp(-16/alpha)/67584+5*exp(-25/alpha)/219648+exp(-36/alpha)/439296+exp(-49/alpha)/6223360+exp(-64/alpha)/171991040+2145.0/8192.0 ;
    break;
  case 18: return -21879*exp(-1/alpha)/81920+1989*exp(-4/alpha)/71680+221*exp(-9/alpha)/122880+51*exp(-16/alpha)/225280+255*exp(-25/alpha)/8200192+17*exp(-36/alpha)/4392960+9*exp(-49/alpha)/23429120+9*exp(-64/alpha)/343982080+exp(-81/alpha)/1111326720+12155.0/49152.0 ;
    break;
  case 20: return -4199*exp(-1/alpha)/16384+12597*exp(-4/alpha)/458752+323*exp(-9/alpha)/172032+323*exp(-16/alpha)/1261568+323*exp(-25/alpha)/8200192+323*exp(-36/alpha)/56229888+57*exp(-49/alpha)/79659008.0+exp(-64/alpha)/14483456.0+exp(-81/alpha)/222265344.0+exp(-100/alpha)/6816137216.0+46189.0/196608.0 ;
    break;
  }
  return 0.0;
}

//for n:0 thru 20 do print("n=",n,": ",string(expand(integrate(integrate(sqrt(alpha/%pi)*exp(-alpha*(t)^2)*cos(t+d)^n*sin(d)*cos(d)^2, t, -inf, inf), d, 0, %pi))));
double int_convoluted_cos_2_sin_cos_n(int n, double sigma)
{
  assert(n >= 0);
  assert(n <= 20);
  double alpha = 1.0/(2.0*sigma*sigma);
  if (n % 2 == 1)
    return 0.0;
  switch (n) {
  case 0: return 2.0/3.0 ;
    break;
  case 2: return exp(-1/alpha)/15+1.0/3.0 ;
    break;
  case 4: return exp(-1/alpha)/15-13*exp(-4/alpha)/420+1.0/4.0 ;
    break;
  case 6: return exp(-1/alpha)/16-13*exp(-4/alpha)/280-11*exp(-9/alpha)/5040+5.0/24.0 ;
    break;
  case 8: return 7*exp(-1/alpha)/120-13*exp(-4/alpha)/240-11*exp(-9/alpha)/2520-61*exp(-16/alpha)/221760+35.0/192.0 ;
    break;
  case 10: return 7*exp(-1/alpha)/128-13*exp(-4/alpha)/224-11*exp(-9/alpha)/1792-61*exp(-16/alpha)/88704-97*exp(-25/alpha)/2306304+21.0/128.0 ;
    break;
  case 12: return 33*exp(-1/alpha)/640-429*exp(-4/alpha)/7168-121*exp(-9/alpha)/16128-61*exp(-16/alpha)/53760-97*exp(-25/alpha)/768768-47*exp(-36/alpha)/6589440+77.0/512.0 ;
    break;
  case 14: return 1001*exp(-1/alpha)/20480-1859*exp(-4/alpha)/30720-1573*exp(-9/alpha)/184320-793*exp(-16/alpha)/506880-97*exp(-25/alpha)/405504-329*exp(-36/alpha)/13178880-193*exp(-49/alpha)/149360640+143.0/1024.0 ;
    break;
  case 16: return 143*exp(-1/alpha)/3072-1859*exp(-4/alpha)/30720-143*exp(-9/alpha)/15360-793*exp(-16/alpha)/405504-485*exp(-25/alpha)/1317888-47*exp(-36/alpha)/878592-193*exp(-49/alpha)/37340160-253*exp(-64/alpha)/1031946240+2145.0/16384.0 ;
    break;
  case 18: return 7293*exp(-1/alpha)/163840-8619*exp(-4/alpha)/143360-2431*exp(-9/alpha)/245760-1037*exp(-16/alpha)/450560-8245*exp(-25/alpha)/16400384-799*exp(-36/alpha)/8785920.0-579*exp(-49/alpha)/46858240.0-759*exp(-64/alpha)/687964160.0-107*exp(-81/alpha)/2222653440.0+12155.0/98304.0 ;
    break;
  case 20: return 4199*exp(-1/alpha)/98304-54587*exp(-4/alpha)/917504-3553*exp(-9/alpha)/344064-19703*exp(-16/alpha)/7569408-31331*exp(-25/alpha)/49201152-15181*exp(-36/alpha)/112459776.0-3667*exp(-49/alpha)/159318016.0-253*exp(-64/alpha)/86900736.0-107*exp(-81/alpha)/444530688.0-397*exp(-100/alpha)/40896823296.0+46189.0/393216.0 ;
    break;
  }
  return 0.0;
}

//for n:0 thru 20 do print("n=",n,": ",string(expand(integrate(integrate(sqrt(alpha/%pi)*exp(-alpha*(t)^2)*cos(t+d)^n*sin(2*d)*sin(d), t, -inf, inf), d, 0, %pi))));
double int_convoluted_2_sin_sin_cos_n(int n, double sigma)
{
  assert(n >= 0);
  assert(n <= 20);
  double pi = TMath::Pi();
  double alpha = 1.0/(2.0*sigma*sigma);
  if (n % 2 == 0)
    return 0.0;
  switch (n) {
  case 1: return pi*exp(-1/(4*alpha))/4 ;
    break;
  case 3: return 3*pi*exp(-1/(4*alpha))/16-pi*exp(-9/(4*alpha))/16 ;
    break;
  case 5: return 5*pi*exp(-1/(4*alpha))/32-5*pi*exp(-9/(4*alpha))/64 ;
    break;
  case 7: return 35*pi*exp(-1/(4*alpha))/256-21*pi*exp(-9/(4*alpha))/256 ;
    break;
  case 9: return 63*pi*exp(-1/(4*alpha))/512-21*pi*exp(-9/(4*alpha))/256 ;
    break;
  case 11: return 231*pi*exp(-1/(4*alpha))/2048-165*pi*exp(-9/(4*alpha))/2048 ;
    break;
  case 13: return 429*pi*exp(-1/(4*alpha))/4096-1287*pi*exp(-9/(4*alpha))/16384 ;
    break;
  case 15: return 6435*pi*exp(-1/(4*alpha))/65536-5005*pi*exp(-9/(4*alpha))/65536 ;
    break;
  case 17: return 12155*pi*exp(-1/(4*alpha))/131072-2431*pi*exp(-9/(4*alpha))/32768 ;
    break;
  case 19: return 46189*pi*exp(-1/(4*alpha))/524288-37791*pi*exp(-9/(4*alpha))/524288 ;
    break;
  }
  return 0.0;
}

//for n:0 thru 20 do print("n=",n,": ",string(expand(integrate(integrate(sqrt(alpha/%pi)*exp(-alpha*(t)^2)*cos(t+d)^n*sin(d)*cos(d), t, -inf, inf), d, 0, %pi))));
double int_convoluted_cos_sin_cos_n(int n, double sigma)
{
  assert(n >= 0);
  assert(n <= 20);
  double alpha = 1.0/(2.0*sigma*sigma);
  if (n % 2 == 0)
    return 0.0;
  switch (n) {
  case 1: return 2*exp(-1/(4*alpha))/3 ;
    break;
  case 3: return exp(-1/(4*alpha))/2-exp(-9/(4*alpha))/10 ;
    break;
  case 5: return 5*exp(-1/(4*alpha))/12-exp(-9/(4*alpha))/8-exp(-25/(4*alpha))/168 ;
    break;
  case 7: return 35*exp(-1/(4*alpha))/96-21*exp(-9/(4*alpha))/160-exp(-25/(4*alpha))/96-exp(-49/(4*alpha))/1440 ;
    break;
  case 9: return 21*exp(-1/(4*alpha))/64-21*exp(-9/(4*alpha))/160-3*exp(-25/(4*alpha))/224-exp(-49/(4*alpha))/640-exp(-81/(4*alpha))/9856 ;
    break;
  case 11: return 77*exp(-1/(4*alpha))/256-33*exp(-9/(4*alpha))/256-55*exp(-25/(4*alpha))/3584-11*exp(-49/(4*alpha))/4608-exp(-81/(4*alpha))/3584-exp(-121/(4*alpha))/59904 ;
    break;
  case 13: return 143*exp(-1/(4*alpha))/512-1287*exp(-9/(4*alpha))/10240-715*exp(-25/(4*alpha))/43008-143*exp(-49/(4*alpha))/46080-39*exp(-81/(4*alpha))/78848-exp(-121/(4*alpha))/18432-exp(-169/(4*alpha))/337920 ;
    break;
  case 15: return 2145*exp(-1/(4*alpha))/8192-1001*exp(-9/(4*alpha))/8192-143*exp(-25/(4*alpha))/8192-91*exp(-49/(4*alpha))/24576-65*exp(-81/(4*alpha))/90112-35*exp(-121/(4*alpha))/319488-exp(-169/(4*alpha))/90112-exp(-225/(4*alpha))/1810432 ;
    break;
  case 17: return 12155*exp(-1/(4*alpha))/49152-2431*exp(-9/(4*alpha))/20480-221*exp(-25/(4*alpha))/12288-1547*exp(-49/(4*alpha))/368640-85*exp(-81/(4*alpha))/90112-85*exp(-121/(4*alpha))/479232-17*exp(-169/(4*alpha))/675840-exp(-225/(4*alpha))/425984-exp(-289/(4*alpha))/9338880 ;
    break;
  case 19: return 46189*exp(-1/(4*alpha))/196608-37791*exp(-9/(4*alpha))/327680-4199*exp(-25/(4*alpha))/229376-2261*exp(-49/(4*alpha))/491520-2907*exp(-81/(4*alpha))/2523136-323*exp(-121/(4*alpha))/1277952-323*exp(-169/(4*alpha))/7208960-171*exp(-225/(4*alpha))/28966912-exp(-289/(4*alpha))/1966080-exp(-361/(4*alpha))/46792704 ;
    break;
  }
  return 0.0;
}

//for n:0 thru 20 do print("n=",n,": ",string(expand(integrate(integrate(sqrt(alpha/%pi)*exp(-alpha*(t)^2)*cos(t+d)^n*sin(d)^2, t, -inf, inf), d, 0, %pi))));
double int_convoluted_sin_2_cos_n(int n, double sigma)
{
  assert(n >= 0);
  assert(n <= 20);
  double pi = TMath::Pi();
  double alpha = 1.0/(2.0*sigma*sigma);
  if (n % 2 == 1)
    return 0.0;
  switch (n) {
  case 0: return pi/2 ;
    break;
  case 2: return pi/4-pi*exp(-1/alpha)/8 ;
    break;
  case 4: return 3*pi/16-pi*exp(-1/alpha)/8 ;
    break;
  case 6: return 5*pi/32-15*pi*exp(-1/alpha)/128 ;
    break;
  case 8: return 35*pi/256-7*pi*exp(-1/alpha)/64 ;
    break;
  case 10: return 63*pi/512-105*pi*exp(-1/alpha)/1024 ;
    break;
  case 12: return 231*pi/2048-99*pi*exp(-1/alpha)/1024 ;
    break;
  case 14: return 429*pi/4096-3003*pi*exp(-1/alpha)/32768 ;
    break;
  case 16: return 6435*pi/65536-715*pi*exp(-1/alpha)/8192 ;
    break;
  case 18: return 12155*pi/131072-21879*pi*exp(-1/alpha)/262144 ;
    break;
  case 20: return 46189*pi/524288-20995*pi*exp(-1/alpha)/262144 ;
    break;
  }
  return 0.0;
}

//for n:0 thru 20 do print("case ",n,": return ",string(float(expand(integrate(integrate(sqrt(alpha/%pi)*(phi+d)^n*exp(-alpha*phi^2), phi, -inf, inf), d, -%pi, %pi)))), "; break;");
double int_inf_convoluted_const_x_n(int n, double sigma)
{
  assert(n >= 0);
  assert(n <= 20);
  double alpha = 1.0/(2.0*sigma*sigma);
  if (n % 2 == 1)
    return 0.0;
  switch (n) {
  case  0 : return  6.283185307179586 ; break; 
  case  2 : return  3.141592653589793/alpha+20.67085112019988 ; break; 
  case  4 : return  62.01255336059963/alpha+4.71238898038469/pow(alpha,2)+122.4078739141126 ; break; 
  case  6 : return  918.0590543558442/alpha+232.5470751022486/pow(alpha,2)+11.78097245096172/pow(alpha,3)+862.9409222219404 ; break; 
  case  8 : return  12081.17291110717/alpha+6426.41338049091/pow(alpha,2)+1085.219683810494/pow(alpha,3)+41.23340357836604/pow(alpha,4)+6624.244296321378 ; break; 
  case  10 : return  149045.496667231/alpha+135913.1952499556/pow(alpha,2)+48198.10035368182/pow(alpha,3)+6104.360721434026/pow(alpha,4)+185.5503161026471/pow(alpha,5)+53491.63963161646 ; break; 
  case  12 : return  
      1765224.107843343/alpha+2459250.695009312/pow(alpha,2)+1495045.147749512/pow(alpha,3)+397634.327917875/pow(alpha,4)+40288.78076146457/pow(alpha,5)+1020.526738564559/pow(alpha,6)+446719.5800943512 
      ; break; 
  case  14 : return  2.0325740894292977E+7/alpha+4.0158848453436054E+7/pow(alpha,2)+3.7298635540974565E+7/pow(alpha,3)+1.7006138555650696E+7/pow(alpha,4)+3618472.384052663/pow(alpha,5)+305523.254107773/pow(alpha,6)+6633.423800669636/pow(alpha,7)+3821086.129251732 ; break; 
  case  16 : return  2.2926516775510389E+8/alpha+6.0977222682878935E+8/pow(alpha,2)+8.0317696906872118E+8/pow(alpha,3)+5.5947953311461842E+8/pow(alpha,4)+2.0407366266780835E+8/pow(alpha,5)+3.6184723840526626E+7/pow(alpha,6)+2618770.749495197/pow(alpha,7)+49750.67850502227/pow(alpha,8)+3.3275831010180339E+7 ; break; 
  case  18 : return  2.5456010722787962E+9/alpha+8.7693926666327229E+9/pow(alpha,2)+1.5549191784134127E+10/pow(alpha,3)+1.5360759533439291E+10/pow(alpha,4)+8.5600368566536617E+9/pow(alpha,5)+2.6019391990145564E+9/pow(alpha,6)+3.9544733911432672E+8/pow(alpha,7)+2.5041995292047825E+7/pow(alpha,8)+422880.7672926893/pow(alpha,9)+2.9384883679977304E+8 ; break; 
  case  20 : return  2.7915639495978439E+10/alpha+1.2091605093324281E+11/pow(alpha,2)+2.7769743444336957E+11/pow(alpha,3)+3.6929330487318555E+11/pow(alpha,4)+2.9185443113534656E+11/pow(alpha,5)+1.3553391689701631E+11/pow(alpha,6)+3.5312031986626122E+10/pow(alpha,7)+4.6959371519826298E+9/pow(alpha,8)+2.6433217252717146E+8/pow(alpha,9)+4017367.289280548/pow(alpha,10)+2.623964937416502E+9 
      ; break; 
  }
  return 0.0;
}

//for n:0 thru 20 do print("case ",n,": return ",string(float(expand(integrate(integrate(sqrt(alpha/%pi)*(phi+d)^n*exp(-alpha*phi^2)*cos(d)^2, phi, -inf, inf), d, -%pi, %pi)))), "; break;");
double int_inf_convoluted_cos_2_x_n(int n, double sigma)
{
  assert(n >= 0);
  assert(n <= 20);
  double alpha = 1.0/(2.0*sigma*sigma);
  if (n % 2 == 1)
    return 0.0;
  switch (n) {
  case  0 : return  3.141592653589793 ; break; 
  case  2 : return  1.570796326794897/alpha+11.90622188689484 ; break; 
  case  4 : return  35.71866566068451/alpha+2.356194490192345/pow(alpha,2)+87.4978246569714 ; break; 
  case  6 : return  656.2336849272856/alpha+133.9449962275669/pow(alpha,2)+5.890486225480862/pow(alpha,3)+693.2958305395289 ; break; 
  case  8 : return  9706.141627553405/alpha+4593.635794490999/pow(alpha,2)+625.0766490619789/pow(alpha,3)+20.61670178918302/pow(alpha,4)+5687.153431714451 ; break; 
  case  10 : return  127960.9522135751/alpha+109194.0933099758/pow(alpha,2)+34452.26845868249/pow(alpha,3)+3516.056150973631/pow(alpha,4)+92.77515805132357/pow(alpha,5)+47830.36426946412 ; break;
  case  12 : return  1578402.020892315/alpha+2111355.71152399/pow(alpha,2)+1201135.026409734/pow(alpha,3)+284231.2147841306/pow(alpha,4)+23205.97059642596/pow(alpha,5)+510.2633692822797/pow(alpha,6)+410181.8769982042 
      ; break; 
  case  14 : return  1.8663275403418235E+7/alpha+3.5908645975300178E+7/pow(alpha,2)+3.2022228291447181E+7/pow(alpha,3)+1.3662910925410721E+7/pow(alpha,4)+2586504.054535588/pow(alpha,5)+175978.6103562302/pow(alpha,6)+3316.711900334818/pow(alpha,7)+3573008.555500609 ; break; 
  case  16 : return  2.1438051333003667E+8/alpha+5.5989826210254693E+8/pow(alpha,2)+7.1817291950600362E+8/pow(alpha,3)+4.8033342437170768E+8/pow(alpha,4)+1.6395493110492867E+8/pow(alpha,5)+2.5865040545355879E+7/pow(alpha,6)+1508388.088767688/pow(alpha,7)+24875.33925251113/pow(alpha,8)+3.1522569930157386E+7 ; break; 
  case  18 : return  2.4114765996571231E+9/alpha+8.2000546348738604E+9/pow(alpha,2)+1.4277405683614956E+10/pow(alpha,3)+1.3735057085552319E+10/pow(alpha,4)+7.3491013928871279E+9/pow(alpha,5)+2.0904253715878406E+9/pow(alpha,6)+2.8266794310281783E+8/pow(alpha,7)+1.4423961098841013E+7/pow(alpha,8)+211440.3836463447/pow(alpha,9)+2.810488910215596E+8 ; break; 
  case  20 : return  2.6699644647049049E+10/alpha+1.1454513848371289E+11/pow(alpha,2)+2.5966839677100656E+11/pow(alpha,3)+3.390883849858551E+11/pow(alpha,4)+2.6096608462549408E+11/pow(alpha,5)+1.1636077205404617E+11/pow(alpha,6)+2.8370058614406403E+10/pow(alpha,7)+3.3566818243459616E+9/pow(alpha,8)+1.5225292270998847E+8/pow(alpha,9)+2008683.644640274/pow(alpha,10)+2.527977317637641E+9 
      ; break; 
  }
  return 0.0;
}

//for n:0 thru 20 do print("case ",n,": return ",string(float(expand(integrate(integrate(sqrt(alpha/%pi)*(phi+d)^n*exp(-alpha*phi^2)*sin(d)^2, phi, -inf, inf), d, -%pi, %pi)))), "; break;");
double int_inf_convoluted_sin_2_x_n(int n, double sigma)
{
  assert(n >= 0);
  assert(n <= 20);
  double alpha = 1.0/(2.0*sigma*sigma);
  if (n % 2 == 1)
    return 0.0;
  switch (n) {
  case  0 : return  3.141592653589793 ; break; 
  case  2 : return  1.570796326794897/alpha+8.764629233305042 ; break; 
  case  4 : return  26.29388769991513/alpha+2.356194490192345/pow(alpha,2)+34.91004925714116 ; break; 
  case  6 : return  261.8253694285586/alpha+98.60207887468172/pow(alpha,2)+5.890486225480862/pow(alpha,3)+169.6450916824115 ; break; 
  case  8 : return  2375.031283553761/alpha+1832.777585999911/pow(alpha,2)+460.1430347485147/pow(alpha,3)+20.61670178918302/pow(alpha,4)+937.0908646069279 ; break; 
  case  10 : return  21084.54445365589/alpha+26719.10193997981/pow(alpha,2)+13745.83189499933/pow(alpha,3)+2588.304570460395/pow(alpha,4)+92.77515805132357/pow(alpha,5)+5661.275362152337 ; break;                               
  case  12 : return  
      186822.0869510286/alpha+347894.9834853215/pow(alpha,2)+293910.121339778/pow(alpha,3)+113403.1131337445/pow(alpha,4)+17082.81016503861/pow(alpha,5)+510.2633692822797/pow(alpha,6)+36537.70309614705 
      ; break; 
  case  14 : return  1662465.490874743/alpha+4250202.478135873/pow(alpha,2)+5276407.249527384/pow(alpha,3)+3343227.630239975/pow(alpha,4)+1031968.329517075/pow(alpha,5)+129544.6437515428/pow(alpha,6)+3316.711900334818/pow(alpha,7)+248077.5737511227 ; break; 
  case  16 : return  1.4884654425067216E+7/alpha+4.9873964726242363E+7/pow(alpha,2)+8.5004049562717617E+7/pow(alpha,3)+7.9146108742910743E+7/pow(alpha,4)+4.0118731562879689E+7/pow(alpha,5)+1.0319683295170747E+7/pow(alpha,6)+1110382.66072751/pow(alpha,7)+24875.33925251113/pow(alpha,8)+1753261.080022954 ; break; 
  case  18 : return  1.3412447262167311E+8/alpha+5.693380317588625E+8/pow(alpha,2)+1.2717861005191698E+9/pow(alpha,3)+1.6257024478869734E+9/pow(alpha,4)+1.2109354637665339E+9/pow(alpha,5)+5.1151382742671597E+8/pow(alpha,6)+1.1277939601150887E+8/pow(alpha,7)+1.0618034193206811E+7/pow(alpha,8)+211440.3836463447/pow(alpha,9)+1.2799945778213412E+7 ; break; 
  case  20 : return  1.21599484892939E+9/alpha+6.3709124495299149E+9/pow(alpha,2)+1.8029037672363007E+10/pow(alpha,3)+3.0204919887330444E+10/pow(alpha,4)+3.0888346509852478E+10/pow(alpha,5)+1.9173144842970131E+10/pow(alpha,6)+6.941973372219717E+9/pow(alpha,7)+1.339255327636668E+9/pow(alpha,8)+1.1207924981718299E+8/pow(alpha,9)+2008683.644640274/pow(alpha,10)+9.5987619778861046E+7 ; break; 
  }
  return 0.0;
}

//for n:0 thru 20 do print("case ",n,": return ",string(float(expand(integrate(integrate(sqrt(alpha/%pi)*(phi+d)^n*exp(-alpha*phi^2)*sin(d), phi, -inf, inf), d, -%pi, %pi)))), "; break;");
double int_inf_convoluted_sin_x_n(int n, double sigma)
{
  assert(n >= 0);
  assert(n <= 20);
  double alpha = 1.0/(2.0*sigma*sigma);
  if (n % 2 == 0)
    return 0.0;
  switch (n) {
  case  1 : return  6.283185307179586 ; break; 
  case  3 : return  9.424777960769379/alpha+24.31344151752212 ; break; 
  case  5 : return  121.5672075876106/alpha+23.56194490192345/pow(alpha,2)+125.7705392201204 ; break; 
  case  7 : return  1320.590661811266/alpha+638.2278398349555/pow(alpha,2)+82.46680715673207/pow(alpha,3)+758.2238083085194 ; break; 
  case  9 : return  13648.02854955342/alpha+11885.31595630138/pow(alpha,2)+3829.367039009733/pow(alpha,3)+371.1006322052943/pow(alpha,4)+5026.084468678731 ; break; 
  case  11 : return  138217.3228886668/alpha+187660.3925563593/pow(alpha,2)+108948.7295994294/pow(alpha,3)+26326.89839319191/pow(alpha,4)+2041.053477129119/pow(alpha,5)+35538.744393114 ; break; 
  case  13 : return  
      1386011.031330336/alpha+2695237.796329141/pow(alpha,2)+2439585.103232659/pow(alpha,3)+1062250.113594437/pow(alpha,4)+205349.807466897/pow(alpha,5)+13266.84760133927/pow(alpha,6)+263310.4159052214 
      ; break; 
  case  15 : return  1.3823796834663689E+7/alpha+3.6382789572466373E+7/pow(alpha,2)+4.7166661435754895E+7/pow(alpha,3)+3.201955447992897E+7/pow(alpha,4)+1.1153626192741588E+7/pow(alpha,5)+1796810.815335348/pow(alpha,6)+99501.35701004454/pow(alpha,7)+2021104.600121215 ; break; 
  case  17 : return  1.3743511271761036E+8/alpha+4.7000909238989449E+8/pow(alpha,2)+8.2467656364292908E+8/pow(alpha,3)+8.018332444078331E+8/pow(alpha,4)+4.3546594092703295E+8/pow(alpha,5)+1.2640776351773798E+8/pow(alpha,6)+1.7454733634686239E+7/pow(alpha,7)+845761.5345853786/pow(alpha,8)+1.5948676302624345E+7 ; break; 
  case  19 : return  1.363611830458992E+9/alpha+5.8753510678547668E+9/pow(alpha,2)+1.3395259133206604E+10/pow(alpha,3)+1.7627461547833862E+10/pow(alpha,4)+1.3711348479374664E+10/pow(alpha,5)+6.2053896582101974E+9/pow(alpha,6)+1.5439805401095133E+9/pow(alpha,7)+1.865474657207092E+8/pow(alpha,8)+8034734.578561097/pow(alpha,9)+1.2868057735972023E+8 ; break; 
  }
  return 0.0;
}

//for n:0 thru 20 do print("case ",n,": return ",string(float(expand(integrate(integrate(sqrt(alpha/%pi)*(phi+d)^n*exp(-alpha*phi^2)*sin(2*d), phi, -inf, inf), d, -%pi, %pi)))), "; break;");
double int_inf_convoluted_2_sin_x_n(int n, double sigma)
{
  assert(n >= 0);
  assert(n <= 20);
  double alpha = 1.0/(2.0*sigma*sigma);
  if (n % 2 == 0)
    return 0.0;
  switch (n) {
  case  1 : return  -3.141592653589793 ; break; 
  case  3 : return  -4.71238898038469/alpha-26.29388769991513 ; break; 
  case  5 : return  -131.4694384995756/alpha-11.78097245096172/pow(alpha,2)-174.5502462857058 ; break; 
  case  7 : return  -1832.777585999911/alpha-690.214552122772/pow(alpha,2)-41.23340357836604/pow(alpha,3)-1187.515641776881 ; break; 
  case  9 : return  -21375.28155198386/alpha-16494.99827399919/pow(alpha,2)-4141.287312736632/pow(alpha,3)-185.5503161026471/pow(alpha,4)-8433.817781462349 ; break; 
  case  11 : return  -231929.9889902145/alpha-293910.121339778/pow(alpha,2)-151204.1508449926/pow(alpha,3)-28471.35027506435/pow(alpha,4)-1020.526738564559/pow(alpha,5)-62274.02898367599 ; break;                                  
  case  13 : return  -2428687.130363355/alpha-4522634.785309188/pow(alpha,2)-3820831.577417114/pow(alpha,3)-1474240.470738678/pow(alpha,4)-222076.5321455019/pow(alpha,5)-6633.423800669636/pow(alpha,6)-474990.140249928 
      ; break; 
  case  15 : return  -2.4936982363121182E+7/alpha-6.3753037172038078E+7/pow(alpha,2)-7.9146108742910743E+7/pow(alpha,3)-5.0148414453599617E+7/pow(alpha,4)-1.5479524942756122E+7/pow(alpha,5)-1943169.656273142/pow(alpha,6)-49750.67850502227/pow(alpha,7)-3721163.606266804 ; break; 
  case  17 : return  -2.5303912522615099E+8/alpha-8.4785740034611607E+8/pow(alpha,2)-1.4450688425662041E+9/pow(alpha,3)-1.3454838486294813E+9/pow(alpha,4)-6.8201843656895471E+8/pow(alpha,5)-1.7543461601790273E+8/pow(alpha,6)-1.8876505232367665E+7/pow(alpha,7)-422880.7672926893/pow(alpha,8)-2.9805438360381901E+7 ; break; 
  case  19 : return  -2.5483649798125954E+9/alpha-1.0817422603417984E+10/pow(alpha,2)-2.4163935909864288E+10/pow(alpha,3)-3.0888346509852478E+10/pow(alpha,4)-2.3007773811564163E+10/pow(alpha,5)-9.7187627211076069E+9/pow(alpha,6)-2.1428085242186687E+9/pow(alpha,7)-2.017426496709294E+8/pow(alpha,8)-4017367.289280548/pow(alpha,9)-2.4319896978524876E+8 ; break; 
  }
  return 0.0;
}

//for n:0 thru 20 do print("case ",n,": return ",string(float(expand(integrate(integrate(sqrt(alpha/%pi)*(phi+d)^n*exp(-alpha*phi^2)*cos(d), phi, -inf, inf), d, -%pi, %pi)))), "; break;");
double int_inf_convoluted_cos_x_n(int n, double sigma)
{
  assert(n >= 0);
  assert(n <= 20);
  double alpha = 1.0/(2.0*sigma*sigma);
  if (n % 2 == 1)
    return 0.0;
  switch (n) {
  case  0 : return  0.0 ; break; 
  case  2 : return  -12.56637061435917 ; break; 
  case  4 : return  -37.69911184307752/alpha-97.25376607008846 ; break; 
  case  6 : return  -729.4032455256634/alpha-141.3716694115407/pow(alpha,2)-754.6232353207233 ; break; 
  case  8 : return  -10564.72529449013/alpha-5105.822718679644/pow(alpha,2)-659.7344572538566/pow(alpha,3)-6065.790466468156 ; break; 
  case  10 : return  -136480.2854955334/alpha-118853.1595630139/pow(alpha,2)-38293.67039009733/pow(alpha,3)-3711.006322052943/pow(alpha,4)-50260.8446867906 ; break; 
  case  12 : return  -1658607.874663904/alpha-2251924.710676324/pow(alpha,2)-1307384.755193153/pow(alpha,3)-315922.780718303/pow(alpha,4)-24492.64172554942/pow(alpha,5)-426464.9327177554 ; break;
  case  14 : return  -1.9404154438628018E+7/alpha-3.7733329148607552E+7/pow(alpha,2)-3.4154191445257187E+7/pow(alpha,3)-1.4871501590322122E+7/pow(alpha,4)-2874897.304536558/pow(alpha,5)-185735.8664187498/pow(alpha,6)-3686345.822659835 ; break; 
  case  16 : return  -2.2118074935461903E+8/alpha-5.8212463315946198E+8/pow(alpha,2)-7.5466658297207832E+8/pow(alpha,3)-5.1231287167886353E+8/pow(alpha,4)-1.784580190838654E+8/pow(alpha,5)-2.8748973045365565E+7/pow(alpha,6)-1592021.712160713/pow(alpha,7)-3.233767360193944E+7 ; break; 
  case  18 : return  -2.4738320296894684E+9/alpha-8.4601636629215393E+9/pow(alpha,2)-1.4844178145556E+10/pow(alpha,3)-1.4432998399342041E+10/pow(alpha,4)-7.8383869366865463E+9/pow(alpha,5)-2.2753397433192844E+9/pow(alpha,6)-3.1418520542435235E+8/pow(alpha,7)-1.5223707622536814E+7/pow(alpha,8)-2.8707617035731125E+8 ; break; 
  case  20 : return  -2.7272236397462036E+10/alpha-1.1750702138356006E+11/pow(alpha,2)-2.6790518266332812E+11/pow(alpha,3)-3.5254923095672754E+11/pow(alpha,4)-2.742269695875127E+11/pow(alpha,5)-1.2410779316420337E+11/pow(alpha,6)-3.0879610802190308E+10/pow(alpha,7)-3.7309493144141827E+9/pow(alpha,8)-1.6069469157122192E+8/pow(alpha,9)-2.5736123940656128E+9 ; break; 
  }
  return 0.0;
}

//void fitter::ct_tagging_part(measurement* meas, double& one, double& two)
//this function calculates the convolution of the product exp(-t/tau) * sin(deltam*t) and exp(-t/tau) * cos(deltam*t) 
void convoluted_exp_sincos(double ct, double sigma, double tau, double deltam, double& one, double& two)
{
  double sigma_2 = sigma*sigma;
  double tau_2 = tau*tau;
  double deltam_2 = deltam*deltam;
  const double sqrt2 = TMath::Sqrt(2.0);

  const bool cdf_way = false;
  if (!cdf_way)
    {
      //safety for large negative decay times
      if (ct < - 6.0*sigma)
	{
	  one = 0.0;
	  two = 0.0;
	  return;
	}
      double c1 = 0.5;
      //      double exp1arg = 0.5*sigma_2/tau_2 - 0.5*deltam_2*sigma_2 - ct/tau;
      double exp1arg = 0.5*sigma_2*(1.0/tau_2 - deltam_2) - ct/tau;
      double exp1 = exp(exp1arg);
      double exp2arg = -deltam*(ct - sigma_2/tau);
      complex<double> exp2(cos(exp2arg), sin(exp2arg)); 
      complex<double> cerfarg(sigma/(tau*sqrt2) - ct/(sigma*sqrt2), +deltam*sigma/sqrt2);
      complex<double> cerf;
      if  (cerfarg.real() < -20.0)
	cerf = complex<double>(2.0,0.0);
      else
	cerf = cErrF_2(cerfarg);//best complex error function
      complex<double> c2(exp2*cerf);
      double im = -c2.imag();//exp*sin
      double re = +c2.real();//exp*cos
      one = c1*exp1*im ;
      two = c1*exp1*re ;
    }
  else //all three methods are analytically the same. There might be numerical differences, though
    {
      if (ct < -6.0*sigma)
	{
	  one = 0.0;
	  two = 0.0;
	  return;
	}

      //current cdf
      std::complex<double> z(deltam*sigma/sqrt2, (sigma/tau-ct/sigma)/sqrt2);
      if (ct<0) {//i do not quite get this part. the calculation should also be corerct for ct < 0, right?
	one= 2.0*nwwerf(z).real()/4.0*exp(-ct*ct/2.0/sigma/sigma);
	two= 2.0*nwwerf(z).imag()/4.0*exp(-ct*ct/2.0/sigma/sigma);
      }
      else {
	one= -2.0*nwwerf(std::conj(z)).real()/tau/4*exp(-ct*ct/2.0/sigma/sigma) +
	  exp(sigma*sigma/2 *(1/tau/tau - deltam*deltam) - ct/tau)*cos(deltam*ct - deltam/tau*sigma*sigma);
	two= +2.0*nwwerf(std::conj(z)).imag()/tau/4*exp(-ct*ct/2.0/sigma/sigma) +
	  exp(sigma*sigma/2 *(1/tau/tau - deltam*deltam) - ct/tau)*sin(deltam*ct - deltam/tau*sigma*sigma);
      } 

    }

  if (std::isnan(one) || std::isnan(two))
    {
      //void convoluted_exp_sincos(double ct, double sigma, double tau, double deltam, double& one, double& two)
      cout << "ct= " << ct << " sigma=" << sigma << " tau=" << tau 
	   << " deltam=" << deltam << " one=" << one << " two=" << two << endl;
      cout << "Error: NaN" << endl;
      assert(0);
      
    }
  return;
};


double integrate_convoluted_exp_cos(double a, double b, double 
sigma, double tau, double deltam)
{
double sigma_2 = sigma*sigma;
double tau_2 = tau*tau;

complex<double> i(0.0, 1.0);
//std::cout<<setprecision(8)<<i.imag()<<std::endl;
//double prefexp = exp(sigma*(2.0+sigma*(-1.0+tau_2*deltam*deltam))/(2.0*tau_2));
complex<double> pref(tau/(-1.0-i*tau*deltam));
complex<double> erfarga (a/(TMath::Sqrt(2)*sigma), 0.0);
complex<double> erfcarga ((-tau*a + sigma_2*(1.0 + 
i*tau*deltam))/(TMath::Sqrt(2)*sigma*tau));
complex<double> exparga(-(tau*deltam - i)*(sigma_2*(tau*deltam - i) + 
2.0*i*tau*a)/(2.0*tau_2));

//std::cout<<setprecision(8)<<prefa.imag()<<std::endl;

complex<double> erfargb  (b/(TMath::Sqrt(2)*sigma), 0.0);
complex<double> erfcargb ((-tau*b + sigma_2*(1.0 + 
i*tau*deltam))/(TMath::Sqrt(2)*sigma*tau));
complex<double> expargb (-(tau*deltam - i)*(sigma_2*(tau*deltam - i) + 
2.0*i*tau*b)/(2.0*tau_2));

//std::cout<<setprecision(8)<<prefb.real()<<std::endl;

complex<double> erfa = ErrF_2(erfarga);
complex<double> erfca = cErrF_2(erfcarga);
complex<double> expa = exp(exparga);

complex<double> erfb = ErrF_2(erfargb);
complex<double> erfcb = cErrF_2(erfcargb);
complex<double> expb = exp(expargb);


complex <double> za = pref * (- erfa + erfca * expa);
complex <double> zb = pref* (- erfb + erfcb * expb);

double two = 0.5*(zb.real() - za.real());

 if (std::isnan(two))
    {
      //void convoluted_exp_sincos(double ct, double sigma, double tau, double deltam, double& one, double& two)
      cout << " sigma=" << sigma << " tau=" << tau
           << " deltam=" << deltam << " two=" << two << endl;
      cout << "Error: NaN" << endl;
      assert(0);

    }
  return two;

}



