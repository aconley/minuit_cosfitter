#include <algorithm>
#include <iterator>
#include <sstream>

#include <utility.h>
#include <cosfitterexcept.h>

using namespace std;

/*!
  \brief Unpacks a parameter vector based on the type of fit being done
*/
void utility::parse_param_values( const std::vector<double>& par,
				  utility::cosmo_fittype fittype,
				  double& om, double& ode, double& w,
				  double& wa, double& alpha, double& beta,
				  double& scriptm1, double& scriptm2) {

  //Figure out how many arguments we expect
  unsigned int argexpected = 0;
  switch (fittype) {
  case utility::omegam_omegade : argexpected = 6; break;
  case utility::flat_omegam : argexpected = 5; break;
  case utility::flat_omegam_w : argexpected = 6; break;
  case utility::omegam_omegade_w : argexpected = 7; break;
  case utility::flat_omegam_w0_wa : argexpected = 7; break;
  default : throw CosFitterExcept("utility","parse_param_values",
				  "Unknown fittype",1);
    break;
  }
  if (par.size() != argexpected) {
    std::stringstream errstrng;
    errstrng << "Expected " << argexpected << " params, got: "
	     << par.size();
    throw CosFitterExcept("utility","parse_param_values",
			  errstrng.str(),2);
  }

  //Set up lumdist
  w=-1; wa=0.0; om=1.0; ode=0.0; alpha=1.5; beta=1.5; 
  scriptm1=23.5; scriptm2=23.5;//Defaults
  switch (fittype) {
  case omegam_omegade :
    om = par[0];
    ode = par[1];
    alpha = par[2];
    beta = par[3];
    scriptm1 = par[4];
    scriptm2 = par[5];
    break;
  case flat_omegam :
    om = par[0];
    ode = 1.0-om;
    alpha = par[1];
    beta = par[2];
    scriptm1 = par[3];
    scriptm2 = par[4];
    break;
  case flat_omegam_w :
    om = par[0];
    ode = 1.0 - om;
    w = par[1];
    alpha = par[2];
    beta = par[3];
    scriptm1 = par[4];
    scriptm2 = par[5];
    break;
  case omegam_omegade_w :
    om = par[0];
    ode = par[1];
    w = par[2];
    alpha = par[3];
    beta = par[4];
    scriptm1 = par[5];
    scriptm2 = par[6];
    break;
  case flat_omegam_w0_wa :
    om = par[0];
    ode = 1.0 - om;
    w = par[1];
    wa = par[2];
    alpha = par[3];
    beta = par[4];
    scriptm1 = par[5];
    scriptm2 = par[6];
    break;
  default :
    throw CosFitterExcept("utility","parse_param_values",
			  "Unknown fittype",4);
  }
}

/*!
  \brief Unpacks a parameter vector based on the type of fit being done,]
   no \f$\mathcal{M}\f$ version
*/
void utility::parse_param_values_nosm( const std::vector<double>& par,
				       utility::cosmo_fittype fittype,
				       double& om, double& ode, double& w,
				       double& wa, double& alpha, 
				       double& beta) {

  //Figure out how many arguments we expect
  unsigned int argexpected = 0;
  switch (fittype) {
  case utility::omegam_omegade : argexpected = 4; break;
  case utility::flat_omegam : argexpected = 3; break;
  case utility::flat_omegam_w : argexpected = 4; break;
  case utility::omegam_omegade_w : argexpected = 5; break;
  case utility::flat_omegam_w0_wa : argexpected = 5; break;
  default : throw CosFitterExcept("utility","parse_param_values_nosm",
				  "Unknown fittype",1);
    break;
  }
  if (par.size() < argexpected) {
    std::stringstream errstrng;
    errstrng << "Expected " << argexpected << " params, got: "
	     << par.size();
    throw CosFitterExcept("utility","parse_param_values_nosm",
			  errstrng.str(),2);
  }

  //Set up lumdist
  w=-1; wa=0.0; om=1.0; ode=0.0; alpha=1.5; beta=1.5;//Defaults
  switch (fittype) {
  case omegam_omegade :
    om = par[0];
    ode = par[1];
    alpha = par[2];
    beta = par[3];
    break;
  case flat_omegam :
    om = par[0];
    ode = 1.0-om;
    alpha = par[1];
    beta = par[2];
    break;
  case flat_omegam_w :
    om = par[0];
    ode = 1.0 - om;
    w = par[1];
    alpha = par[2];
    beta = par[3];
    break;
  case omegam_omegade_w :
    om = par[0];
    ode = par[1];
    w = par[2];
    alpha = par[3];
    beta = par[4];
    break;
  case flat_omegam_w0_wa :
    om = par[0];
    w = par[1];
    wa = par[2];
    alpha = par[3];
    beta = par[4];
    break;
  default :
    throw CosFitterExcept("utility","parse_param_values_nosm",
			  "Unknown fittype",4);
  }
}

/*! 
  Breaks an input string up into a vector of string, which correspond
  to the input string split on spaces.  Cheerfully stolen from Rob
  Knop.
*/
void utility::stringwords(const std::string &ins,
			  std::vector<std::string> &words) {
  string s,tmp;
  unsigned int i,p;
  int first,last;

  s = ins;

  // Trim spaces from beginning and end

  first=s.find_first_not_of(" ");
  if (first==-1) {
    s="";
  } else {
    last=s.find_last_not_of(" ");
    s=s.substr(first,last-first+1);
  }

  words.clear();

  p=s.find_first_not_of(" \t\r\n");
  if (p>=s.size()) return;
  while ((i=s.find_first_of(" \t\r\n",p))) {
    tmp=s.substr(p,i-p);
    words.push_back(tmp);
    p=s.find_first_not_of(" \t\r\n",i);
    if (p>=s.size()) return;
  }
  tmp=s.substr(p);
  words.push_back(tmp);
}
