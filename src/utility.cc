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
				  double& scriptm) {

  //Figure out how many arguments we expect
  unsigned int argexpected = 0;
  switch (fittype) {
  case utility::omegam_omegade : argexpected = 5; break;
  case utility::flat_omegam : argexpected = 4; break;
  case utility::flat_omegam_w : argexpected = 5; break;
  case utility::omegam_omegade_w : argexpected = 6; break;
  case utility::flat_omegam_w0_wa : argexpected = 6; break;
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
  w=-1; wa=0.0; om=1.0; ode=0.0; alpha=1.5; beta=1.5; scriptm=23.5;//Defaults
  switch (fittype) {
  case omegam_omegade :
    om = par[0];
    ode = par[1];
    alpha = par[2];
    beta = par[3];
    scriptm = par[4];
    break;
  case flat_omegam :
    om = par[0];
    ode = 1.0-om;
    alpha = par[1];
    beta = par[2];
    scriptm = par[3];
    break;
  case flat_omegam_w :
    om = par[0];
    ode = 1.0 - om;
    w = par[1];
    alpha = par[2];
    beta = par[3];
    scriptm = par[4];
    break;
  case omegam_omegade_w :
    om = par[0];
    ode = par[1];
    w = par[2];
    alpha = par[3];
    beta = par[4];
    scriptm = par[5];
    break;
  case flat_omegam_w0_wa :
    om = par[0];
    ode = 1.0 - om;
    w = par[1];
    wa = par[2];
    alpha = par[3];
    beta = par[4];
    scriptm = par[5];
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
  \brief Unpacks a parameter vector based on the type of fit being done,
    including multilinear alpha/beta
*/
void utility::parse_param_values_multi( const std::vector<double>& par,
					utility::cosmo_fittype fittype,
					double& om, double& ode, double& w,
					double& wa, double& alpha1, 
					double& alpha2, double& beta1,
					double& beta2, double& scriptm) {

  //Figure out how many arguments we expect
  unsigned int argexpected = 0;
  switch (fittype) {
  case utility::omegam_omegade : argexpected = 7; break;
  case utility::flat_omegam : argexpected = 6; break;
  case utility::flat_omegam_w : argexpected = 7; break;
  case utility::omegam_omegade_w : argexpected = 8; break;
  case utility::flat_omegam_w0_wa : argexpected = 8; break;
  default : throw CosFitterExcept("utility","parse_param_values_multi",
				  "Unknown fittype",1);
    break;
  }
  if (par.size() != argexpected) {
    std::stringstream errstrng;
    errstrng << "Expected " << argexpected << " params, got: "
	     << par.size();
    throw CosFitterExcept("utility","parse_param_values_multi",
			  errstrng.str(),2);
  }

  //Set up lumdist
  w=-1; wa=0.0; om=1.0; ode=0.0; alpha1=1.5; alpha2=1.5;
  beta1=2.0; beta2=2.0; scriptm=23.5;//Defaults
  switch (fittype) {
  case omegam_omegade :
    om = par[0];
    ode = par[1];
    alpha1 = par[2];
    alpha2 = par[3];
    beta1 = par[4];
    beta2 = par[5];
    scriptm = par[6];
    break;
  case flat_omegam :
    om = par[0];
    ode = 1.0-om;
    alpha1 = par[1];
    alpha2 = par[2];
    beta1 = par[3];
    beta2 = par[4];
    scriptm = par[5];
    break;
  case flat_omegam_w :
    om = par[0];
    ode = 1.0 - om;
    w = par[1];
    alpha1 = par[2];
    alpha2 = par[3];
    beta1 = par[4];
    beta2 = par[5];
    scriptm = par[6];
    break;
  case omegam_omegade_w :
    om = par[0];
    ode = par[1];
    w = par[2];
    alpha1 = par[3];
    alpha2 = par[4];
    beta1 = par[5];
    beta2 = par[6];
    scriptm = par[7];
    break;
  case flat_omegam_w0_wa :
    om = par[0];
    w = par[1];
    wa = par[2];
    alpha1 = par[3];
    alpha2 = par[4];
    beta1 = par[5];
    beta2 = par[6];
    scriptm = par[7];
    break;
  default :
    throw CosFitterExcept("utility","parse_param_values_multi",
			  "Unknown fittype",4);
  }
}

/*!
  \brief Unpacks a parameter vector based on the type of fit being done,]
   no \f$\mathcal{M}\f$ version
*/
void 
utility::parse_param_values_multi_nosm( const std::vector<double>& par,
					utility::cosmo_fittype fittype,
					double& om, double& ode, double& w,
					double& wa, double& alpha1, 
					double& alpha2, double& beta1,
					double& beta2) {

  //Figure out how many arguments we expect
  unsigned int argexpected = 0;
  switch (fittype) {
  case utility::omegam_omegade : argexpected = 6; break;
  case utility::flat_omegam : argexpected = 5; break;
  case utility::flat_omegam_w : argexpected = 6; break;
  case utility::omegam_omegade_w : argexpected = 7; break;
  case utility::flat_omegam_w0_wa : argexpected = 7; break;
  default : throw CosFitterExcept("utility","parse_param_values_multi_nosm",
				  "Unknown fittype",1);
    break;
  }
  if (par.size() != argexpected) {
    std::stringstream errstrng;
    errstrng << "Expected " << argexpected << " params, got: "
	     << par.size();
    throw CosFitterExcept("utility","parse_param_values_multi_nosm",
			  errstrng.str(),2);
  }

  //Set up lumdist
  w=-1; wa=0.0; om=1.0; ode=0.0; alpha1=1.5; alpha2=1.5;
  beta1=2.0; beta2=2.0;//Defaults
  switch (fittype) {
  case omegam_omegade :
    om = par[0];
    ode = par[1];
    alpha1 = par[2];
    alpha2 = par[3];
    beta1 = par[4];
    beta2 = par[5];
    break;
  case flat_omegam :
    om = par[0];
    ode = 1.0-om;
    alpha1 = par[1];
    alpha2 = par[2];
    beta1 = par[3];
    beta2 = par[4];
    break;
  case flat_omegam_w :
    om = par[0];
    ode = 1.0 - om;
    w = par[1];
    alpha1 = par[2];
    alpha2 = par[3];
    beta1 = par[4];
    beta2 = par[5];
    break;
  case omegam_omegade_w :
    om = par[0];
    ode = par[1];
    w = par[2];
    alpha1 = par[3];
    alpha2 = par[4];
    beta1 = par[5];
    beta2 = par[6];
    break;
  case flat_omegam_w0_wa :
    om = par[0];
    w = par[1];
    wa = par[2];
    alpha1 = par[3];
    alpha2 = par[4];
    beta1 = par[5];
    beta2 = par[6];
    break;
  default :
    throw CosFitterExcept("utility","parse_param_values_multi_nosm",
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
