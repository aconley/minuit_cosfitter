//utility.h

#ifndef __cosfitter_utility__
#define __cosfitter_utility__

#include <string>
#include <vector>

/*!
  \brief Utility functions
*/
namespace utility {

  enum cosmo_fittype { flat_omegam, omegam_omegade, flat_omegam_w,
		       omegam_omegade_w, flat_omegam_w0_wa }; //!<

  void parse_param_values( const std::vector<double>&, utility::cosmo_fittype,
			   double&, double&, double&,
			   double&, double&, double&, double&); 
  void parse_param_values_nosm( const std::vector<double>&, 
				utility::cosmo_fittype,
				double&, double&, double&,
				double&, double&, double& ); 

  void parse_param_values_multi( const std::vector<double>&, utility::cosmo_fittype,
				 double&, double&, double&, double&, double&,
				 double&, double&, double&, double&); 
  void parse_param_values_multi_nosm( const std::vector<double>&, 
				      utility::cosmo_fittype,
				      double&, double&, double&, double&,
				      double&, double&, double&, double& ); 

  void stringwords(const std::string& ins, std::vector<std::string>& words); //!< Breaks up an input string into a vector


}

#endif
