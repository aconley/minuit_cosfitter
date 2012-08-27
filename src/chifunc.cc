#include<cassert>
#include<cmath>

#include<iostream>
#include<iomanip>
#include<sstream>
#include<cmath>

#include <chifunc.h>
#include <auxconstraint.h>
#include <cosfitterexcept.h>

/////////////////////////////////////////////////////////////////////
//                           chifunc_base                          //
/////////////////////////////////////////////////////////////////////
chifunc_base::chifunc_base( const SNeData& data, 
			    const std::map<unsigned int,double>& idisp, 
			    double pz, utility::cosmo_fittype fit ) : 

  sne( data ), fittype(fit), intrinsicdisp(idisp), pecz(pz), errdef(1.0), 
  isprepped(false), includebao(false), includewmap(false),
  useEisenstein(false) {

  nsn = sne.size();
  dl.resize( nsn );

  //Calculate the luminosity distance using fiducial cosmo
  lmdist.getLumDist( sne, dl, 0.3, 0.7 );
}

void chifunc_base::print(std::ostream& os) const {
  os << "Summary of SN data:" << std::endl;
  if (isprepped) os << " (errors do not include width and colour corrections)"
		    << std::endl;
  for (unsigned int i = 0; i < nsn; ++i) {
    os << std::left << std::setw(8) << sne[i].name << " z: " << std::right <<
       std::setw(6) << std::fixed << std::setprecision(4) << sne[i].zhel << 
      " s: " << std::right << std::setw(5) << std::fixed << 
      std::setprecision(3) << sne[i].widthpar << " c: " << 
      std::fixed << std::setprecision(3) << std::right << 
      std::setw(6) << sne[i].colourpar;
    if (isprepped) os << " var: " << std::setw(6) << std::fixed <<
      std::right << pre_vars[i] << " sigma: " << std::fixed << 
      std::setw(5) << std::right << sqrt(pre_vars[i]);
    os << std::endl;
  }
}

/////////////////////////////////////////////////////////////////////
//                              chifunc                            //
/////////////////////////////////////////////////////////////////////
chifunc::chifunc( const SNeData& data, 
		  const std::map<unsigned int,double>& idisp, 
		  double pz, utility::cosmo_fittype fit ) : 
  chifunc_base( data, idisp, pz, fit ), fixalphaset(false), fixbetaset(false),
  errsfixed(false) { }

/*!
  Updates the inverse covariance matrix

  \param[in] alpha \f$\alpha\f$
  \param[in] beta \f$\beta\f$

  If the requisite \f$\alpha, \beta\f$ covariance matrix
  isn't available, alpha and/or beta are ignored.
*/
void chifunc::computeInvCovMatrix( float alpha, float beta ) const {

  int status;
  float working_alpha, working_beta;

  if (sne.areErrorsDiagonal())
    throw CosFitterExcept("chifunc","computeInvCovMatrix",
			  "No cov matrix present -- errors are diagonal!",1);

  //Now update the diagonal errors
  //We could check to see if they are non-zero, but in reality the
  // cost of this function is so dominated by the inversion it
  // just doesn't matter
  working_vec.resize( nsn );

  if (fixalphaset) working_alpha=fixalphaval; else working_alpha=alpha;
  if (fixbetaset) working_beta=fixbetaval; else working_beta=beta;

  float asq = working_alpha * working_alpha;
  float bsq = working_beta * working_beta;
  float ab = working_alpha * working_beta;

  if (std::isnan(working_alpha)) 
    throw CosFitterExcept("chifunc","computeInvCovMatrix",
			  "Non finite alpha",4);
  if (std::isnan(working_beta)) 
    throw CosFitterExcept("chifunc","computeInvCovMatrix",
			  "Non finite beta",8);

  if ( errsfixed ) {
    //Then the alpha and beta terms were already included
    // in pre_vars by calcPreErrs
    for (unsigned int i = 0; i < nsn; ++i)
      working_vec[i] = pre_vars[i]; //Needed because they are diff types
  } else if ( fixalphaset ) {
    //The alpha term was included, but the beta and cross term wasn't
    for (unsigned int i = 0; i < nsn; ++i)
      working_vec[i] = pre_vars[i] + bsq * sne[i].var_colourpar
        - 2.0 * working_beta  * sne[i].cov_mag_colourpar
        - 2.0 * ab * sne[i].cov_widthpar_colourpar;
  } else if ( fixbetaset ) {
    //Likewise, but for beta
    for (unsigned int i = 0; i < nsn; ++i)
      working_vec[i] = pre_vars[i] + asq * sne[i].var_widthpar
        + 2.0 * working_alpha  * sne[i].cov_mag_widthpar
        - 2.0 * ab * sne[i].cov_widthpar_colourpar;
  } else {
    //Neither included
    for (unsigned int i = 0; i < nsn; ++i)
      working_vec[i] = pre_vars[i] + asq * sne[i].var_widthpar
	+ bsq * sne[i].var_colourpar
        + 2.0 * working_alpha * sne[i].cov_mag_widthpar
        - 2.0 * working_beta  * sne[i].cov_mag_colourpar
        - 2.0 * ab * sne[i].cov_widthpar_colourpar;
  }

  //Build the combined cov matrix and add in the current diagonal
  invcovmatrix = sne.getCombinedCovMatrix(working_alpha,working_beta);
  for (unsigned int i = 0; i < nsn; ++i)
    invcovmatrix[i][i] += working_vec[i];

  //The actual inversion
  status = 0;
  invcovmatrix.invert(status);
  if (status) {
    std::stringstream errstr;
    errstr << "Error inverting cov matrix with alpha: " 
	   << working_alpha << " beta: " << working_beta;
    throw CosFitterExcept("chifunc","computeInvCovMatrix",
			  errstr.str(),4);
  }
}

void chifunc::SetAlphaErrFix( double val ) {
  fixalphaset = true; fixalphaval = val;
}

void chifunc::SetBetaErrFix( double val ) {
  fixbetaset = true; fixbetaval = val;
}

void chifunc::CalcPreErrs( ) {

  bool fixedintrinsic; //!< True if there is only one value for all SNe
  unsigned int nintrinsicdisp;

  double currintrinsicsq;
  nintrinsicdisp = intrinsicdisp.size();
  if ( nintrinsicdisp == 0) {
    //Nothing, this is bad
    throw CosFitterExcept("chifunc","CalcPreErrs","No intrinsic disp set",1);
  } else if ( nintrinsicdisp == 1 ) {
    //Use the same for all even if it doesn't match the data set number
    fixedintrinsic = true;
    std::map<unsigned int,double>::const_iterator it = intrinsicdisp.begin();
    currintrinsicsq = (it->second)*(it->second);
  } else {
    fixedintrinsic = false;
    currintrinsicsq = 0.0;
  }

  isprepped = false;
  errsfixed = false;
  pre_vars.resize( nsn );
  const double zfacsq = (5.0/log(10.0))*(5.0/log(10.0));
  double dzerrsq, emptyfac;
  //This is the part that can be done even if alpha/beta aren't
  // fixed for error propagation
  if (fixedintrinsic) {
    for (unsigned int i = 0; i < nsn; ++i) {
      emptyfac = (1.0+sne[i].zcmb)/(sne[i].zcmb*(1+0.5*sne[i].zcmb));
      dzerrsq = pecz*pecz + sne[i].var_z;
      dzerrsq *= zfacsq*emptyfac*emptyfac;
      pre_vars[i] = sne[i].var_mag + dzerrsq + currintrinsicsq;
    }
  } else {
    std::map< unsigned int, double >::const_iterator it;
    for (unsigned int i = 0; i < nsn; ++i) {
      it = intrinsicdisp.find( sne[i].dataset );
      if (it == intrinsicdisp.end()) {
	std::stringstream errstr;
	errstr << "Unknown dataset number " << sne[i].dataset;
	errstr << " for sn: " << sne[i].name;
	throw CosFitterExcept("chifunc","calcPreErrs",
			      errstr.str(),2);
      }
      currintrinsicsq = it->second * it->second;
      emptyfac = (1.0+sne[i].zcmb)/(sne[i].zcmb*(1+0.5*sne[i].zcmb));
      dzerrsq = pecz*pecz + sne[i].var_z;
      dzerrsq *= zfacsq*emptyfac*emptyfac;
      pre_vars[i] = sne[i].var_mag + dzerrsq + currintrinsicsq;
    }
  }

  //Handle the fixed parts
  errsfixed = fixalphaset && fixbetaset;
  if ( fixalphaset || fixbetaset ) {
    if (errsfixed) {
      double aerrsq = fixalphaval*fixalphaval;    
      double berrsq = fixbetaval*fixbetaval;  
      double aberr = fixalphaval * fixbetaval;
      for (unsigned int i = 0; i < nsn; ++i)
	pre_vars[i] += aerrsq * sne[i].var_widthpar
	  + berrsq * sne[i].var_colourpar
	  + (2.0 * fixalphaval * sne[i].cov_mag_widthpar)
	  - (2.0 * fixbetaval * sne[i].cov_mag_colourpar)
	  - (2.0 * aberr * sne[i].cov_widthpar_colourpar);
    } else if (fixalphaset) {
      double aerrsq = fixalphaval*fixalphaval;    
      for (unsigned int i = 0; i < nsn; ++i)
	pre_vars[i] += aerrsq * sne[i].var_widthpar
	  + (2.0 * fixalphaval * sne[i].cov_mag_widthpar);
    } else if (fixbetaset) {
      double berrsq = fixbetaval*fixbetaval;    
      for (unsigned int i = 0; i < nsn; ++i)
	pre_vars[i] += berrsq * sne[i].var_colourpar
	  - (2.0 * fixbetaval * sne[i].cov_mag_colourpar);
    }
    
    //Pre-invert cov matrix if we can
    if ( errsfixed && ! sne.areErrorsDiagonal() ) 
      computeInvCovMatrix(fixalphaval,fixbetaval); 
  }
  isprepped = true;
}



double chifunc::EstimateScriptm( const std::vector<double>& par ) const {
  if (isprepped != true ) throw CosFitterExcept("chifunc","EstimateScriptm",
						"Not prepped",1);

  double w, wa, om, ode, alpha, beta;
  utility::parse_param_values_nosm( par, fittype, om, ode, w, wa, alpha, 
				    beta );

  int st = lmdist.getLumDist( sne, dl, om, ode, w, wa );
  if (st != 0) return 1e20; //Failure

  double dmag, scriptm, totweight;
  scriptm = 0.0;
  totweight = 0.0;
  if (errsfixed) {
    //Pre-vars has all error info already
    for (unsigned int i = 0; i < nsn; ++i) {
      dmag = sne[i].mag - (dl[i] - alpha*(sne[i].widthpar-1.0)
			   + beta * sne[i].colourpar );
      scriptm += dmag / pre_vars[i];
      totweight += 1.0/pre_vars[i];
    }
  } else {
    double invvar, alphaval, betaval, alphabetaval;
    double alphasq, betasq;

    //Both aren't set, but one could be
    if (fixalphaset) alphaval = fixalphaval; else alphaval=alpha;
    if (fixbetaset) betaval = fixbetaval; else betaval=beta;

    alphabetaval = alphaval * betaval;
    alphasq = alpha*alpha;
    betasq = beta*beta;

    for (unsigned int i = 0; i < nsn; ++i) {
      invvar = 1.0 / (pre_vars[i] 
		      + alphasq * sne[i].var_widthpar
		      + betasq * sne[i].var_colourpar
		      + 2 * alphaval * sne[i].cov_mag_widthpar
		      - 2 * betaval * sne[i].cov_mag_colourpar
		      - 2 * alphabetaval * sne[i].cov_widthpar_colourpar ); 
      dmag = sne[i].mag - (dl[i] - alpha*(sne[i].widthpar-1.0)
			   + beta * sne[i].colourpar );
      scriptm += dmag * invvar;
      totweight += invvar;
    }
  }

  return scriptm/totweight;
}

void chifunc_base::SetUseKomatsuForm() { 
  lmdist.setUseKomatsuForm();
  E05.SetUseKomatsuForm();
  P09.SetUseKomatsuForm();
  WMAP7.SetUseKomatsuForm();
  WMAP7_P09.SetUseKomatsuForm();
} 
void chifunc_base::UnsetUseKomatsuForm() { 
  lmdist.unsetUseKomatsuForm();
  E05.UnsetUseKomatsuForm();
  P09.UnsetUseKomatsuForm();
  WMAP7.UnsetUseKomatsuForm();
  WMAP7_P09.UnsetUseKomatsuForm();
} 

void chifunc_base::SetAtrans(double val) { 
  lmdist.setAtrans(val); 
  E05.SetAtrans(val);
  P09.SetAtrans(val);
  WMAP7.SetAtrans(val);
  WMAP7_P09.SetAtrans(val);
} 

/*!
  \param[in] par Cosmological and nuisance parameters.  Not all
    may be present, but the order should be \f$\Omega_m , \Omega_{DE},
    w , \alpha , \beta\f$.
  \returns The \f$\chi^2\f$
*/
double chifunc::GetRMS(const std::vector<double>& par) const {
  if (isprepped != true ) throw CosFitterExcept("GetRMS","GetRMS",
						"Not prepped",1);

  double w, wa, om, ode, alpha, beta, scriptm;
  utility::parse_param_values( par, fittype, om, ode, w, wa, alpha, beta,
			       scriptm );

  int st = lmdist.getLumDist( sne, dl, om, ode, w, wa );
  if (st != 0) return 1e20; //Failure

  double diffmag, rms;
  rms = 0.0;
  for (unsigned int i = 0; i < nsn; ++i) {
    diffmag = sne[i].mag - (dl[i] - alpha * (sne[i].widthpar-1.0) +
			    + beta * sne[i].colourpar + scriptm);
    rms += diffmag * diffmag;
  }

  return sqrt( rms / static_cast<double>(nsn) );
}

/*!
  \param[in] om \f$\Omega_m
  \param[in] ode \f$\Omega_{DE}\f$
  \param[in] w   \f$w\f$
  \param[in] wa  \f$w_a\f$
  \param[in] alpha \f$\alpha\f$
  \param[in] beta \f$\beta\f$
  \param[in] scriptm \f$\mathcal{M}\f$
  \returns The \f$\chi^2\f$ of the SN fit

  This handles only the case where the SN are uncorrelated
*/
double chifunc::GetDiagonalChisq(double om, double ode, double w, double wa,
				 double alpha, double beta, double scriptm) 
  const {

  double diffmag, chisq = 0.0;
  if (errsfixed) {
    //Errors completely pre-calculated
    for (unsigned int i = 0; i < nsn; ++i) {
      diffmag = sne[i].mag - (dl[i] - alpha * (sne[i].widthpar-1.0) +
			      + beta * sne[i].colourpar + scriptm);
      chisq += diffmag*diffmag/pre_vars[i];
    }
  } else {
    //Here we are updating the errors on every step
    double alphaval, betaval, alphabetaval, alphasq, betasq;

    //Both aren't set, but one could be
    if (fixalphaset) alphaval = fixalphaval; else alphaval=alpha;
    if (fixbetaset) betaval = fixbetaval; else betaval=beta;

    alphabetaval = alphaval * betaval;
    alphasq = alphaval * alphaval;
    betasq = betaval * betaval;

    double invvar;
    if (fixalphaset) {
      //alpha included, beta+cross terms not in pre-vars
      for (unsigned int i = 0; i < nsn; ++i) {
	diffmag = sne[i].mag - (dl[i] - alpha * (sne[i].widthpar-1.0) +
				+ beta * sne[i].colourpar + scriptm);
	invvar = 1.0 / ( pre_vars[i] 
			 + betasq * sne[i].var_colourpar
			 - 2.0 * betaval * sne[i].cov_mag_colourpar
			 - 2.0 * alphabetaval * sne[i].cov_widthpar_colourpar);
	chisq += diffmag*diffmag * invvar;
      }
    } else if (fixbetaset) {
      //Beta included, alpha+cross terms not
      for (unsigned int i = 0; i < nsn; ++i) {
	diffmag = sne[i].mag - (dl[i] - alpha * (sne[i].widthpar-1.0) +
				+ beta * sne[i].colourpar + scriptm);
	invvar = 1.0 / ( pre_vars[i] 
			 + alphasq * sne[i].var_widthpar
			 + 2.0 * alphaval * sne[i].cov_mag_widthpar
			 - 2.0 * alphabetaval * sne[i].cov_widthpar_colourpar);

	chisq += diffmag*diffmag * invvar;
      }
    } else {
      //Neither fixed, none of alpha/beta already in pre-vars
      for (unsigned int i = 0; i < nsn; ++i) {
	diffmag = sne[i].mag - (dl[i] - alpha * (sne[i].widthpar-1.0) +
				+ beta * sne[i].colourpar + scriptm);
	invvar = 1.0 / ( pre_vars[i] 
			 + alphasq * sne[i].var_widthpar
			 + betasq * sne[i].var_colourpar
			 + 2.0 * alphaval * sne[i].cov_mag_widthpar
			 - 2.0 * betaval * sne[i].cov_mag_colourpar
			 - 2.0 * alphabetaval * sne[i].cov_widthpar_colourpar);

	chisq += diffmag*diffmag * invvar;
      }
    }
  }

  /*
  printf("Om: %5.3f Ode: %5.3f W: %6.3f WA: %6.3f alpha: %5.3f beta: %5.3f sm: %6.3f\n",
	 om,ode,w,wa,alpha,beta,scriptm);
  printf("Chisq: %11.3f\n",chisq);
  */

  return chisq;
}

/*!
  \param[in] om \f$\Omega_m
  \param[in] ode \f$\Omega_{DE}\f$
  \param[in] w   \f$w\f$
  \param[in] wa  \f$w_a\f$
  \param[in] alpha \f$\alpha\f$
  \param[in] beta \f$\beta\f$
  \param[in] scriptm \f$\mathcal{M}\f$
  \returns The \f$\chi^2\f$ of the SN fit

  This handles the case where the SN are correlated.
*/
double chifunc::GetChisq(double om, double ode, double w, double wa,
			 double alpha, double beta, double scriptm) const {

  double chisq;
  static std::vector<double> diffmag(nsn);

  for (unsigned int i = 0; i < nsn; ++i)
    diffmag[i] = sne[i].mag - (dl[i] - alpha * (sne[i].widthpar-1.0) +
			       + beta * sne[i].colourpar + scriptm);

  if ( ! errsfixed ) {
    //If errsfixed, this was computed in CalcPreErrs
    computeInvCovMatrix( alpha, beta );
  }

  //\chi^2 = \Delta \vec{m} \mathbf{V}^{-1} \cdot \Delta \vec{m}
  //This requires that the cov matrix product be correct -- i.e.,
  // the symmetry is not used, so the reflection better have happened
  int status = 0;
  working_vec = invcovmatrix.mult( diffmag, status );
  chisq = diffmag[0] * working_vec[0];
  for (unsigned int i = 1; i < nsn; ++i)
    chisq += diffmag[i]*working_vec[i];

  if (status != 0)
    throw CosFitterExcept("chifunc","GetChisq","Error inverting cov matrix",1);

  return chisq;
}

/*!
  \param[in] par Cosmological and nuisance parameters.  Not all
    may be present, but the order should be \f$\Omega_m , \Omega_{DE},
    w , wa, \alpha , \beta\f$.
  \returns The \f$\chi^2\f$
*/
double chifunc::operator()(const std::vector<double>& par) const {
  if (isprepped != true ) throw CosFitterExcept("chifunc","operator()",
						"Not prepped",1);

  double w, wa, om, ode, alpha, beta, scriptm;
  utility::parse_param_values( par, fittype, om, ode, w, wa, alpha, beta,
			       scriptm );

  int st = lmdist.getLumDist( sne, dl, om, ode, w, wa );
  if (st != 0) return 1e20; //Failure

  if (std::isnan(alpha)) {
    std::cerr << "Warning -- non finite alpha.  Results probably unreliable"
	      << std::endl;
    return 1e20;
  }
  if (std::isnan(beta)) {
    std::cerr << "Warning -- non finite beta.  Results probably unreliable"
	      << std::endl;
    return 1e20;
  }

  double chisq;
  if (sne.areErrorsDiagonal()) 
    chisq = GetDiagonalChisq(om,ode,w,wa,alpha,beta,scriptm); 
  else 
    chisq = GetChisq(om,ode,w,wa,alpha,beta,scriptm);
  
  if ( (includebao && includewmap) && (!useEisenstein) ) {
    //Altogether
    WMAP7_P09.setFittype(fittype);
    chisq += WMAP7_P09( om, ode, w, wa );
  } else {
    if (includebao) {
      if (useEisenstein) {
	E05.setFittype(fittype);
	chisq += E05( om, ode, w, wa );
      } else {
	P09.setFittype(fittype);
	chisq += P09( om, ode, w, wa );
      }
    }
    if (includewmap) {
      WMAP7.setFittype(fittype);
      chisq += WMAP7( om, ode, w, wa );
    }
  }

  return chisq;

}


void chifunc::print(std::ostream& os, const std::vector<double>& par) const {

  double w, wa, om, ode, alpha, beta, scriptm;
  utility::parse_param_values( par, fittype, om, ode, w, wa, alpha, beta,
			       scriptm );

  int st = lmdist.getLumDist( sne, dl, om, ode, w, wa );
  if (st != 0) {
    os << "Unable to compute cosmology" << std::endl;
    return;
  }

  if (!isprepped) {
    os << "Data not prepared for full output" << std::endl;
    return;
  }

  double diffmag, err, corrmag;
  os << "Summary of SN data:" << std::endl;
  if (! sne.areErrorsDiagonal() )
    os << " Output ignores correlations between SN" << std::endl;
  os << std::left << std::setw(8) << "Name" << "  " << std::setw(6) << "zcmb"
     << "  " << std::setw(7) << "mag" << "  " << std::setw(8) << "corr mag"
     << "  " << std::setw(7) << "lumdist" << "  " << std::setw(7) << "err"
     << "  " << std::setw(7) << "diff" << std::endl;
  for (unsigned int i = 0; i < nsn; ++i) {
    if (errsfixed) {
      //Errors completely pre-calculated
      corrmag = sne[i].mag + alpha*(sne[i].widthpar-1.0)
	- beta*sne[i].colourpar;
      diffmag = corrmag - dl[i] - scriptm;
      err = pre_vars[i];
    } else {
      //Here we are updating the errors on every step
      double alphaval, betaval, alphabetaval, alphasq, betasq;
      
      //Both aren't set, but one could be
      if (fixalphaset) alphaval = fixalphaval; else alphaval=alpha;
      if (fixbetaset) betaval = fixbetaval; else betaval=beta;
      
      alphabetaval = alphaval * betaval;
      alphasq = alphaval * alphaval;
      betasq = betaval * betaval;

      if (fixalphaset) {
	//alpha included, beta+cross terms not in pre-vars
	corrmag = sne[i].mag + alpha*(sne[i].widthpar-1.0)
	  - beta*sne[i].colourpar;
	diffmag = corrmag - dl[i] - scriptm;
	err = sqrt( pre_vars[i] 
		    + betasq * sne[i].var_colourpar
		    - 2.0 * betaval * sne[i].cov_mag_colourpar
		    - 2.0 * alphabetaval * sne[i].cov_widthpar_colourpar);
      } else if (fixbetaset) {
	//Beta included, alpha+cross terms not
	corrmag = sne[i].mag + alpha*(sne[i].widthpar-1.0)
	  - beta*sne[i].colourpar;
	diffmag = corrmag - dl[i] - scriptm;
	err = sqrt( pre_vars[i] 
		    + alphasq * sne[i].var_widthpar
		    + 2.0 * alphaval * sne[i].cov_mag_widthpar
		    - 2.0 * alphabetaval * sne[i].cov_widthpar_colourpar);
      } else {
	//Neither fixed, none of alpha/beta already in pre-vars
	corrmag = sne[i].mag + alpha*(sne[i].widthpar-1.0)
	  - beta*sne[i].colourpar;
	diffmag = corrmag - dl[i] - scriptm;
	err = sqrt( pre_vars[i] 
		    + alphasq * sne[i].var_widthpar
		    + betasq * sne[i].var_colourpar
		    + 2.0 * alphaval * sne[i].cov_mag_widthpar
		    - 2.0 * betaval * sne[i].cov_mag_colourpar
		    - 2.0 * alphabetaval * sne[i].cov_widthpar_colourpar);
      }
    }
    
    os << std::left << std::setw(8) << sne[i].name << std::fixed << std::right
       << "  " << std::setprecision(4) << std::setw(6) << sne[i].zcmb
       << "  " << std::setprecision(4) << std::setw(7) << sne[i].mag
       << "  " << std::setprecision(4) << std::setw(8) << corrmag
       << "  " << std::setprecision(4) << std::setw(7) << dl[i]+scriptm
       << "  " << std::setprecision(4) << std::setw(7) << diffmag
       << "  " << std::setprecision(4) << std::setw(7) << err
       << std::endl;
  }
}


////////////////////////////////////////////////////////////////////////
//                      chifunc_multilin                              //
////////////////////////////////////////////////////////////////////////
chifunc_multilin::chifunc_multilin( const SNeData& data, 
				    const std::map<unsigned int,double>& idisp, 
				    double pz, utility::cosmo_fittype fit ) : 
  chifunc_base( data, idisp, pz, fit ), widthcut(0.8), colourcut(1.0) {} 

void chifunc_multilin::CalcPreErrs( ) {

  bool fixedintrinsic; //!< True if there is only one value for all SNe
  unsigned int nintrinsicdisp;

  double currintrinsicsq;
  nintrinsicdisp = intrinsicdisp.size();
  if ( nintrinsicdisp == 0) {
    //Nothing, this is bad
    throw CosFitterExcept("chifunc","CalcPreErrs","No intrinsic disp set",1);
  } else if ( nintrinsicdisp == 1 ) {
    //Use the same for all even if it doesn't match the data set number
    fixedintrinsic = true;
    std::map<unsigned int,double>::const_iterator it = intrinsicdisp.begin();
    currintrinsicsq = (it->second)*(it->second);
  } else {
    fixedintrinsic = false;
    currintrinsicsq = 0.0;
  }

  isprepped = false;
  pre_vars.resize( nsn );
  const double zfacsq = (5.0/log(10.0))*(5.0/log(10.0));
  double dzerrsq, emptyfac;
  if (fixedintrinsic) {
    for (unsigned int i = 0; i < nsn; ++i) {
      emptyfac = (1.0+sne[i].zcmb)/(sne[i].zcmb*(1+0.5*sne[i].zcmb));
      dzerrsq = pecz*pecz + sne[i].var_z;
      dzerrsq *= zfacsq*emptyfac*emptyfac;
      pre_vars[i] = sne[i].var_mag + dzerrsq + currintrinsicsq;
    }
  } else {
    std::map< unsigned int, double >::const_iterator it;
    for (unsigned int i = 0; i < nsn; ++i) {
      it = intrinsicdisp.find( sne[i].dataset );
      if (it == intrinsicdisp.end()) {
	std::stringstream errstr;
	errstr << "Unknown dataset number " << sne[i].dataset;
	errstr << " for sn: " << sne[i].name;
	throw CosFitterExcept("chifunc","calcPreErrs",
			      errstr.str(),2);
      }
      currintrinsicsq = it->second * it->second;
      emptyfac = (1.0+sne[i].zcmb)/(sne[i].zcmb*(1+0.5*sne[i].zcmb));
      dzerrsq = pecz*pecz + sne[i].var_z;
      dzerrsq *= zfacsq*emptyfac*emptyfac;
      pre_vars[i] = sne[i].var_mag + dzerrsq + currintrinsicsq;
    }
  }

  isprepped = true;
}

double chifunc_multilin::EstimateScriptm( const std::vector<double>& par ) 
  const {
  if (isprepped != true ) throw CosFitterExcept("chifunc_multilin",
						"EstimateScriptm",
						"Not prepped",1);
  double w, wa, om, ode, alpha1, alpha2, beta1, beta2;
  utility::parse_param_values_multi_nosm( par, fittype, om, ode, w, wa, alpha1,
					  alpha2, beta1, beta2 );

  int st = lmdist.getLumDist( sne, dl, om, ode, w, wa );
  if (st != 0) return 1e20; //Failure

  double dmag, scriptm, totweight;
  scriptm = 0.0;
  totweight = 0.0;

  double invvar, alphaval, betaval, alphabetaval;
  double widthconst, colourconst, alphasq, betasq;

  widthconst = (alpha1 - alpha2)*(widthcut - 1);
  colourconst = (beta1 - beta2)*colourcut;

  for (unsigned int i = 0; i < nsn; ++i) {
    if ( sne[i].widthpar < widthcut ) alphaval = alpha1; else
      alphaval = alpha2;
    if ( sne[i].colourpar < colourcut ) betaval = beta1; else
      betaval = beta2;
    alphabetaval = alphaval * betaval;
    alphasq = alphaval*alphaval;
    betasq = betaval*betaval;
    invvar = 1.0 / (pre_vars[i] 
		    + alphasq * sne[i].var_widthpar
		    + betasq * sne[i].var_colourpar
		    + 2 * alphaval * sne[i].cov_mag_widthpar
		    - 2 * betaval * sne[i].cov_mag_colourpar
		    - 2 * alphabetaval * sne[i].cov_widthpar_colourpar ); 
    dmag = sne[i].mag - (dl[i] - alphaval*(sne[i].widthpar-1.0)
			 + betaval * sne[i].colourpar );
    if (sne[i].widthpar >= widthcut ) dmag += widthconst;
    if (sne[i].colourpar >= colourcut ) dmag += colourconst;
    scriptm += dmag * invvar;
    totweight += invvar;
  }

  return scriptm/totweight;
}


/*!
  \param[in] par Cosmological and nuisance parameters.  Not all
    may be present, but the order should be \f$\Omega_m , \Omega_{DE},
    w , \alpha_1 , \alpha_2, \beta_1, \beta_2\f$.
  \returns The \f$\chi^2\f$
*/
double chifunc_multilin::GetRMS(const std::vector<double>& par) const {
  if (isprepped != true ) throw CosFitterExcept("GetRMS","GetRMS",
						"Not prepped",1);

  double w, wa, om, ode, alpha1, alpha2, beta1, beta2, scriptm;
  utility::parse_param_values_multi( par, fittype, om, ode, w, wa, alpha1, 
				     alpha2, beta1, beta2, scriptm );

  int st = lmdist.getLumDist( sne, dl, om, ode, w, wa );
  if (st != 0) return 1e20; //Failure

  double widthconst, colourconst;
  widthconst = (alpha1 - alpha2)*(widthcut - 1);
  colourconst = (beta1 - beta2)*colourcut;

  double diffmag, rms;
  rms = 0.0;
  for (unsigned int i = 0; i < nsn; ++i) {
    diffmag = sne[i].mag - dl[i] - scriptm;
    if (sne[i].widthpar < widthcut) 
      diffmag += alpha1*(sne[i].widthpar-1.0); else 
      diffmag += alpha2*(sne[i].widthpar-1.0) + widthconst;
    if (sne[i].colourpar < colourcut) 
      diffmag -= beta1*sne[i].colourpar; else
      diffmag -= beta2*sne[i].colourpar + colourconst;
    rms += diffmag * diffmag;
  }
  return sqrt( rms / static_cast<double>(nsn) );
}

/*!
  \param[in] om \f$\Omega_m
  \param[in] ode \f$\Omega_{DE}\f$
  \param[in] w   \f$w\f$
  \param[in] wa  \f$w_a\f$
  \param[in] alpha \f$\alpha\f$
  \param[in] beta \f$\beta\f$
  \param[in] scriptm \f$\mathcal{M}\f$
  \returns The \f$\chi^2\f$ of the SN fit

*/
double chifunc_multilin::GetChisq(double om, double ode, double w, double wa,
				  double alpha1, double alpha2, 
				  double beta1, double beta2, double scriptm) 
  const {

  double diffmag, chisq = 0.0;

  //Here we are updating the errors on every step
  double alphaval, betaval, alphabetaval, alphasq, betasq, invvar;
  double widthconst, colourconst;

  widthconst = (alpha1 - alpha2)*(widthcut - 1);
  colourconst = (beta1 - beta2)*colourcut;

  for (unsigned int i = 0; i < nsn; ++i) {
    diffmag = sne[i].mag - dl[i] - scriptm;
    if (sne[i].widthpar < widthcut) 
      diffmag += alpha1*(sne[i].widthpar-1.0); else 
      diffmag += alpha2*(sne[i].widthpar-1.0) + widthconst;
    if (sne[i].colourpar < colourcut) 
      diffmag -= beta1*sne[i].colourpar; else
      diffmag -= beta2*sne[i].colourpar + colourconst;
    
    if ( sne[i].widthpar < widthcut ) alphaval = alpha1; else
      alphaval = alpha2;
    if ( sne[i].colourpar < colourcut ) betaval = beta1; else
      betaval = beta2;

    alphabetaval = alphaval * betaval;
    alphasq = alphaval * alphaval;
    betasq = betaval * betaval;
    
    invvar = 1.0 / ( pre_vars[i] 
		     + alphasq * sne[i].var_widthpar
		     + betasq * sne[i].var_colourpar
		     + 2.0 * alphaval * sne[i].cov_mag_widthpar
		     - 2.0 * betaval * sne[i].cov_mag_colourpar
		     - 2.0 * alphabetaval * sne[i].cov_widthpar_colourpar);
    
    chisq += diffmag*diffmag * invvar;
  }
  return chisq;
}

/*!
  \param[in] par Cosmological and nuisance parameters.  Not all
    may be present, but the order should be \f$\Omega_m , \Omega_{DE},
    w , \alpha_1 , \alpha_2 , \beta_1 , \beta_2 \f$.
  \returns The \f$\chi^2\f$
*/
double chifunc_multilin::operator()(const std::vector<double>& par) const {
  if (isprepped != true ) throw CosFitterExcept("chifunc_multilin",
						"operator()",
						"Not prepped",1);

  double w, wa, om, ode, alpha1, alpha2, beta1, beta2, scriptm;
  utility::parse_param_values_multi( par, fittype, om, ode, w, wa, alpha1, 
				     alpha2, beta1, beta2, scriptm );
  int st = lmdist.getLumDist( sne, dl, om, ode, w, wa );
  if (st != 0) return 1e20; //Failure

  if (std::isnan(alpha1) || std::isnan(alpha2)) {
    std::cerr << "Warning -- non finite alpha.  Results probably unreliable"
	      << std::endl;
    return 1e20;
  }
  if (std::isnan(beta1) || std::isnan(beta2) ) {
    std::cerr << "Warning -- non finite beta.  Results probably unreliable"
	      << std::endl;
    return 1e20;
  }

  double chisq;
  chisq = GetChisq(om,ode,w,wa,alpha1,alpha2,beta1,beta2,scriptm);


  if ( (includebao && includewmap) && (!useEisenstein) ) {
    //Altogether
    WMAP7_P09.setFittype(fittype);
    chisq += WMAP7_P09( om, ode, w, wa );
  } else {
    if (includebao) {
      if (useEisenstein) {
	E05.setFittype(fittype);
	chisq += E05( om, ode, w, wa );
      } else {
	P09.setFittype(fittype);
	chisq += P09( om, ode, w, wa );
      }
    }
    if (includewmap) {
	WMAP7.setFittype(fittype);
	chisq += WMAP7( om, ode, w, wa );
    }
  }

  return chisq;
}

void chifunc_multilin::print(std::ostream& os, 
			     const std::vector<double>& par) const {
  os << "Extended output not supported for chifunc_multilin" << std::endl;
}
